#include "numerics/schemes.hpp"

#include <algorithm>
#include <cmath>

#include "numerics/face_data.hpp"  // kNConv
#include "numerics/llf.hpp"        // minmod_flux_limiter

namespace ns::numerics {
namespace {
using physics::dvi::P;
using physics::dvi::RHO;
using physics::dvi::U;
using physics::dvi::V;

// Legacy Max3/Min3 semantics (>= / <= chains).
constexpr double max3(double a, double b, double c) {
    if (a >= b && a >= c) return a;
    if (b >= c && b >= a) return b;
    return c;
}
constexpr double min3(double a, double b, double c) {
    if (a <= b && a <= c) return a;
    if (b <= c && b <= a) return b;
    return c;
}

}  // namespace

double van_albada_limiter(double del_plus, double del_minus) {
    const double r = del_plus / (del_minus + 1e-15);
    return (r * r + r) / (r * r + 1.0);
}

double max_eigenvalue_nkfds(double ur, double ul, double vr, double vl, double ar,
                            double al, double nx, double ny) {
    const double vn_r = std::fabs(ur * nx + vr * ny);
    const double vn_l = std::fabs(ul * nx + vl * ny);
    // Legacy L1 == L2 (both plain |vn|); L3 subtracts a.
    const double l3r = vn_r - ar, l3l = vn_l - al;
    return std::max(max3(vn_l, vn_l, l3l), max3(vn_r, vn_r, l3r));
}

double min_eigenvalue_nkfds(double ur, double ul, double vr, double vl, double ar,
                            double al, double nx, double ny) {
    const double vn_r = std::fabs(ur * nx + vr * ny);
    const double vn_l = std::fabs(ul * nx + vl * ny);
    const double l3r = vn_r - ar, l3l = vn_l - al;
    return std::max(min3(vn_l, vn_l, l3l), min3(vn_r, vn_r, l3r));
}

double round_to(double dec, double n) {
    const double sp = std::round(std::pow(10.0, n) * dec);
    return sp / std::pow(10.0, n);
}

double movers_wave_speed(double fr, double fl, double ur, double ul, double l_max,
                         double l_min) {
    constexpr double epsilon = 1e-16;
    const double delta_f = round_to(fr - fl, 3.0);
    const double delta_u = round_to(ur - ul, 3.0);

    double s = 0.0;
    if (std::fabs(delta_f) < epsilon)
        return 0.0;
    else if (std::fabs(delta_u) < epsilon)
        s = l_min;
    else if (std::fabs(delta_u) > epsilon && std::fabs(delta_f) > epsilon)
        s = std::fabs(delta_f / delta_u);
    else
        s = l_min;

    if (std::fabs(s) >= l_max) return l_max;
    if (std::fabs(s) <= l_min) return l_min;
    return s;
}

double nkfds_entropy_switch(double rho_l, double rho_r, double un_l, double un_r,
                            double p_l, double p_r, double ric, double e_max,
                            double cp) {
    const double gamma = 1.4;  // legacy reads global Gamma here
    const double cv = cp / gamma;
    const double rr = cp - cv;
    const double t_l = p_l / rho_l / rr;
    const double t_r = p_r / rho_r / rr;

    const double dd =
        (rho_r - rho_l) * std::log((rho_r / rho_l) * std::pow((t_l / t_r), 2.5)) +
        rho_l / (2 * rr * t_l) * ((rho_r / rho_l) + (t_l / t_r)) * (un_l - un_r) *
            (un_l - un_r) +
        2.5 * (rho_l * (t_l - t_r) / t_r + rho_r * (t_r - t_l) / t_l);
    const double s_l = cv * std::log(p_l / (rho_l * (gamma - 1))) - rr * std::log(rho_l);
    const double s_r = cv * std::log(p_r / (rho_r * (gamma - 1))) - rr * std::log(rho_r);
    const double del_s = s_l - s_r;

    if ((dd > 0) && (std::fabs(del_s) <= 1 * e_max)) return 1.0 * (-ric);
    return 0;
}

namespace {

/// emax window and per-face shared prologue data for the NKFDS-MOVERS family.
struct MoversCommon {
    double e_max = 0;

    explicit MoversCommon(const StructuredMesh& mesh, const MultiField<double>& dv,
                          const physics::IdealGas& eos) {
        const int nx = mesh.nx(), ny = mesh.ny();
        const double gam1 = eos.gamma - 1.0;
        // legacy window: i in [2,id1] (= [0,nx]) x j in [2,jb] (= [0,ny))
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i <= nx; ++i) {
                const double temp =
                    dv(P, i, j) / gam1 +
                    0.5 * dv(RHO, i, j) *
                        (dv(U, i, j) * dv(U, i, j) + dv(V, i, j) * dv(V, i, j));
                if (temp > e_max) e_max = temp;
            }
    }
};

}  // namespace

void nkfds_movers_dissipation_first_order(const StructuredMesh& mesh,
                                          const MeshMetrics& metrics,
                                          const MultiField<double>& dv,
                                          const physics::IdealGas& eos, double cp,
                                          MultiField<double>& diss) {
    const int nx = mesh.nx(), ny = mesh.ny();
    const double gamma = eos.gamma;
    const double gam1 = gamma - 1.0;
    const double ggm1 = gamma / gam1;
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;
    const MoversCommon common(mesh, dv, eos);

    double fr[4], fl[4], ur[4], ul[4];

    // ---- i-direction -------------------------------------------------------
    for (int j = 0; j < ny; ++j) {
        for (int f = 0; f <= nx; ++f) {
            const int im1 = f - 1;
            const double ds = std::sqrt(si(f, j).x * si(f, j).x + si(f, j).y * si(f, j).y);
            const double nzn = si(f, j).x / ds;
            const double nyn = si(f, j).y / ds;

            const double rhol = dv(RHO, im1, j), ulv = dv(U, im1, j),
                         vlv = dv(V, im1, j), pl = dv(P, im1, j);
            const double al = std::sqrt(gamma * pl / rhol);
            const double hl = ggm1 * pl / rhol + 0.5 * (ulv * ulv + vlv * vlv);
            const double el = pl / gam1 + 0.5 * rhol * (ulv * ulv + vlv * vlv);

            const double rhor = dv(RHO, f, j), urv = dv(U, f, j),
                         vrv = dv(V, f, j), pr = dv(P, f, j);
            const double ar = std::sqrt(gamma * pr / rhor);
            const double hr = ggm1 * pr / rhor + 0.5 * (urv * urv + vrv * vrv);
            const double er = pr / gam1 + 0.5 * rhor * (urv * urv + vrv * vrv);

            const double max_eig =
                max_eigenvalue_nkfds(urv, ulv, vrv, vlv, ar, al, nzn, nyn);
            const double min_eig =
                min_eigenvalue_nkfds(urv, ulv, vrv, vlv, ar, al, nzn, nyn);

            fl[0] = (rhol * ulv) * nzn + (rhol * vlv) * nyn;
            fl[1] = (rhol * ulv * ulv + pl) * nzn + (rhol * ulv * vlv) * nyn;
            fl[2] = (rhol * ulv * vlv) * nzn + (rhol * vlv * vlv + pl) * nyn;
            fl[3] = (rhol * ulv * hl) * nzn + (rhol * vlv * hl) * nyn;

            fr[0] = (rhor * urv) * nzn + (rhor * vrv) * nyn;
            fr[1] = (rhor * urv * urv + pr) * nzn + (rhor * urv * vrv) * nyn;
            fr[2] = (rhor * urv * vrv) * nzn + (rhor * vrv * vrv + pr) * nyn;
            fr[3] = (rhor * urv * hr) * nzn + (rhor * vrv * hr) * nyn;

            ur[0] = rhor; ur[1] = rhor * urv; ur[2] = rhor * vrv; ur[3] = er;
            ul[0] = rhol; ul[1] = rhol * ulv; ul[2] = rhol * vlv; ul[3] = el;

            const double fed =
                nkfds_entropy_switch(rhol, rhor, ulv, urv, pl, pr, max_eig,
                                     common.e_max, cp);
            for (int k = 0; k < kNConv; ++k) {
                double fd = 0.5 * movers_wave_speed(fr[k], fl[k], ur[k], ul[k], max_eig,
                                                    min_eig);
                if (fed > 0) fd = fed;
                fd *= (ur[k] - ul[k]);

                diss(k, f, j) -= fd * ds;
                diss(k, im1, j) += fd * ds;
            }
        }
    }

    // ---- j-direction -------------------------------------------------------
    for (int i = 0; i < nx; ++i) {
        for (int f = 0; f <= ny; ++f) {
            const int jm1 = f - 1;
            const double ds = std::sqrt(sj(i, f).x * sj(i, f).x + sj(i, f).y * sj(i, f).y);
            const double nzn = sj(i, f).x / ds;
            const double nyn = sj(i, f).y / ds;

            const double rhol = dv(RHO, i, jm1), ulv = dv(U, i, jm1),
                         vlv = dv(V, i, jm1), pl = dv(P, i, jm1);
            const double al = std::sqrt(gamma * pl / rhol);
            const double hl = ggm1 * pl / rhol + 0.5 * (ulv * ulv + vlv * vlv);
            const double el = pl / gam1 + 0.5 * rhol * (ulv * ulv + vlv * vlv);

            const double rhor = dv(RHO, i, f), urv = dv(U, i, f),
                         vrv = dv(V, i, f), pr = dv(P, i, f);
            const double ar = std::sqrt(gamma * pr / rhor);
            const double hr = ggm1 * pr / rhor + 0.5 * (urv * urv + vrv * vrv);
            const double er = pr / gam1 + 0.5 * rhor * (urv * urv + vrv * vrv);

            const double max_eig =
                max_eigenvalue_nkfds(urv, ulv, vrv, vlv, ar, al, nzn, nyn);
            const double min_eig =
                min_eigenvalue_nkfds(urv, ulv, vrv, vlv, ar, al, nzn, nyn);

            fl[0] = (rhol * ulv) * nzn + (rhol * vlv) * nyn;
            fl[1] = (rhol * ulv * ulv + pl) * nzn + (rhol * ulv * vlv) * nyn;
            fl[2] = (rhol * ulv * vlv) * nzn + (rhol * vlv * vlv + pl) * nyn;
            fl[3] = (rhol * ulv * hl) * nzn + (rhol * vlv * hl) * nyn;

            fr[0] = (rhor * urv) * nzn + (rhor * vrv) * nyn;
            fr[1] = (rhor * urv * urv + pr) * nzn + (rhor * urv * vrv) * nyn;
            fr[2] = (rhor * urv * vrv) * nzn + (rhor * vrv * vrv + pr) * nyn;
            fr[3] = (rhor * urv * hr) * nzn + (rhor * vrv * hr) * nyn;

            ur[0] = rhor; ur[1] = rhor * urv; ur[2] = rhor * vrv; ur[3] = er;
            ul[0] = rhol; ul[1] = rhol * ulv; ul[2] = rhol * vlv; ul[3] = el;

            const double fed =
                nkfds_entropy_switch(rhol, rhor, ulv, urv, pl, pr, max_eig,
                                     common.e_max, cp);
            for (int k = 0; k < kNConv; ++k) {
                double fd = 0.5 * movers_wave_speed(fr[k], fl[k], ur[k], ul[k], max_eig,
                                                    min_eig);
                if (fed > 0) fd = fed;
                fd *= (ur[k] - ul[k]);

                diss(k, i, f) -= fd * ds;
                diss(k, i, jm1) += fd * ds;
            }
        }
    }
}

void nkfds_movers_dissipation_second_order(const StructuredMesh& mesh,
                                           const MeshMetrics& metrics,
                                           const MultiField<double>& dv,
                                           const MultiField<double>& dui,
                                           const MultiField<double>& duj,
                                           const physics::IdealGas& eos, double cp,
                                           MultiField<double>& diss) {
    const int nx = mesh.nx(), ny = mesh.ny();
    const double gamma = eos.gamma;
    const double gam1 = gamma - 1.0;
    const double ggm1 = gamma / gam1;
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;
    const MoversCommon common(mesh, dv, eos);

    double fr[4], fl[4], ur[4], ul[4];

    // ---- i-direction -------------------------------------------------------
    for (int j = 0; j < ny; ++j) {
        for (int f = 0; f <= nx; ++f) {
            const int im1 = f - 1, ip1 = f + 1;
            const double ds = std::sqrt(si(f, j).x * si(f, j).x + si(f, j).y * si(f, j).y);
            const double nzn = si(f, j).x / ds;
            const double nyn = si(f, j).y / ds;

            double deltl[4], deltr[4];
            for (int k = 0; k < 4; ++k) {
                deltl[k] =
                    0.5 * dui(k, im1, j) * minmod_flux_limiter(dui(k, f, j), dui(k, im1, j));
                deltr[k] =
                    0.5 * dui(k, f, j) * minmod_flux_limiter(dui(k, ip1, j), dui(k, f, j));
            }

            const double rhol = dv(RHO, im1, j) + deltl[0];
            const double ulv = dv(U, im1, j) + deltl[1];
            const double vlv = dv(V, im1, j) + deltl[2];
            const double pl = dv(P, im1, j) + deltl[3];
            const double al = std::sqrt(gamma * pl / rhol);
            const double hl = ggm1 * pl / rhol + 0.5 * (ulv * ulv + vlv * vlv);
            const double el = pl / gam1 + 0.5 * rhol * (ulv * ulv + vlv * vlv);

            const double rhor = dv(RHO, f, j) - deltr[0];
            const double urv = dv(U, f, j) - deltr[1];
            const double vrv = dv(V, f, j) - deltr[2];
            const double pr = dv(P, f, j) - deltr[3];
            const double ar = std::sqrt(gamma * pr / rhor);
            const double hr = ggm1 * pr / rhor + 0.5 * (urv * urv + vrv * vrv);
            const double er = pr / gam1 + 0.5 * rhor * (urv * urv + vrv * vrv);

            const double max_eig =
                max_eigenvalue_nkfds(urv, ulv, vrv, vlv, ar, al, nzn, nyn);
            const double min_eig =
                min_eigenvalue_nkfds(urv, ulv, vrv, vlv, ar, al, nzn, nyn);

            fl[0] = (rhol * ulv) * nzn + (rhol * vlv) * nyn;
            fl[1] = (rhol * ulv * ulv + pl) * nzn + (rhol * ulv * vlv) * nyn;
            fl[2] = (rhol * ulv * vlv) * nzn + (rhol * vlv * vlv + pl) * nyn;
            fl[3] = (rhol * ulv * hl) * nzn + (rhol * vlv * hl) * nyn;

            fr[0] = (rhor * urv) * nzn + (rhor * vrv) * nyn;
            fr[1] = (rhor * urv * urv + pr) * nzn + (rhor * urv * vrv) * nyn;
            fr[2] = (rhor * urv * vrv) * nzn + (rhor * vrv * vrv + pr) * nyn;
            fr[3] = (rhor * urv * hr) * nzn + (rhor * vrv * hr) * nyn;

            ur[0] = rhor; ur[1] = rhor * urv; ur[2] = rhor * vrv; ur[3] = er;
            ul[0] = rhol; ul[1] = rhol * ulv; ul[2] = rhol * vlv; ul[3] = el;

            const double fed =
                nkfds_entropy_switch(rhol, rhor, ulv, urv, pl, pr, max_eig,
                                     common.e_max, cp);
            for (int k = 0; k < kNConv; ++k) {
                double fd;
                if (fed > 0)
                    fd = 0.5 * fed * (ur[k] - ul[k]);
                else
                    fd = 0.5 * movers_wave_speed(fr[k], fl[k], ur[k], ul[k], max_eig,
                                                 min_eig) * (ur[k] - ul[k]);

                diss(k, f, j) -= fd * ds;
                diss(k, im1, j) += fd * ds;
            }
        }
    }

    // ---- j-direction -------------------------------------------------------
    for (int i = 0; i < nx; ++i) {
        for (int f = 0; f <= ny; ++f) {
            const int jm1 = f - 1, jp1 = f + 1;
            const double ds = std::sqrt(sj(i, f).x * sj(i, f).x + sj(i, f).y * sj(i, f).y);
            const double nzn = sj(i, f).x / ds;
            const double nyn = sj(i, f).y / ds;

            double deltl[4], deltr[4];
            for (int k = 0; k < 4; ++k) {
                deltl[k] =
                    0.5 * duj(k, i, jm1) * minmod_flux_limiter(duj(k, i, f), duj(k, i, jm1));
                deltr[k] =
                    0.5 * duj(k, i, f) * minmod_flux_limiter(duj(k, i, jp1), duj(k, i, f));
            }

            const double rhol = dv(RHO, i, jm1) + deltl[0];
            const double ulv = dv(U, i, jm1) + deltl[1];
            const double vlv = dv(V, i, jm1) + deltl[2];
            const double pl = dv(P, i, jm1) + deltl[3];
            const double al = std::sqrt(gamma * pl / rhol);
            const double hl = ggm1 * pl / rhol + 0.5 * (ulv * ulv + vlv * vlv);
            const double el = pl / gam1 + 0.5 * rhol * (ulv * ulv + vlv * vlv);

            const double rhor = dv(RHO, i, f) - deltr[0];
            const double urv = dv(U, i, f) - deltr[1];
            const double vrv = dv(V, i, f) - deltr[2];
            const double pr = dv(P, i, f) - deltr[3];
            const double ar = std::sqrt(gamma * pr / rhor);
            const double hr = ggm1 * pr / rhor + 0.5 * (urv * urv + vrv * vrv);
            const double er = pr / gam1 + 0.5 * rhor * (urv * urv + vrv * vrv);

            const double max_eig =
                max_eigenvalue_nkfds(urv, ulv, vrv, vlv, ar, al, nzn, nyn);
            const double min_eig =
                min_eigenvalue_nkfds(urv, ulv, vrv, vlv, ar, al, nzn, nyn);

            fl[0] = (rhol * ulv) * nzn + (rhol * vlv) * nyn;
            fl[1] = (rhol * ulv * ulv + pl) * nzn + (rhol * ulv * vlv) * nyn;
            fl[2] = (rhol * ulv * vlv) * nzn + (rhol * vlv * vlv + pl) * nyn;
            fl[3] = (rhol * ulv * hl) * nzn + (rhol * vlv * hl) * nyn;

            fr[0] = (rhor * urv) * nzn + (rhor * vrv) * nyn;
            fr[1] = (rhor * urv * urv + pr) * nzn + (rhor * urv * vrv) * nyn;
            fr[2] = (rhor * urv * vrv) * nzn + (rhor * vrv * vrv + pr) * nyn;
            fr[3] = (rhor * urv * hr) * nzn + (rhor * vrv * hr) * nyn;

            ur[0] = rhor; ur[1] = rhor * urv; ur[2] = rhor * vrv; ur[3] = er;
            ul[0] = rhol; ul[1] = rhol * ulv; ul[2] = rhol * vlv; ul[3] = el;

            const double fed =
                nkfds_entropy_switch(rhol, rhor, ulv, urv, pl, pr, max_eig,
                                     common.e_max, cp);
            for (int k = 0; k < kNConv; ++k) {
                double fd;
                if (fed > 0)
                    fd = 0.5 * fed * (ur[k] - ul[k]);
                else
                    fd = 0.5 * movers_wave_speed(fr[k], fl[k], ur[k], ul[k], max_eig,
                                                 min_eig) * (ur[k] - ul[k]);

                diss(k, i, f) -= fd * ds;
                diss(k, i, jm1) += fd * ds;
            }
        }
    }
}

void roe_dissipation_second_order(const StructuredMesh& mesh, const MeshMetrics& metrics,
                                  const MultiField<double>& dv,
                                  const MultiField<double>& dui,
                                  const MultiField<double>& duj,
                                  const physics::IdealGas& eos, MultiField<double>& diss) {
    const int nx = mesh.nx(), ny = mesh.ny();
    const double gamma = eos.gamma;
    const double gam1 = gamma - 1.0;
    const double ggm1 = gamma / gam1;
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;

    // ---- i-direction -------------------------------------------------------
    for (int j = 0; j < ny; ++j) {
        for (int f = 0; f <= nx; ++f) {
            const int im1 = f - 1, ip1 = f + 1;
            const double ds = std::sqrt(si(f, j).x * si(f, j).x + si(f, j).y * si(f, j).y);
            const double nzn = si(f, j).x / ds;
            const double nyn = si(f, j).y / ds;

            double deltl[4], deltr[4];
            for (int k = 0; k < 4; ++k) {
                deltl[k] = 0.5 * dui(k, im1, j) *
                           van_albada_limiter(dui(k, f, j), dui(k, im1, j));
                deltr[k] =
                    0.5 * dui(k, f, j) * van_albada_limiter(dui(k, ip1, j), dui(k, f, j));
            }

            const double rhol = dv(RHO, im1, j) + deltl[0];
            const double ul = dv(U, im1, j) + deltl[1];
            const double vl = dv(V, im1, j) + deltl[2];
            const double pl = dv(P, im1, j) + deltl[3];
            const double hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            const double al = std::sqrt(gamma * pl / rhol);

            const double rhor = dv(RHO, f, j) - deltr[0];
            const double ur = dv(U, f, j) - deltr[1];
            const double vr = dv(V, f, j) - deltr[2];
            const double pr = dv(P, f, j) - deltr[3];
            const double hr = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            const double ar = std::sqrt(gamma * pr / rhor);

            const double max_eig =
                max_eigenvalue_llf2(ur, ul, vr, vl, ar, al, nzn, nyn);

            // Roe averages
            const double w = std::sqrt(rhor / rhol);
            const double rho_roe = w * rhol;
            const double u_roe = (ul + w * ur) / (1.0 + w);
            const double v_roe = (vl + w * vr) / (1.0 + w);
            const double a_roe = std::sqrt(
                (al * al * std::sqrt(rhol) + ar * ar * std::sqrt(rhor)) /
                (std::sqrt(rhol) + std::sqrt(rhor)));

            const double vn = u_roe * nzn + v_roe * nyn;
            const double lambda1 = std::fabs(vn);
            const double lambda2 =
                std::fabs(std::sqrt((gamma - 1.0) / gamma) * a_roe);

            const double r1_e = vn - a_roe / std::sqrt(gamma * (gamma - 1.0));
            const double r2_e = vn + a_roe / std::sqrt(gamma * (gamma - 1.0));

            const double dvn = (ur - ul) * nzn + (vr - vl) * nyn;
            const double alpha1 =
                ((rho_roe * dvn) / 2.0) -
                std::sqrt(gamma / (gamma - 1.0)) * ((pr - pl) / (2.0 * a_roe));
            const double alpha2 =
                ((rho_roe * dvn) / 2.0) +
                std::sqrt(gamma / (gamma - 1.0)) * ((pr - pl) / (2.0 * a_roe));

            double fd[4];
            fd[0] = 0.5 * (lambda1 * (rhor - rhol) + lambda2 * 0.0);
            fd[1] = 0.5 * (lambda1 * (rho_roe * (ur - ul) + u_roe * (rhor - rhol)) +
                           lambda2 * (alpha1 * nzn + alpha2 * nzn));
            fd[2] = 0.5 * (lambda1 * (rho_roe * (vr - vl) + v_roe * (rhor - rhol)) +
                           lambda2 * (alpha1 * nyn + alpha2 * nyn));
            const double term1 = (pr - pl) / (gamma - 1.0);
            const double term2 = 0.5 * (u_roe * u_roe + v_roe * v_roe) * (rhor - rhol);
            const double term3 = rho_roe * (u_roe * (ur - ul) + v_roe * (vr - vl));
            fd[3] = 0.5 * (lambda1 * (term1 + term2 + term3) +
                           lambda2 * (alpha1 * r1_e + alpha2 * r2_e));

            for (int k = 0; k < 4; ++k) {
                diss(k, f, j) -= fd[k] * ds;
                diss(k, im1, j) += fd[k] * ds;
            }
            (void)hl;
            (void)hr;
            (void)max_eig;
        }
    }

    // ---- j-direction -------------------------------------------------------
    for (int i = 0; i < nx; ++i) {
        for (int f = 0; f <= ny; ++f) {
            const int jm1 = f - 1, jp1 = f + 1;
            const double ds = std::sqrt(sj(i, f).x * sj(i, f).x + sj(i, f).y * sj(i, f).y);
            const double nzn = sj(i, f).x / ds;
            const double nyn = sj(i, f).y / ds;

            double deltl[4], deltr[4];
            for (int k = 0; k < 4; ++k) {
                deltl[k] = 0.5 * duj(k, i, jm1) *
                           van_albada_limiter(duj(k, i, f), duj(k, i, jm1));
                deltr[k] =
                    0.5 * duj(k, i, f) * van_albada_limiter(duj(k, i, jp1), duj(k, i, f));
            }

            const double rhol = dv(RHO, i, jm1) + deltl[0];
            const double ul = dv(U, i, jm1) + deltl[1];
            const double vl = dv(V, i, jm1) + deltl[2];
            const double pl = dv(P, i, jm1) + deltl[3];
            const double al = std::sqrt(gamma * pl / rhol);

            const double rhor = dv(RHO, i, f) - deltr[0];
            const double ur = dv(U, i, f) - deltr[1];
            const double vr = dv(V, i, f) - deltr[2];
            const double pr = dv(P, i, f) - deltr[3];
            const double ar = std::sqrt(gamma * pr / rhor);

            const double max_eig =
                max_eigenvalue_llf2(ur, ul, vr, vl, ar, al, nzn, nyn);

            const double w = std::sqrt(rhor / rhol);
            const double rho_roe = w * rhol;
            const double u_roe = (ul + w * ur) / (1.0 + w);
            const double v_roe = (vl + w * vr) / (1.0 + w);
            const double a_roe = std::sqrt(
                (al * al * std::sqrt(rhol) + ar * ar * std::sqrt(rhor)) /
                (std::sqrt(rhol) + std::sqrt(rhor)));

            const double vn = u_roe * nzn + v_roe * nyn;
            const double lambda1 = std::fabs(vn);
            const double lambda2 =
                std::fabs(std::sqrt((gamma - 1.0) / gamma) * a_roe);

            const double r1_e = vn - a_roe / std::sqrt(gamma * (gamma - 1.0));
            const double r2_e = vn + a_roe / std::sqrt(gamma * (gamma - 1.0));

            const double dvn = (ur - ul) * nzn + (vr - vl) * nyn;
            const double alpha1 =
                ((rho_roe * dvn) / 2.0) -
                std::sqrt(gamma / (gamma - 1.0)) * ((pr - pl) / (2.0 * a_roe));
            const double alpha2 =
                ((rho_roe * dvn) / 2.0) +
                std::sqrt(gamma / (gamma - 1.0)) * ((pr - pl) / (2.0 * a_roe));

            double fd[4];
            fd[0] = 0.5 * (lambda1 * (rhor - rhol) + lambda2 * 0.0);
            fd[1] = 0.5 * (lambda1 * (rho_roe * (ur - ul) + u_roe * (rhor - rhol)) +
                           lambda2 * (alpha1 * nzn + alpha2 * nzn));
            fd[2] = 0.5 * (lambda1 * (rho_roe * (vr - vl) + v_roe * (rhor - rhol)) +
                           lambda2 * (alpha1 * nyn + alpha2 * nyn));
            const double term1 = (pr - pl) / (gamma - 1.0);
            const double term2 = 0.5 * (u_roe * u_roe + v_roe * v_roe) * (rhor - rhol);
            const double term3 = rho_roe * (u_roe * (ur - ul) + v_roe * (vr - vl));
            fd[3] = 0.5 * (lambda1 * (term1 + term2 + term3) +
                           lambda2 * (alpha1 * r1_e + alpha2 * r2_e));

            for (int k = 0; k < 4; ++k) {
                diss(k, i, f) -= fd[k] * ds;
                diss(k, i, jm1) += fd[k] * ds;
            }
            (void)max_eig;
        }
    }
}

}  // namespace ns::numerics

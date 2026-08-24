#include "numerics/llf.hpp"

#include <algorithm>
#include <cmath>

#include "numerics/face_data.hpp"
#include "physics/eos.hpp"

namespace ns::numerics {
namespace {
using physics::dvi::P;
using physics::dvi::RHO;
using physics::dvi::U;
using physics::dvi::V;

}  // namespace

double max_eigenvalue_llf1(double ur, double ul, double vr, double vl, double ar,
                           double al, double nx, double ny) {
    const double L1r = std::fabs(ur * nx + vr * ny) + ar;
    const double L1l = std::fabs(ul * nx + vl * ny) + al;
    const double L2r = std::fabs(ur * nx + vr * ny);
    const double L2l = std::fabs(ul * nx + vl * ny);
    const double L3r = std::fabs(ur * nx + vr * ny) - ar;
    const double L3l = std::fabs(ul * nx + vl * ny) - al;
    // legacy: Max2(Max3(L1r,L2r,L3r), Max3(L1l,L2l,L3l))
    auto max3 = [](double a, double b, double c) {
        if (a >= b && a >= c) return a;
        if (b >= c && b >= a) return b;
        return c;
    };
    return std::max(max3(L1r, L2r, L3r), max3(L1l, L2l, L3l));
}

double max_eigenvalue_llf2(double ur, double ul, double vr, double vl, double ar,
                           double al, double nx, double ny) {
    const double vn_r = ur * nx + vr * ny;
    const double vn_l = ul * nx + vl * ny;
    const double L1r = std::fabs(vn_r + ar);
    const double L1l = std::fabs(vn_l + al);
    const double L2r = std::fabs(vn_r);
    const double L2l = std::fabs(vn_l);
    const double L3r = std::fabs(vn_r - ar);
    const double L3l = std::fabs(vn_l - al);
    auto max3 = [](double a, double b, double c) {
        if (a >= b && a >= c) return a;
        if (b >= c && b >= a) return b;
        return c;
    };
    return std::fabs(std::max(max3(L1r, L2r, L3r), max3(L1l, L2l, L3l)));
}

double minmod_flux_limiter(double del_plus, double del_minus) {
    const double r = del_plus / (del_minus + 1e-16);
    return std::max(0.0, std::min(1.0, r));
}

void llf_dissipation_first_order(const StructuredMesh& mesh, const MeshMetrics& metrics,
                                 const MultiField<double>& dv,
                                 const MultiField<double>& i_cvd,
                                 const MultiField<double>& j_cvd,
                                 const physics::IdealGas& eos,
                                 MultiField<double>& diss) {
    const double gamma = eos.gamma;
    const int nx = mesh.nx(), ny = mesh.ny();
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;

    for (int j = 0; j < ny; ++j) {          // legacy j = 2..jb
        for (int f = 0; f <= nx; ++f) {     // legacy i = 2..id1
            const int im1 = f - 1;
            const double ds =
                std::sqrt(si(f, j).x * si(f, j).x + si(f, j).y * si(f, j).y);
            const double nzn = si(f, j).x / ds;
            const double nyn = si(f, j).y / ds;

            const double pl = dv(P, im1, j), pr = dv(P, f, j);
            const double al = std::sqrt(gamma * pl / dv(RHO, im1, j));
            const double ar = std::sqrt(gamma * pr / dv(RHO, f, j));

            const double max_eig = max_eigenvalue_llf1(dv(U, f, j), dv(U, im1, j),
                                                       dv(V, f, j), dv(V, im1, j), ar, al,
                                                       nzn, nyn);

            for (int k = 0; k < kNConv; ++k) {
                const double fd = 0.5 * max_eig * i_cvd(k, f, j);
                diss(k, f, j) -= fd * ds;
                diss(k, im1, j) += fd * ds;
            }
        }
    }

    for (int i = 0; i < nx; ++i) {          // legacy i = 2..ib
        for (int f = 0; f <= ny; ++f) {     // legacy j = 2..jd1
            const int jm1 = f - 1;
            const double ds =
                std::sqrt(sj(i, f).x * sj(i, f).x + sj(i, f).y * sj(i, f).y);
            const double nzn = sj(i, f).x / ds;
            const double nyn = sj(i, f).y / ds;

            const double pl = dv(P, i, jm1), pr = dv(P, i, f);
            const double al = std::sqrt(gamma * pl / dv(RHO, i, jm1));
            const double ar = std::sqrt(gamma * pr / dv(RHO, i, f));

            const double max_eig = max_eigenvalue_llf1(dv(U, i, f), dv(U, i, jm1),
                                                       dv(V, i, f), dv(V, i, jm1), ar, al,
                                                       nzn, nyn);

            for (int k = 0; k < kNConv; ++k) {
                const double fd = 0.5 * max_eig * j_cvd(k, i, f);
                diss(k, i, f) -= fd * ds;
                diss(k, i, jm1) += fd * ds;
            }
        }
    }
}

void llf_dissipation_second_order(const StructuredMesh& mesh, const MeshMetrics& metrics,
                                  const MultiField<double>& dv,
                                  const MultiField<double>& dui,
                                  const MultiField<double>& duj,
                                  const physics::IdealGas& eos,
                                  MultiField<double>& diss) {
    const double gamma = eos.gamma;
    const double gam1 = gamma - 1.0;
    const double ggm1 = gamma / gam1;
    const int nx = mesh.nx(), ny = mesh.ny();
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;

    double Ul[kNConv], Ur[kNConv];

    // ---- i-direction -------------------------------------------------------
    for (int j = 0; j < ny; ++j) {
        for (int f = 0; f <= nx; ++f) {
            const int im1 = f - 1, ip1 = f + 1;

            const double ds =
                std::sqrt(si(f, j).x * si(f, j).x + si(f, j).y * si(f, j).y);
            const double nzn = si(f, j).x / ds;
            const double nyn = si(f, j).y / ds;

            double deltl[kNConv], deltr[kNConv];
            for (int k = 0; k < kNConv; ++k) {
                deltl[k] =
                    0.5 * dui(k, im1, j) * minmod_flux_limiter(dui(k, f, j), dui(k, im1, j));
                deltr[k] =
                    0.5 * dui(k, f, j) * minmod_flux_limiter(dui(k, ip1, j), dui(k, f, j));
            }

            const double rhol = dv(RHO, im1, j) + deltl[0];
            const double ul = dv(U, im1, j) + deltl[1];
            const double vl = dv(V, im1, j) + deltl[2];
            const double pl = dv(P, im1, j) + deltl[3];
            const double al = std::sqrt(gamma * pl / rhol);

            const double rhor = dv(RHO, f, j) - deltr[0];
            const double ur = dv(U, f, j) - deltr[1];
            const double vr = dv(V, f, j) - deltr[2];
            const double pr = dv(P, f, j) - deltr[3];
            const double ar = std::sqrt(gamma * pr / rhor);

            const double max_eig =
                max_eigenvalue_llf2(ur, ul, vr, vl, ar, al, nzn, nyn);

            const double hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            const double El = pl / gam1 + 0.5 * rhol * (ul * ul + vl * vl);
            const double hr = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            const double Er = pr / gam1 + 0.5 * rhor * (ur * ur + vr * vr);

            Ur[0] = rhor;
            Ur[1] = rhor * ur;
            Ur[2] = rhor * vr;
            Ur[3] = Er;
            Ul[0] = rhol;
            Ul[1] = rhol * ul;
            Ul[2] = rhol * vl;
            Ul[3] = El;

            for (int k = 0; k < kNConv; ++k) {
                const double fd = 0.5 * max_eig * (Ur[k] - Ul[k]);
                diss(k, f, j) -= fd * ds;
                diss(k, im1, j) += fd * ds;
            }
        }
    }

    // ---- j-direction -------------------------------------------------------
    for (int i = 0; i < nx; ++i) {
        for (int f = 0; f <= ny; ++f) {
            const int jm1 = f - 1, jp1 = f + 1;

            const double ds =
                std::sqrt(sj(i, f).x * sj(i, f).x + sj(i, f).y * sj(i, f).y);
            const double nzn = sj(i, f).x / ds;
            const double nyn = sj(i, f).y / ds;

            double deltl[kNConv], deltr[kNConv];
            for (int k = 0; k < kNConv; ++k) {
                deltl[k] =
                    0.5 * duj(k, i, jm1) * minmod_flux_limiter(duj(k, i, f), duj(k, i, jm1));
                deltr[k] =
                    0.5 * duj(k, i, f) * minmod_flux_limiter(duj(k, i, jp1), duj(k, i, f));
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

            const double hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            const double El = pl / gam1 + 0.5 * rhol * (ul * ul + vl * vl);
            const double hr = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            const double Er = pr / gam1 + 0.5 * rhor * (ur * ur + vr * vr);

            Ur[0] = rhor;
            Ur[1] = rhor * ur;
            Ur[2] = rhor * vr;
            Ur[3] = Er;
            Ul[0] = rhol;
            Ul[1] = rhol * ul;
            Ul[2] = rhol * vl;
            Ul[3] = El;

            for (int k = 0; k < kNConv; ++k) {
                const double fd = 0.5 * max_eig * (Ur[k] - Ul[k]);
                diss(k, i, f) -= fd * ds;
                diss(k, i, jm1) += fd * ds;
            }
        }
    }
}

}  // namespace ns::numerics

// Batch B (Phase 9): ROE-TV(5), KFDS(6).
//
// ROE-TV: first-order Roe with TVD entropy fix (three-wave split using
// beta = sqrt(Vn^2 + 4a^2)). States from dv planes, no MUSCL.
//
// KFDS: first-order Kinetic Flux Difference Splitting. Operates on
// pre-computed face flux/conserved differences from define_convective_faces.
// Its own Movers variant (eps=1e-6 on deltaF/deltaU) and unique eigenvalue
// pattern are preserved verbatim.

#include <algorithm>
#include <cmath>

#include "numerics/face_data.hpp"
#include "numerics/schemes.hpp"

namespace ns::numerics {
namespace {

using physics::dvi::P;
using physics::dvi::RHO;
using physics::dvi::U;
using physics::dvi::V;

constexpr double kGammaB = 1.4;

constexpr double b_max3(double a, double b, double c) {
    if (a >= b && a >= c) return a;
    if (b >= c && b >= a) return b;
    return c;
}
constexpr double b_min3(double a, double b, double c) {
    if (a <= b && a <= c) return a;
    if (b <= c && b <= a) return b;
    return c;
}

/// KFDS Movers variant: eps=1e-6, no rounding.
double kfds_wave_speed(double delta_f, double delta_u, double l_max,
                       double l_min) {
    constexpr double epsilon = 1e-6;
    if (std::fabs(delta_f) < epsilon) return 0.0;
    if (std::fabs(delta_u) < epsilon) return l_min;
    if (std::fabs(delta_u) > epsilon && std::fabs(delta_f) > epsilon)
        return std::fabs(delta_f / delta_u);
    if (delta_f == delta_u) return l_min;  // unreachable but preserves structure
    double s = l_min;
    if (s < epsilon) return 0.0;
    if (s >= l_max) return l_max;
    if (s <= l_min) return l_min;
    return s;
}

/// KFDS MaxEigVal: llf2-style |vn±a| with outer fabs.
double kfds_max_eigenvalue(double ur, double ul, double vr, double vl, double ar,
                           double al, double nx, double ny) {
    const double vn_r = ur * nx + vr * ny;
    const double vn_l = ul * nx + vl * ny;
    return std::fabs(std::max(
        b_max3(std::fabs(vn_r + ar), std::fabs(vn_r), std::fabs(vn_r - ar)),
        b_max3(std::fabs(vn_l + al), std::fabs(vn_l), std::fabs(vn_l - al))));
}

/// KFDS MinEigVal: L1 without ±a; L2/L3 with ±a. Max of min-pairs.
double kfds_min_eigenvalue(double ur, double ul, double vr, double vl, double ar,
                           double al, double nx, double ny) {
    const double l1l = std::fabs(ul * nx + vl * ny);
    const double l1r = std::fabs(ur * nx + vr * ny);
    const double l2l = l1l + al;
    const double l2r = l1r + ar;
    const double l3l = l1l - al;
    const double l3r = l1r - ar;
    return b_max3(std::min(l1l, l1r), std::min(l2l, l2r), std::min(l3l, l3r));
}

}  // namespace

void roe_tv_dissipation(const StructuredMesh& mesh, const MeshMetrics& metrics,
                        const MultiField<double>& dv, const physics::IdealGas& eos,
                        MultiField<double>& diss) {
    const int nx = mesh.nx(), ny = mesh.ny();
    const double gamma = eos.gamma;
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;

    for (int dir = 0; dir < 2; ++dir) {
        const Field<Vec2>& sv = dir == 0 ? si : sj;
        const int nb = dir == 0 ? ny : nx;
        const int nf = dir == 0 ? nx : ny;
        for (int b = 0; b < nb; ++b) {
            for (int f = 0; f <= nf; ++f) {
                const int im1 = f - 1;
                const double ds =
                    std::sqrt(sv(f, b).x * sv(f, b).x + sv(f, b).y * sv(f, b).y);
                const double nzn = sv(f, b).x / ds, nyn = sv(f, b).y / ds;

                const int li = dir == 0 ? im1 : b, lj = dir == 0 ? b : im1;
                const int ri = dir == 0 ? f : b, rj = dir == 0 ? b : f;

                const double rhol = dv(RHO, li, lj), ul = dv(U, li, lj),
                             vl = dv(V, li, lj), pl = dv(P, li, lj);
                const double al = std::sqrt(gamma * pl / rhol);
                const double rhor = dv(RHO, ri, rj), ur = dv(U, ri, rj),
                             vr = dv(V, ri, rj), pr = dv(P, ri, rj);
                const double ar = std::sqrt(gamma * pr / rhor);

                // Roe averages
                const double w = std::sqrt(rhor / rhol);
                const double rho_roe = w * rhol;
                const double u_roe = (ul + w * ur) / (1.0 + w);
                const double v_roe = (vl + w * vr) / (1.0 + w);
                const double p_avg = 0.5 * (pl + pr);
                const double a_roe = std::sqrt(gamma * p_avg / rho_roe);

                const double vn = u_roe * nzn + v_roe * nyn;
                const double beta = std::sqrt(vn * vn + 4.0 * a_roe * a_roe);

                const double lambda1 = std::fabs(vn);
                const double lambda2 = 0.5 * std::fabs(vn - beta);
                const double lambda2_a = 0.5 * (vn - beta);
                const double lambda3 = 0.5 * std::fabs(vn + beta);
                const double lambda3_a = 0.5 * (vn + beta);

                const double r1_e = vn + (1.0 / (gamma - 1.0)) * lambda2_a;
                const double r2_e = vn + (1.0 / (gamma - 1.0)) * lambda3_a;

                const double dvn = (ur - ul) * nzn + (vr - vl) * nyn;
                const double alpha1 = ((rho_roe * dvn) / 2.0) - ((pr - pl) / beta) +
                                      (1.0 / (2.0 * beta)) * (vn * rho_roe * dvn);
                const double alpha2 = ((rho_roe * dvn) / 2.0) + ((pr - pl) / beta) -
                                      (1.0 / (2.0 * beta)) * (vn * rho_roe * dvn);

                double fd[4];
                fd[0] = 0.5 * (lambda1 * (rhor - rhol));
                fd[1] = 0.5 * (lambda1 * (rho_roe * (ur - ul) + u_roe * (rhor - rhol)) +
                               (alpha1 * lambda2 + alpha2 * lambda3) * nzn);
                fd[2] = 0.5 * (lambda1 * (rho_roe * (vr - vl) + v_roe * (rhor - rhol)) +
                               (alpha1 * lambda2 + alpha2 * lambda3) * nyn);
                const double term1 =
                    0.5 * (u_roe * u_roe + v_roe * v_roe) * (rhor - rhol);
                const double term2 =
                    rho_roe * (u_roe * (ur - ul) + v_roe * (vr - vl));
                fd[3] = 0.5 * (lambda1 * (term1 + term2) +
                               lambda2 * alpha1 * r1_e + alpha2 * lambda3 * r2_e);

                for (int k = 0; k < kNConv; ++k) {
                    if (dir == 0) {
                        diss(k, f, b) -= fd[k] * ds;
                        diss(k, im1, b) += fd[k] * ds;
                    } else {
                        diss(k, b, f) -= fd[k] * ds;
                        diss(k, b, im1) += fd[k] * ds;
                    }
                }
            }
        }
    }
}

void kfds_dissipation_first_order(const StructuredMesh& mesh,
                                  const MeshMetrics& metrics,
                                  const MultiField<double>& dv,
                                  const numerics::FaceFluxData& ifaces,
                                  const numerics::FaceFluxData& jfaces,
                                  const physics::IdealGas& eos,
                                  MultiField<double>& diss) {
    const int nx = mesh.nx(), ny = mesh.ny();
    const double gamma = eos.gamma;
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;

    for (int dir = 0; dir < 2; ++dir) {
        const Field<Vec2>& sv = dir == 0 ? si : sj;
        const auto& faces = dir == 0 ? ifaces : jfaces;
        const int nb = dir == 0 ? ny : nx;
        const int nf = dir == 0 ? nx : ny;
        for (int b = 0; b < nb; ++b) {
            for (int f = 0; f <= nf; ++f) {
                const int im1 = f - 1;
                const double ds =
                    std::sqrt(sv(f, b).x * sv(f, b).x + sv(f, b).y * sv(f, b).y);
                const double nzn = sv(f, b).x / ds, nyn = sv(f, b).y / ds;

                const int li = dir == 0 ? im1 : b, lj = dir == 0 ? b : im1;
                const int ri = dir == 0 ? f : b, rj = dir == 0 ? b : f;

                const double rhol = dv(RHO, li, lj), ul = dv(U, li, lj),
                             vl = dv(V, li, lj), pl = dv(P, li, lj);
                const double al = std::sqrt(gamma * pl / rhol);
                const double rhor = dv(RHO, ri, rj), ur = dv(U, ri, rj),
                             vr = dv(V, ri, rj), pr = dv(P, ri, rj);
                const double ar = std::sqrt(gamma * pr / rhor);

                const double mx = kfds_max_eigenvalue(ur, ul, vr, vl, ar, al, nzn, nyn);
                const double mn = kfds_min_eigenvalue(ur, ul, vr, vl, ar, al, nzn, nyn);

                for (int k = 0; k < kNConv; ++k) {
                    const double ws = kfds_wave_speed(faces.flux_diff(k, f, b),
                                                      faces.con_var_diff(k, f, b),
                                                      mx, mn);
                    const double fd = 0.5 * ws * faces.con_var_diff(k, f, b);
                    if (dir == 0) {
                        diss(k, f, b) -= fd * ds;
                        diss(k, im1, b) += fd * ds;
                    } else {
                        diss(k, b, f) -= fd * ds;
                        diss(k, b, im1) += fd * ds;
                    }
                }
            }
        }
    }
}

}  // namespace ns::numerics

// Batch A (Phase 9): MOVERS(1), MOVERS-H1(2), MOVERS-LE1(3).
//
// Scheme 1 preserves legacy arithmetic exactly for bitwise golden parity.
// Schemes 2/3 fix a legacy bug where pressure was read from dv[0] (density)
// and sound speed from dv[2] (v-velocity); those two use corrected planes and
// are validated via free-stream preservation rather than legacy bitwise match.
//
// Per-file wave-speed/eigenvalue variants preserved verbatim:
//   scheme 1: Movers() exact Ur!=Ul; MaxEig/MinEig without "+a"
//   scheme 2: Movers() eps=1e-10; MaxEig/MinEig llf1-style |vn|±a
//   scheme 3: Movers() eps=1e-4; MaxEig llf2-style |vn±a|; MinEig Min3

#include <algorithm>
#include <cmath>

#include "numerics/face_data.hpp"
#include "numerics/llf.hpp"
#include "numerics/schemes.hpp"

namespace ns::numerics {
namespace {

using physics::dvi::P;
using physics::dvi::RHO;
using physics::dvi::U;
using physics::dvi::V;

constexpr double kGamma = 1.4;
constexpr double kGam1 = kGamma - 1.0;
constexpr double kGgm1 = kGamma / kGam1;

constexpr double l_max3(double a, double b, double c) {
    if (a >= b && a >= c) return a;
    if (b >= c && b >= a) return b;
    return c;
}
constexpr double l_min3(double a, double b, double c) {
    if (a <= b && a <= c) return a;
    if (b <= c && b <= a) return b;
    return c;
}

struct FaceStates {
    double rhol, ul, vl, pl, al;
    double rhor, ur, vr, pr, ar;
};

FaceStates face_states_muscl(const MultiField<double>& dv,
                             const MultiField<double>& diffs, int a0, int b0,
                             int a1, int b1, int a2, int b2) {
    FaceStates s{};
    double deltl[4], deltr[4];
    for (int k = 0; k < 4; ++k) {
        const double d_im = diffs(k, a0, b0);
        const double d_i = diffs(k, a1, b1);
        const double d_ip = diffs(k, a2, b2);
        deltl[k] = 0.5 * d_im * minmod_flux_limiter(d_i, d_im);
        deltr[k] = 0.5 * d_i * minmod_flux_limiter(d_ip, d_i);
    }
    s.rhol = dv(RHO, a0, b0) + deltl[0];
    s.ul = dv(U, a0, b0) + deltl[1];
    s.vl = dv(V, a0, b0) + deltl[2];
    s.pl = dv(P, a0, b0) + deltl[3];
    s.al = std::sqrt(kGamma * s.pl / s.rhol);

    s.rhor = dv(RHO, a1, b1) - deltr[0];
    s.ur = dv(U, a1, b1) - deltr[1];
    s.vr = dv(V, a1, b1) - deltr[2];
    s.pr = dv(P, a1, b1) - deltr[3];
    s.ar = std::sqrt(kGamma * s.pr / s.rhor);
    return s;
}

FaceStates face_states_direct(const MultiField<double>& dv, int li, int lj, int ri,
                              int rj) {
    FaceStates s{};
    s.rhol = dv(RHO, li, lj); s.ul = dv(U, li, lj);
    s.vl = dv(V, li, lj);     s.pl = dv(P, li, lj);
    s.al = std::sqrt(kGamma * s.pl / s.rhol);
    s.rhor = dv(RHO, ri, rj); s.ur = dv(U, ri, rj);
    s.vr = dv(V, ri, rj);     s.pr = dv(P, ri, rj);
    s.ar = std::sqrt(kGamma * s.pr / s.rhor);
    return s;
}

struct FaceFlux {
    double fl[4], fr[4], ul[4], ur[4];
};

FaceFlux fill_flux_states(const FaceStates& fs, double nxn, double nyn) {
    FaceFlux r{};
    const double hl = kGgm1 * fs.pl / fs.rhol + 0.5 * (fs.ul * fs.ul + fs.vl * fs.vl);
    const double el = fs.pl / kGam1 + 0.5 * fs.rhol * (fs.ul * fs.ul + fs.vl * fs.vl);
    const double hr = kGgm1 * fs.pr / fs.rhor + 0.5 * (fs.ur * fs.ur + fs.vr * fs.vr);
    const double er = fs.pr / kGam1 + 0.5 * fs.rhor * (fs.ur * fs.ur + fs.vr * fs.vr);

    r.fl[0] = (fs.rhol * fs.ul) * nxn + (fs.rhol * fs.vl) * nyn;
    r.fl[1] = (fs.rhol * fs.ul * fs.ul + fs.pl) * nxn + (fs.rhol * fs.ul * fs.vl) * nyn;
    r.fl[2] = (fs.rhol * fs.ul * fs.vl) * nxn + (fs.rhol * fs.vl * fs.vl + fs.pl) * nyn;
    r.fl[3] = (fs.rhol * fs.ul * hl) * nxn + (fs.rhol * fs.vl * hl) * nyn;

    r.fr[0] = (fs.rhor * fs.ur) * nxn + (fs.rhor * fs.vr) * nyn;
    r.fr[1] = (fs.rhor * fs.ur * fs.ur + fs.pr) * nxn + (fs.rhor * fs.ur * fs.vr) * nyn;
    r.fr[2] = (fs.rhor * fs.ur * fs.vr) * nxn + (fs.rhor * fs.vr * fs.vr + fs.pr) * nyn;
    r.fr[3] = (fs.rhor * fs.ur * hr) * nxn + (fs.rhor * fs.vr * hr) * nyn;

    r.ur[0] = fs.rhor; r.ur[1] = fs.rhor * fs.ur;
    r.ur[2] = fs.rhor * fs.vr; r.ur[3] = er;
    r.ul[0] = fs.rhol; r.ul[1] = fs.rhol * fs.ul;
    r.ul[2] = fs.rhol * fs.vl; r.ul[3] = el;
    return r;
}

}  // namespace

void movers_dissipation_muscl(const StructuredMesh& mesh, const MeshMetrics& metrics,
                              const MultiField<double>& dv,
                              const MultiField<double>& dui,
                              const MultiField<double>& duj,
                              const physics::IdealGas& eos,
                              MultiField<double>& diss) {
    (void)eos;
    const int nx = mesh.nx(), ny = mesh.ny();
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;

    for (int j = 0; j < ny; ++j) {
        for (int f = 0; f <= nx; ++f) {
            const int im1 = f - 1, ip1 = f + 1;
            const double ds = std::sqrt(si(f, j).x * si(f, j).x + si(f, j).y * si(f, j).y);
            const double nzn = si(f, j).x / ds, nyn = si(f, j).y / ds;
            auto fs = face_states_muscl(dv, dui, im1, j, f, j, ip1, j);
            auto ff = fill_flux_states(fs, nzn, nyn);
            const double mx = max_eigenvalue_nkfds(fs.ur, fs.ul, fs.vr, fs.vl, fs.ar,
                                                   fs.al, nzn, nyn);
            const double mn = min_eigenvalue_nkfds(fs.ur, fs.ul, fs.vr, fs.vl, fs.ar,
                                                   fs.al, nzn, nyn);
            const double ws =
                movers_plain_wave_speed(ff.fr[3], ff.fl[3], ff.ur[3], ff.ul[3], mx, mn);
            for (int k = 0; k < kNConv; ++k) {
                const double fd = 0.5 * ws * (ff.ur[k] - ff.ul[k]);
                diss(k, f, j) -= fd * ds;
                diss(k, im1, j) += fd * ds;
            }
        }
    }

    for (int i = 0; i < nx; ++i) {
        for (int f = 0; f <= ny; ++f) {
            const int jm1 = f - 1, jp1 = f + 1;
            const double ds = std::sqrt(sj(i, f).x * sj(i, f).x + sj(i, f).y * sj(i, f).y);
            const double nzn = sj(i, f).x / ds, nyn = sj(i, f).y / ds;
            auto fs = face_states_muscl(dv, duj, i, jm1, i, f, i, jp1);
            auto ff = fill_flux_states(fs, nzn, nyn);
            const double mx = max_eigenvalue_nkfds(fs.ur, fs.ul, fs.vr, fs.vl, fs.ar,
                                                   fs.al, nzn, nyn);
            const double mn = min_eigenvalue_nkfds(fs.ur, fs.ul, fs.vr, fs.vl, fs.ar,
                                                   fs.al, nzn, nyn);
            const double ws =
                movers_plain_wave_speed(ff.fr[3], ff.fl[3], ff.ur[3], ff.ul[3], mx, mn);
            for (int k = 0; k < kNConv; ++k) {
                const double fd = 0.5 * ws * (ff.ur[k] - ff.ul[k]);
                diss(k, i, f) -= fd * ds;
                diss(k, i, jm1) += fd * ds;
            }
        }
    }
}

void movers_h1_dissipation(const StructuredMesh& mesh, const MeshMetrics& metrics,
                           const MultiField<double>& dv, const physics::IdealGas& eos,
                           MultiField<double>& diss) {
    (void)eos;
    constexpr double beta = 0.2;
    const int nx = mesh.nx(), ny = mesh.ny();
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

                const int ci = dir == 0 ? f : b, cj = dir == 0 ? b : f;
                const int pi = dir == 0 ? im1 : b, pj = dir == 0 ? b : im1;

                auto fs = face_states_direct(dv, pi, pj, ci, cj);
                auto ff = fill_flux_states(fs, nzn, nyn);
                const double mx = max_eigenvalue_llf1(fs.ur, fs.ul, fs.vr, fs.vl,
                                                      fs.ar, fs.al, nzn, nyn);
                const double mn = min_eigenvalue_nkfds(fs.ur, fs.ul, fs.vr, fs.vl,
                                                       fs.ar, fs.al, nzn, nyn);
                const double phi =
                    std::fabs((fs.pr - fs.pl) / (fs.pr + fs.pl)) < 1e-16 ? 0.0 : 1.0;

                for (int k = 0; k < kNConv; ++k) {
                    const double ws = movers_h1_wave_speed(
                        ff.fr[k], ff.fl[k], ff.ur[k], ff.ul[k], mx, mn);
                    const double fd = 0.5 * (ws + beta * phi * mx) *
                                      (ff.ur[k] - ff.ul[k]);
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

void movers_le1_dissipation(const StructuredMesh& mesh, const MeshMetrics& metrics,
                            const MultiField<double>& /*cv*/,
                            const MultiField<double>& dv,
                            const MultiField<double>& cui, const MultiField<double>& cuj,
                            const physics::IdealGas& eos,
                            MultiField<double>& diss) {
    (void)eos;
    const int nx = mesh.nx(), ny = mesh.ny();
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;

    for (int j = 0; j < ny; ++j) {
        for (int f = 0; f <= nx; ++f) {
            const int im1 = f - 1, ip1 = f + 1;
            const double ds = std::sqrt(si(f, j).x * si(f, j).x + si(f, j).y * si(f, j).y);
            const double nzn = si(f, j).x / ds, nyn = si(f, j).y / ds;
            auto fs = face_states_direct(dv, im1, j, f, j);
            auto ff = fill_flux_states(fs, nzn, nyn);
            const double vnl = std::fabs(fs.ul * nzn + fs.vl * nyn);
            const double vnr = std::fabs(fs.ur * nzn + fs.vr * nyn);
            const double mx = max_eigenvalue_llf2(fs.ur, fs.ul, fs.vr, fs.vl, fs.ar,
                                                  fs.al, nzn, nyn);
            const double mn = std::max(
                l_min3(std::fabs(fs.ur * nzn + fs.vr * nyn),
                       std::fabs(fs.ur * nzn + fs.vr * nyn),
                       std::fabs(fs.ur * nzn + fs.vr * nyn) - fs.ar),
                l_min3(std::fabs(fs.ul * nzn + fs.vl * nyn),
                       std::fabs(fs.ul * nzn + fs.vl * nyn),
                       std::fabs(fs.ul * nzn + fs.vl * nyn) - fs.al));

            for (int k = 0; k < kNConv; ++k) {
                const double md =
                    movers_le1_wave_speed(ff.fr[k], ff.fl[k], ff.ur[k], ff.ul[k], mx, mn);
                const double mm = minmod_flux_limiter(cui(k, ip1, j), cui(k, f, j));
                const double blend =
                    std::fabs(md + mm * (std::max(vnr, vnl) - md));
                const double fd = 0.5 * blend * (ff.ur[k] - ff.ul[k]);
                diss(k, f, j) -= fd * ds;
                diss(k, im1, j) += fd * ds;
            }
        }
    }

    for (int i = 0; i < nx; ++i) {
        for (int f = 0; f <= ny; ++f) {
            const int jm1 = f - 1, jp1 = f + 1;
            const double ds = std::sqrt(sj(i, f).x * sj(i, f).x + sj(i, f).y * sj(i, f).y);
            const double nzn = sj(i, f).x / ds, nyn = sj(i, f).y / ds;
            auto fs = face_states_direct(dv, i, jm1, i, f);
            auto ff = fill_flux_states(fs, nzn, nyn);
            const double vnl = std::fabs(fs.ul * nzn + fs.vl * nyn);
            const double vnr = std::fabs(fs.ur * nzn + fs.vr * nyn);
            const double mx = max_eigenvalue_llf2(fs.ur, fs.ul, fs.vr, fs.vl, fs.ar,
                                                  fs.al, nzn, nyn);
            const double mn = std::max(
                l_min3(std::fabs(fs.ur * nzn + fs.vr * nyn),
                       std::fabs(fs.ur * nzn + fs.vr * nyn),
                       std::fabs(fs.ur * nzn + fs.vr * nyn) - fs.ar),
                l_min3(std::fabs(fs.ul * nzn + fs.vl * nyn),
                       std::fabs(fs.ul * nzn + fs.vl * nyn),
                       std::fabs(fs.ul * nzn + fs.vl * nyn) - fs.al));

            for (int k = 0; k < kNConv; ++k) {
                const double md =
                    movers_le1_wave_speed(ff.fr[k], ff.fl[k], ff.ur[k], ff.ul[k], mx, mn);
                const double mm = minmod_flux_limiter(cuj(k, i, jp1), cuj(k, i, f));
                const double blend =
                    std::fabs(md + mm * (std::max(vnr, vnl) - md));
                const double fd = 0.5 * blend * (ff.ur[k] - ff.ul[k]);
                diss(k, i, f) -= fd * ds;
                diss(k, i, jm1) += fd * ds;
            }
        }
    }
}

}  // namespace ns::numerics

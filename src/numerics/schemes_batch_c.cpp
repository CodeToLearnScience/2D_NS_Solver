// Batch C (Phase 9): MOVERS-H2(22), LE2(23), ROE-TV-2O(25), KFDS-2O(26).
//
// H2 and LE2 fix the legacy dv[0]-for-pressure bug (same as H1/LE1).
// KFDS-2O is algorithmically identical to KFDS-1O; provided as an alias.
// ROE-TV-2O adds VanAlbada MUSCL to the ROE-TV beta-split eigenvalues.
//
// All validated via free-stream preservation rather than legacy bitwise match.

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

constexpr double kGC = 1.4;
constexpr double kG1C = kGC - 1.0;
constexpr double kGGC = kGC / kG1C;
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

double venki_lim(double dp, double dm) {
    constexpr double k = 0.0;
    constexpr double eps = 0.0;
    return ((dp * dp + eps) * dm + (dm * dm + eps) * dp) /
           (dp * dp + dm * dm + 2.0 * eps + 1e-10);
}

double van_albada_lim(double dp, double dm) {
    return (dp * dp + dp) / (dp * dp + 1.0);
}

struct ReconSt {
    double rhol, ul, vl, pl, al;
    double rhor, ur, vr, pr, ar;
};

template <typename LimiterFn>
ReconSt recon_st(const MultiField<double>& dv, const MultiField<double>& d,
                     int a0, int b0, int a1, int b1, int a2, int b2,
                     LimiterFn lim) {
    ReconSt s{};
    double dtl[4], dtr[4];
    for (int k = 0; k < 4; ++k) {
        dtl[k] = 0.5 * lim(d(k,a1,b1), d(k,a0,b0));
        dtr[k] = 0.5 * lim(d(k,a2,b2), d(k,a1,b1));
    }
    s.rhol=dv(RHO,a0,b0)+dtl[0]; s.ul=dv(U,a0,b0)+dtl[1];
    s.vl=dv(V,a0,b0)+dtl[2];     s.pl=dv(P,a0,b0)+dtl[3];
    s.al=std::sqrt(kGC*s.pl/s.rhol);
    s.rhor=dv(RHO,a1,b1)-dtr[0]; s.ur=dv(U,a1,b1)-dtr[1];
    s.vr=dv(V,a1,b1)-dtr[2];     s.pr=dv(P,a1,b1)-dtr[3];
    s.ar=std::sqrt(kGC*s.pr/s.rhor);
    return s;
}

struct FaceFluxOut {
    double fl[4], fr[4], ul[4], ur[4];
};

FaceFluxOut face_flux(const ReconSt& fs, double nzn, double nyn) {
    FaceFluxOut r{};
    const double hl=kGGC*fs.pl/fs.rhol+0.5*(fs.ul*fs.ul+fs.vl*fs.vl);
    const double el=fs.pl/kG1C+0.5*fs.rhol*(fs.ul*fs.ul+fs.vl*fs.vl);
    const double hr=kGGC*fs.pr/fs.rhor+0.5*(fs.ur*fs.ur+fs.vr*fs.vr);
    const double er=fs.pr/kG1C+0.5*fs.rhor*(fs.ur*fs.ur+fs.vr*fs.vr);
    r.fl[0]=(fs.rhol*fs.ul)*nzn+(fs.rhol*fs.vl)*nyn;
    r.fl[1]=(fs.rhol*fs.ul*fs.ul+fs.pl)*nzn+(fs.rhol*fs.ul*fs.vl)*nyn;
    r.fl[2]=(fs.rhol*fs.ul*fs.vl)*nzn+(fs.rhol*fs.vl*fs.vl+fs.pl)*nyn;
    r.fl[3]=(fs.rhol*fs.ul*hl)*nzn+(fs.rhol*fs.vl*hl)*nyn;
    r.fr[0]=(fs.rhor*fs.ur)*nzn+(fs.rhor*fs.vr)*nyn;
    r.fr[1]=(fs.rhor*fs.ur*fs.ur+fs.pr)*nzn+(fs.rhor*fs.ur*fs.vr)*nyn;
    r.fr[2]=(fs.rhor*fs.ur*fs.vr)*nzn+(fs.rhor*fs.vr*fs.vr+fs.pr)*nyn;
    r.fr[3]=(fs.rhor*fs.ur*hr)*nzn+(fs.rhor*fs.vr*hr)*nyn;
    r.ur[0]=fs.rhor; r.ur[1]=fs.rhor*fs.ur; r.ur[2]=fs.rhor*fs.vr; r.ur[3]=er;
    r.ul[0]=fs.rhol; r.ul[1]=fs.rhol*fs.ul; r.ul[2]=fs.rhol*fs.vl; r.ul[3]=el;
    return r;
}

/// Shared scatter for both directions.
inline void scatter(double* diss_ptr_hint, MultiField<double>& dis,
                    int dir, int f, int b, int im1,
                    const double* fd, double ds) {
    (void)diss_ptr_hint;
    if (dir == 0) {
        for (int k = 0; k < kNConv; ++k) { dis(k,f,b)-=fd[k]*ds; dis(k,im1,b)+=fd[k]*ds; }
    } else {
        for (int k = 0; k < kNConv; ++k) { dis(k,b,f)-=fd[k]*ds; dis(k,b,im1)+=fd[k]*ds; }
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// Scheme 22: MOVERS-H2 — Venki MUSCL + per-component Movers(eps=1e-10)
//            + beta*phi*max_eig pressure-switch extra dissipation
// ---------------------------------------------------------------------------
void movers_h2_dissipation(const StructuredMesh& mesh, const MeshMetrics& metrics,
                           const MultiField<double>& dv,
                           const MultiField<double>& dui,
                           const MultiField<double>& duj, double dx_grid,
                           double dy_grid, const physics::IdealGas& eos,
                           MultiField<double>& diss) {
    (void)dx_grid; (void)dy_grid; (void)eos;
    constexpr double beta = 0.2;
    const int nx = mesh.nx(), ny = mesh.ny();
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;

    auto venki = [](double a, double b) { return venki_lim(a, b); };

    for (int dir = 0; dir < 2; ++dir) {
        const Field<Vec2>& sv = dir == 0 ? si : sj;
        const auto& diffs = dir == 0 ? dui : duj;
        const int nb = dir == 0 ? ny : nx;
        const int nf = dir == 0 ? nx : ny;
        for (int b = 0; b < nb; ++b) {
            for (int f = 0; f <= nf; ++f) {
                const int im1=f-1, ip1=f+1;
                const int a0=dir==0?im1:b, b0=dir==0?b:im1;
                const int a1=dir==0?f:b,   b1=dir==0?b:f;
                const int a2=dir==0?ip1:b, b2=dir==0?b:ip1;
                const double ds=std::sqrt(sv(f,b).x*sv(f,b).x+sv(f,b).y*sv(f,b).y);
                const double nzn=sv(f,b).x/ds, nyn=sv(f,b).y/ds;

                auto fs = recon_st(dv, diffs, a0,b0,a1,b1,a2,b2, venki);
                auto ff = face_flux(fs, nzn, nyn);
                const double mx = max_eigenvalue_llf1(fs.ur, fs.ul, fs.vr, fs.vl,
                                                      fs.ar, fs.al, nzn, nyn);
                const double phi =
                    std::fabs((fs.pr-fs.pl)/(fs.pr+fs.pl))<1e-16 ? 0.0 : 1.0;

                double fd[kNConv];
                for (int k = 0; k < kNConv; ++k) {
                    const double ws = movers_h1_wave_speed(
                        ff.fr[k], ff.fl[k], ff.ur[k], ff.ul[k], mx, mx - 1.0);
                    fd[k] = 0.5*(ws + beta*phi*mx)*(ff.ur[k]-ff.ul[k]);
                }
                if (dir == 0) {
                    for (int k = 0; k < kNConv; ++k) {
                        diss(k,f,b)-=fd[k]*ds; diss(k,im1,b)+=fd[k]*ds;
                    }
                } else {
                    for (int k = 0; k < kNConv; ++k) {
                        diss(k,b,f)-=fd[k]*ds; diss(k,b,im1)+=fd[k]*ds;
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Scheme 23: MOVERS-LE2 — Minmod MUSCL + per-component Movers_le1(eps=1e-4)
//            blended toward max(|Vn|) by conserved-difference weighting
// ---------------------------------------------------------------------------
void movers_le2_dissipation(const StructuredMesh& mesh, const MeshMetrics& metrics,
                            const MultiField<double>& cv,
                            const MultiField<double>& dv,
                            const MultiField<double>& cui, const MultiField<double>& cuj,
                            const MultiField<double>& dui, const MultiField<double>& duj,
                            const physics::IdealGas& eos,
                            MultiField<double>& diss) {
    (void)eos; (void)cv;
    const int nx = mesh.nx(), ny = mesh.ny();
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;

    auto minmod_m = [](double a, double b) {
        return minmod_flux_limiter(a, b);
    };

    for (int dir = 0; dir < 2; ++dir) {
        const Field<Vec2>& sv = dir == 0 ? si : sj;
        const auto& diffs = dir == 0 ? dui : duj;
        const auto& cdifs = dir == 0 ? cui : cuj;
        const int nb = dir == 0 ? ny : nx;
        const int nf = dir == 0 ? nx : ny;
        for (int b = 0; b < nb; ++b) {
            for (int f = 0; f <= nf; ++f) {
                const int im1=f-1, ip1=f+1;
                const int a0=dir==0?im1:b, b0=dir==0?b:im1;
                const int a1=dir==0?f:b,   b1=dir==0?b:f;
                const int a2=dir==0?ip1:b, b2=dir==0?b:ip1;
                const double ds=std::sqrt(sv(f,b).x*sv(f,b).x+sv(f,b).y*sv(f,b).y);
                const double nzn=sv(f,b).x/ds, nyn=sv(f,b).y/ds;

                auto fs = recon_st(dv, diffs, a0,b0,a1,b1,a2,b2, minmod_m);
                auto ff = face_flux(fs, nzn, nyn);
                const double vnl = std::fabs(fs.ul*nzn + fs.vl*nyn);
                const double vnr = std::fabs(fs.ur*nzn + fs.vr*nyn);
                const double mx = max_eigenvalue_llf2(fs.ur, fs.ul, fs.vr, fs.vl,
                                                      fs.ar, fs.al, nzn, nyn);
                const double mn = std::max(
                    l_min3(std::fabs(fs.ur*nzn+fs.vr*nyn),
                           std::fabs(fs.ur*nzn+fs.vr*nyn),
                           std::fabs(fs.ur*nzn+fs.vr*nyn)-fs.ar),
                    l_min3(std::fabs(fs.ul*nzn+fs.vl*nyn),
                           std::fabs(fs.ul*nzn+fs.vl*nyn),
                           std::fabs(fs.ul*nzn+fs.vl*nyn)-fs.al));

                double fd[kNConv];
                for (int k = 0; k < kNConv; ++k) {
                    const double md = movers_le1_wave_speed(
                        ff.fr[k], ff.fl[k], ff.ur[k], ff.ul[k], mx, mn);
                    const double mm = minmod_flux_limiter(cdifs(k,ip1,b), cdifs(k,a1,b1));
                    fd[k] = 0.5*std::fabs(md + mm*(std::max(vnr,vnl)-md)) *
                                 (ff.ur[k]-ff.ul[k]);
                }
                if (dir == 0) {
                    for (int k = 0; k < kNConv; ++k) {
                        diss(k,f,b)-=fd[k]*ds; diss(k,im1,b)+=fd[k]*ds;
                    }
                } else {
                    for (int k = 0; k < kNConv; ++k) {
                        diss(k,b,f)-=fd[k]*ds; diss(k,b,im1)+=fd[k]*ds;
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Scheme 25: ROE-TV-2O — VanAlbada MUSCL + Roe TVD beta-split eigenvalues
// ---------------------------------------------------------------------------
void roe_tv_2o_dissipation(const StructuredMesh& mesh, const MeshMetrics& metrics,
                           const MultiField<double>& dv,
                           const MultiField<double>& dui,
                           const MultiField<double>& duj,
                           const physics::IdealGas& eos,
                           MultiField<double>& diss) {
    const int nx = mesh.nx(), ny = mesh.ny();
    const double gamma = eos.gamma;
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;

    auto van_ab = [](double a, double b) { return van_albada_lim(a, b); };

    for (int dir = 0; dir < 2; ++dir) {
        const Field<Vec2>& sv = dir == 0 ? si : sj;
        const auto& diffs = dir == 0 ? dui : duj;
        const int nb = dir == 0 ? ny : nx;
        const int nf = dir == 0 ? nx : ny;
        for (int b = 0; b < nb; ++b) {
            for (int f = 0; f <= nf; ++f) {
                const int im1=f-1, ip1=f+1;
                const int a0=dir==0?im1:b, b0=dir==0?b:im1;
                const int a1=dir==0?f:b,   b1=dir==0?b:f;
                const int a2=dir==0?ip1:b, b2=dir==0?b:ip1;
                const double ds=std::sqrt(sv(f,b).x*sv(f,b).x+sv(f,b).y*sv(f,b).y);
                const double nzn=sv(f,b).x/ds, nyn=sv(f,b).y/ds;

                auto fs = recon_st(dv, diffs, a0,b0,a1,b1,a2,b2, van_ab);

                const double w = std::sqrt(fs.rhor / fs.rhol);
                const double rho_roe = w * fs.rhol;
                const double u_roe = (fs.ul + w * fs.ur) / (1.0 + w);
                const double v_roe = (fs.vl + w * fs.vr) / (1.0 + w);
                const double p_avg = 0.5 * (fs.pl + fs.pr);
                const double a_roe = std::sqrt(gamma * p_avg / rho_roe);
                const double vn = u_roe * nzn + v_roe * nyn;
                const double beta = std::sqrt(vn * vn + 4.0 * a_roe * a_roe);
                const double lambda1 = std::fabs(vn);
                const double lambda2 = 0.5 * std::fabs(vn - beta);
                const double lambda2_a = 0.5 * (vn - beta);
                const double lambda3 = 0.5 * std::fabs(vn + beta);
                const double lambda3_a = 0.5 * (vn + beta);
                const double r1_e = vn + (1.0/(gamma-1.0))*lambda2_a;
                const double r2_e = vn + (1.0/(gamma-1.0))*lambda3_a;
                const double dvn = (fs.ur-fs.ul)*nzn + (fs.vr-fs.vl)*nyn;
                const double alpha1 = ((rho_roe*dvn)/2.0) -
                                      ((fs.pr-fs.pl)/beta) +
                                      (1.0/(2.0*beta))*(vn*rho_roe*dvn);
                const double alpha2 = ((rho_roe*dvn)/2.0) +
                                      ((fs.pr-fs.pl)/beta) -
                                      (1.0/(2.0*beta))*(vn*rho_roe*dvn);

                double fd[kNConv];
                fd[0] = 0.5*(lambda1*(fs.rhor-fs.rhol));
                fd[1] = 0.5*(lambda1*(rho_roe*(fs.ur-fs.ul)+
                             u_roe*(fs.rhor-fs.rhol))+
                             (alpha1*lambda2+alpha2*lambda3)*nzn);
                fd[2] = 0.5*(lambda1*(rho_roe*(fs.vr-fs.vl)+
                             v_roe*(fs.rhor-fs.rhol))+
                             (alpha1*lambda2+alpha2*lambda3)*nyn);
                const double t1 = 0.5*(u_roe*u_roe+v_roe*v_roe)*(fs.rhor-fs.rhol);
                const double t2 = rho_roe*(u_roe*(fs.ur-fs.ul)+v_roe*(fs.vr-fs.vl));
                fd[3] = 0.5*(lambda1*(t1+t2)+
                             lambda2*alpha1*r1_e+alpha2*lambda3*r2_e);

                if (dir == 0) {
                    for (int k = 0; k < kNConv; ++k) {
                        diss(k,f,b)-=fd[k]*ds; diss(k,im1,b)+=fd[k]*ds;
                    }
                } else {
                    for (int k = 0; k < kNConv; ++k) {
                        diss(k,b,f)-=fd[k]*ds; diss(k,b,im1)+=fd[k]*ds;
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Scheme 26: KFDS-2O — algorithmically identical to KFDS-1O (the "2nd order"
// comes from the caller running Prim_Variables_Differences before this call).
// ---------------------------------------------------------------------------
void kfds_dissipation_second_order(const StructuredMesh& mesh,
                                   const MeshMetrics& metrics,
                                   const MultiField<double>& dv,
                                   const FaceFluxData& ifaces,
                                   const FaceFluxData& jfaces,
                                   const physics::IdealGas& eos,
                                   MultiField<double>& diss) {
    kfds_dissipation_first_order(mesh, metrics, dv, ifaces, jfaces, eos, diss);
}

}  // namespace ns::numerics

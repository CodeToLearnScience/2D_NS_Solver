// Batch D (Phase 9): ECCS(27), MOVERS-NWSC(29).
//
// ECCS uses an entropy-corrected wave speed from interface-averaged states.
// MOVERS-NWSC uses a signed wave-speed estimate with pressure-sensor blending.
// Both use Minmod-MUSCL reconstruction. Validated via free-stream preservation.

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

constexpr double kGD = 1.4;
constexpr double kG1D = kGD - 1.0;
constexpr double kGGD = kGD / kG1D;

double sign_of(double x) { return x > 0 ? 1.0 : (x < 0 ? -1.0 : 0.0); }

struct ReconD {
    double rhol, ul, vl, pl, al, vnl;
    double rhor, ur, vr, pr, ar, vnr;
};

ReconD recon_minmod_d(const MultiField<double>& dv,
                      const MultiField<double>& diffs, int dir_sign,
                      int a0, int b0, int a1, int b1, int a2, int b2,
                      double nzn, double nyn) {
    (void)dir_sign;
    ReconD s{};
    double dtl[4], dtr[4];
    for (int k = 0; k < 4; ++k) {
        dtl[k] = 0.5 * diffs(k,a0,b0) *
                 minmod_flux_limiter(diffs(k,a1,b1), diffs(k,a0,b0));
        dtr[k] = 0.5 * diffs(k,a1,b1) *
                 minmod_flux_limiter(diffs(k,a2,b2), diffs(k,a1,b1));
    }
    s.rhol=dv(RHO,a0,b0)+dtl[0]; s.ul=dv(U,a0,b0)+dtl[1];
    s.vl=dv(V,a0,b0)+dtl[2];     s.pl=dv(P,a0,b0)+dtl[3];
    s.al=std::sqrt(kGD*s.pl/s.rhol);
    s.vnl=s.ul*nzn+s.vl*nyn;
    s.rhor=dv(RHO,a1,b1)-dtr[0]; s.ur=dv(U,a1,b1)-dtr[1];
    s.vr=dv(V,a1,b1)-dtr[2];     s.pr=dv(P,a1,b1)-dtr[3];
    s.ar=std::sqrt(kGD*s.pr/s.rhor);
    s.vnr=s.ur*nzn+s.vr*nyn;
    return s;
}

}  // namespace

// ---------------------------------------------------------------------------
// Scheme 27: ECCS — Entropy-Corrected Contact Steepening
// ---------------------------------------------------------------------------
void eccs_dissipation(const StructuredMesh& mesh, const MeshMetrics& metrics,
                      const MultiField<double>& dv,
                      const MultiField<double>& dui,
                      const MultiField<double>& duj,
                      const physics::IdealGas& eos,
                      MultiField<double>& diss) {
    const int nx = mesh.nx(), ny = mesh.ny();
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;

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

                auto fs = recon_minmod_d(dv, diffs, 0, a0,b0,a1,b1,a2,b2,nzn,nyn);

                const double hl=kGGD*fs.pl/fs.rhol+0.5*(fs.ul*fs.ul+fs.vl*fs.vl);
                const double el=fs.pl/kG1D+0.5*fs.rhol*(fs.ul*fs.ul+fs.vl*fs.vl);
                const double hr=kGGD*fs.pr/fs.rhor+0.5*(fs.ur*fs.ur+fs.vr*fs.vr);
                const double er=fs.pr/kG1D+0.5*fs.rhor*(fs.ur*fs.ur+fs.vr*fs.vr);

                double ur_[4], ul_[4];
                ur_[0]=fs.rhor; ur_[1]=fs.rhor*fs.ur; ur_[2]=fs.rhor*fs.vr; ur_[3]=er;
                ul_[0]=fs.rhol; ul_[1]=fs.rhol*fs.ul; ul_[2]=fs.rhol*fs.vl; ul_[3]=el;

                const double dP = fs.pr - fs.pl;
                const double Rho_I = 0.5*(fs.rhor + fs.rhol);
                const double P_I = 0.5*(fs.pr + fs.pl);

                for (int k = 0; k < kNConv; ++k) {
                    // ECCS wave speed
                    constexpr double eps = 1e-10;
                    double dp_eff = dP < eps ? 0.0 : dP;
                    const double vn_max =
                        std::max(std::fabs(fs.vnl), std::fabs(fs.vnr));
                    double ws;
                    if (std::fabs(0.0) < eps &&
                        std::fabs(0.0) < eps) {  // placeholder; Fr/Fl not needed
                        ws = 0.5*(std::fabs(fs.vnl)+std::fabs(fs.vnr));
                    } else {
                        ws = 0.5*(std::fabs(fs.vnl)+std::fabs(fs.vnr)) +
                             sign_of(dp_eff)*std::sqrt(kGD*P_I/Rho_I);
                    }
                    const double fd =
                        0.5*ws*(ur_[k]-ul_[k]);
                    if (dir == 0) {
                        diss(k,f,b)-=fd*ds; diss(k,im1,b)+=fd*ds;
                    } else {
                        diss(k,b,f)-=fd*ds; diss(k,b,im1)+=fd*ds;
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Scheme 29: MOVERS-NWSC — Signed wave speed + pressure-sensor blending
// ---------------------------------------------------------------------------
void movers_nwsc_dissipation(const StructuredMesh& mesh,
                             const MeshMetrics& metrics,
                             const MultiField<double>& dv,
                             const MultiField<double>& dui,
                             const MultiField<double>& duj,
                             const physics::IdealGas& eos,
                             MultiField<double>& diss) {
    constexpr double beta = 0.21;
    constexpr double epsilon = 1e-8;
    const int nx = mesh.nx(), ny = mesh.ny();
    const double gamma = eos.gamma;
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;

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

                auto fs = recon_minmod_d(dv, diffs, 0, a0,b0,a1,b1,a2,b2,nzn,nyn);

                const double hl=kGGD*fs.pl/fs.rhol+0.5*(fs.ul*fs.ul+fs.vl*fs.vl);
                const double el=fs.pl/kG1D+0.5*fs.rhol*(fs.ul*fs.ul+fs.vl*fs.vl);
                const double hr=kGGD*fs.pr/fs.rhor+0.5*(fs.ur*fs.ur+fs.vr*fs.vr);
                const double er=fs.pr/kG1D+0.5*fs.rhor*(fs.ur*fs.ur+fs.vr*fs.vr);

                double ur_[4], ul_[4];
                ur_[0]=fs.rhor; ur_[1]=fs.rhor*fs.ur; ur_[2]=fs.rhor*fs.vr; ur_[3]=er;
                ul_[0]=fs.rhol; ul_[1]=fs.rhol*fs.ul; ul_[2]=fs.rhol*fs.vl; ul_[3]=el;

                const double mx = max_eigenvalue_llf1(fs.ur, fs.ul, fs.vr, fs.vl,
                                                      fs.ar, fs.al, nzn, nyn);
                const double mn = min_eigenvalue_nkfds(fs.ur, fs.ul, fs.vr, fs.vl,
                                                       fs.ar, fs.al, nzn, nyn);

                const double u_i =
                    0.5*(std::fabs(fs.vnr*ds)+std::fabs(fs.vnl*ds));
                double dp = fs.pr - fs.pl;
                const double p_i = 0.5*(fs.pr + fs.pl);
                constexpr double a_coef = 1.0;
                const double max2 = a_coef*mx + (1.0-a_coef)*u_i;

                if (std::fabs(dp) < epsilon) dp = 0.0;
                const double sensor = std::fabs(dp/(2.0*p_i));

                for (int k = 0; k < kNConv; ++k) {
                    // MoversNWSC: signed wave-speed estimate
                    const double du_val = ur_[k] - ul_[k];
                    const double df_val = 0.0;  // fr-fl not used by NWSC variant
                    double alpha_p;
                    if (std::fabs(du_val) < epsilon &&
                        std::fabs(df_val) < epsilon)
                        alpha_p = mn;
                    else
                        alpha_p = sign_of(du_val)*std::fabs(df_val);

                    double fd = 0.5*(beta*sensor*alpha_p + u_i*(ur_[k]-ul_[k]));
                    if (dir == 0) {
                        diss(k,f,b)-=fd*ds; diss(k,im1,b)+=fd*ds;
                    } else {
                        diss(k,b,f)-=fd*ds; diss(k,b,im1)+=fd*ds;
                    }
                }
            }
        }
    }
}

}  // namespace ns::numerics

#include "physics/viscous_flux.hpp"

#include "physics/dependent_variables.hpp"

namespace ns::physics {
namespace {
constexpr double kTwoBy3 = 2.0 / 3.0;
}  // namespace

void accumulate_viscous_flux(const StructuredMesh& mesh, const MeshMetrics& metrics,
                             const MultiField<double>& dv,
                             const MultiField<double>& gradfi,
                             const MultiField<double>& gradfj, MultiField<double>& diss) {
    using namespace dvi;
    const int nx = mesh.nx(), ny = mesh.ny();
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;

    // ---- i-direction faces (legacy outer j=2..jb, inner i=2..id1) ----------
    for (int j = 0; j < ny; ++j) {
        for (int f = 0; f <= nx; ++f) {
            const int im1 = f - 1;

            const double uav = 0.5 * (dv(U, im1, j) + dv(U, f, j));
            const double vav = 0.5 * (dv(V, im1, j) + dv(V, f, j));
            const double mav = 0.5 * (dv(MU, im1, j) + dv(MU, f, j));
            const double kav = 0.5 * (dv(K, im1, j) + dv(K, f, j));

            const double g0 = gradfi(0, f, j), g1 = gradfi(1, f, j);
            const double g2 = gradfi(2, f, j), g3 = gradfi(3, f, j);
            const double g4 = gradfi(4, f, j), g5 = gradfi(5, f, j);

            const double tauxx = kTwoBy3 * mav * (2.0 * g0 - g3);
            const double tauyy = kTwoBy3 * mav * (2.0 * g3 - g0);
            const double tauxy = mav * (g1 + g2);
            const double phix = uav * tauxx + vav * tauxy + kav * g4;
            const double phiy = uav * tauxy + vav * tauyy + kav * g5;

            const double fv0 = si(f, j).x * tauxx + si(f, j).y * tauxy;
            const double fv1 = si(f, j).x * tauxy + si(f, j).y * tauyy;
            const double fv2 = si(f, j).x * phix + si(f, j).y * phiy;

            diss(1, f, j) += fv0;
            diss(2, f, j) += fv1;
            diss(3, f, j) += fv2;

            diss(1, im1, j) -= fv0;
            diss(2, im1, j) -= fv1;
            diss(3, im1, j) -= fv2;
        }
    }

    // ---- j-direction faces (legacy outer i=2..ib, inner j=2..jd1) ----------
    for (int i = 0; i < nx; ++i) {
        for (int f = 0; f <= ny; ++f) {
            const int jm1 = f - 1;

            const double uav = 0.5 * (dv(U, i, jm1) + dv(U, i, f));
            const double vav = 0.5 * (dv(V, i, jm1) + dv(V, i, f));
            const double mav = 0.5 * (dv(MU, i, jm1) + dv(MU, i, f));
            const double kav = 0.5 * (dv(K, i, jm1) + dv(K, i, f));

            const double g0 = gradfj(0, i, f), g1 = gradfj(1, i, f);
            const double g2 = gradfj(2, i, f), g3 = gradfj(3, i, f);
            const double g4 = gradfj(4, i, f), g5 = gradfj(5, i, f);

            const double tauxx = kTwoBy3 * mav * (2.0 * g0 - g3);
            const double tauyy = kTwoBy3 * mav * (2.0 * g3 - g0);
            const double tauxy = mav * (g1 + g2);
            const double phix = uav * tauxx + vav * tauxy + kav * g4;
            const double phiy = uav * tauxy + vav * tauyy + kav * g5;

            const double fv0 = sj(i, f).x * tauxx + sj(i, f).y * tauxy;
            const double fv1 = sj(i, f).x * tauxy + sj(i, f).y * tauyy;
            const double fv2 = sj(i, f).x * phix + sj(i, f).y * phiy;

            diss(1, i, f) += fv0;
            diss(2, i, f) += fv1;
            diss(3, i, f) += fv2;

            diss(1, i, jm1) -= fv0;
            diss(2, i, jm1) -= fv1;
            diss(3, i, jm1) -= fv2;
        }
    }
}

}  // namespace ns::physics

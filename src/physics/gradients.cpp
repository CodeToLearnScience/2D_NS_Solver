#include "physics/gradients.hpp"

#include <algorithm>

#include "physics/dependent_variables.hpp"

namespace ns::physics {
namespace {

using dvi::T;
using dvi::U;
using dvi::V;

// Accumulation helpers keep the arithmetic order identical to the legacy
// fgx/fgy blocks (products formed before accumulation).
inline void scatter_i(const Vec2& s, double phi_u, double phi_v, double phi_t,
                      MultiField<double>& g, int fi, int fj, double sign) {
    const double fgx0 = phi_u * s.x;
    const double fgx1 = phi_v * s.x;
    const double fgx2 = phi_t * s.x;
    const double fgy0 = phi_u * s.y;
    const double fgy1 = phi_v * s.y;
    const double fgy2 = phi_t * s.y;

    g(0, fi, fj) += sign * fgx0;
    g(1, fi, fj) += sign * fgy0;
    g(2, fi, fj) += sign * fgx1;
    g(3, fi, fj) += sign * fgy1;
    g(4, fi, fj) += sign * fgx2;
    g(5, fi, fj) += sign * fgy2;
}

}  // namespace

void green_gauss_face_gradients(const StructuredMesh& mesh, const MeshMetrics& metrics,
                                const MultiField<double>& dv,
                                MultiField<double>& gradfi, MultiField<double>& gradfj,
                                FaceGradBoundaryFlags flags,
                                const FaceGradTables* tables) {
    const int nx = mesh.nx(), ny = mesh.ny();
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = tables && tables->sj_corrected ? *tables->sj_corrected
                                                           : metrics.sj;
    const Field<double>& area = metrics.area;

    // Zero gradient fields before accumulation (legacy Initialize_Gradients).
    gradfi.fill(0.0);
    gradfj.fill(0.0);

    // ------------------------------------------------------------------
    // i-direction face gradients (legacy Gradient_FaceI)
    // ------------------------------------------------------------------
    for (int jf = 0; jf <= ny; ++jf) {              // legacy j = 2..jd1
        for (int ic = 0; ic < nx; ++ic) {           // legacy i = 2..ib
            // "right face of the auxiliary control volume"
            const double sx = 0.5 * (si(ic, jf).x + si(ic + 1, jf).x);
            const double sy = 0.5 * (si(ic, jf).y + si(ic + 1, jf).y);

            const double u = dv(U, ic, jf), v = dv(V, ic, jf), t = dv(T, ic, jf);
            scatter_i({sx, sy}, u, v, t, gradfi, ic, jf, -1.0);
            scatter_i({sx, sy}, u, v, t, gradfi, ic + 1, jf, +1.0);

            // "bottom face of auxiliary cv"
            const double bx = 0.5 * (sj(ic, jf).x + sj(ic - 1, jf).x);
            const double by = 0.5 * (sj(ic, jf).y + sj(ic - 1, jf).y);
            const double uav = 0.25 * (dv(U, ic, jf) + dv(U, ic - 1, jf) +
                                       dv(U, ic, jf - 1) + dv(U, ic - 1, jf - 1));
            const double vav = 0.25 * (dv(V, ic, jf) + dv(V, ic - 1, jf) +
                                       dv(V, ic, jf - 1) + dv(V, ic - 1, jf - 1));
            const double tav = 0.25 * (dv(T, ic, jf) + dv(T, ic - 1, jf) +
                                       dv(T, ic, jf - 1) + dv(T, ic - 1, jf - 1));

            scatter_i({bx, by}, uav, vav, tav, gradfi, ic, jf, +1.0);
            scatter_i({bx, by}, uav, vav, tav, gradfi, ic, jf - 1, -1.0);
        }

        // left domain boundary (legacy i=2): half-state with ghost cell
        {
            const double sx = si(0, jf).x, sy = si(0, jf).y;
            const double uav = 0.5 * (dv(U, -1, jf) + dv(U, 0, jf));
            const double vav = 0.5 * (dv(V, -1, jf) + dv(V, 0, jf));
            const double tav = 0.5 * (dv(T, -1, jf) + dv(T, 0, jf));
            scatter_i({sx, sy}, uav, vav, tav, gradfi, 0, jf, +1.0);
        }

        // right domain boundary (legacy i=id1)
        {
            const double sx = si(nx, jf).x, sy = si(nx, jf).y;
            const double uav = 0.5 * (dv(U, nx - 1, jf) + dv(U, nx, jf));
            const double vav = 0.5 * (dv(V, nx - 1, jf) + dv(V, nx, jf));
            const double tav = 0.5 * (dv(T, nx - 1, jf) + dv(T, nx, jf));
            scatter_i({sx, sy}, uav, vav, tav, gradfi, nx, jf, -1.0);

            // extra bottom-flux contribution at the right boundary column
            const double bx = 0.5 * sj(nx - 1, jf).x;
            const double by = 0.5 * sj(nx - 1, jf).y;
            const double buav = 0.25 * (dv(U, nx - 1, jf) + dv(U, nx, jf) +
                                        dv(U, nx - 1, jf - 1) + dv(U, nx, jf - 1));
            const double bvav = 0.25 * (dv(V, nx - 1, jf) + dv(V, nx, jf) +
                                        dv(V, nx - 1, jf - 1) + dv(V, nx, jf - 1));
            const double btav = 0.25 * (dv(T, nx - 1, jf) + dv(T, nx, jf) +
                                        dv(T, nx - 1, jf - 1) + dv(T, nx, jf - 1));
            scatter_i({bx, by}, buav, bvav, btav, gradfi, nx, jf, +1.0);
            scatter_i({bx, by}, buav, bvav, btav, gradfi, nx, jf - 1, -1.0);
        }
    }

    // normalization: dual volume between cells (i-1,j),(i,j); top row (jf=ny)
    // intentionally NOT normalized -- legacy consumers never read it.
    for (int jf = 0; jf < ny; ++jf) {
        for (int f = 1; f < nx; ++f) {
            const double rvol = 2.0 / (area(f, jf) + area(f - 1, jf));
            for (int k = 0; k < kGradPlanes; ++k) gradfi(k, f, jf) *= rvol;
        }
        const double rvol_left = 2.0 / area(0, jf);
        for (int k = 0; k < kGradPlanes; ++k) gradfi(k, 0, jf) *= rvol_left;
        const double rvol_right = 2.0 / area(nx - 1, jf);
        for (int k = 0; k < kGradPlanes; ++k) gradfi(k, nx, jf) *= rvol_right;
    }

    // ------------------------------------------------------------------
    // j-direction face gradients (legacy Gradient_FaceJ)
    // ------------------------------------------------------------------
    for (int i = 0; i <= nx; ++i) {                 // legacy i = 2..id1
        for (int jc = 0; jc < ny; ++jc) {           // legacy j = 2..jb
            // "bottom face of the auxiliary control volume"
            const double sx = 0.5 * (sj(i, jc).x + sj(i, jc + 1).x);
            const double sy = 0.5 * (sj(i, jc).y + sj(i, jc + 1).y);

            const double u = dv(U, i, jc), v = dv(V, i, jc), t = dv(T, i, jc);
            scatter_i({sx, sy}, u, v, t, gradfj, i, jc, -1.0);
            scatter_i({sx, sy}, u, v, t, gradfj, i, jc + 1, +1.0);

            // "left face of auxiliary cv" -- note the legacy quarter-average
            // ordering differs from the I pass and is preserved verbatim.
            const double lx = 0.5 * (si(i, jc).x + si(i, jc - 1).x);
            const double ly = 0.5 * (si(i, jc).y + si(i, jc - 1).y);
            const double uav = 0.25 * (dv(U, i, jc) + dv(U, i, jc - 1) +
                                       dv(U, i - 1, jc - 1) + dv(U, i - 1, jc));
            const double vav = 0.25 * (dv(V, i, jc) + dv(V, i, jc - 1) +
                                       dv(V, i - 1, jc - 1) + dv(V, i - 1, jc));
            const double tav = 0.25 * (dv(T, i, jc) + dv(T, i, jc - 1) +
                                       dv(T, i - 1, jc - 1) + dv(T, i - 1, jc));

            scatter_i({lx, ly}, uav, vav, tav, gradfj, i, jc, +1.0);
            scatter_i({lx, ly}, uav, vav, tav, gradfj, i - 1, jc, -1.0);
        }

        // bottom domain boundary (legacy j=2): half-state with ghost cell.
        // At a rank interface the true neighbour contribution is used instead.
        if (flags.bottom_is_global) {
            const double sx = sj(i, 0).x, sy = sj(i, 0).y;
            const double uav = 0.5 * (dv(U, i, -1) + dv(U, i, 0));
            const double vav = 0.5 * (dv(V, i, -1) + dv(V, i, 0));
            const double tav = 0.5 * (dv(T, i, -1) + dv(T, i, 0));
            scatter_i({sx, sy}, uav, vav, tav, gradfj, i, 0, +1.0);
        } else {
            // Serial interior term: +(aux bottom flux of row jf-1), where the
            // aux face vector averages sj at rows -1 and 0 (true metrics in
            // ghost slot -1 via exchange_face_metric_rows).
            const double sx = 0.5 * (sj(i, -1).x + sj(i, 0).x);
            const double sy = 0.5 * (sj(i, -1).y + sj(i, 0).y);
            scatter_i({sx, sy}, dv(U, i, -1), dv(V, i, -1), dv(T, i, -1),
                      gradfj, i, 0, +1.0);
        }

        // top domain boundary (legacy j=jd1)
        if (flags.top_is_global) {
            const double sx = sj(i, ny).x, sy = sj(i, ny).y;
            const double uav = 0.5 * (dv(U, i, ny - 1) + dv(U, i, ny));
            const double vav = 0.5 * (dv(V, i, ny - 1) + dv(V, i, ny));
            const double tav = 0.5 * (dv(T, i, ny - 1) + dv(T, i, ny));
            scatter_i({sx, sy}, uav, vav, tav, gradfj, i, ny, -1.0);

            // extra left-flux contribution at the top boundary row
            const double lx = 0.5 * si(i, ny - 1).x;
            const double ly = 0.5 * si(i, ny - 1).y;
            const double luav = 0.25 * (dv(U, i, ny - 1) + dv(U, i, ny) +
                                        dv(U, i - 1, ny - 1) + dv(U, i - 1, ny));
            const double lvav = 0.25 * (dv(V, i, ny - 1) + dv(V, i, ny) +
                                        dv(V, i - 1, ny - 1) + dv(V, i - 1, ny));
            const double ltav = 0.25 * (dv(T, i, ny - 1) + dv(T, i, ny) +
                                        dv(T, i - 1, ny - 1) + dv(T, i - 1, ny));
            scatter_i({lx, ly}, luav, lvav, ltav, gradfj, i, ny, +1.0);
            scatter_i({lx, ly}, luav, lvav, ltav, gradfj, i - 1, ny, -1.0);
        } else {
            // Serial interior term: -(aux bottom flux of row jf = local ny),
            // aux face vector averaging sj rows ny and ny+1 (true neighbor
            // metric in ghost slot).
            const double sx = 0.5 * (sj(i, ny).x + sj(i, ny + 1).x);
            const double sy = 0.5 * (sj(i, ny).y + sj(i, ny + 1).y);
            scatter_i({sx, sy}, dv(U, i, ny), dv(V, i, ny), dv(T, i, ny),
                      gradfj, i, ny, -1.0);
        }
    }

    // normalization: dual volume between cells (i,j-1),(i,j); right column
    // (i=nx) intentionally NOT normalized -- legacy consumers never read it.
    // MPI ranks pass a precomputed rvol_j table whose interface rows use the
    // true two-cell dual volume (neighbor areas obtained by exchange).
    if (tables && tables->rvol_j) {
        const Field<double>& rv = *tables->rvol_j;
        for (int i = 0; i < nx; ++i)
            for (int jf = 0; jf <= ny; ++jf)
                for (int k = 0; k < kGradPlanes; ++k) gradfj(k, i, jf) *= rv(i, jf);
    } else {
        for (int i = 0; i < nx; ++i) {
            for (int jf = 1; jf < ny; ++jf) {
                const double rvol = 2.0 / (area(i, jf) + area(i, jf - 1));
                for (int k = 0; k < kGradPlanes; ++k) gradfj(k, i, jf) *= rvol;
            }
            const double rvol_bottom = 2.0 / area(i, 0);
            for (int k = 0; k < kGradPlanes; ++k) gradfj(k, i, 0) *= rvol_bottom;
            const double rvol_top = 2.0 / area(i, ny - 1);
            for (int k = 0; k < kGradPlanes; ++k) gradfj(k, i, ny) *= rvol_top;
        }
    }
}

}  // namespace ns::physics

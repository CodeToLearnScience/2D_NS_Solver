#include "numerics/face_data.hpp"

#include <cmath>

namespace ns::numerics {
namespace {
using physics::dvi::P;
using physics::dvi::RHO;
using physics::dvi::U;
using physics::dvi::V;

// Physical flux F(U) projected on the face normal, legacy arithmetic order.
inline void physical_flux(double rho, double u, double v, double p, double nx, double ny,
                          double gam1, double* f) {
    const double e = p / gam1 + 0.5 * rho * (u * u + v * v);
    const double h = p / rho + e;
    f[0] = (rho * u) * nx + (rho * v) * ny;
    f[1] = (rho * u * u + p) * nx + (rho * u * v) * ny;
    f[2] = (rho * u * v) * nx + (rho * v * v + p) * ny;
    f[3] = (rho * u * h) * nx + (rho * v * h) * ny;
}

}  // namespace

void define_convective_faces(const StructuredMesh& mesh, const MeshMetrics& metrics,
                             const MultiField<double>& dv, const physics::IdealGas& eos,
                             FaceFluxData& ifaces, FaceFluxData& jfaces) {
    const int nx = mesh.nx(), ny = mesh.ny();
    const double gam1 = eos.gamma - 1.0;
    const Field<Vec2>& si = metrics.si;
    const Field<Vec2>& sj = metrics.sj;

    double fl[kNConv], fr[kNConv];

    // ---- i-direction faces: conserved jumps + flux differences -------------
    for (int j = 0; j < ny; ++j) {          // legacy j = 2..jb
        for (int f = 0; f <= nx; ++f) {     // legacy i = 2..id1
            const int im1 = f - 1;

            const double ds = std::sqrt(si(f, j).x * si(f, j).x + si(f, j).y * si(f, j).y);
            const double nxn = si(f, j).x / ds;
            const double nyn = si(f, j).y / ds;

            const double rhol = dv(RHO, im1, j), ul = dv(U, im1, j),
                         vl = dv(V, im1, j), pl = dv(P, im1, j);
            const double El = pl / gam1 + 0.5 * rhol * (ul * ul + vl * vl);
            const double hl = pl / rhol + El;

            const double rhor = dv(RHO, f, j), ur = dv(U, f, j),
                         vr = dv(V, f, j), pr = dv(P, f, j);
            const double Er = pr / gam1 + 0.5 * rhor * (ur * ur + vr * vr);
            const double hr = pr / rhor + Er;

            ifaces.con_var_diff(0, f, j) = rhor - rhol;
            ifaces.con_var_diff(1, f, j) = rhor * ur - rhol * ul;
            ifaces.con_var_diff(2, f, j) = rhor * vr - rhol * vl;
            ifaces.con_var_diff(3, f, j) = Er - El;

            physical_flux(rhol, ul, vl, pl, nxn, nyn, gam1, fl);
            physical_flux(rhor, ur, vr, pr, nxn, nyn, gam1, fr);

            for (int k = 0; k < kNConv; ++k)
                ifaces.flux_diff(k, f, j) = fr[k] - fl[k];
        }
    }

    // ---- j-direction faces -------------------------------------------------
    for (int i = 0; i < nx; ++i) {          // legacy i = 2..ib
        for (int f = 0; f <= ny; ++f) {     // legacy j = 2..jd1
            const int jm1 = f - 1;

            const double ds = std::sqrt(sj(i, f).x * sj(i, f).x + sj(i, f).y * sj(i, f).y);
            const double nxn = sj(i, f).x / ds;
            const double nyn = sj(i, f).y / ds;

            const double rhol = dv(RHO, i, jm1), ul = dv(U, i, jm1),
                         vl = dv(V, i, jm1), pl = dv(P, i, jm1);
            const double El = pl / gam1 + 0.5 * rhol * (ul * ul + vl * vl);
            const double hl = pl / rhol + El;

            const double rhor = dv(RHO, i, f), ur = dv(U, i, f),
                         vr = dv(V, i, f), pr = dv(P, i, f);
            const double Er = pr / gam1 + 0.5 * rhor * (ur * ur + vr * vr);
            const double hr = pr / rhor + Er;

            jfaces.con_var_diff(0, i, f) = rhor - rhol;
            jfaces.con_var_diff(1, i, f) = rhor * ur - rhol * ul;
            jfaces.con_var_diff(2, i, f) = rhor * vr - rhol * vl;
            jfaces.con_var_diff(3, i, f) = Er - El;

            physical_flux(rhol, ul, vl, pl, nxn, nyn, gam1, fl);
            physical_flux(rhor, ur, vr, pr, nxn, nyn, gam1, fr);

            for (int k = 0; k < kNConv; ++k)
                jfaces.flux_diff(k, i, f) = fr[k] - fl[k];
        }
    }

    // ---- average-state flux (the surviving vcont definition) ---------------
    for (int j = 0; j < ny; ++j) {
        for (int f = 0; f <= nx; ++f) {
            const double rhoa = 0.5 * (dv(RHO, f - 1, j) + dv(RHO, f, j));
            const double rhoua =
                0.5 * (dv(U, f - 1, j) * dv(RHO, f - 1, j) + dv(U, f, j) * dv(RHO, f, j));
            const double rhova =
                0.5 * (dv(V, f - 1, j) * dv(RHO, f - 1, j) + dv(V, f, j) * dv(RHO, f, j));
            const double rhoea =
                0.5 * (dv(P, f - 1, j) / (eos.gamma - 1.0) +
                       0.5 * dv(RHO, f - 1, j) *
                           (dv(U, f - 1, j) * dv(U, f - 1, j) +
                            dv(V, f - 1, j) * dv(V, f - 1, j)) +
                       dv(P, f, j) / (eos.gamma - 1.0) +
                       0.5 * dv(RHO, f, j) *
                           (dv(U, f, j) * dv(U, f, j) + dv(V, f, j) * dv(V, f, j)));
            const double pa = 0.5 * (dv(P, f - 1, j) + dv(P, f, j));

            const double vcont = (rhoua * si(f, j).x + rhova * si(f, j).y) / rhoa;

            ifaces.avg_flux(0, f, j) = vcont * rhoa;
            ifaces.avg_flux(1, f, j) = vcont * rhoua + pa * si(f, j).x;
            ifaces.avg_flux(2, f, j) = vcont * rhova + pa * si(f, j).y;
            ifaces.avg_flux(3, f, j) = vcont * (rhoea + pa);
        }
    }

    for (int i = 0; i < nx; ++i) {
        for (int f = 0; f <= ny; ++f) {
            const double rhoa = 0.5 * (dv(RHO, i, f - 1) + dv(RHO, i, f));
            const double rhoua =
                0.5 * (dv(U, i, f - 1) * dv(RHO, i, f - 1) + dv(U, i, f) * dv(RHO, i, f));
            const double rhova =
                0.5 * (dv(V, i, f - 1) * dv(RHO, i, f - 1) + dv(V, i, f) * dv(RHO, i, f));
            const double rhoea =
                0.5 * (dv(P, i, f - 1) / (eos.gamma - 1.0) +
                       0.5 * dv(RHO, i, f - 1) *
                           (dv(U, i, f - 1) * dv(U, i, f - 1) +
                            dv(V, i, f - 1) * dv(V, i, f - 1)) +
                       dv(P, i, f) / (eos.gamma - 1.0) +
                       0.5 * dv(RHO, i, f) *
                           (dv(U, i, f) * dv(U, i, f) + dv(V, i, f) * dv(V, i, f)));
            const double pa = 0.5 * (dv(P, i, f - 1) + dv(P, i, f));

            const double vcont = (rhoua * sj(i, f).x + rhova * sj(i, f).y) / rhoa;

            jfaces.avg_flux(0, i, f) = vcont * rhoa;
            jfaces.avg_flux(1, i, f) = vcont * rhoua + pa * sj(i, f).x;
            jfaces.avg_flux(2, i, f) = vcont * rhova + pa * sj(i, f).y;
            jfaces.avg_flux(3, i, f) = vcont * (rhoea + pa);
        }
    }
}

void assemble_rhs_from_average_fluxes(const StructuredMesh& mesh,
                                      const MeshMetrics& metrics,
                                      const MultiField<double>& i_avg,
                                      const MultiField<double>& j_avg,
                                      const MultiField<double>& diss,
                                      MultiField<double>& rhs) {
    (void)metrics;
    const int nx = mesh.nx(), ny = mesh.ny();

    for (int k = 0; k < kNConv; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) rhs(k, i, j) = -diss(k, i, j);

    for (int j = 0; j < ny; ++j) {
        for (int f = 0; f <= nx; ++f) {
            for (int k = 0; k < kNConv; ++k) {
                rhs(k, f, j) += i_avg(k, f, j);
                rhs(k, f - 1, j) -= i_avg(k, f, j);
            }
        }
    }

    for (int i = 0; i < nx; ++i) {
        for (int f = 0; f <= ny; ++f) {
            for (int k = 0; k < kNConv; ++k) {
                rhs(k, i, f) += j_avg(k, i, f);
                rhs(k, i, f - 1) -= j_avg(k, i, f);
            }
        }
    }
}

void primitive_differences(const MultiField<double>& dv, MultiField<double>& dui,
                           MultiField<double>& duj) {
    constexpr int planes = 4;
    const int ng = dv.ng();

    // legacy Prim_Variables_Differences: dui over i=1..id2 per j=2..jb row,
    // then the leftmost ghost column copies its neighbour.
    for (int j = 0; j < dv.nj(); ++j) {
        for (int i = 1 - ng; i < dv.ni() + ng; ++i) {
            for (int k = 0; k < planes; ++k)
                dui(k, i, j) = dv(k, i, j) - dv(k, i - 1, j);
        }
        for (int k = 0; k < planes; ++k) dui(k, -ng, j) = dui(k, 1 - ng, j);
    }

    // duj over j=1..jd2 per i=2..ib column, then bottom ghost row copies.
    for (int i = 0; i < dv.ni(); ++i) {
        for (int j = 1 - ng; j < dv.nj() + ng; ++j) {
            for (int k = 0; k < planes; ++k)
                duj(k, i, j) = dv(k, i, j) - dv(k, i, j - 1);
        }
        for (int k = 0; k < planes; ++k) duj(k, i, -ng) = duj(k, i, 1 - ng);
    }
}

}  // namespace ns::numerics

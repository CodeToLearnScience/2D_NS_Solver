#pragma once

// Face flux definitions shared by all inviscid schemes (legacy FluxConvDefns):
// conserved-variable jumps, flux differences, and the surviving average-state
// (vcont-form) flux per face. Canvases are face-indexed ((ni+1)x(nj+1)).
//
// Divergence from legacy: the first average-flux definition (simple 0.5(Fl+Fr))
// is NOT computed -- FluxConvDefns unconditionally overwrites it later in the
// same call, so it was dead work.

#include "fields/field.hpp"
#include "mesh/mesh.hpp"
#include "physics/eos.hpp"

namespace ns::numerics {

inline constexpr int kNConv = 4;

struct FaceFluxData {
    FaceFluxData() = default;
    FaceFluxData(int nv, int nifaces, int njfaces, int ng)
        : avg_flux(nv, nifaces, njfaces, ng),
          flux_diff(nv, nifaces, njfaces, ng),
          con_var_diff(nv, nifaces, njfaces, ng) {}

    /// Canvas is (nx+1) x (ny+1), one entry per face.
    MultiField<double> avg_flux;      // average-state physical flux * face vector
    MultiField<double> flux_diff;     // F(Ur)-F(Ul)
    MultiField<double> con_var_diff;  // Ur-Ul
};

/// Computes both direction datasets (i-faces and j-faces).
void define_convective_faces(const StructuredMesh& mesh, const MeshMetrics& metrics,
                             const MultiField<double>& dv, const physics::IdealGas& eos,
                             FaceFluxData& ifaces, FaceFluxData& jfaces);

/// Legacy Avg_Conv_Flux1: rhs = -diss + scattered average fluxes.
void assemble_rhs_from_average_fluxes(const StructuredMesh& mesh,
                                      const MeshMetrics& metrics,
                                      const MultiField<double>& i_avg,
                                      const MultiField<double>& j_avg,
                                      const MultiField<double>& diss,
                                      MultiField<double>& rhs);

/// Primitive differences (legacy Prim_Variables_Differences): dui[k][i][j] =
/// dv[k][i][j]-dv[k][i-1][j] on the cell canvas; boundary ghost rows are only
/// partially filled, exactly as in the legacy kernel (stale entries were never
/// read there either -- documented, preserved).
void primitive_differences(const MultiField<double>& dv, MultiField<double>& dui,
                           MultiField<double>& duj);

}  // namespace ns::numerics

#pragma once

// Face gradients of (u, v, T) via the legacy dual-control-volume Green-Gauss
// scheme, transcribed index-for-index for bit parity:
//
//   gradfi[plane][i][j] : gradient at the i-face between cells (i-1,j),(i,j)
//   gradfj[plane][i][j] : gradient at the j-face between cells (i,j-1),(i,j)
//
// Planes per direction: [0]=du/dx [1]=du/dy [2]=dv/dx [3]=dv/dy [4]=dT/dx
// [5]=dT/dy. Face canvases are (ni+1)x(nj+1). NOTE (legacy behaviour kept):
// the normalization pass leaves the top row of gradfi and the right column of
// gradfj un-divided by the dual volume -- the legacy consumers never read
// those entries.

#include "fields/field.hpp"
#include "mesh/mesh.hpp"

namespace ns::physics {

/// Number of gradient planes stored per face direction.
inline constexpr int kGradPlanes = 6;

/// Which local mesh edges coincide with global domain edges.  Serial runs
/// set both true; j-slab MPI ranks set them from HaloExchange neighbors.
struct FaceGradBoundaryFlags {
    bool bottom_is_global = true;
    bool top_is_global = true;
};

/// Optional corrected-ghost sj metrics and per-row normalization table for
/// j-face gradients (MPI ranks supply these for exact serial parity).
struct FaceGradTables {
    const Field<Vec2>* sj_corrected = nullptr;  ///< sj whose ghost rows -1/nf+1 hold true neighbor face metrics
    const Field<double>* rvol_j = nullptr;      ///< per-(i,jface) normalization factors, sized nx x (ny+1)
};

void green_gauss_face_gradients(const StructuredMesh& mesh, const MeshMetrics& metrics,
                                const MultiField<double>& dv,
                                MultiField<double>& gradfi, MultiField<double>& gradfj,
                                FaceGradBoundaryFlags flags = {},
                                const FaceGradTables* tables = nullptr);

}  // namespace ns::physics

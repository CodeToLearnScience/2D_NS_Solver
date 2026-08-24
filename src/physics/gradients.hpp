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

void green_gauss_face_gradients(const StructuredMesh& mesh, const MeshMetrics& metrics,
                                const MultiField<double>& dv,
                                MultiField<double>& gradfi, MultiField<double>& gradfj);

}  // namespace ns::physics

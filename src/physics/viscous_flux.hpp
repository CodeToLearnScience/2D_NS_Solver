#pragma once

// Viscous flux accumulation onto the RHS dissipation field (legacy
// Flux_Viscous), transcribed for bit parity. Face stresses are built from the
// Green-Gauss face gradients; contributions are scattered +/- across each
// face into diss planes [1..3] (momentum-x, momentum-y, energy).

#include "fields/field.hpp"
#include "mesh/mesh.hpp"
#include "physics/gradients.hpp"

namespace ns::physics {

void accumulate_viscous_flux(const StructuredMesh& mesh, const MeshMetrics& metrics,
                             const MultiField<double>& dv,
                             const MultiField<double>& gradfi,
                             const MultiField<double>& gradfj, MultiField<double>& diss);

}  // namespace ns::physics

#pragma once

// Pointwise dependent variables (primitive + transport planes) computed from
// conserved variables, filling the legacy dv layout:
//   [0] rho, [1] u, [2] v, [3] p, [4] T, [5] a, [6] mu, [7] k
//
// The whole ghosted canvas is filled (matching Dependent_Variables, which the
// driver calls over [0..id2]); positivity checks are the caller's policy.

#include "fields/field.hpp"
#include "physics/eos.hpp"
#include "physics/transport.hpp"

namespace ns::physics {

void compute_dependent_variables(const MultiField<double>& cv, MultiField<double>& dv,
                                 const IdealGas& eos, Formulation formulation,
                                 const TransportScaling& scaling);

/// Single-point variant (legacy Dependent_Variables_One), for boundary kernels.
void compute_dependent_variables_at(const MultiField<double>& cv, MultiField<double>& dv,
                                    int i, int j, const IdealGas& eos,
                                    Formulation formulation,
                                    const TransportScaling& scaling);

}  // namespace ns::physics

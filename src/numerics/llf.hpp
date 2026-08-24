#pragma once

// Lax-Friedrichs (LLF) dissipation ports.
//
// NOTE: the legacy code has TWO different MaxEigVal definitions -- one in
// diss_llf1.cpp (|vn| +/- a) and a different one in diss_llf2_prim.cpp
// (|vn +/- a|, outer fabs). Both are preserved faithfully as separate
// functions because scheme parity depends on it.

#include "fields/field.hpp"
#include "mesh/mesh.hpp"
#include "physics/eos.hpp"

namespace ns::numerics {

/// Legacy MaxEigVal from diss_llf1.cpp.
[[nodiscard]] double max_eigenvalue_llf1(double ur, double ul, double vr, double vl,
                                         double ar, double al, double nx, double ny);

/// Legacy MaxEigVal from diss_llf2_prim.cpp.
[[nodiscard]] double max_eigenvalue_llf2(double ur, double ul, double vr, double vl,
                                         double ar, double al, double nx, double ny);

/// Legacy Minmod_Flux limiter: r = dp/(dm + 1e-16), clamped to [0, 1].
[[nodiscard]] double minmod_flux_limiter(double del_plus, double del_minus);

/// Diss_LLF1: first-order LLF dissipation accumulated onto `diss` using
/// precomputed face conserved-variable jumps.
void llf_dissipation_first_order(const StructuredMesh& mesh, const MeshMetrics& metrics,
                                 const MultiField<double>& dv,
                                 const MultiField<double>& i_con_var_diff,
                                 const MultiField<double>& j_con_var_diff,
                                 const physics::IdealGas& eos, MultiField<double>& diss);

/// Diss_LLF2_Prim: second-order MUSCL LLF dissipation from primitive
/// differences (dui/duj from primitive_differences).
void llf_dissipation_second_order(const StructuredMesh& mesh, const MeshMetrics& metrics,
                                  const MultiField<double>& dv,
                                  const MultiField<double>& dui,
                                  const MultiField<double>& duj,
                                  const physics::IdealGas& eos,
                                  MultiField<double>& diss);

}  // namespace ns::numerics

#pragma once

// Inviscid scheme dissipation ports for the baseline-critical schemes:
//   - MOVERS/NKFDS family (legacy schemes 7 and 30)
//   - two-wave Roe with Van Albada MUSCL (legacy scheme 24)
//
// Faithfully preserved legacy quirks (see REFACTORING_PLAN.md):
//   * a third MaxEigVal variant without the "+a" term, plus MinEigVal
//   * Movers() rounds flux/state jumps to 3 decimals before comparing
//   * the Mahalanobis entropy switch takes cp as a parameter -- legacy read
//     the global Cp (= 1/((gamma-1)*M_inf^2) in ND runs)

#include "fields/field.hpp"
#include "mesh/mesh.hpp"
#include "physics/eos.hpp"

namespace ns::numerics {

[[nodiscard]] double van_albada_limiter(double del_plus, double del_minus);

/// Legacy MaxEigVal from the NKFDS/MOVERS files (no "+a" term).
[[nodiscard]] double max_eigenvalue_nkfds(double ur, double ul, double vr, double vl,
                                          double ar, double al, double nx, double ny);

/// Legacy MinEigVal.
[[nodiscard]] double min_eigenvalue_nkfds(double ur, double ul, double vr, double vl,
                                          double ar, double al, double nx, double ny);

/// Legacy rounding(dec, 3) used inside Movers().
[[nodiscard]] double round_to(double dec, double n);

/// Legacy Movers() wave-speed estimate on 3-decimal-rounded jumps.
[[nodiscard]] double movers_wave_speed(double fr, double fl, double ur, double ul,
                                       double l_max, double l_min);

/// Legacy Mahalanobis() entropy switch; returns -ric or 0 (NaN-safe: with
/// degenerate inputs the legacy comparisons evaluate false and return 0).
[[nodiscard]] double nkfds_entropy_switch(double rho_l, double rho_r, double un_l,
                                          double un_r, double p_l, double p_r,
                                          double ric, double e_max, double cp);

/// Scheme 7: first-order NKFDS-MOVERS dissipation.
void nkfds_movers_dissipation_first_order(const StructuredMesh& mesh,
                                          const MeshMetrics& metrics,
                                          const MultiField<double>& dv,
                                          const physics::IdealGas& eos, double cp,
                                          MultiField<double>& diss);

/// Scheme 30: second-order NKFDS-MOVERS dissipation (Minmod MUSCL).
void nkfds_movers_dissipation_second_order(const StructuredMesh& mesh,
                                           const MeshMetrics& metrics,
                                           const MultiField<double>& dv,
                                           const MultiField<double>& dui,
                                           const MultiField<double>& duj,
                                           const physics::IdealGas& eos, double cp,
                                           MultiField<double>& diss);

/// Scheme 24: two-wave Roe dissipation with Van Albada MUSCL reconstruction.
void roe_dissipation_second_order(const StructuredMesh& mesh, const MeshMetrics& metrics,
                                  const MultiField<double>& dv,
                                  const MultiField<double>& dui,
                                  const MultiField<double>& duj,
                                  const physics::IdealGas& eos, MultiField<double>& diss);

}  // namespace ns::numerics

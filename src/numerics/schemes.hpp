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

#include <limits>

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

/// Scans the legacy emax window ([0..nx] x [0..ny)) of total energy.
double compute_emax_window(const StructuredMesh& mesh, const MultiField<double>& dv,
                           const physics::IdealGas& eos);

/// Scheme 7: first-order NKFDS-MOVERS dissipation. `emax_override` (NaN =
/// scan locally) lets an MPI driver supply the global window maximum.
void nkfds_movers_dissipation_first_order(const StructuredMesh& mesh,
                                          const MeshMetrics& metrics,
                                          const MultiField<double>& dv,
                                          const physics::IdealGas& eos, double cp,
                                          MultiField<double>& diss,
                                          double emax_override = std::numeric_limits<double>::quiet_NaN());

/// Scheme 30: second-order NKFDS-MOVERS dissipation (Minmod MUSCL).
void nkfds_movers_dissipation_second_order(const StructuredMesh& mesh,
                                           const MeshMetrics& metrics,
                                           const MultiField<double>& dv,
                                           const MultiField<double>& dui,
                                           const MultiField<double>& duj,
                                           const physics::IdealGas& eos, double cp,
                                           MultiField<double>& diss,
                                           double emax_override = std::numeric_limits<double>::quiet_NaN());

/// Scheme 24: two-wave Roe dissipation with Van Albada MUSCL reconstruction.
///
/// Batch A additions (Phase 9). Each preserves its own legacy wave-speed /
/// eigenvalue variant exactly. H1 and LE1 fix a legacy bug where pressure was
/// read from dv[0] (density) and sound speed from dv[2] (v-velocity);
/// bitwise parity with legacy is therefore NOT attempted for those two —
/// they are validated via free-stream preservation and analytic-face tests.

/// Legacy Movers() variant from diss_movers.cpp (scheme 1): exact Ur!=Ul
/// comparison, no epsilon, no rounding.
[[nodiscard]] double movers_plain_wave_speed(double fr, double fl, double ur,
                                             double ul, double l_max, double l_min);

/// Legacy Movers() variant from diss_movers_h1.cpp (eps = 1e-10).
[[nodiscard]] double movers_h1_wave_speed(double fr, double fl, double ur,
                                          double ul, double l_max, double l_min);

/// Legacy Movers() variant from diss_movers_le1.cpp (eps = 1e-4).
[[nodiscard]] double movers_le1_wave_speed(double fr, double fl, double ur,
                                           double ul, double l_max, double l_min);

/// Scheme 1: MOVERS with Minmod-MUSCL reconstruction; single wave speed
/// evaluated from the energy-component flux/state jumps, applied to all k.
void movers_dissipation_muscl(const StructuredMesh& mesh,
                              const MeshMetrics& metrics,
                              const MultiField<double>& dv,
                              const MultiField<double>& dui,
                              const MultiField<double>& duj,
                              const physics::IdealGas& eos,
                              MultiField<double>& diss);

/// Scheme 2: MOVERS-H1 — first-order, per-component Movers wave speed plus
/// beta*phi*max_eig extra dissipation (beta=0.2 pressure-switch).
void movers_h1_dissipation(const StructuredMesh& mesh, const MeshMetrics& metrics,
                           const MultiField<double>& dv, const physics::IdealGas& eos,
                           MultiField<double>& diss);

/// Scheme 3: MOVERS-LE1 — first-order, per-component Movers blended toward
/// max(|Vn|) by conserved-difference MinMod weighting.
void movers_le1_dissipation(const StructuredMesh& mesh, const MeshMetrics& metrics,
                            const MultiField<double>& cv, const MultiField<double>& dv,
                            const MultiField<double>& cui, const MultiField<double>& cuj,
                            const physics::IdealGas& eos, MultiField<double>& diss);

/// Conserved-variable differences (legacy Conv_Variables_Differences): same
/// stencil pattern as primitive_differences but on cv planes.
void conserved_differences(const MultiField<double>& cv, MultiField<double>& cui,
                           MultiField<double>& cuj);
void roe_dissipation_second_order(const StructuredMesh& mesh, const MeshMetrics& metrics,
                                  const MultiField<double>& dv,
                                  const MultiField<double>& dui,
                                  const MultiField<double>& duj,
                                  const physics::IdealGas& eos, MultiField<double>& diss);

}  // namespace ns::numerics

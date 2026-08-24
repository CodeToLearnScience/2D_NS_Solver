#pragma once

// Sutherland transport coefficients, transcribed verbatim from the legacy
// dependent-variables kernel (both nondimensional and dimensional forms).
//
// Legacy constants: C0 = 288.15 K, C1 = 110.4 K (Sutherland), and
//   ND:   mu(T) = [(1 + C1/C0) / (C1/C0 + T)] * T^{3/2} / Re_inf
//         k(T)  = mu / ((gamma-1) * Mach_inf^2 * Prandtl)
//   DIM:  mu(T) = [(C0 + C1) / (C1 + T)] * (T/C0)^{3/2} * ref_visc
//         k(T)  = mu * cp / Prandtl

namespace ns::physics {

inline constexpr double kSutherlandT0 = 288.15;  // legacy C0
inline constexpr double kSutherlandS = 110.4;    // legacy C1

enum class Formulation { Nondimensional, Dimensional };

struct TransportScaling {
    double re_inf = 1.0;      // used by the ND form
    double mach_inf = 1.0;    // used by the ND form
    double ref_visc = 0.0;    // used by the DIM form (mu at reference temperature)
    double cp = 1005.0;       // used by the DIM form
};

[[nodiscard]] double sutherland_viscosity_nd(double t_nd, const TransportScaling& s);

[[nodiscard]] double sutherland_viscosity_dim(double t, const TransportScaling& s);

[[nodiscard]] double conductivity(double mu, Formulation f,
                                  const TransportScaling& s, double gamma, double prandtl);

}  // namespace ns::physics

#pragma once

// Calorically perfect ideal-gas EOS. All legacy hard-coded gamma/R/Pr values
// become members with the same defaults, so ported results are unchanged.

#include <cmath>

namespace ns::physics {

/// Index helpers matching the legacy plane conventions (parity porting aid).
namespace cvi {
inline constexpr int RHO = 0, RHOU = 1, RHOV = 2, E = 3;
}
namespace dvi {
inline constexpr int RHO = 0, U = 1, V = 2, P = 3, T = 4, A = 5, MU = 6, K = 7;
}

struct PrimitiveState {
    double rho = 0.0, u = 0.0, v = 0.0, p = 0.0;
};

struct ConservedState {
    double rho = 0.0, rhou = 0.0, rhov = 0.0, e = 0.0;
};

class IdealGas {
public:
    double gamma = 1.4;
    double gas_constant = 286.92;
    double prandtl = 0.72;

    [[nodiscard]] constexpr double gamma_minus_1() const noexcept { return gamma - 1.0; }

    [[nodiscard]] constexpr double pressure(double rho, double rhou, double rhov,
                                            double energy) const noexcept {
        return gamma_minus_1() *
               (energy - 0.5 * (rhou * rhou + rhov * rhov) / rho);
    }

    [[nodiscard]] constexpr double sound_speed(double rho, double p) const noexcept {
        return std::sqrt(gamma * p / rho);
    }

    [[nodiscard]] constexpr double temperature(double rho, double p) const noexcept {
        return p / (rho * gas_constant);
    }

    [[nodiscard]] constexpr ConservedState conserved(const PrimitiveState& w) const noexcept {
        return {w.rho, w.rho * w.u, w.rho * w.v,
                w.p / gamma_minus_1() +
                    0.5 * w.rho * (w.u * w.u + w.v * w.v)};
    }

    [[nodiscard]] constexpr PrimitiveState primitive(const ConservedState& c) const noexcept {
        const double p = pressure(c.rho, c.rhou, c.rhov, c.e);
        return {c.rho, c.rhou / c.rho, c.rhov / c.rho, p};
    }
};

}  // namespace ns::physics

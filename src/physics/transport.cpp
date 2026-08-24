#include "physics/transport.hpp"

#include <cmath>

namespace ns::physics {

double sutherland_viscosity_nd(double t_nd, const TransportScaling& s) {
    const double c_ratio = kSutherlandS / kSutherlandT0;
    return ((1.0 + c_ratio) / (c_ratio + t_nd)) * std::pow(t_nd, 1.5) / s.re_inf;
}

double sutherland_viscosity_dim(double t, const TransportScaling& s) {
    return ((kSutherlandT0 + kSutherlandS) / (kSutherlandS + t)) *
           std::pow((t / kSutherlandT0), 1.5) * s.ref_visc;
}

double conductivity(double mu, Formulation f, const TransportScaling& s, double gamma,
                    double prandtl) {
    if (f == Formulation::Nondimensional)
        return mu / ((gamma - 1.0) * s.mach_inf * s.mach_inf * prandtl);
    return mu * s.cp / prandtl;
}

}  // namespace ns::physics

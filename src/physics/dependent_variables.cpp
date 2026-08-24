#include "physics/dependent_variables.hpp"

namespace ns::physics {

void compute_dependent_variables(const MultiField<double>& cv, MultiField<double>& dv,
                                 const IdealGas& eos, Formulation formulation,
                                 const TransportScaling& scaling) {
    for (int i = -cv.ng(); i < cv.ni() + cv.ng(); ++i) {
        for (int j = -cv.ng(); j < cv.nj() + cv.ng(); ++j) {
            compute_dependent_variables_at(cv, dv, i, j, eos, formulation, scaling);
        }
    }
}

void compute_dependent_variables_at(const MultiField<double>& cv, MultiField<double>& dv,
                                    int i, int j, const IdealGas& eos,
                                    Formulation formulation,
                                    const TransportScaling& scaling) {
    using cvi::E;
    using cvi::RHO;
    using cvi::RHOU;
    using cvi::RHOV;
    using dvi::A;
    using dvi::K;
    using dvi::MU;
    using dvi::P;
    using dvi::T;
    using dvi::U;
    using dvi::V;

    const double rho = cv(RHO, i, j);
    const double rhou = cv(RHOU, i, j);
    const double rhov = cv(RHOV, i, j);
    const double e = cv(E, i, j);

    const double p = eos.pressure(rho, rhou, rhov, e);

    dv(RHO, i, j) = rho;
    dv(U, i, j) = rhou / rho;
    dv(V, i, j) = rhov / rho;
    dv(P, i, j) = p;
    dv(T, i, j) = eos.temperature(rho, p);
    dv(A, i, j) = eos.sound_speed(rho, p);

    if (formulation == Formulation::Nondimensional) {
        dv(MU, i, j) = sutherland_viscosity_nd(dv(T, i, j), scaling);
    } else {
        dv(MU, i, j) = sutherland_viscosity_dim(dv(T, i, j), scaling);
    }
    dv(K, i, j) = conductivity(dv(MU, i, j), formulation, scaling, eos.gamma,
                               eos.prandtl);
}

}  // namespace ns::physics

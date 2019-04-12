#include "global_declarations.h"
#include "basic_functions.h"

void Limiter_Reference(double *&limref) {

    double temp, rho, cs, mach;

    flo_intr_extr = "external";

    if (flo_intr_extr == "external") {

        limref[0] = rhoinf;
        limref[1] = sqrt(uinf * uinf + vinf * vinf);
        limref[2] = limref[1];
        limref[3] = pinf;
    } else {

        temp = ttinl * pow(pout / ptinl, (Gamma - 1.0) / Gamma);
        rho = pout / (Rgas * temp);
        cs = sqrt(Gamma * pout / rho);
        mach = sqrt(2.0 * ((ttinl / temp) - 1) / (Gamma - 1.0));
        limref[0] = rho;
        limref[1] = mach * cs;
        limref[2] = limref[1];
        limref[3] = pout;
    }

    //std::cout<<"limiter reference values"<<"\t"<<limref[3]<<std::endl;
}

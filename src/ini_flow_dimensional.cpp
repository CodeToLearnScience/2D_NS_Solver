//
// Created by Maruthi NH on 13-07-2018.
//
#include "../inc/global_declarations.h"

void InitFlowDimensional(int id2, int jd2, double Re_inf, double Machinf, double Lref,  double alpha, double pinf,
                         double ***&cv, double ***&dv) {

    Cp      = Rgas*GamaByGamaMin1;
    rhoinf  = pinf/(Rgas*tinf);
    qinf    = Machinf*sqrt((Gamma-1.0)*Cp*tinf);
    alpha   = alpha/rad;
    uinf    = qinf*cos(alpha);
    vinf    = qinf*sin(alpha);
    rhoref  = rhoinf;
    velref  = qinf;
    refvisc = rhoinf*qinf*Lref/Re_inf;

    std::cout<<"the values are at Ic_cond:"<<rhoinf<<"\t"<<Rgas<<"\t"<<qinf<<"\t"<<uinf<<"\t"<<pinf<<"\t"<<refvisc<<
             std::endl;
    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            cv[0][i][j] = 0.0;
            cv[1][i][j] = 0.0;
            cv[2][i][j] = 0.0;
            cv[3][i][j] = 0.0;

            dv[0][i][j] = 0.0; //rho
            dv[1][i][j] = 0.0; //u
            dv[2][i][j] = 0.0; //v
            dv[3][i][j] = 0.0; //p
            dv[4][i][j] = 0.0; //T
            dv[5][i][j] = 0.0; //a
            dv[6][i][j] = 0.0; //Mu
            dv[7][i][j] = 0.0; //K
        }
    }

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            dv[0][i][j] = rhoinf; //rho
            dv[1][i][j] = uinf; //u
            dv[2][i][j] = vinf; //v
            dv[3][i][j] = pinf; //p
            dv[4][i][j] = tinf; //T
            dv[5][i][j] = sqrt(Gamma*pinf/rhoinf); //a
            dv[6][i][j] = ((C0+C1)/(C1+dv[4][i][j]))*pow((dv[4][i][j]/C0),1.5)*refvisc; //Mu
            dv[7][i][j] = dv[6][i][j]*Cp/Pr; //K

            /*std::cout <<"I.Cs: Dep Var"<<"\t"<<i<<"\t"<<j<<"\t"<<dv[0][i][j]<<"\t"<<dv[1][i][j]<<"\t"<< dv[2][i][j]
                      <<"\t"<<dv[3][i][j]<<"\t"<<dv[4][i][j]<<"\t"<<dv[5][i][j]<<"\t"<<dv[6][i][j]<<"\t"<<dv[7][i][j]
                      << std::endl;*/
        }
    }

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            cv[0][i][j] = dv[0][i][j];
            cv[1][i][j] = dv[0][i][j] * dv[1][i][j];
            cv[2][i][j] = dv[0][i][j] * dv[2][i][j];
            cv[3][i][j] = dv[3][i][j]/(Gamma-1.0)+0.5*dv[0][i][j]*(dv[1][i][j]*dv[1][i][j]+dv[2][i][j]*dv[2][i][j]);

            /*std::cout<<"I.Cs: Con Var" << "\t" << i<<"\t"<<j<<"\t"<<cv[0][i][j] << "\t" << cv[1][i][j]<< "\t"
                     << cv[2][i][j]<< "\t" << cv[3][i][j] << std::endl;*/
        }
    }
}

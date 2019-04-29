//
// Created by Maruthi NH on 13-07-2018.
//

#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"
void InitFlowNonDimensional(int id2, int jd2, double Re_inf, double Machinf, double Lref,  double alpha, double ***&cv,
                            double ***&dv) {

    if(eqn_type==0){
        std::cout<<"solving Euler equations in Non-dimensional form"<<std::endl;
        //std::cout<<"I am  here in initialization"<<std::endl;
        //rhoinf = 1.0, tinf=1.0, qinf = 1.0;
        rhoinf  = 1.0, tinf=1.0/1.4, qinf = Machinf;
        Rgas    = 1.0;
        Cp      = 1.0/((Gamma-1.0)*Machinf*Machinf);
        pinf    = rhoinf*Rgas*tinf;
        alpha   = alpha/rad;
        uinf    = qinf*cos(alpha);
        vinf    = qinf*sin(alpha);
    }
    else if(eqn_type==1){
        std::cout<<"solving Navier-Stokes equations in Non-dimensional form"<<std::endl;
        //std::cout<<"I am  here in initialization"<<std::endl;
        rhoinf  = 1.0, tinf=1.0, qinf = 1.0;
        refvisc = Lref/Re_inf;
        Rgas    = 1.0/(Gamma*Machinf*Machinf);
        Cp      = 1.0/((Gamma-1.0)*Machinf*Machinf);
        pinf    = rhoinf*Rgas*tinf;
        alpha   = alpha/rad;
        uinf    = qinf*cos(alpha);
        vinf    = qinf*sin(alpha);
    }
    else{
        std::cout<<"wrong equation type: please enter 0 for Euler and 1 for NSE"<<std::endl;
        exit(0);
    }

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
            dv[6][i][j] = ((1.0+C1/C0)/(C1/C0+dv[4][i][j]))*pow(dv[4][i][j],1.5)/Re_inf; //Mu
            dv[7][i][j] = dv[6][i][j]/(GammaMinus1*Machinf*Machinf*Pr); //K

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

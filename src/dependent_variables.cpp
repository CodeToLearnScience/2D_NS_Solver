#include "../inc/global_declarations.h"

void Dependent_Variables(int id2, int jd2, double ***&cv, double ***&dv) {

    double rhoq;

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {

            dv[0][i][j] = cv[0][i][j];
            dv[1][i][j] = cv[1][i][j]/cv[0][i][j];
            dv[2][i][j] = cv[2][i][j]/cv[0][i][j];
            rhoq = cv[1][i][j] * cv[1][i][j] + cv[2][i][j] * cv[2][i][j];
            dv[3][i][j] = GammaMinus1 * (cv[3][i][j] - 0.5 * rhoq / cv[0][i][j]);  //p
            dv[4][i][j] = dv[3][i][j] / (dv[0][i][j] * Rgas);                      //temp
            dv[5][i][j] = sqrt(Gamma * dv[3][i][j] / dv[0][i][j]);                 //sound speed
            if(dimen==0){
                dv[6][i][j] = ((1.0+C1/C0)/(C1/C0+dv[4][i][j]))*pow(dv[4][i][j],1.5)/Re_inf; //Mu
                dv[7][i][j] = dv[6][i][j]/(GammaMinus1*Machinf*Machinf*Pr); //K
            }
            else if(dimen==1){
                dv[6][i][j] = ((C0+C1)/(C1+dv[4][i][j]))*pow((dv[4][i][j]/C0),1.5)*refvisc; //Mu
                dv[7][i][j] = dv[6][i][j]*Cp/Pr; //K
            }
            /*std::cout <<"I.Cs: Dep Var"<<"\t"<<i<<"\t"<<j<<"\t"<<dv[0][i][j]<<"\t"<<dv[1][i][j]<<"\t"<< dv[2][i][j]
                      <<"\t"<<dv[3][i][j]<<"\t"<<dv[4][i][j]<<"\t"<<dv[5][i][j]<<"\t"<<dv[6][i][j]<<"\t"<<dv[7][i][j]
                      << std::endl;*/
        }
    }

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {

            if (dv[0][i][j] < 0) {
                std::cout << "in dependent var fun: density became negative at i=" << i << "\t" << "and at j=" << j
                << std::endl;
                std::cout << dv[0][i][j] << std::endl;
                exit(0);
            }
            else if (dv[3][i][j] < 0) {
                std::cout << "in dependent var fun: pressure became negative at i=" << i << "\t" << "and at j=" << j
                << std::endl;
                std::cout << dv[3][i][j] << std::endl;
                exit(0);
            }
        }
    }

}

void Dependent_Variables_One(int i, int j, double ***&cv, double ***&dv) {

    double rhoq=0.0;

    dv[0][i][j] = cv[0][i][j];
    dv[1][i][j] = cv[1][i][j]/cv[0][i][j];
    dv[2][i][j] = cv[2][i][j]/cv[0][i][j];
    rhoq = cv[1][i][j] * cv[1][i][j] + cv[2][i][j] * cv[2][i][j];
    dv[3][i][j] = GammaMinus1 * (cv[3][i][j] - 0.5 * rhoq / cv[0][i][j]);  //p
    dv[4][i][j] = dv[3][i][j] / (dv[0][i][j] * Rgas);                      //temp
    dv[5][i][j] = sqrt(Gamma * dv[3][i][j] / dv[0][i][j]);                 //sound speed
    if(dimen==0){
        dv[6][i][j] = ((1.0+C1/C0)/(C1/C0+dv[4][i][j]))*pow(dv[4][i][j],1.5)/Re_inf; //Mu
        dv[7][i][j] = dv[6][i][j]/(GammaMinus1*Machinf*Machinf*Pr); //K
    }
    else if(dimen==1){
        dv[6][i][j] = ((C0+C1)/(C1+dv[4][i][j]))*pow((dv[4][i][j]/C0),1.5)*refvisc; //Mu
        dv[7][i][j] = dv[6][i][j]*Cp/Pr; //K
    }

    if ((dv[0][i][j] != dv[0][i][j]) or (dv[1][i][j] != dv[1][i][j]) or (dv[2][i][j] != dv[2][i][j])) {
        std::cout << "The values at the boundary from Dependent_Variable_One fun: i=" << i << "\t" <<"j="<< j << "\t"
        << "rho="<<dv[0][i][j]<<"\t"<<"u="<<dv[1][i][j]<<"\t"<<"v="<<dv[2][i][j]<<"\t"<<"p="<<dv[3][i][j]<<std::endl;
        exit(0);
    }
}

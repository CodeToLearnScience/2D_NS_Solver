#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"


void BC_Eulerwall(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv) {

    double rho = 0.0, rhoe = 0.0;
    int dum1=0, dum2=0, ins1=0, ins2=0;
    maxwchg = 1.0;

    if (bind == 1 or bind == 4) {
        dum1 = bnode;
        dum2 = bnode - 1;
        ins1 = bnode + 1;
        ins2 = bnode + 2;
    } else if (bind == 2 or bind == 3) {
        dum1 = bnode;
        dum2 = bnode + 1;
        ins1 = bnode - 1;
        ins2 = bnode - 2;
    } else {
        std::cout << "Euler b.c., please check the boundary index, 1=bottom, 2=right, 3=top, 4=left and bind is =" <<bind<<std::endl;
        exit(0);
    }

    //std::cout<<"the dummy and inside nodes are for side: bind="<<bind<<"\t"<<"dum1="<<dum1<<"\t"<<"dum2="<<dum2<<"\t"<<"ins1="<<ins1<<"\t"<<"ins2="
    //<<ins2<<std::endl;

    if (bind == 1 or bind == 3) {
        for (int i = sbind; i <= ebind; i++) {
            /*rho = 2.0 * cv[0][i][ins1] - cv[0][i][ins2];
            rhoe = 2.0 * cv[3][i][ins1] - cv[3][i][ins2];
            if ((rho <= 0) or (rhoe <= 0) or ((fabs(cv[0][i][dum1]) - rho) > maxwchg * cv[0][i][dum1])
                or ((fabs(cv[3][i][dum1]) - rhoe) > maxwchg * cv[3][i][dum1])) {*/
            dv[0][i][dum1] = dv[0][i][ins1];
            dv[1][i][dum1] = dv[1][i][ins1];
            dv[2][i][dum1] = -dv[2][i][ins1];
            dv[3][i][dum1] = dv[3][i][ins1];
            dv[4][i][dum1] = dv[3][i][dum1] / (dv[0][i][dum1] * Rgas);                      //temp
            dv[5][i][dum1] = sqrt(Gamma * dv[3][i][dum1] / dv[0][i][dum1]);                 //sound speed
            if(dimen==0){
                dv[6][i][dum1] = ((1.0+C1/C0)/(C1/C0+dv[4][i][dum1]))*pow(dv[4][i][dum1],1.5)/Re_inf; //Mu
                dv[7][i][dum1] = dv[6][i][dum1]/(GammaMinus1*Machinf*Machinf*Pr); //K
            }
            else if(dimen==1){
                dv[6][i][dum1] = ((C0+C1)/(C1+dv[4][i][dum1]))*pow((dv[4][i][dum1]/C0),1.5)*refvisc; //Mu
                dv[7][i][dum1] = dv[6][i][dum1]*Cp/Pr; //K
            }
                /*std::cout << "I am here in this euler wall" << std::endl;
                exit(0);
            } else {
                cv[0][i][dum1] = rho;
                cv[1][i][dum1] = 2.0 * cv[1][i][ins1] - cv[1][i][ins2];
                cv[2][i][dum1] = 2.0 * cv[2][i][ins1] - cv[2][i][ins2];
                cv[3][i][dum1] = rhoe;
            }*/
            dv[0][i][dum2] = dv[0][i][dum1];
            dv[1][i][dum2] = dv[1][i][dum1];
            dv[2][i][dum2] = dv[2][i][dum1];
            dv[3][i][dum2] = dv[3][i][dum1];
            dv[4][i][dum2] = dv[4][i][dum1];
            dv[5][i][dum2] = dv[5][i][dum1];
            dv[6][i][dum2] = dv[6][i][dum1];
            dv[7][i][dum2] = dv[7][i][dum1];
        }
    }

    if (bind == 2 or bind == 4) {
        for (int j = sbind; j <= ebind; j++) {
            /*rho = 2.0 * cv[0][ins1][j] - cv[0][ins1][j];
            rhoe = 2.0 * cv[3][ins1][j] - cv[3][ins1][j];
            if ((rho <= 0) or (rhoe <= 0) or ((fabs(cv[0][dum1][j]) - rho) > maxwchg * cv[0][dum1][j])
                or ((fabs(cv[3][dum1][j]) - rhoe) > maxwchg * cv[3][dum1][j])) {*/
            dv[0][dum1][j] = dv[0][ins1][j]; //rho
            dv[1][dum1][j] = -dv[0][ins1][j]; //u
            dv[2][dum1][j] = dv[0][ins1][j]; //v
            dv[3][dum1][j] = dv[3][ins1][j]; //p
            dv[4][dum1][j] = dv[3][dum1][j] / (dv[0][dum1][j] * Rgas); //T
            dv[5][dum1][j] = sqrt(Gamma*dv[3][dum1][j]/dv[0][dum1][j]); //a
            if(dimen==0){
                dv[6][dum1][j] = ((1.0+C1/C0)/(C1/C0+dv[4][dum1][j]))*pow(dv[4][dum1][j],1.5)/Re_inf; //Mu
                dv[7][dum1][j] = dv[6][dum1][j]/(GammaMinus1*Machinf*Machinf*Pr); //K
            }
            else if(dimen==1){
                dv[6][dum1][j] = ((C0+C1)/(C1+dv[4][dum1][j]))*pow((dv[4][dum1][j]/C0),1.5)*refvisc; //Mu
                dv[7][dum1][j] = dv[6][dum1][j]*Cp/Pr; //K
            }
            /*} else {
                cv[0][dum1][j] = rho;
                cv[1][dum1][j] = 2.0 * cv[1][ins1][j] - cv[1][ins2][j];
                cv[2][dum1][j] = 2.0 * cv[2][ins1][j] - cv[2][ins2][j];
                cv[3][dum1][j] = rhoe;
            }*/
            dv[0][dum2][j] = dv[0][dum1][j];
            dv[1][dum2][j] = dv[1][dum1][j];
            dv[2][dum2][j] = dv[2][dum1][j];
            dv[3][dum2][j] = dv[3][dum1][j];
            dv[4][dum2][j] = dv[4][dum1][j];
            dv[5][dum2][j] = dv[5][dum1][j];
            dv[6][dum2][j] = dv[6][dum1][j];
            dv[7][dum2][j] = dv[7][dum1][j];
        }
    }

    /*for(int j=sbind;j<=ebind;j++){
        std::cout<<"the EULER_wall boundary values are="<<cv[3][idum1][j]<<std::endl;
    }*/
}






#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

void Bc_Transmitive(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv) {

    int dum1, dum2, ins1, ins2;

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
        std::cout << "Transmissive b.c., please check the boundary index, 1=bottom, 2=right, 3=top, 4=left"
                  << std::endl;
        exit(0);
    }
//std::cout<<"dum1 and dum2: inflow="<<dum1<<"\t"<<dum2<<std::endl;
    if (bind == 1 or bind == 3) {
        for (int i = sbind; i <= ebind; i++) {

            dv[0][i][dum1] = dv[0][i][ins1]; //rho
            dv[1][i][dum1] = dv[1][i][ins1]; //u
            dv[2][i][dum1] = dv[2][i][ins1]; //v
            dv[3][i][dum1] = dv[3][i][ins1]; //p
            dv[4][i][dum1] = dv[3][i][dum1] / (dv[0][i][dum1] * Rgas);  //T
            dv[5][i][dum1] = sqrt(Gamma*dv[3][i][dum1]/dv[0][i][dum1]); //a
            if(dimen==0){
                dv[6][i][dum1] = ((1.0+C1/C0)/(C1/C0+dv[4][i][dum1]))*pow(dv[4][i][dum1],1.5)/Re_inf; //Mu
                dv[7][i][dum1] = dv[6][i][dum1]/(GammaMinus1*Machinf*Machinf*Pr); //K
            }
            else if(dimen==1){
                dv[6][i][dum1] = ((C0+C1)/(C1+dv[4][i][dum1]))*pow((dv[4][i][dum1]/C0),1.5)*refvisc; //Mu
                dv[7][i][dum1] = dv[6][i][dum1]*Cp/Pr; //K
            }

            dv[0][i][dum2] = dv[0][i][dum1];
            dv[1][i][dum2] = dv[1][i][dum1];
            dv[2][i][dum2] = dv[2][i][dum1];
            dv[3][i][dum2] = dv[3][i][dum1];
            dv[4][i][dum2] = dv[4][i][dum1];
            dv[5][i][dum2] = dv[5][i][dum1];
            dv[6][i][dum2] = dv[6][i][dum1];
            dv[7][i][dum2] = dv[7][i][dum1];

        }
    } else if (bind == 2 or bind == 4) {
        for (int j = sbind; j <= ebind; j++) {

            dv[0][dum1][j] = dv[0][ins1][j]; //rho
            dv[1][dum1][j] = dv[1][ins1][j]; //u
            dv[2][dum1][j] = dv[2][ins1][j]; //v
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

            dv[0][dum2][j] = dv[0][dum1][j];
            dv[1][dum2][j] = dv[1][dum1][j];
            dv[2][dum2][j] = dv[2][dum1][j];
            dv[3][dum2][j] = dv[3][dum1][j];
            dv[4][dum2][j] = dv[4][dum1][j];
            dv[5][dum2][j] = dv[5][dum1][j];
            dv[6][dum2][j] = dv[6][dum1][j];
            dv[7][dum2][j] = dv[7][dum1][j];

  //          std::cout<<"values of cv[dum1] and cv[dum2]="<<cv[0][dum1][j]<<"\t"<<cv[0][dum2][j]<<std::endl;
            if ((dv[0][dum1][j] != dv[0][dum1][j]) or (dv[1][dum1][j] != dv[1][dum1][j]) or
                (dv[2][dum1][j] != dv[2][dum1][j]) or (dv[3][dum1][j] != dv[3][dum1][j])) {
                std::cout << "The values of conserved var in trans b. c. at i, j=" << dum1 << "\t" << j << "\t"
                          << dv[0][dum1][j]
                          << "\t" << dv[1][dum1][j] << "\t" << dv[2][dum1][j] << "\t" << dv[3][dum1][j]
                          << std::endl;
            }
        }
    }
}

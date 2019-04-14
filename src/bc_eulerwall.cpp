#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"


void BC_Eulerwall(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv) {

    double rho = 0.0, rhoe = 0.0;
    int dum1, dum2, ins1, ins2;
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
    if (bind == 1 or bind == 3) {
        for (int i = sbind; i <= ebind; i++) {
            /*rho = 2.0 * cv[0][i][ins1] - cv[0][i][ins2];
            rhoe = 2.0 * cv[3][i][ins1] - cv[3][i][ins2];
            if ((rho <= 0) or (rhoe <= 0) or ((fabs(cv[0][i][dum1]) - rho) > maxwchg * cv[0][i][dum1])
                or ((fabs(cv[3][i][dum1]) - rhoe) > maxwchg * cv[3][i][dum1])) {*/
                cv[0][i][dum1] = cv[0][i][ins1];
                cv[1][i][dum1] = cv[1][i][ins1];
                cv[2][i][dum1] = cv[2][i][ins1];
                cv[3][i][dum1] = cv[3][i][ins1];
                /*std::cout << "I am here in this euler wall" << std::endl;
                exit(0);
            } else {
                cv[0][i][dum1] = rho;
                cv[1][i][dum1] = 2.0 * cv[1][i][ins1] - cv[1][i][ins2];
                cv[2][i][dum1] = 2.0 * cv[2][i][ins1] - cv[2][i][ins2];
                cv[3][i][dum1] = rhoe;
            }*/
            cv[0][i][dum2] = cv[0][i][dum1];
            cv[1][i][dum2] = cv[1][i][dum1];
            cv[2][i][dum2] = cv[2][i][dum1];
            cv[3][i][dum2] = cv[3][i][dum1];

            /*cv[0][i][dum1] = 2.0*cv[0][i][ins1] - cv[0][i][ins2];
            cv[1][i][dum1] = 2.0*cv[1][i][ins1] - cv[1][i][ins2];
            cv[2][i][dum1] = 2.0*cv[2][i][ins1] - cv[2][i][ins2];
            cv[3][i][dum1] = 2.0*cv[3][i][ins1] - cv[3][i][ins2];*/

            /*cv[0][i][dum2] = 2.0 * cv[0][i][dum1] - cv[0][i][ins1];
            cv[1][i][dum2] = 2.0 * cv[1][i][dum1] - cv[1][i][ins1];
            cv[2][i][dum2] = 2.0 * cv[2][i][dum1] - cv[2][i][ins1];
            cv[3][i][dum2] = 2.0 * cv[3][i][dum1] - cv[3][i][ins1];*/

            Dependent_Variables_One(i, dum1, cv, dv);
            Dependent_Variables_One(i, dum2, cv, dv);
        }
    }

    if (bind == 2 or bind == 4) {
        for (int j = sbind; j <= ebind; j++) {
            /*rho = 2.0 * cv[0][ins1][j] - cv[0][ins1][j];
            rhoe = 2.0 * cv[3][ins1][j] - cv[3][ins1][j];
            if ((rho <= 0) or (rhoe <= 0) or ((fabs(cv[0][dum1][j]) - rho) > maxwchg * cv[0][dum1][j])
                or ((fabs(cv[3][dum1][j]) - rhoe) > maxwchg * cv[3][dum1][j])) {*/
                cv[0][dum1][j] = cv[0][ins1][j];
                cv[1][dum1][j] = cv[1][ins1][j];
                cv[2][dum1][j] = cv[2][ins1][j];
                cv[3][dum1][j] = cv[3][ins1][j];
            /*} else {
                cv[0][dum1][j] = rho;
                cv[1][dum1][j] = 2.0 * cv[1][ins1][j] - cv[1][ins2][j];
                cv[2][dum1][j] = 2.0 * cv[2][ins1][j] - cv[2][ins2][j];
                cv[3][dum1][j] = rhoe;
            }*/
            cv[0][dum2][j] = cv[0][dum1][j];
            cv[1][dum2][j] = cv[1][dum1][j];
            cv[2][dum2][j] = cv[2][dum1][j];
            cv[3][dum2][j] = cv[3][dum1][j];
            /*cv[0][dum1][j] = 2.0*cv[0][ins1][j] - cv[0][ins2][j];
            cv[1][dum1][j] = 2.0*cv[1][ins1][j] - cv[1][ins2][j];
            cv[2][dum1][j] = 2.0*cv[2][ins1][j] - cv[2][ins2][j];
            cv[3][dum1][j] = 2.0*cv[3][ins1][j] - cv[3][ins2][j];*/

            /*cv[0][dum2][j] = 2.0 * cv[0][dum1][j] - cv[0][ins1][j];
            cv[1][dum2][j] = 2.0 * cv[1][dum1][j] - cv[1][ins1][j];
            cv[2][dum2][j] = 2.0 * cv[2][dum1][j] - cv[2][ins1][j];
            cv[3][dum2][j] = 2.0 * cv[3][dum1][j] - cv[3][ins1][j];*/

            Dependent_Variables_One(dum1, j, cv, dv);
            Dependent_Variables_One(dum2, j, cv, dv);
        }
    }

    /*for(int j=sbind;j<=ebind;j++){
        std::cout<<"the EULER_wall boundary values are="<<cv[3][idum1][j]<<std::endl;
    }*/
}






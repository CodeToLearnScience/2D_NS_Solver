#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

void Bc_Prescribed_Inflow(int bind, int bnode, int sbind, int ebind, double rhoinf, double uinf, double vinf,
                          double pinf, double ***&cv, double ***&dv) {

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
        std::cout << "Inflo b.c., please check the boundary index, 1=bottom, 2=right, 3=top, 4=left and bind is =" <<bind<<std::endl;
        exit(0);
    }
        //std::cout<<"dum1 and dum2: inflow="<<dum1<<"\t"<<dum2<<"\t"<<rhoinf<<"\t"<<uinf<<"\t"<<vinf<<"\t"<<E_inf<<std::endl;
    //std::cout<<"values of cv[dum1] and cv[dum2]="<<cv[dum1][2]<<"\t"<<cv[dum2][2]<<std::endl;
    if (bind == 1 or bind == 3) {
        for (int i = sbind; i <= ebind; i++) {

            dv[0][i][dum1] = rhoinf; //rho
            dv[1][i][dum1] = uinf; //u
            dv[2][i][dum1] = vinf; //v
            dv[3][i][dum1] = pinf; //p
            dv[4][i][dum1] = dv[3][i][dum1] / (dv[0][i][dum1] * Rgas); //T
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

            dv[0][dum1][j] = rhoinf; //rho
            dv[1][dum1][j] = uinf; //u
            dv[2][dum1][j] = vinf; //v
            dv[3][dum1][j] = pinf; //p
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

            if ((dv[0][dum1][j] != dv[0][dum1][j]) or (dv[1][dum1][j] != dv[1][dum1][j]) or
                (dv[2][dum1][j] != dv[2][dum1][j]) or (dv[3][dum1][j] != dv[3][dum1][j])) {
                std::cout << "The values of conserved var at inflow at i, j=" << dum1 << "\t" << j << "\t"
                          << cv[0][dum1][j] << "\t" << cv[1][dum1][j] << "\t" << cv[2][dum1][j] << "\t"
                          << cv[3][dum1][j] << std::endl;
                exit(0);
            }
        }
    }
}


void Bc_Prescribed_Inflow2(int ib, int id1, int id2, int jb, int jd1, int jd2, double loc, double rhol, double rhor,
                           double ul, double ur, double vl, double vr, double pl, double pr, int bind, int sbind,
                           int ebind, double **&x, double **&y, double ***&cv) {

    double El, Er;
    int idum1, idum2, jdum1, jdum2, iins1, iins2, jins1, jins2;

    El = pl / (Gamma - 1.0) + 0.5 * rhol * (ul * ul + vl * vl);
    Er = pr / (Gamma - 1.0) + 0.5 * rhor * (ur * ur + vr * vr);

    if (bind == 1) {

        jdum1 = 1;
        jdum2 = 0;
        jins1 = 2;
        jins2 = 3;

    } else if (bind == 2) {

        idum1 = id1;
        idum2 = id2;
        iins1 = ib;
        iins2 = ib - 1;

    } else if (bind == 3) {

        jdum1 = jd1;
        jdum2 = jd2;
        jins1 = jb;
        jins2 = jb - 1;

    } else if (bind == 4) {

        idum1 = 1;
        idum2 = 0;
        iins1 = 2;
        iins2 = 3;
    }

    if (bind == 1 or bind == 3) {
        for (int i = sbind; i <= ebind; i++) {

            if (y[i][jdum1] > 0.5) {

                cv[0][i][jdum1] = rhol;
                cv[1][i][jdum1] = rhol * ul;
                cv[2][i][jdum1] = rhol * vl;
                cv[3][i][jdum1] = El;

                cv[0][i][jdum2] = 2.0 * cv[0][i][jdum1] - cv[0][i][jins1];
                cv[1][i][jdum2] = 2.0 * cv[1][i][jdum1] - cv[1][i][jins1];
                cv[2][i][jdum2] = 2.0 * cv[2][i][jdum1] - cv[2][i][jins1];
                cv[3][i][jdum2] = 2.0 * cv[3][i][jdum1] - cv[3][i][jins1];

                Dependent_Variables_One(i, jdum1, cv, dv);
                Dependent_Variables_One(i, jdum2, cv, dv);
            } else {

                cv[0][i][jdum1] = rhor;
                cv[1][i][jdum1] = rhor * ur;
                cv[2][i][jdum1] = rhor * vr;
                cv[3][i][jdum1] = Er;

                cv[0][i][jdum2] = 2.0 * cv[0][i][jdum1] - cv[0][i][jins1];
                cv[1][i][jdum2] = 2.0 * cv[1][i][jdum1] - cv[1][i][jins1];
                cv[2][i][jdum2] = 2.0 * cv[2][i][jdum1] - cv[2][i][jins1];
                cv[3][i][jdum2] = 2.0 * cv[3][i][jdum1] - cv[3][i][jins1];

                Dependent_Variables_One(i, jdum1, cv, dv);
                Dependent_Variables_One(i, jdum2, cv, dv);

            }
        }
    } else if (bind == 2 or bind == 4) {
        for (int j = sbind; j <= ebind; j++) {

            if (y[idum1][j] > 0.5) {
                cv[0][idum1][j] = rhol;
                cv[1][idum1][j] = rhol * ul;
                cv[2][idum1][j] = rhol * vl;
                cv[3][idum1][j] = El;

                cv[0][idum2][j] = 2.0 * cv[0][idum1][j] - cv[0][iins1][j];
                cv[1][idum2][j] = 2.0 * cv[1][idum1][j] - cv[1][iins1][j];
                cv[2][idum2][j] = 2.0 * cv[2][idum1][j] - cv[2][iins1][j];
                cv[3][idum2][j] = 2.0 * cv[3][idum1][j] - cv[3][iins1][j];

                Dependent_Variables_One(idum1, j, cv, dv);
                Dependent_Variables_One(idum1, j, cv, dv);
            } else {
                cv[0][idum1][j] = rhor;
                cv[1][idum1][j] = rhor * ur;
                cv[2][idum1][j] = rhor * vr;
                cv[3][idum1][j] = Er;

                cv[0][idum2][j] = 2.0 * cv[0][idum1][j] - cv[0][iins1][j];
                cv[1][idum2][j] = 2.0 * cv[1][idum1][j] - cv[1][iins1][j];
                cv[2][idum2][j] = 2.0 * cv[2][idum1][j] - cv[2][iins1][j];
                cv[3][idum2][j] = 2.0 * cv[3][idum1][j] - cv[3][iins1][j];

                Dependent_Variables_One(idum1, j, cv, dv);
                Dependent_Variables_One(idum1, j, cv, dv);
            }

        }
    }
}
//
// Created by DNS_Lab on 26-03-2018.
//
#include "global_declarations.h"

void BC_InFlow(int bind, int bnode, int sbind, int ebind, double ***&cv,
                double ***&dv) {

    double nx, ny, ds, vnor, vtan, rho, u, v, p, E, uabs, unorm, cosa, c02, cc02, rinv, dis, cb, tb, rhob, pb, ub, vb,
           uabsb;
    int dum1, dum2, ins1, ins2;
    maxwchg = 0.0;

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
        std::cout << "Inflo b.c., please check the boundary index, 1=bottom, 2=right, 3=top, 4=left" << std::endl;
        exit(0);
    }
    if (bind == 1 or bind == 3) {

        // if(sbind==2) sbind = 1;
        // if(ebind==ib) ebind = id1;

        for (int i = sbind; i <= ebind; i++) {
            if (bind == 1) {
                ds = sqrt(sj[0][i][ins1] * sj[0][i][ins1] + sj[1][i][ins1] * sj[1][i][ins1]);
                nx = sj[0][i][ins1] / ds;
                ny = sj[1][i][ins1] / ds;
            } else {
                ds = sqrt(sj[0][i][dum1] * sj[0][i][dum1] + sj[1][i][dum1] * sj[1][i][dum1]);
                nx = -sj[0][i][dum1] / ds;
                ny = -sj[1][i][dum1] / ds;

                /* ds = sqrt(sj[0][i][ins1] * sj[0][i][ins1] + sj[1][i][ins1] * sj[1][i][ins1]);
                 nx = -sj[0][i][ins1] / ds;
                 ny = -sj[1][i][ins1] / ds;*/
            }

            rho = dv[0][i][ins1];
            u   = dv[1][i][ins1];
            v   = dv[2][i][ins1];
            p   = dv[3][i][ins1];

            uabs  = sqrt(u*u + v*v);
            unorm = u*nx + v*ny;

            if (uabs < 1e-20)cosa = 1.0;
            else cosa = -unorm/uabs;

            c02  = dv[5][i][ins1]*dv[5][i][ins1] + GammaMinus1*uabs*uabs;
            rinv = unorm - 2.0*dv[5][i][ins1]/GammaMinus1;
            dis  = (GammaMinus1*cosa*cosa+2.0)*c02/(GammaMinus1*rinv*rinv)-0.5*GammaMinus1;

            if (dis<0.0){
                std::cout<<"dis value in Inflow B.C="<<dis<<std::endl;
                dis = 1e-20;
            }

            cb    = -rinv*(GammaMinus1/(GammaMinus1*cosa*cosa+2.0))*(1.0+cosa*sqrt(dis));
            cc02  = std::min(cb*cb/c02, 1.0);
            tb    = ttinl * cc02;
            pb    = ptinl * (tb/ttinl)*GamaByGamaMin1;
            rhob  = pb/(Rgas*tb);
            uabsb = 2.0*Cp*(ttinl-tb);
            uabsb = sqrt(uabsb);
            ub    = uabsb*cos(betainl);
            vb    = uabsb*sin(betainl);

            dv[0][i][dum1] = rhob;
            dv[1][i][dum1] = ub;
            dv[2][i][dum1] = vb;
            dv[3][i][dum1] = pb;
            dv[4][i][dum1] = tb;                      //temp
            dv[5][i][dum1] = sqrt(Gamma * dv[3][i][dum1] / dv[0][i][dum1]);                 //sound speed
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

            /*std::cout<<"NS wall boundary values are at bind="<<bind<<"\t"<<"dum1="<<dum1<<"\t"<<"i="<<i<<"\t"<<
                     cv[0][i][dum1]<<"\t"<<cv[1][i][dum1]<<"\t"<<cv[2][i][dum1]<<"\t"<<cv[3][i][dum1]<<"\t"<<
                     cv[0][i][ins1]<<"\t"<<cv[1][i][ins1]<<"\t"<<cv[2][i][ins1]<<"\t"<<cv[3][i][ins1]<<std::endl;*/
        }
    }

    if (bind == 2 or bind == 4) {

        // if(sbind==2) sbind = 1;
        // if(ebind==jb) ebind = jd1;

        for (int j = sbind; j <= ebind; j++) {

            if (bind == 4) {
                ds = sqrt(si[0][ins1][j] * si[0][ins1][j] + si[1][ins1][j] * si[1][ins1][j]);
                nx = si[0][ins1][j] / ds;
                ny = si[1][ins1][j] / ds;
            } else {
                ds = sqrt(si[0][dum1][j] * si[0][dum1][j] + si[1][dum1][j] * si[1][dum1][j]);
                nx = -si[0][dum1][j] / ds;
                ny = -si[1][dum1][j] / ds;
            }

            rho = dv[0][ins1][j];
            u   = dv[1][ins1][j];
            v   = dv[2][ins1][j];
            p   = dv[3][ins1][j];

            uabs  = sqrt(u*u + v*v);
            unorm = u*nx + v*ny;

            if (uabs < 1e-20)cosa = 1.0;
            else cosa = -unorm/uabs;

            c02  = dv[5][ins1][j]*dv[5][ins1][j] + GammaMinus1*uabs*uabs;
            rinv = unorm - 2.0*dv[5][ins1][j]/GammaMinus1;
            dis  = (GammaMinus1*cosa*cosa+2.0)*c02/(GammaMinus1*rinv*rinv)-0.5*GammaMinus1;

            if (dis<0.0){
                std::cout<<"dis value in Inflow B.C="<<dis<<std::endl;
                dis = 1e-20;
            }

            cb    = -rinv*(GammaMinus1/(GammaMinus1*cosa*cosa+2.0))*(1.0+cosa*sqrt(dis));
            cc02  = std::min(cb*cb/c02, 1.0);
            tb    = ttinl * cc02;
            pb    = ptinl * (tb/ttinl)*GamaByGamaMin1;
            rhob  = pb/(Rgas*tb);
            uabsb = 2.0*Cp*(ttinl-tb);
            uabsb = sqrt(uabsb);
            ub    = uabsb*cos(betainl);
            vb    = uabsb*sin(betainl);

            dv[0][dum1][j] = rhob; //rho
            dv[1][dum1][j] = ub; //u
            dv[2][dum1][j] = vb; //v
            dv[3][dum1][j] = pb; //p
            dv[4][dum1][j] = tb; //T
            dv[5][dum1][j] = sqrt(Gamma*pinf/rhoinf); //a
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

            /*std::cout<<"NS wall boundary values are at bind="<<bind<<"\t"<<"dum1="<<dum1<<"\t"<<"j="<<j<<"\t"<<
                     cv[0][dum1][j]<<"\t"<<cv[1][dum1][j]<<"\t"<<cv[2][dum1][j]<<"\t"<<cv[3][dum1][j]<<std::endl;*/
        }
    }
}
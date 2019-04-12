#include "global_declarations.h"
#include "basic_functions.h"

void BC_wall(int bind, int bnode, int sbind, int ebind, double ***&cv,
             double ***&dv) {
    double nx, ny, ds, vnor, vtan, u, v, p, E;
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
            vnor = (dv[1][i][ins1]*nx + dv[2][i][ins1]*ny);
            vtan = (dv[1][i][ins1]*ny - dv[2][i][ins1]*nx);
            //p    = 0.5*(3.0*dv[0][i][ins1]-dv[0][i][ins2]);

            vnor = -vnor;
            u    = vnor*nx + vtan*ny;
            v    = vnor*ny - vtan*nx;

            dv[0][i][dum1] = dv[0][i][ins1];
            dv[1][i][dum1] = u;
            dv[2][i][dum1] = v;
            //dv[1][i][dum1] = -dv[1][i][ins1];
            //dv[2][i][dum1] = -dv[2][i][ins1];
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

            vnor = (dv[1][ins1][j]*nx + dv[2][ins1][j]*ny);
            vtan = (dv[1][ins1][j]*ny - dv[2][ins1][j]*nx);
            //p    = 0.5*(3.0*dv[0][ins1][j]-dv[0][ins2][j]);
            vnor = -vnor;
            u    = vnor*nx + vtan*ny;
            v    = vnor*ny - vtan*nx;

            dv[0][dum1][j] = dv[0][ins1][j]; //rho
            dv[1][dum1][j] = u; //u
            dv[2][dum1][j] = v; //v
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

            /*std::cout<<"NS wall boundary values are at bind="<<bind<<"\t"<<"dum1="<<dum1<<"\t"<<"j="<<j<<"\t"<<
                     cv[0][dum1][j]<<"\t"<<cv[1][dum1][j]<<"\t"<<cv[2][dum1][j]<<"\t"<<cv[3][dum1][j]<<std::endl;*/
        }
    }
}

void BC_NS_wall(int bind, int bnode, int sbind, int ebind, double ***&cv,
             double ***&dv) {

    double nx, ny, ds, vnor, vtan, u, v, p, E;
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
            vnor = (dv[1][i][ins1]*nx + dv[2][i][ins1]*ny);
            vtan = (dv[1][i][ins1]*ny - dv[2][i][ins1]*nx);
            //p    = 0.5*(3.0*dv[0][i][ins1]-dv[0][i][ins2]);

            vnor = -vnor;
            vtan = -vtan;
            u    = vnor*nx + vtan*ny;
            v    = vnor*ny - vtan*nx;

            dv[0][i][dum1] = dv[0][i][ins1];
            dv[1][i][dum1] = u;
            dv[2][i][dum1] = v;
            //dv[1][i][dum1] = -dv[1][i][ins1];
            //dv[2][i][dum1] = -dv[2][i][ins1];
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

            vnor = (dv[1][ins1][j]*nx + dv[2][ins1][j]*ny);
            vtan = (dv[1][ins1][j]*ny - dv[2][ins1][j]*nx);
            //p    = 0.5*(3.0*dv[0][ins1][j]-dv[0][ins2][j]);
            vnor = -vnor;
            vtan = -vtan;
            u    = vnor*nx + vtan*ny;
            v    = vnor*ny - vtan*nx;

            dv[0][dum1][j] = dv[0][ins1][j]; //rho
            dv[1][dum1][j] = u; //u
            dv[2][dum1][j] = v; //v
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

            /*std::cout<<"NS wall boundary values are at bind="<<bind<<"\t"<<"dum1="<<dum1<<"\t"<<"j="<<j<<"\t"<<
                     cv[0][dum1][j]<<"\t"<<cv[1][dum1][j]<<"\t"<<cv[2][dum1][j]<<"\t"<<cv[3][dum1][j]<<std::endl;*/
        }
    }
}

void BC_NS_wallIsothermal(int bind, int bnode, int sbind, int ebind, double ***&cv,
                double ***&dv, double twall) {

    double nx, ny, ds, vnor, vtan, u, v, p, E;
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
            vnor = (dv[1][i][ins1]*nx + dv[2][i][ins1]*ny);
            vtan = (dv[1][i][ins1]*ny - dv[2][i][ins1]*nx);
            //p    = 0.5*(3.0*dv[0][i][ins1]-dv[0][i][ins2]);

            vnor = -vnor;
            vtan = -vtan;
            u    = vnor*nx + vtan*ny;
            v    = vnor*ny - vtan*nx;


            dv[1][i][dum1] = u;
            dv[2][i][dum1] = v;
            //dv[1][i][dum1] = -dv[1][i][ins1];
            //dv[2][i][dum1] = -dv[2][i][ins1];
            dv[3][i][dum1] = dv[3][i][ins1];
            dv[4][i][dum1] = 2.0*twall-dv[4][i][ins1];                      //temp
            dv[0][i][dum1] = dv[3][i][dum1]/(Rgas*dv[4][i][dum1]);  //density--cal frm T and p
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

            vnor = (dv[1][ins1][j]*nx + dv[2][ins1][j]*ny);
            vtan = (dv[1][ins1][j]*ny - dv[2][ins1][j]*nx);
            //p    = 0.5*(3.0*dv[0][ins1][j]-dv[0][ins2][j]);
            //p    = dv[3][ins1][j];
            vnor = -vnor;
            vtan = -vtan;
            u    = vnor*nx + vtan*ny;
            v    = vnor*ny - vtan*nx;

            dv[1][dum1][j] = u; //u
            dv[2][dum1][j] = v; //v
            dv[3][dum1][j] = dv[3][ins1][j]; //p
            //dv[4][dum1][j] = (2.0*twall-dv[4][ins1][j]);//2.0*twall-dv[4][ins1][j]; //T
            dv[4][dum1][j] = twall; //T
            dv[0][dum1][j] = dv[3][dum1][j]/(Rgas*dv[4][dum1][j]); //rho
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

            /*std::cout<<"NS wall boundary values are at bind="<<bind<<"\t"<<"dum1="<<dum1<<"\t"<<"j="<<j<<"\t"<<
                     cv[0][dum1][j]<<"\t"<<cv[1][dum1][j]<<"\t"<<cv[2][dum1][j]<<"\t"<<cv[3][dum1][j]<<std::endl;*/
        }
    }
}
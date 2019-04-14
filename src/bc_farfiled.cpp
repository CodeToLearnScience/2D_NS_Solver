//
// Created by maruthinh on 17/2/18.
//
#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"
void BC_FarField(int bind, int bnode, int sbind, int ebind, double ***&si, double ***&sj, double ***&cv,
                double ***&dv, double rhoinf, double uinf, double vinf, double pinf) {

    double nx=0.0, ny=0.0, ds, vnor, vtan, u, v, p, E, *rhof, *uf, *vf, *pf, bet, cir;
    int dum1, dum2, ins1, ins2, nbel;
    nbel = ebind-sbind+10;
    rhof = new double[nbel]; uf = new double[nbel]; vf = new double[nbel]; pf = new double[nbel];
    double gam1, ggm1, gmr, gmg, xa, ya, dist, angle, sn, dn, qv2, vc;
    double rhoe, ue, ve, pe, qqe, qn, crho0, rhoa, ua, va, pa, pb, ul, vl, pl, sgn;
    maxwchg = 0.0;

    if (bind == 1 or bind == 4){
        dum1 = bnode;
        dum2 = bnode - 1;
        ins1 = bnode + 1;
        ins2 = bnode + 2;
    } else if (bind == 2 or bind == 3){
        dum1 = bnode;
        dum2 = bnode + 1;
        ins1 = bnode - 1;
        ins2 = bnode - 2;
    } else {
        std::cout << "Inflo b.c., please check the boundary index, 1=bottom, 2=right, 3=top, 4=left" << std::endl;
        exit(0);
    }

    for(int i=sbind; i<=ebind; i++){
        rhof[i] = rhoinf;
        uf[i]   = uinf;
        vf[i]   = vinf;
        pf[i]   = pinf;

        /*std::cout<<"FarField B. C. values of rhof="<<rhof[i]<<"\t"<<"uf="<<uf[i]<<"\t"<<"vf="<<vf[i]<<"\t"
                 <<"pf="<<pf[i]<<std::endl;*/
    }

    /*Forces(bind, bnode, sbind, ebind, ndvar, x, y, dv, si, sj);

    bet = sqrt(1.0-Machinf*Machinf);

    cir = 0.25*Lref*cl*qinf/pi;
    if (bind == 1 or bind == 3) {
        for (int i = sbind; i <= ebind; i++){
            gam1 = Gamma-1.0;
            ggm1 = Gamma/gam1;
            gmr = 1.0/Gamma;
            gmg = gam1/Gamma;
            xa = x[i][ins1]-Xref;
            ya = y[i][ins1]-Yref;
            dist = sqrt(xa*xa+ya*ya);
            angle = atan2(xa, ya);
            sn = sin(angle - alpha);
            dn = 1-Machinf*Machinf*sn*sn;
            vc = cir*bet/(dn*dist);
            uf[i] = uinf + vc*sin(angle);
            vf[i] = vinf + vc*cos(angle);
            qv2 = uf[i]*uf[i] + vf[i]*vf[i];
            pf[i] = (pow(pinf, gmg)+gmg*rhoinf*(qinf*qinf-qv2))/(pow(pow(2.0*pinf,gmr),ggm1));
            rhof[i] = pow(rhoinf*(pf[i]/pinf), gmr);
        }
    }
    else if (bind == 2 or bind == 4) {
        for (int j = sbind; j <= ebind; j++){
                gam1 = Gamma-1.0;
                ggm1 = Gamma/gam1;
                gmr = 1.0/Gamma;
                gmg = gam1/Gamma;
                xa = x[ins1][j]-Xref;
                ya = y[ins1][j]-Yref;
                dist = sqrt(xa*xa+ya*ya);
                angle = atan2(xa, ya);
                sn = sin(angle - alpha);
                dn = 1-Machinf*Machinf*sn*sn;
                vc = cir*bet/(dn*dist);
                uf[j] = uinf + vc*sin(angle);
                vf[j] = vinf + vc*cos(angle);
                qv2 = uf[j]*uf[j] + vf[j]*vf[j];
                pf[j] = (pow(pinf, gmg)+gmg*rhoinf*(qinf*qinf-qv2))/(pow(pow(2.0*pinf,gmr),ggm1));
                rhof[j] = pow(rhoinf*(pf[j]/pinf), gmr);
        }
    }*/

/**computation of values at the boundary**/

    if(bind == 1 or bind == 3){
        for (int i = sbind; i <= ebind; i++){
            if (bind == 1) {
                ds = sqrt(sj[0][i][ins1] * sj[0][i][ins1] + sj[1][i][ins1] * sj[1][i][ins1]);
                nx = sj[0][i][ins1] / ds;
                ny = sj[1][i][ins1] / ds;
            } else {
                ds = sqrt(sj[0][i][dum1] * sj[0][i][dum1] + sj[1][i][dum1] * sj[1][i][dum1]);
                nx = -sj[0][i][dum1] / ds;
                ny = -sj[1][i][dum1] / ds;
            }
            gam1 = Gamma-1.0;
            rhoe = dv[0][i][ins1];
            ue = dv[1][i][ins1];
            ve = dv[2][i][ins1];
            qqe = ue*ue + ve*ve;
            pe = dv[3][i][ins1];

            if(Machinf<1.0){ //subsonic flow
                qn = ue*nx+ve*ny;
                crho0 = dv[5][i][ins1]*dv[0][i][ins1];
                if(qn<0.0){
                    rhoa = rhof[i];
                    ua = uf[i];
                    va = vf[i];
                    pa = pf[i];
                    ul = ue;
                    vl = ve;
                    pl = pe;
                    sgn =-1.0;
                    pb = 0.5*(pa+pl-crho0*(nx*(ua-ul)+ny*(va-vl)));
                }
                else{
                    rhoa = rhoe;
                    ua   = ue;
                    va   = ve;
                    pa   = pe;
                    ul   = uf[i];
                    vl   = vf[i];
                    pl   = pf[i];
                    sgn  = +1.0;
                    pb   = pf[i];
                }

                dv[0][i][dum1] = rhoa + (pb-pa)/pow(dv[5][i][ins1], 2.0);
                dv[1][i][dum1] = (ua+sgn*nx*(pa-pb)/crho0);
                dv[2][i][dum1] = (va+sgn*ny*(pa-pb)/crho0);
                dv[3][i][dum1] = pb;
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
            }
            else{ //supersonic flow
                qn = ue*nx + ve*ny;
                if(qn<0.0){
                    dv[0][i][dum1] = rhoinf;
                    dv[1][i][dum1] = uinf;
                    dv[2][i][dum1] = vinf;
                    dv[3][i][dum1] = pinf;
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
                }
                else{

                    dv[0][i][dum1] = rhoe;
                    dv[1][i][dum1] = ue;
                    dv[2][i][dum1] = ve;
                    dv[3][i][dum1] = pe;
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
                }
            }

            dv[0][i][dum2] = dv[0][i][dum1];
            dv[1][i][dum2] = dv[1][i][dum1];
            dv[2][i][dum2] = dv[2][i][dum1];
            dv[3][i][dum2] = dv[3][i][dum1];
            dv[4][i][dum2] = dv[4][i][dum1];
            dv[5][i][dum2] = dv[5][i][dum1];
            dv[6][i][dum2] = dv[6][i][dum1];
            dv[7][i][dum2] = dv[7][i][dum1];

            /*std::cout<<"the FarField boundary values are at bind="<<bind<<"\t"<<"dum1="<<dum1<<"\t"<<"i="<<i<<"\t"<<
                     cv[0][i][dum1]<<"\t"<<cv[1][i][dum1]<<"\t"<<cv[2][i][dum1]<<"\t"<<cv[3][i][dum1]<<std::endl;*/
        }
    }

    if (bind == 2 or bind == 4) {
        for (int j = sbind; j <= ebind; j++){
            if (bind == 4) {
                ds = sqrt(si[0][ins1][j] * si[0][ins1][j] + si[1][ins1][j] * si[1][ins1][j]);
                nx = si[0][ins1][j] / ds;
                ny = si[1][ins1][j] / ds;
            } else {
                ds = sqrt(si[0][dum1][j] * si[0][dum1][j] + si[1][dum1][j] * si[1][dum1][j]);
                nx = -si[0][dum1][j] / ds;
                ny = -si[1][dum1][j] / ds;
            }
            gam1  = Gamma-1.0;
            rhoe  = dv[0][ins1][j];
            ue    = dv[1][ins1][j];
            ve    = dv[2][ins1][j];
            pe    = dv[3][ins1][j];

            if(Machinf<1.0){ //subsonic flow
                qn = ue*nx+ve*ny;
                crho0 = dv[5][ins1][j]*dv[0][ins1][j];
                if(qn<0.0){
                    rhoa = rhof[j];
                    ua = uf[j];
                    va = vf[j];
                    pa = pf[j];
                    ul = ue;
                    vl = ve;
                    pl = pe;
                    sgn =-1.0;
                    pb = 0.5*(pa+pl-crho0*(nx*(ua-ul)+ny*(va-vl)));
                    //std::cout<<"In Far field, at bind="<<bind<<"\t"<<"j="<<j<<"\t"<<"pb="<<pb<<std::endl;
                }
                else{
                    rhoa = rhoe;
                    ua   = ue;
                    va   = ve;
                    pa   = pe;
                    ul   = uf[j];
                    vl   = vf[j];
                    pl   = pf[j];
                    sgn  = +1.0;
                    pb   = pl;
                }
                dv[0][dum1][j] = rhoa + (pb-pa)/pow(dv[5][ins1][j], 2.0);
                dv[1][dum1][j] = (ua+sgn*nx*(pa-pb)/crho0);
                dv[2][dum1][j] = (va+sgn*ny*(pa-pb)/crho0);
                dv[3][dum1][j] = pb;
                dv[4][dum1][j] = dv[3][dum1][j] / (dv[0][dum1][j] * Rgas);                      //temp
                dv[5][dum1][j] = sqrt(Gamma * dv[3][dum1][j] / dv[0][dum1][j]);                 //sound speed
                if(dimen==0){
                    dv[6][dum1][j] = ((1.0+C1/C0)/(C1/C0+dv[4][dum1][j]))*pow(dv[4][dum1][j],1.5)/Re_inf; //Mu
                    dv[7][dum1][j] = dv[6][dum1][j]/(GammaMinus1*Machinf*Machinf*Pr); //K
                }
                else if(dimen==1){
                    dv[6][dum1][j] = ((C0+C1)/(C1+dv[4][dum1][j]))*pow((dv[4][dum1][j]/C0),1.5)*refvisc; //Mu
                    dv[7][dum1][j] = dv[6][dum1][j]*Cp/Pr; //K
                }
            }
            else{ //supersonic flow
                qn = ue*nx + ve*ny;
                if(qn<0.0){
                    dv[0][dum1][j] = rhoinf;
                    dv[1][dum1][j] = uinf;
                    dv[2][dum1][j] = vinf;
                    dv[3][dum1][j] = pinf;
                    dv[4][dum1][j] = dv[3][dum1][j] / (dv[0][dum1][j] * Rgas);                      //temp
                    dv[5][dum1][j] = sqrt(Gamma * dv[3][dum1][j] / dv[0][dum1][j]);                 //sound speed
                    if(dimen==0){
                        dv[6][dum1][j] = ((1.0+C1/C0)/(C1/C0+dv[4][dum1][j]))*pow(dv[4][dum1][j],1.5)/Re_inf; //Mu
                        dv[7][dum1][j] = dv[6][dum1][j]/(GammaMinus1*Machinf*Machinf*Pr); //K
                    }
                    else if(dimen==1){
                        dv[6][dum1][j] = ((C0+C1)/(C1+dv[4][dum1][j]))*pow((dv[4][dum1][j]/C0),1.5)*refvisc; //Mu
                        dv[7][dum1][j] = dv[6][dum1][j]*Cp/Pr; //K
                    }
                }
                else{
                    dv[0][dum1][j] = rhoe;
                    dv[1][dum1][j] = ue;
                    dv[2][dum1][j] = ve;
                    dv[3][dum1][j] = pe;
                    dv[4][dum1][j] = dv[3][dum1][j] / (dv[0][dum1][j] * Rgas);                      //temp
                    dv[5][dum1][j] = sqrt(Gamma * dv[3][dum1][j] / dv[0][dum1][j]);                 //sound speed
                    if(dimen==0){
                        dv[6][dum1][j] = ((1.0+C1/C0)/(C1/C0+dv[4][dum1][j]))*pow(dv[4][dum1][j],1.5)/Re_inf; //Mu
                        dv[7][dum1][j] = dv[6][dum1][j]/(GammaMinus1*Machinf*Machinf*Pr); //K
                    }
                    else if(dimen==1){
                        dv[6][dum1][j] = ((C0+C1)/(C1+dv[4][dum1][j]))*pow((dv[4][dum1][j]/C0),1.5)*refvisc; //Mu
                        dv[7][dum1][j] = dv[6][dum1][j]*Cp/Pr; //K
                    }
                }

            }
            dv[0][dum2][j] = dv[0][dum1][j];
            dv[1][dum2][j] = dv[1][dum1][j];
            dv[2][dum2][j] = dv[2][dum1][j];
            dv[3][dum2][j] = dv[3][dum1][j];
            dv[4][dum2][j] = dv[4][dum1][j];
            dv[5][dum2][j] = dv[5][dum1][j];
            dv[6][dum2][j] = dv[6][dum1][j];
            dv[7][dum2][j] = dv[7][dum1][j];

            /*std::cout<<"FarField B. C. values are at bind="<<bind<<"\t"<<"dum1="<<dum1<<"\t"<<"j="<<j<<"\t"<<
                     cv[0][dum1][j]<<"\t"<<cv[1][dum1][j]<<"\t"<<cv[2][dum1][j]<<"\t"<<cv[3][dum1][j]<<std::endl;*/
        }
    }
delete rhof, uf, vf, pf;
}






/*
void BC_FarField(int bind, int bnode, int sbind, int ebind, double ***&si, double ***&sj, double ***&cv,
                 double ***&dv, double rhoinf, double uinf, double vinf, double pinf) {

    double nx=0.0, ny=0.0, ds, vnor, vtan, u, v, p, E, *rhof, *uf, *vf, *pf, bet, cir;
    int dum1, dum2, ins1, ins2, nbel;
    nbel = ebind-sbind+10;
    rhof = new double[nbel]; uf = new double[nbel]; vf = new double[nbel]; pf = new double[nbel];
    double gam1, ggm1, gmr, gmg, xa, ya, dist, angle, sn, dn, qv2, vc;
    double rhoe, ue, ve, pe, qqe, qn, crho0, rhoa, ua, va, pa, pb, ul, vl, pl, sgn;
    maxwchg = 0.0;

    if (bind == 1 or bind == 4){
        dum1 = bnode;
        dum2 = bnode - 1;
        ins1 = bnode + 1;
        ins2 = bnode + 2;
    } else if (bind == 2 or bind == 3){
        dum1 = bnode;
        dum2 = bnode + 1;
        ins1 = bnode - 1;
        ins2 = bnode - 2;
    } else {
        std::cout << "Inflo b.c., please check the boundary index, 1=bottom, 2=right, 3=top, 4=left" << std::endl;
        exit(0);
    }

    for(int i=sbind; i<=ebind; i++){
        rhof[i] = rhoinf;
        uf[i] = uinf;
        vf[i] = vinf;
        pf[i] = pinf;
    }

    Forces(bind, bnode, sbind, ebind, ndvar, x, y, dv, si, sj);

    bet = sqrt(1.0-Machinf*Machinf);

    cir = 0.25*Lref*cl*qinf/pi;
    if (bind == 1 or bind == 3) {
        for (int i = sbind; i <= ebind; i++){
            gam1 = Gamma-1.0;
            ggm1 = Gamma/gam1;
            gmr = 1.0/Gamma;
            gmg = gam1/Gamma;
            xa = x[i][ins1]-Xref;
            ya = y[i][ins1]-Yref;
            dist = sqrt(xa*xa+ya*ya);
            angle = atan2(xa, ya);
            sn = sin(angle - alpha);
            dn = 1-Machinf*Machinf*sn*sn;
            vc = cir*bet/(dn*dist);
            uf[i] = uinf + vc*sin(angle);
            vf[i] = vinf + vc*cos(angle);
            qv2 = uf[i]*uf[i] + vf[i]*vf[i];
            pf[i] = (pow(pinf, gmg)+gmg*rhoinf*(qinf*qinf-qv2))/(pow(pow(2.0*pinf,gmr),ggm1));
            rhof[i] = pow(rhoinf*(pf[i]/pinf), gmr);
        }
    }
    else if (bind == 2 or bind == 4) {
        for (int j = sbind; j <= ebind; j++){
                gam1 = Gamma-1.0;
                ggm1 = Gamma/gam1;
                gmr = 1.0/Gamma;
                gmg = gam1/Gamma;
                xa = x[ins1][j]-Xref;
                ya = y[ins1][j]-Yref;
                dist = sqrt(xa*xa+ya*ya);
                angle = atan2(xa, ya);
                sn = sin(angle - alpha);
                dn = 1-Machinf*Machinf*sn*sn;
                vc = cir*bet/(dn*dist);
                uf[j] = uinf + vc*sin(angle);
                vf[j] = vinf + vc*cos(angle);
                qv2 = uf[j]*uf[j] + vf[j]*vf[j];
                pf[j] = (pow(pinf, gmg)+gmg*rhoinf*(qinf*qinf-qv2))/(pow(pow(2.0*pinf,gmr),ggm1));
                rhof[j] = pow(rhoinf*(pf[j]/pinf), gmr);
        }
    }

//computation of values at the boundary

    if(bind == 1 or bind == 3){
        for (int i = sbind; i <= ebind; i++){
            if (bind == 1) {
                ds = sqrt(sj[0][i][ins1] * sj[0][i][ins1] + sj[1][i][ins1] * sj[1][i][ins1]);
                nx = sj[0][i][ins1] / ds;
                ny = sj[1][i][ins1] / ds;
            } else {
                ds = sqrt(sj[0][i][dum1] * sj[0][i][dum1] + sj[1][i][dum1] * sj[1][i][dum1]);
                nx = -sj[0][i][dum1] / ds;
                ny = -sj[1][i][dum1] / ds;
            }
            gam1 = Gamma-1.0;
            rhoe = cv[0][i][ins1];
            ue = cv[1][i][ins1]/rhoe;
            ve = cv[2][i][ins1]/rhoe;
            qqe = ue*ue + ve*ve;
            pe = dv[0][i][ins1];

            if(Machinf<1.0){ //subsonic flow
                qn = ue*nx+ve*ny;
                crho0 = dv[2][i][ins1]*cv[0][i][ins1];
                if(qn<0.0){
                    rhoa = rhof[i];
                    ua = uf[i];
                    va = vf[i];
                    pa = pf[i];
                    ul = ue;
                    vl = ve;
                    pl = pe;
                    sgn =-1.0;
                    pb = 0.5*(pa+pl-crho0*(nx*(ua-ul)+ny*(va-vl)));
                }
                else{
                    rhoa = rhoe;
                    ua   = ue;
                    va   = ve;
                    pa   = pe;
                    ul   = uf[i];
                    vl   = vf[i];
                    pl   = pf[i];
                    sgn  = +1.0;
                    pb   = pf[i];
                }
                cv[0][i][dum1] = rhoa + (pb-pa)/pow(dv[2][i][ins1], 2.0);
                cv[1][i][dum1] = cv[0][i][dum1]*(ua+sgn*nx*(pa-pb)/crho0);
                cv[2][i][dum1] = cv[0][i][dum1]*(va+sgn*ny*(pa-pb)/crho0);
                cv[3][i][dum1] = pb/gam1 + 0.5*(cv[1][i][dum1]*cv[1][i][dum1] + cv[2][i][dum1]*cv[2][i][dum1])/cv[0][i][dum1];
            }
            else{ //supersonic flow
                qn = ue*nx + ve*ny;
                if(qn<0.0){
                    cv[0][i][dum1] = rhoinf;
                    cv[1][i][dum1] = rhoinf*uinf;
                    cv[2][i][dum1] = rhoinf*vinf;
                    cv[3][i][dum1] = pinf/gam1 + 0.5*rhoinf*qinf*qinf;
                }
                else{
                    cv[0][i][dum1] = rhoe;
                    cv[1][i][dum1] = rhoe*ue;
                    cv[2][i][dum1] = rhoe*ve;
                    cv[3][i][dum1] = pe/gam1 + 0.5*rhoe*qqe;
                }
            }
            cv[0][i][dum2] = cv[0][i][dum1];
            cv[1][i][dum2] = cv[1][i][dum1];
            cv[2][i][dum2] = cv[2][i][dum1];
            cv[3][i][dum2] = cv[3][i][dum1];

            Dependent_Variables_One(i, dum1, cv, dv);
            Dependent_Variables_One(i, dum2, cv, dv);

        }
    }

    if (bind == 2 or bind == 4) {
        for (int j = sbind; j <= ebind; j++){
            if (bind == 4) {
                ds = sqrt(si[0][ins1][j] * si[0][ins1][j] + si[1][ins1][j] * si[1][ins1][j]);
                nx = si[0][ins1][j] / ds;
                ny = si[1][ins1][j] / ds;
            } else {
                ds = sqrt(si[0][dum1][j] * si[0][dum1][j] + si[1][dum1][j] * si[1][dum1][j]);
                nx = -si[0][dum1][j] / ds;
                ny = -si[1][dum1][j] / ds;
            }
            gam1 = Gamma-1.0;
            rhoe = cv[0][ins1][j];
            ue = cv[1][ins1][j]/rhoe;
            ve = cv[2][ins1][j]/rhoe;
            pe = dv[0][ins1][j];

            if(Machinf<1.0){ //subsonic flow
                qn = ue*nx+ve*ny;
                crho0 = dv[2][ins1][j]*cv[0][ins1][j];
                if(qn<0.0){
                    rhoa = rhof[j];
                    ua = uf[j];
                    va = vf[j];
                    pa = pf[j];
                    ul = ue;
                    vl = ve;
                    pl = pe;
                    sgn =-1.0;
                    pb = 0.5*(pa+pl-crho0*(nx*(ua-ul)+ny*(va-vl)));
                    //std::cout<<"In Far field, at bind="<<bind<<"\t"<<"j="<<j<<"\t"<<"pb="<<pb<<std::endl;
                }
                else{
                    rhoa = rhoe;
                    ua   = ue;
                    va   = ve;
                    pa   = pe;
                    ul   = uf[j];
                    vl   = vf[j];
                    pl   = pf[j];
                    sgn  = +1.0;
                    pb   = pl;
                }
                cv[0][dum1][j] = rhoa + (pb-pa)/pow(dv[2][ins1][j], 2.0);
                cv[1][dum1][j] = cv[0][dum1][j]*(ua+sgn*nx*(pa-pb)/crho0);
                cv[2][dum1][j] = cv[0][dum1][j]*(va+sgn*ny*(pa-pb)/crho0);
                cv[3][dum1][j] = pb/gam1 + 0.5*(cv[1][dum1][j]*cv[1][dum1][j] + cv[2][dum1][j]*cv[2][dum1][j])/cv[0][dum1][j];
            }
            else{ //supersonic flow
                qn = ue*nx + ve*ny;
                if(qn<0.0){
                    cv[0][dum1][j] = rhoinf;
                    cv[1][dum1][j] = rhoinf*uinf;
                    cv[2][dum1][j] = rhoinf*vinf;
                    cv[3][dum1][j] = pinf/gam1 + 0.5*rhoinf*qinf*qinf;
                }
                else{
                    cv[0][dum1][j] = rhoe;
                    cv[1][dum1][j] = rhoe*ue;
                    cv[2][dum1][j] = rhoe*ve;
                    cv[3][dum1][j] = pe/gam1 + 0.5*rhoe*qqe;
                }

            }
            cv[0][dum2][j] = cv[0][dum1][j];
            cv[1][dum2][j] = cv[1][dum1][j];
            cv[2][dum2][j] = cv[2][dum1][j];
            cv[3][dum2][j] = cv[3][dum1][j];

            Dependent_Variables_One(j, dum1, cv, dv);
            Dependent_Variables_One(j, dum2, cv, dv);

        }
    }
    delete rhof, uf, vf, pf;
}*/
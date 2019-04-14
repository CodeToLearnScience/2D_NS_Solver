//
// Created by Maruthi NH on 14-06-2018.
//

#include "../inc/global_declarations.h"

void FluxConvDefns(int id2, int jd2, int id1, int jd1, int ib, int jb, double ***&cv, double ***&dv,
                   double ***&iAvgFlux, double ***&jAvgFlux, double ***&iFluxDiff, double ***&jFluxDiff,
                   double ***&iConVarDiff, double ***&jConVarDiff){

    double rhoa, rhoua, rhova, rhoea, pa, vcont;
    int im1, jm1;
    double nx, ny, ds, rhoinv, rhol, ul, vl, rhor, ur, vr, pl, pr, El, Er, ar, al, hl, hr, gam1, ggm1;
    double *fr, *fl;

    fr = new double[nconv];
    fl = new double[nconv];

    gam1 = Gamma - 1.0;
    ggm1 = Gamma / gam1;

    //Conserved variable differences
        for (int j = 2; j <= jb; j++) {
            for (int i = 2; i <= id1; i++) {
                im1 = i-1;

                ds = sqrt(si[0][i][j] * si[0][i][j] + si[1][i][j] * si[1][i][j]);
                nx = si[0][i][j] / ds;
                ny = si[1][i][j] / ds;

                rhol = dv[0][im1][j];
                ul = dv[1][im1][j];
                vl = dv[2][im1][j];
                pl = dv[3][im1][j];
                El = pl / (gam1) + 0.5 * rhol * (ul * ul + vl * vl);
                hl = pl / rhol + El;
                al = sqrt(Gamma * pl / rhol);

                rhor = dv[0][i][j];
                ur = dv[1][i][j];
                vr = dv[2][i][j];
                pr = dv[3][i][j];
                Er = pr / (gam1) + 0.5 * rhor * (ur * ur + vr * vr);
                hr = pr / rhor + Er;
                ar = sqrt(Gamma * pr / rhor);

                iConVarDiff[0][i][j] = rhor - rhol;
                iConVarDiff[1][i][j] = rhor*ur - rhol*ul;
                iConVarDiff[2][i][j] = rhor*vr - rhol*vl;
                iConVarDiff[3][i][j] = Er - El;

                fl[0] = (rhol * ul) * nx + (rhol * vl) * ny;
                fl[1] = (rhol * ul * ul + pl) * nx + (rhol * ul * vl) * ny;
                fl[2] = (rhol * ul * vl) * nx + (rhol * vl * vl + pl) * ny;
                fl[3] = (rhol * ul * hl) * nx + (rhol * vl * hl) * ny;

                fr[0] = (rhor * ur) * nx + (rhor * vr) * ny;
                fr[1] = (rhor * ur * ur + pr) * nx + (rhor * ur * vr) * ny;
                fr[2] = (rhor * ur * vr) * nx + (rhor * vr * vr + pr) * ny;
                fr[3] = (rhor * ur * hr) * nx + (rhor * vr * hr) * ny;

                for(int k=0; k<nconv; k++){
                    iAvgFlux[k][i][j] = 0.5*(fl[k]+fr[k]);
                    iFluxDiff[k][i][j] = fr[k] - fl[k];
                }
            }
        }

    //Conserved variable differences
    for (int j = 2; j <= jd1; j++) {
        for (int i = 2; i <= ib; i++) {
            jm1 = j-1;

            ds = sqrt(sj[0][i][j] * sj[0][i][j] + sj[1][i][j] * sj[1][i][j]);
            nx = sj[0][i][j] / ds;
            ny = sj[1][i][j] / ds;

            rhol = dv[0][i][jm1];
            ul = dv[1][i][jm1];
            vl = dv[2][i][jm1];
            pl = dv[3][i][jm1];
            El = pl / (gam1) + 0.5 * rhol * (ul * ul + vl * vl);
            hl = pl / rhol + El;
            al = sqrt(Gamma * pl / rhol);

            rhor = dv[0][i][j];
            ur = dv[1][i][j];
            vr = dv[2][i][j];
            pr = dv[3][i][j];
            Er = pr / (gam1) + 0.5 * rhor * (ur * ur + vr * vr);
            hr = pr / rhor + Er;
            ar = sqrt(Gamma * pr / rhor);

            jConVarDiff[0][i][j] = rhor - rhol;
            jConVarDiff[1][i][j] = rhor*ur - rhol*ul;
            jConVarDiff[2][i][j] = rhor*vr - rhol*vl;
            jConVarDiff[3][i][j] = Er - El;

            fl[0] = (rhol * ul) * nx + (rhol * vl) * ny;
            fl[1] = (rhol * ul * ul + pl) * nx + (rhol * ul * vl) * ny;
            fl[2] = (rhol * ul * vl) * nx + (rhol * vl * vl + pl) * ny;
            fl[3] = (rhol * ul * hl) * nx + (rhol * vl * hl) * ny;

            fr[0] = (rhor * ur) * nx + (rhor * vr) * ny;
            fr[1] = (rhor * ur * ur + pr) * nx + (rhor * ur * vr) * ny;
            fr[2] = (rhor * ur * vr) * nx + (rhor * vr * vr + pr) * ny;
            fr[3] = (rhor * ur * hr) * nx + (rhor * vr * hr) * ny;

            for(int k=0; k<nconv; k++){
                jAvgFlux[k][i][j] = 0.5*(fl[k]+fr[k]);
                jFluxDiff[k][i][j] = fr[k] - fl[k];
            }
        }
    }

    /*for(int k=0;k<nconv; k++){
        for (int i = 2; i <= id1; i++) {
            iConVarDiff[k][i][1] = iConVarDiff[k][i][2];
            iConVarDiff[k][i][0] = iConVarDiff[k][i][1];

            jConVarDiff[k][i][1] = jConVarDiff[k][i][2];
            jConVarDiff[k][i][0] = jConVarDiff[k][i][1];

            iConVarDiff[k][i][id2] = iConVarDiff[k][i][id1];

            jConVarDiff[k][i][id2] = jConVarDiff[k][i][id1];
        }
    }

    for(int k=0;k<nconv; k++){
        for (int j = 0; j <= jd2; j++) {
            iConVarDiff[k][1][j] = iConVarDiff[k][2][j];
            iConVarDiff[k][0][j] = iConVarDiff[k][1][j];

            jConVarDiff[k][1][j] = jConVarDiff[k][2][j];
            jConVarDiff[k][0][j] = jConVarDiff[k][1][j];

            iConVarDiff[k][id2][j] = iConVarDiff[k][id1][j];

            jConVarDiff[k][id2][j] = jConVarDiff[k][id1][j];

        }
    }*/

    /*for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {

            std::cout<<"I.Cs: Con Var" << "\t" << i<<"\t"<<j<<"\t"<<iConVarDiff[0][i][j] << "\t"
                     << iConVarDiff[1][i][j]<< "\t" << iConVarDiff[2][i][j]<< "\t" << iConVarDiff[3][i][j] << std::endl;
        }
    }
    exit(0);*/

    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= id1; i++) {

            rhoa = 0.5 * (dv[0][i - 1][j] + dv[0][i][j]);
            rhoua = 0.5 * (dv[1][i - 1][j] * dv[0][i - 1][j] + dv[1][i][j] * dv[0][i][j]);
            rhova = 0.5 * (dv[2][i - 1][j] * dv[0][i - 1][j] + dv[2][i][j] * dv[0][i][j]);
            rhoea = 0.5 * (dv[3][i - 1][j] / (Gamma - 1.0) +
                           0.5 * dv[0][i - 1][j] * (dv[1][i - 1][j] * dv[1][i - 1][j] + dv[2][i - 1][j] *
                                                                                        dv[2][i - 1][j]) +
                           dv[3][i][j] / (Gamma - 1.0) + 0.5 * dv[0][i][j] * (dv[1][i][j] * dv[1][i][j] + dv[2][i][j] *
                                                                                                          dv[2][i][j]));
            pa = 0.5 * (dv[3][i - 1][j] + dv[3][i][j]);

            vcont = (rhoua * si[0][i][j] + rhova * si[1][i][j]) / rhoa;

            iAvgFlux[0][i][j] = vcont * rhoa;
            iAvgFlux[1][i][j] = vcont * rhoua + pa * si[0][i][j];
            iAvgFlux[2][i][j] = vcont * rhova + pa * si[1][i][j];
            iAvgFlux[3][i][j] = vcont * (rhoea + pa);
        }
    }
    for (int j = 2; j <= jd1; j++) {
        for (int i = 2; i <= ib; i++) {

           rhoa =  0.5 * (dv[0][i][j-1] + dv[0][i][j]);
           rhoua = 0.5 * (dv[1][i][j-1]*dv[0][i][j-1] + dv[1][i][j]*dv[0][i][j]);
           rhova = 0.5 * (dv[2][i][j-1]*dv[0][i][j-1] + dv[2][i][j]*dv[0][i][j]);
           rhoea = 0.5 * (dv[3][i][j-1]/(Gamma-1.0)+0.5*dv[0][i][j-1]*(dv[1][i][j-1]*dv[1][i][j-1]+dv[2][i][j-1] *
                   dv[2][i][j-1]) + dv[3][i][j]/(Gamma-1.0)+0.5*dv[0][i][j]*(dv[1][i][j]*dv[1][i][j]+dv[2][i][j] * dv[2][i][j]));
           pa =    0.5 * (dv[3][i][j-1] + dv[3][i][j]);

           vcont = (rhoua * sj[0][i][j] + rhova * sj[1][i][j]) / rhoa;

           jAvgFlux[0][i][j] = vcont * rhoa;
           jAvgFlux[1][i][j] = vcont * rhoua + pa * sj[0][i][j];
           jAvgFlux[2][i][j] = vcont * rhova + pa * sj[1][i][j];
           jAvgFlux[3][i][j] = vcont * (rhoea + pa);
       }
    }

    /*for(int k=0;k<nconv; k++){
        for (int i = 2; i <= ib; i++) {
            iAvgFlux[k][i][1] = iAvgFlux[k][i][2];
            iAvgFlux[k][i][0] = iAvgFlux[k][i][1];

            jAvgFlux[k][i][1] = jAvgFlux[k][i][2];
            jAvgFlux[k][i][0] = jAvgFlux[k][i][1];

            iAvgFlux[k][i][id1] = iAvgFlux[k][i][ib];
            iAvgFlux[k][i][id2] = iAvgFlux[k][i][id1];

            jAvgFlux[k][i][id1] = jAvgFlux[k][i][ib];
            jAvgFlux[k][i][id2] = jAvgFlux[k][i][id1];
        }
    }*/

    /*for(int k=0;k<nconv; k++){
        for (int j = 0; j <= jd2; j++) {
            iAvgFlux[k][1][j] = iAvgFlux[k][2][j];
            iAvgFlux[k][0][j] = iAvgFlux[k][1][j];

            jAvgFlux[k][1][j] = jAvgFlux[k][2][j];
            jAvgFlux[k][0][j] = jAvgFlux[k][1][j];

            iAvgFlux[k][id1][j] = iAvgFlux[k][ib][j];
            iAvgFlux[k][id2][j] = iAvgFlux[k][id1][j];

            jAvgFlux[k][id1][j] = jAvgFlux[k][ib][j];
            jAvgFlux[k][id2][j] = jAvgFlux[k][id1][j];

        }
    }

    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= ib; i++) {
            im1 = i - 1;

            ds = sqrt(si[0][i][j] * si[0][i][j] + si[1][i][j] * si[1][i][j]);
            nx = si[0][i][j] / ds;
            ny = si[1][i][j] / ds;

            //left and right states

            rhoinv = 1.0 / cv[0][im1][j];
            rhol = cv[0][im1][j];
            ul = cv[1][im1][j] * rhoinv;
            vl = cv[2][im1][j] * rhoinv;
            pl = dv[0][im1][j];
            hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);

            rhoinv = 1.0 / cv[0][i][j];
            rhor = cv[0][i][j];
            ur = cv[1][i][j] * rhoinv;
            vr = cv[2][i][j] * rhoinv;
            pr = dv[0][i][j];
            hr = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);

            fl[0] = (rhol * ul) * nx + (rhol * vl) * ny;
            fl[1] = (rhol * ul * ul + pl) * nx + (rhol * ul * vl) * ny;
            fl[2] = (rhol * ul * vl) * nx + (rhol * vl * vl + pl) * ny;
            fl[3] = (rhol * ul * hl) * nx + (rhol * vl * hl) * ny;

            fr[0] = (rhor * ur) * nx + (rhor * vr) * ny;
            fr[1] = (rhor * ur * ur + pr) * nx + (rhor * ur * vr) * ny;
            fr[2] = (rhor * ur * vr) * nx + (rhor * vr * vr + pr) * ny;
            fr[3] = (rhor * ur * hr) * nx + (rhor * vr * hr) * ny;

            for(int k=0;k<nconv;k++){
                iFluxDiff[k][i][j] = fr[k] - fl[k];
            }

            jm1 = j - 1;

            ds = sqrt(sj[0][i][j] * sj[0][i][j] + sj[1][i][j] * sj[1][i][j]);
            nx = sj[0][i][j] / ds;
            ny = sj[1][i][j] / ds;

            //left and right states

            rhoinv = 1.0 / cv[0][i][jm1];
            rhol = cv[0][i][jm1];
            ul = cv[1][i][jm1] * rhoinv;
            vl = cv[2][i][jm1] * rhoinv;
            pl = dv[0][i][jm1];
            hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);

            rhoinv = 1.0 / cv[0][i][j];
            rhor = cv[0][i][j];
            ur = cv[1][i][j] * rhoinv;
            vr = cv[2][i][j] * rhoinv;
            pr = dv[0][i][j];
            hr = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);

            fl[0] = (rhol * ul) * nx + (rhol * vl) * ny;
            fl[1] = (rhol * ul * ul + pl) * nx + (rhol * ul * vl) * ny;
            fl[2] = (rhol * ul * vl) * nx + (rhol * vl * vl + pl) * ny;
            fl[3] = (rhol * ul * hl) * nx + (rhol * vl * hl) * ny;

            fr[0] = (rhor * ur) * nx + (rhor * vr) * ny;
            fr[1] = (rhor * ur * ur + pr) * nx + (rhor * ur * vr) * ny;
            fr[2] = (rhor * ur * vr) * nx + (rhor * vr * vr + pr) * ny;
            fr[3] = (rhor * ur * hr) * nx + (rhor * vr * hr) * ny;

            for(int k=0;k<nconv;k++){
                jFluxDiff[k][i][j] = fr[k] - fl[k];
            }
        }
    }

    /*for(int k=0;k<nconv; k++){
        for (int i = 2; i <= ib; i++) {
            iFluxDiff[k][i][1] = iFluxDiff[k][i][2];
            iFluxDiff[k][i][0] = iFluxDiff[k][i][1];

            jFluxDiff[k][i][1] = jFluxDiff[k][i][2];
            jFluxDiff[k][i][0] = jFluxDiff[k][i][1];

            iFluxDiff[k][i][id1] = iFluxDiff[k][i][ib];
            iFluxDiff[k][i][id2] = iFluxDiff[k][i][id1];

            jFluxDiff[k][i][id1] = jFluxDiff[k][i][ib];
            jFluxDiff[k][i][id2] = jFluxDiff[k][i][id1];
        }
    }

    for(int k=0;k<nconv; k++){
        for (int j = 0; j <= jd2; j++) {
            iFluxDiff[k][1][j] = iFluxDiff[k][2][j];
            iFluxDiff[k][0][j] = iFluxDiff[k][1][j];

            jFluxDiff[k][1][j] = jFluxDiff[k][2][j];
            jFluxDiff[k][0][j] = jFluxDiff[k][1][j];

            iFluxDiff[k][id1][j] = iFluxDiff[k][ib][j];
            iFluxDiff[k][id2][j] = iFluxDiff[k][id1][j];

            jFluxDiff[k][id1][j] = jFluxDiff[k][ib][j];
            jFluxDiff[k][id2][j] = jFluxDiff[k][id1][j];
        }
    }*/
}


#include "global_declarations.h"
#include "basic_functions.h"

void
Flux_Boundary(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv, double ***&si, double ***&sj,
              double ***&rhs) {

    //computes convective fluxes at for boundaries  
    double sx, sy, rhol, ul, vl, rhor, ur, vr, pl, pr, hl, hr, gam1, ggm1, qsll, qslr, pavg;
    double *fc;
    int dum1, dum2, ins1, ins2;

    fc = new double[nconv];

    gam1 = Gamma - 1.0;
    ggm1 = Gamma / gam1;

    if (bind == 1 or bind == 4) {
        dum1 = bnode;
        ins1 = bnode + 1;
    }
    else if (bind == 2 or bind == 3) {
        dum1 = bnode;
        ins1 = bnode - 1;
    }
    else {
        std::cout << "Flux at wall:boundaries please check the boundary index, 1=bottom, 2=right, 3=top, 4=left and bind is =" <<bind<<std::endl;
        exit(0);
    }

    if (bind == 1 or bind == 3) {
        for (int i = sbind; i <= ebind; i++) {
            if (bind == 1) {
                sx = sj[0][i][ins1];
                sy = sj[1][i][ins1];
            } else {
                sx = -sj[0][i][dum1];
                sy = -sj[1][i][dum1];
            }

            rhol = dv[0][i][dum1];
            ul = dv[1][i][dum1];
            vl = dv[2][i][dum1];
            pl = dv[3][i][dum1];
            hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            qsll = (ul * sx + vl * sy) * rhol;

            rhor = dv[0][i][ins1];
            ur = dv[1][i][ins1];
            vr = dv[2][i][ins1];
            pr = dv[3][i][ins1];
            hr = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            qslr = (ur * sx + vr * sy) * rhor;

            pavg =  0.5 * (pl + pr);
            fc[0] = 0.5 * (qsll + qslr);
            fc[1] = 0.5 * (qsll * ul + qslr * ur) + pavg * sx;
            fc[2] = 0.5 * (qsll * vl + qslr * vr) + pavg * sy;
            fc[3] = 0.5 * (qsll * hl + qslr * hr);

            rhs[0][i][ins1] = rhs[0][i][ins1] + fc[0];
            rhs[1][i][ins1] = rhs[1][i][ins1] + fc[1];
            rhs[2][i][ins1] = rhs[2][i][ins1] + fc[2];
            rhs[3][i][ins1] = rhs[3][i][ins1] + fc[3];

        }

    } else if (bind == 2 or bind == 4) {
        for (int j = sbind; j <= ebind; j++) {
            if (bind == 4) {
                sx = si[0][ins1][j];
                sy = si[1][ins1][j];
            } else {
                sx = -si[0][dum1][j];
                sy = -si[1][dum1][j];
            }

            rhol = dv[0][dum1][j];
            ul = dv[1][dum1][j];
            vl = dv[2][dum1][j];
            pl = dv[3][dum1][j];
            hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            qsll = (ul * sx + vl * sy) * rhol;

            rhor = dv[0][ins1][j];
            ur = dv[1][ins1][j];
            vr = dv[2][ins1][j];
            pr = dv[3][ins1][j];
            hr = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            qslr = (ur * sx + vr * sy) * rhor;

            pavg =  0.5 * (pl + pr);
            fc[0] = 0.5 * (qsll + qslr);
            fc[1] = 0.5 * (qsll * ul + qslr * ur) + pavg * sx;
            fc[2] = 0.5 * (qsll * vl + qslr * vr) + pavg * sy;
            fc[3] = 0.5 * (qsll * hl + qslr * hr);

            rhs[0][ins1][j] = rhs[0][ins1][j] + fc[0];
            rhs[1][ins1][j] = rhs[1][ins1][j] + fc[1];
            rhs[2][ins1][j] = rhs[2][ins1][j] + fc[2];
            rhs[3][ins1][j] = rhs[3][ins1][j] + fc[3];

        }

    }

    /* for(int j=sbind;j<=ebind;j++){
             std::cout<<"Flux boundary values are="<<rhs[0][iins1][j]<<std::endl;
     }*/

    delete[] fc;
}

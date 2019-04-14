//
// Created by maruthinh on 16/3/18.
//

#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

void Avg_Flux1(int ib, int id1, int jb, int jd1, double ***&cv, double ***&dv, double ***&si, double ***&sj, double ***&diss, double ***&rhs) {

    double *fc;
    double pavg, vcont, rhol, qsl, rhor, qsr;

    fc = new double[nconv];

    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= ib; i++) {
            rhs[0][i][j] = -diss[0][i][j];
            rhs[1][i][j] = -diss[1][i][j];
            rhs[2][i][j] = -diss[2][i][j];
            rhs[3][i][j] = -diss[3][i][j];
        }
    }

    for (int j = 2; j <= jb; j++) {
        //flux in i direction except at the boundaries
        for (int i = 3; i <= ib; i++) {
            rhol = dv[0][i - 1][j] + cv[3][i - 1][j];
            qsl = (cv[1][i - 1][j] * si[0][i][j] + cv[2][i - 1][j] * si[1][i][j]) / cv[0][i - 1][j];

            rhor = dv[0][i][j] + cv[3][i][j];
            qsr = (cv[1][i][j] * si[0][i][j] + cv[2][i][j] * si[1][i][j]) / cv[0][i][j];

            pavg = 0.5 * (dv[0][i - 1][j] + dv[0][i][j]);

            fc[0] = 0.5 * (qsl * cv[0][i - 1][j] + qsr * cv[0][i][j]);
            fc[1] = 0.5 * (qsl * cv[1][i - 1][j] + qsr * cv[1][i][j]) + pavg * si[0][i][j];
            fc[2] = 0.5 * (qsl * cv[2][i - 1][j] + qsr * cv[2][i][j]) + pavg * si[1][i][j];
            fc[3] = 0.5 * (qsl * rhol + qsr * rhor);

            rhs[0][i][j] = rhs[0][i][j] + fc[0];
            rhs[1][i][j] = rhs[1][i][j] + fc[1];
            rhs[2][i][j] = rhs[2][i][j] + fc[2];
            rhs[3][i][j] = rhs[3][i][j] + fc[3];

            rhs[0][i - 1][j] = rhs[0][i - 1][j] - fc[0];
            rhs[1][i - 1][j] = rhs[1][i - 1][j] - fc[1];
            rhs[2][i - 1][j] = rhs[2][i - 1][j] - fc[2];
            rhs[3][i - 1][j] = rhs[3][i - 1][j] - fc[3];
        }
    }
    for (int j = 3; j <= jb; j++) {
        //flux in j-direction except at boundaries
        for (int i = 2; i <= ib; i++) {
            rhol = dv[0][i][j-1] + cv[3][i][j-1];
            qsl = (cv[1][i][j-1] * sj[0][i][j] + cv[2][i][j-1] * sj[1][i][j]) / cv[0][i][j-1];

            rhor = dv[0][i][j] + cv[3][i][j];
            qsr = (cv[1][i][j] * sj[0][i][j] + cv[2][i][j] * sj[1][i][j]) / cv[0][i][j];

            pavg = 0.5 * (dv[0][i][j-1] + dv[0][i][j]);

            fc[0] = 0.5 * (qsl * cv[0][i][j-1] + qsr * cv[0][i][j]);
            fc[1] = 0.5 * (qsl * cv[1][i][j-1] + qsr * cv[1][i][j]) + pavg * sj[0][i][j];
            fc[2] = 0.5 * (qsl * cv[2][i][j-1] + qsr * cv[2][i][j]) + pavg * sj[1][i][j];
            fc[3] = 0.5 * (qsl * rhol + qsr * rhor);

            rhs[0][i][j] = rhs[0][i][j] + fc[0];
            rhs[1][i][j] = rhs[1][i][j] + fc[1];
            rhs[2][i][j] = rhs[2][i][j] + fc[2];
            rhs[3][i][j] = rhs[3][i][j] + fc[3];

            rhs[0][i][j - 1] = rhs[0][i][j-1] - fc[0];
            rhs[1][i][j - 1] = rhs[1][i][j-1] - fc[1];
            rhs[2][i][j - 1] = rhs[2][i][j-1] - fc[2];
            rhs[3][i][j - 1] = rhs[3][i][j-1] - fc[3];
            }
        }


    Flux_Boundary(2, id1, 2, 75, cv, dv, si, sj, rhs);
    Flux_Boundary(2, id1, 76, jb, cv, dv, si, sj, rhs);
    Flux_Boundary(1, 1, 2, ib, cv, dv, si, sj, rhs);
    Flux_Boundary(3, id1, 2, ib, cv, dv, si, sj, rhs);
    Flux_ViscWall(4, 1, 2, jb, cv, dv, si, sj, rhs);

    delete[] fc;
}


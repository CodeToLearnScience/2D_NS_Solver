//
// Created by maruthinh on 24/2/18.
//

#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

void Imp_Res_Smoo(double **&dirs, double ***&epsij, double ***&rhs){

    int im1, ip1, jm1, jp1;
    double t;
    //soultion of the tridiagonal system in i-direction

    for (int j = 2; j <= jb; j++) {
        dirs[1][j] = 0.0;
        rhs[0][1][j] = 0.0;
        rhs[1][1][j] = 0.0;
        rhs[2][1][j] = 0.0;
        rhs[3][1][j] = 0.0;
    }


    for (int i = 2; i <= ib; i++) {
        im1 = i-1;
        for (int j = 2; j <= jb; j++) {
            t = 1.0/(1.0+2.0*epsij[0][i][j]-epsij[0][i][j]*dirs[im1][j]);
            dirs[i][j] = t*epsij[0][i][j];
            rhs[0][i][j] = t*(rhs[0][i][j] + epsij[0][i][j]*rhs[0][im1][j]);
            rhs[1][i][j] = t*(rhs[1][i][j] + epsij[0][i][j]*rhs[1][im1][j]);
            rhs[2][i][j] = t*(rhs[2][i][j] + epsij[0][i][j]*rhs[2][im1][j]);
            rhs[3][i][j] = t*(rhs[3][i][j] + epsij[0][i][j]*rhs[3][im1][j]);
        }
    }

    for (int i = ib-1; i <= 2; i--) {
        ip1 = i+1;
        for (int j = 2; j <= jb; j++) {
            rhs[0][i][j] = rhs[0][i][j] + dirs[i][j]*rhs[0][ip1][j];
            rhs[1][i][j] = rhs[1][i][j] + dirs[i][j]*rhs[1][ip1][j];
            rhs[2][i][j] = rhs[2][i][j] + dirs[i][j]*rhs[2][ip1][j];
            rhs[3][i][j] = rhs[3][i][j] + dirs[i][j]*rhs[3][ip1][j];
        }
    }

    //solution of tridiagonal system in j-direction
    for (int i = 2; i <= ib; i++) {
        dirs[i][1] = 0.0;
        rhs[0][i][1] = 0.0;
        rhs[1][i][1] = 0.0;
        rhs[2][i][1] = 0.0;
        rhs[3][i][1] = 0.0;
    }


    for (int j = 2; j <= jb; j++) {
        jm1 = j-1;
        for (int i = 2; i <= ib; i++) {
            t = 1.0/(1.0+2.0*epsij[1][i][j]-epsij[1][i][j]*dirs[i][jm1]);
            dirs[i][j] = t*epsij[1][i][j];
            rhs[0][i][j] = t*(rhs[0][i][j] + epsij[1][i][j]*rhs[0][i][jm1]);
            rhs[1][i][j] = t*(rhs[1][i][j] + epsij[1][i][j]*rhs[1][i][jm1]);
            rhs[2][i][j] = t*(rhs[2][i][j] + epsij[1][i][j]*rhs[2][i][jm1]);
            rhs[3][i][j] = t*(rhs[3][i][j] + epsij[1][i][j]*rhs[3][i][jm1]);
        }
    }

    for (int j = jb-1; j <= 2; j--) {
        jp1 = j+1;
        for (int i = 2; i <= ib; i++) {
            rhs[0][i][j] = rhs[0][i][j] + dirs[i][j]*rhs[0][i][jp1];
            rhs[1][i][j] = rhs[1][i][j] + dirs[i][j]*rhs[1][i][jp1];
            rhs[2][i][j] = rhs[2][i][j] + dirs[i][j]*rhs[2][i][jp1];
            rhs[3][i][j] = rhs[3][i][j] + dirs[i][j]*rhs[3][i][jp1];
        }
    }
}

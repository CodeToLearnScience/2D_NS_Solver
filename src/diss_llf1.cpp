#include "global_declarations.h"
#include "basic_functions.h"

template<typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny);

void Diss_LLF1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
               double ***&sj, double ***&diss) {

    int im1, jm1;
    double nx, ny, ds, rhol, ul, vl, rhor, ur, vr, pl, pr;
    double ar, al, max_eig;
    double *fd, *fr, *fl, *Ur, *Ul;

    fd = new double[nconv];
    fr = new double[nconv];
    fl = new double[nconv];
    Ur = new double[nconv];
    Ul = new double[nconv];

    for(int k=0; k<nconv; k++){
        for (int j = 2; j <= jb; j++) {
            for (int i = 2; i <= ib; i++) {
                rhs[k][i][j] = -diss[k][i][j];
            }
        }
    }

    //compute in i-direction
    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= id1; i++) {

            im1 = i - 1;

            ds = sqrt(si[0][i][j] * si[0][i][j] + si[1][i][j] * si[1][i][j]);
            nx = si[0][i][j] / ds;
            ny = si[1][i][j] / ds;

            //left and right states

            rhol = dv[0][im1][j];
            ul = dv[1][im1][j];
            vl = dv[2][im1][j];
            pl = dv[3][im1][j];
            al = sqrt(Gamma * pl / rhol);

            rhor = dv[0][i][j];
            ur = dv[1][i][j];
            vr = dv[2][i][j];
            pr = dv[3][i][j];
            ar = sqrt(Gamma * pr / rhor);

            max_eig = MaxEigVal(ur, ul, vr, vl, ar, al, nx, ny);

           //final dissipation terms
            for(int k=0; k<nconv; k++){
                //fd[k] = 0.5 * (max_eig) * (Ur[k] - Ul[k]);
                fd[k] = 0.5 * (max_eig) * (iConVarDiff[k][i][j]);
                diss[k][i][j] = diss[k][i][j] - fd[k] * ds;
                //rhs[k][i][j] = -diss[k][i][j] - iAvgFlux[k][i][j];
                //rhs[k][i][j] = rhs[k][i][j] + iAvgFlux[k][i][j];
                //rhs[k][i][j] = iAvgFlux[k][i][j] - fd[k]*ds;
                diss[k][im1][j] = diss[k][im1][j] + fd[k] * ds;
                //rhs[k][im1][j] = -diss[k][im1][j] + iAvgFlux[k][i][j];
                //rhs[k][im1][j] = - iAvgFlux[k][i][j] + fd[k]*ds;
                //rhs[k][im1][j] = rhs[k][im1][j] - iAvgFlux[k][i][j];
            }
        }
    }

    //compute in j-direction
    for (int i = 2; i <= ib; i++) {
        for (int j = 2; j <= jd1; j++) {

            jm1 = j - 1;

            ds = sqrt(sj[0][i][j] * sj[0][i][j] + sj[1][i][j] * sj[1][i][j]);
            nx = sj[0][i][j] / ds;
            ny = sj[1][i][j] / ds;

            //left and right states

            rhol = dv[0][i][jm1];
            ul = dv[1][i][jm1];
            vl = dv[2][i][jm1];
            pl = dv[3][i][jm1];
            al = sqrt(Gamma * pl / rhol);


            rhor = dv[0][i][j];
            ur = dv[1][i][j];
            vr = dv[2][i][j];
            pr = dv[3][i][j];
            ar = sqrt(Gamma * pr / rhor);

            max_eig = MaxEigVal(ur, ul, vr, vl, ar, al, nx, ny);

            for(int k=0; k<nconv; k++){
                //fd[k] = 0.5 * (max_eig) * (Ur[k] - Ul[k]);
                fd[k] = 0.5 * (max_eig) * (jConVarDiff[k][i][j]);
                diss[k][i][j] = diss[k][i][j] - fd[k] * ds;
                //rhs[k][i][j] = rhs[k][i][j] + jAvgFlux[k][i][j];
                //rhs[k][i][j] = -diss[k][i][j] - jAvgFlux[k][i][j];
                //rhs[k][i][j] = jAvgFlux[k][i][j] - fd[k]*ds;
                diss[k][i][jm1] = diss[k][i][jm1] + fd[k] * ds;
                //rhs[k][i][jm1] = -diss[k][i][jm1] + jAvgFlux[k][i][j];
                //rhs[k][i][jm1] = - jAvgFlux[k][i][j] + fd[k]*ds;
                //rhs[k][i][j - 1] = rhs[k][i][j - 1] - jAvgFlux[k][i][j];
            }
        }
    }

    /*for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            std::cout<<"the computed dissipation values="<<diss[3][1][j]<<std::endl;
        }
    }*/

    delete[] fd;
    delete[] fr;
    delete[] fl;
    delete[] Ur;
    delete[] Ul;
}


template<typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny) {
    double L1r, L2r, L3r, L1l, L2l, L3l;
    L1r = fabs(ur * nx + vr * ny) + ar;
    L1l = fabs(ul * nx + vl * ny) + al;
    L2r = fabs(ur * nx + vr * ny);
    L2l = fabs(ul * nx + vl * ny);
    L3r = fabs(ur * nx + vr * ny) - ar;
    L3l = fabs(ul * nx + vl * ny) - al;

    //return Max3(max(L1l,L1r), max(L2l,L2r), max(L3l,L3r));
    return Max2(Max3(L1r, L2r, L3r), Max3(L1l, L2l, L3l));
}



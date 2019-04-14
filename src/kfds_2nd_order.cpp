//
// Created by Maruthi NH on 14-06-2018.
//

#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

template <typename T>
T Movers(T deltaF, T deltaU, T L_max, T L_min);

template<typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny);

template<typename T>
T MinEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny);

void KFDS_2ndOrder(int id2, int jd2, int id1, int jd1, int ib, int jb, double ***&cv, double ***&dv,
                   double ***&iAvgFlux, double ***&jAvgFlux, double ***&iFluxDiff, double ***&jFluxDiff,
                   double ***&iConVarDiff, double ***&jConVarDiff){


    int im1, jm1;
    double nx, ny, ds, rhol, ul, vl, rhor, ur, vr, pl, pr, gam1, ggm1;
    double ar, al, max_eig, min_eig;
    double *fd, *fr, *fl, *Ur, *Ul;

    fd = new double[nconv];
    fr = new double[nconv];
    fl = new double[nconv];
    Ur = new double[nconv];
    Ul = new double[nconv];

    gam1 = Gamma - 1.0;
    ggm1 = Gamma / gam1;

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
            min_eig = MinEigVal(ur, ul, vr, vl, ar, al, nx, ny);

            //diffusive flux
            for(int k=0;k<nconv; k++){
                fd[k] = 0.5 * Movers(iFluxDiff[k][i][j], iConVarDiff[k][i][j], max_eig, min_eig) * (iConVarDiff[k][i][j]);
                //fd[k] = 0.5 * max_eig * (iConVarDiff[k][i][j]);
                diss[k][i][j] = diss[k][i][j] - fd[k] * ds;
                diss[k][im1][j] = diss[k][im1][j] + fd[k] * ds;

                //rhs[k][i][j] = -diss[k][i][j] + iAvgFlux[k][i][j];
                /*rhs[1][i][j] = rhs[1][i][j] + fc[1];
                rhs[2][i][j] = rhs[2][i][j] + fc[2];
                rhs[3][i][j] = rhs[3][i][j] + fc[3];*/

                //rhs[k][im1][j] = rhs[k][im1][j] - iAvgFlux[k][i][j];
                /*rhs[1][i - 1][j] = rhs[1][i - 1][j] - fc[1];
                rhs[2][i - 1][j] = rhs[2][i - 1][j] - fc[2];
                rhs[3][i - 1][j] = rhs[3][i - 1][j] - fc[3];*/
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
            min_eig = MinEigVal(ur, ul, vr, vl, ar, al, nx, ny);

            //diffusive flux
            for(int k=0;k<nconv; k++){
                fd[k] = 0.5 * Movers(jFluxDiff[k][i][j], jConVarDiff[k][i][j], max_eig, min_eig) * (jConVarDiff[k][i][j]);
                //fd[k] = 0.5 * max_eig * (jConVarDiff[k][i][j]);
                diss[k][i][j] = diss[k][i][j] - fd[k] * ds;
                diss[k][i][jm1] = diss[k][i][jm1] + fd[k] * ds;

                //rhs[k][i][j] = -diss[k][i][j] + jAvgFlux[k][i][j];
                /*rhs[1][i][j] = rhs[1][i][j] + fc[1];
                rhs[2][i][j] = rhs[2][i][j] + fc[2];
                rhs[3][i][j] = rhs[3][i][j] + fc[3];*/

                //rhs[k][i][jm1] = rhs[k][i][jm1] - jAvgFlux[k][i][j];
                /*rhs[1][i][j - 1] = rhs[1][i][j - 1] - fc[1];
                rhs[2][i][j - 1] = rhs[2][i][j - 1] - fc[2];
                rhs[3][i][j - 1] = rhs[3][i][j - 1] - fc[3];*/
            }
        }
    }

    /*for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            std::cout<<"the computed dissipation values="<<diss[2][i][j]<<std::endl;
        }
    }*/

    delete[] fd;
    delete[] fr;
    delete[] fl;
    delete[] Ur;
    delete[] Ul;

}

template <typename T>
T Movers(T deltaF, T deltaU, T L_max, T L_min){

    double S;
    const double epsilon=1e-6;

    if (fabs(deltaF)<epsilon) return 0.0;
    else if (fabs(deltaU)<epsilon) return L_min;
    else if (fabs(deltaU)>epsilon and fabs(deltaF)>epsilon) S=fabs(deltaF/deltaU);
    else S=L_min;

    if (S<epsilon)	return 0;
    else if (S>=L_max) return L_max;
    else if (S<=L_min) return L_min;
    else return S;
}

//template <typename T>
//T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){
//
//	double S;
//
//	if (Ur!=Ul) S=fabs(((Fr-Fl)/(Ur-Ul)));
//	else S=L_min;
//
//	if (fabs(S)>=L_max) return L_max;
//	else if(fabs(S)<=L_min) return L_min;
//        else return S;
//}
/*
template<typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny) {
    double L1r, L2r, L3r, L1l, L2l, L3l;
    L1r = fabs(ur * nx + vr * ny);
    L1l = fabs(ul * nx + vl * ny);
    L2r=fabs(ur*nx+vr*ny)+ar; 	   L2l=fabs(ul*nx+vl*ny)+al;
    L3r=fabs(ur*nx+vr*ny)-ar;  L3l=fabs(ul*nx+vl*ny)-al;

    return Max3(std::max(L1l,L1r), std::max(L2l,L2r), std::max(L3l,L3r));
    //return Max2(L1r, L1l) + Max2(ar, al);
}*/


template<typename T>
T MinEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny) {
    double L1r, L2r, L3r, L1l, L2l, L3l;
    L1r = fabs(ur * nx + vr * ny);
    L1l = fabs(ul * nx + vl * ny);
    L2r=fabs(ur*nx+vr*ny)+ar; 	   L2l=fabs(ul*nx+vl*ny)+al;
    L3r=fabs(ur*nx+vr*ny)-ar;  L3l=fabs(ul*nx+vl*ny)-al;

    return Max3(std::min(L1l,L1r), std::min(L2l,L2r), std::min(L3l,L3r));
    //return Min2(L1r, L1l);
}

template<typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny) {
    double L1r, L2r, L3r, L1l, L2l, L3l;
    L1r = fabs(ur * nx + vr * ny + ar);
    L1l = fabs(ul * nx + vl * ny + al);
    L2r = fabs(ur * nx + vr * ny);
    L2l = fabs(ul * nx + vl * ny);
    L3r = fabs(ur * nx + vr * ny - ar);
    L3l = fabs(ul * nx + vl * ny - al);

    //return Max3(max(L1l,L1r), max(L2l,L2r), max(L3l,L3r));
    return fabs(Max2(Max3(L1r, L2r, L3r), Max3(L1l, L2l, L3l)));
}
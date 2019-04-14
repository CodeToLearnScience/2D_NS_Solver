#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

template<typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny);

template<typename T>
T Minmod_Flux(T del_p, T del_m);

template<typename T>
T VanAlbada_Limiter(T del_plus, T del_minus);

template<typename T>
T Venki_Limiter(T del_plus, T del_minus, T dx);


void Diss_LLF2_Prim(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                    double ***&sj, double ***&diss) {

    int im1, jm1, ip1, jp1;
    double nx, ny, ds, rhol, ul, vl, rhor, ur, vr, pl, pr, hl, hr, El, Er, gam1, ggm1;
    double ar, al, max_eig;
    double *fd, *fr, *fl, *Ur, *Ul, *deltl, *deltr;
    //double area_avg;

    fd = new double[nconv];
    fr = new double[nconv];
    fl = new double[nconv];
    Ur = new double[nconv];
    Ul = new double[nconv];
    deltl = new double[nconv];
    deltr = new double[nconv];

    gam1 = Gamma - 1.0;
    ggm1 = Gamma / gam1;

    //compute in i-direction
    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= id1; i++) {

            im1 = i - 1;
            ip1 = i + 1;

            ds = sqrt(si[0][i][j] * si[0][i][j] + si[1][i][j] * si[1][i][j]);
            nx = si[0][i][j] / ds;
            ny = si[1][i][j] / ds;


            //area_avg = 0.5*(area[im1][j]+area[i][j]);
            //left and right states
            for(int k=0; k<4; k++){
                deltl[k] = 0.5*dui[k][im1][j]*Minmod_Flux(dui[k][i][j], dui[k][im1][j]);
                deltr[k] = 0.5*dui[k][i][j]*Minmod_Flux(dui[k][ip1][j], dui[k][i][j]);

                /*deltl[k] = 0.5*dui[k][im1][j]*VanAlbada_Limiter(dui[k][i][j], dui[k][im1][j]);
                deltr[k] = 0.5*dui[k][i][j]*VanAlbada_Limiter(dui[k][ip1][j], dui[k][i][j]);*/

                /*deltl[k] = 0.5 * Venki_Limiter(dui[k][i][j], dui[k][im1][j], dx);
                deltr[k] = 0.5 * Venki_Limiter(dui[k][ip1][j], dui[k][i][j], dx);*/
            }

            rhol = dv[0][im1][j] + deltl[0];
            ul = dv[1][im1][j] + deltl[1];
            vl = dv[2][im1][j] + deltl[2];
            pl = dv[3][im1][j] + deltl[3];
            hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            El = pl / (gam1) + 0.5 * rhol * (ul * ul + vl * vl);
            al = sqrt(Gamma * pl / rhol);

            rhor = dv[0][i][j] - deltr[0];
            ur = dv[1][i][j]  - deltr[1];
            vr = dv[2][i][j]  - deltr[2];
            pr = dv[3][i][j] - deltr[3];
            hr = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            Er = pr / (gam1) + 0.5 * rhor * (ur * ur + vr * vr);
            ar = sqrt(Gamma * pr / rhor);

            max_eig = MaxEigVal(ur, ul, vr, vl, ar, al, nx, ny);

            fl[0] = (rhol * ul) * nx + (rhol * vl) * ny;
            fl[1] = (rhol * ul * ul + pl) * nx + (rhol * ul * vl) * ny;
            fl[2] = (rhol * ul * vl) * nx + (rhol * vl * vl + pl) * ny;
            fl[3] = (rhol * ul * hl) * nx + (rhol * vl * hl) * ny;

            fr[0] = (rhor * ur) * nx + (rhor * vr) * ny;
            fr[1] = (rhor * ur * ur + pr) * nx + (rhor * ur * vr) * ny;
            fr[2] = (rhor * ur * vr) * nx + (rhor * vr * vr + pr) * ny;
            fr[3] = (rhor * ur * hr) * nx + (rhor * vr * hr) * ny;

            Ur[0] = rhor;
            Ur[1] = rhor * ur;
            Ur[2] = rhor * vr;
            Ur[3] = Er;

            Ul[0] = rhol;
            Ul[1] = rhol * ul;
            Ul[2] = rhol * vl;
            Ul[3] = El;

            for(int k=0; k<4; k++){
                fd[k] = 0.5 * (max_eig) * (Ur[k] - Ul[k]);
                //fd[k] = 0.5 * (max_eig) * (iConVarDiff[k][i][j]);
                //final dissipation terms
                diss[k][i][j] = diss[k][i][j] - fd[k] * ds;
                diss[k][im1][j] = diss[k][im1][j] + fd[k] * ds;
            }



        }
    }

    //compute in j-direction
    for (int i = 2; i <= ib; i++) {
        for (int j = 2; j <= jd1; j++) {

            jm1 = j - 1;
            jp1 = j + 1;

            ds = sqrt(sj[0][i][j] * sj[0][i][j] + sj[1][i][j] * sj[1][i][j]);
            nx = sj[0][i][j] / ds;
            ny = sj[1][i][j] / ds;


            // area_avg = 0.5*(area[i][jm1]+area[i][j]);

            //left and right states
            for(int k=0; k<4; k++){
                deltl[k] = 0.5*duj[k][i][jm1]*Minmod_Flux(duj[k][i][j], duj[k][i][jm1]);
                deltr[k] = 0.5*duj[k][i][j]*Minmod_Flux(duj[k][i][jp1], duj[k][i][j]);

                /*deltl[k] = 0.5*duj[k][i][jm1]*VanAlbada_Limiter(duj[k][i][j], duj[k][i][jm1]);
                deltr[k] = 0.5*duj[k][i][j]*VanAlbada_Limiter(duj[k][i][jp1], duj[k][i][j]);*/

                /*deltl[k] = 0.5 * Venki_Limiter(duj[k][i][j], duj[k][i][jm1], dy);
                deltr[k] = 0.5 * Venki_Limiter(duj[k][i][jp1], duj[k][i][j], dy);*/
            }

            rhol = dv[0][i][jm1] + deltl[0];
            ul = dv[1][i][jm1]  + deltl[1];
            vl = dv[2][i][jm1]  + deltl[2];
            pl = dv[3][i][jm1] + deltl[3];
            hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            El = pl / (gam1) + 0.5 * rhol * (ul * ul + vl * vl);
            al = sqrt(Gamma * pl / rhol);


            rhor = dv[0][i][j] - deltr[0];
            ur = dv[1][i][j]  - deltr[1];
            vr = dv[2][i][j]  - deltr[2];
            pr = dv[3][i][j] - deltr[3];
            hr = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            Er = pr / (gam1) + 0.5 * rhor * (ur * ur + vr * vr);
            ar = sqrt(Gamma * pr / rhor);

            max_eig = MaxEigVal(ur, ul, vr, vl, ar, al, nx, ny);

            fl[0] = (rhol * ul) * nx + (rhol * vl) * ny;
            fl[1] = (rhol * ul * ul + pl) * nx + (rhol * ul * vl) * ny;
            fl[2] = (rhol * ul * vl) * nx + (rhol * vl * vl + pl) * ny;
            fl[3] = (rhol * ul * hl) * nx + (rhol * vl * hl) * ny;

            fr[0] = (rhor * ur) * nx + (rhor * vr) * ny;
            fr[1] = (rhor * ur * ur + pr) * nx + (rhor * ur * vr) * ny;
            fr[2] = (rhor * ur * vr) * nx + (rhor * vr * vr + pr) * ny;
            fr[3] = (rhor * ur * hr) * nx + (rhor * vr * hr) * ny;

            Ur[0] = rhor;
            Ur[1] = rhor * ur;
            Ur[2] = rhor * vr;
            Ur[3] = Er;

            Ul[0] = rhol;
            Ul[1] = rhol * ul;
            Ul[2] = rhol * vl;
            Ul[3] = El;

            for(int k=0; k<4; k++){
                fd[k] = 0.5 * (max_eig) * (Ur[k] - Ul[k]);
                //fd[k] = 0.5 * (max_eig) * (jConVarDiff[k][i][j]);
                //final dissipation terms
                diss[k][i][j]   = diss[k][i][j] - fd[k] * ds;
                diss[k][i][jm1] = diss[k][i][jm1] + fd[k] * ds;
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
    delete[] deltl;
    delete[] deltr;
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

/*
template<typename T>
T Minmod_Flux(T del_p, T del_m) {

    double r, eps = 1e-6;

    r = 0.0;

    r = del_p / (del_m + 1e-16);

    //if (fabs(del_p)<eps) {
        
      //  r = Sign(del_p)*del_m/eps;

   // }

    return std::max(0.0, std::min(1.0, r));

}

template<typename T>
T VanAlbada_Limiter(T del_plus, T del_minus) {

    const double omega = 1e-12;
    double r;

    //del_minus = Sign(del_minus)*(fabs(del_minus)+omega);

    r = del_plus / (del_minus + 1e-15);

    return ((r * r + r) / (r * r + 1.0));
}

template<typename T>
T Venki_Limiter(T del_plus, T del_minus, T dx) {

    double const k = 0.0;
    double eps;

    eps = pow(k * dx, 3);

    return ((del_plus * del_plus + eps) * del_minus + (del_minus * del_minus + eps) * del_plus) /
           (del_plus * del_plus + del_minus * del_minus + 2.0 * eps + 1e-10);
}*/

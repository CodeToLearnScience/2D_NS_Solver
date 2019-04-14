#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"


void Avg_Conv_Flux2(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double **&area,
                    double ***&si, double ***&sj, double ***&diss, double ***&rhs) {

    double *fc, *deltl, *deltr;
    int im1, jm1, ip1, jp1;
    double rhoinv, rhol, ul, vl, rhor, ur, vr, pl, pr, hl, hr, qslr, qsll, gam1, ggm1, pavg, nx, ny, ds, dx, dy;

    fc = new double[nconv];
    deltl = new double[nconv];
    deltr = new double[nconv];

    gam1 = Gamma - 1.0;
    ggm1 = Gamma / gam1;

    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= id1; i++) {

            rhs[0][i][j] = -diss[0][i][j];
            rhs[1][i][j] = -diss[1][i][j];
            rhs[2][i][j] = -diss[2][i][j];
            rhs[3][i][j] = -diss[3][i][j];
        }
    }

    for (int j = 2; j <= jb; j++) {
        //flux in i direction
        for (int i = 2; i <= id1; i++) {

            im1 = i - 1;
            ip1 = i + 1;

            ds = sqrt(si[0][i][j] * si[0][i][j] + si[1][i][j] * si[1][i][j]);
            nx = si[0][i][j] / ds;
            ny = si[1][i][j] / ds;
            dx = 0.5 * (x[i][j] - x[i - 1][j]);

            deltl[0] = 0.5 * dui[0][im1][j] * Minmod_Flux(dui[0][i][j], dui[0][im1][j]);
            deltr[0] = 0.5 * dui[0][i][j] * Minmod_Flux(dui[0][ip1][j], dui[0][i][j]);
            deltl[1] = 0.5 * dui[1][im1][j] * Minmod_Flux(dui[1][i][j], dui[1][im1][j]);
            deltr[1] = 0.5 * dui[1][i][j] * Minmod_Flux(dui[1][ip1][j], dui[1][i][j]);
            deltl[2] = 0.5 * dui[2][im1][j] * Minmod_Flux(dui[2][i][j], dui[2][im1][j]);
            deltr[2] = 0.5 * dui[2][i][j] * Minmod_Flux(dui[2][ip1][j], dui[2][i][j]);
            deltl[3] = 0.5 * dui[3][im1][j] * Minmod_Flux(dui[3][i][j], dui[3][im1][j]);
            deltr[3] = 0.5 * dui[3][i][j] * Minmod_Flux(dui[3][ip1][j], dui[3][i][j]);

            /*deltl[0] = 0.5*dui[0][im1][j]*VanAlbada_Limiter(dui[0][i][j], dui[0][im1][j]);
            deltr[0] = 0.5*dui[0][i][j]*VanAlbada_Limiter(dui[0][ip1][j], dui[0][i][j]);
            deltl[1] = 0.5*dui[1][im1][j]*VanAlbada_Limiter(dui[1][i][j], dui[1][im1][j]);
            deltr[1] = 0.5*dui[1][i][j]*VanAlbada_Limiter(dui[1][ip1][j], dui[1][i][j]);
            deltl[2] = 0.5*dui[2][im1][j]*VanAlbada_Limiter(dui[2][i][j], dui[2][im1][j]);
            deltr[2] = 0.5*dui[2][i][j]*VanAlbada_Limiter(dui[2][ip1][j], dui[2][i][j]);
            deltl[3] = 0.5*dui[3][im1][j]*VanAlbada_Limiter(dui[3][i][j], dui[3][im1][j]);
            deltr[3] = 0.5*dui[3][i][j]*VanAlbada_Limiter(dui[3][ip1][j], dui[3][i][j]);*/

            /*deltl[0] = 0.5 * Venki_Limiter(dui[0][i][j], dui[0][im1][j], dx);
            deltr[0] = 0.5 * Venki_Limiter(dui[0][ip1][j], dui[0][i][j], dx);
            deltl[1] = 0.5 * Venki_Limiter(dui[1][i][j], dui[1][im1][j], dx);
            deltr[1] = 0.5 * Venki_Limiter(dui[1][ip1][j], dui[1][i][j], dx);
            deltl[2] = 0.5 * Venki_Limiter(dui[2][i][j], dui[2][im1][j], dx);
            deltr[2] = 0.5 * Venki_Limiter(dui[2][ip1][j], dui[2][i][j], dx);
            deltl[3] = 0.5 * Venki_Limiter(dui[3][i][j], dui[3][im1][j], dx);
            deltr[3] = 0.5 * Venki_Limiter(dui[3][ip1][j], dui[3][i][j], dx);*/

            rhol = dv[0][im1][j] + deltl[0];
            ul   = dv[1][im1][j] + deltl[1];
            vl   = dv[2][im1][j] + deltl[2];
            pl   = dv[3][im1][j] + deltl[3];
            hl   = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            qsll = (ul * si[0][i][j] + vl * si[1][i][j]) * rhol;

            rhor = dv[0][i][j] - deltr[0];
            ur   = dv[1][i][j] - deltr[1];
            vr   = dv[2][i][j] - deltr[2];
            pr   = dv[3][i][j] - deltr[3];
            hr   = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            qslr = (ur * si[0][i][j] + vr * si[1][i][j]) * rhor;

            pavg  = 0.5 * (pl + pr);
            fc[0] = 0.5 * (qsll + qslr);
            fc[1] = 0.5 * (qsll * ul + qslr * ur) + pavg * si[0][i][j];
            fc[2] = 0.5 * (qsll * vl + qslr * vr) + pavg * si[1][i][j];
            fc[3] = 0.5 * (qsll * hl + qslr * hr);

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
    for (int j = 2; j <= jd1; j++) {
        //flux in j-direction except at boundaries

        for (int i = 2; i <= ib; i++) {

            jm1 = j - 1;
            jp1 = j + 1;

            ds = sqrt(sj[0][i][j] * sj[0][i][j] + sj[1][i][j] * sj[1][i][j]);
            nx = sj[0][i][j] / ds;
            ny = sj[1][i][j] / ds;
            dy = 0.5 * (y[i][j] - y[i][j - 1]);

            deltl[0] = 0.5 * duj[0][i][jm1] * Minmod_Flux(duj[0][i][j], duj[0][i][jm1]);
            deltr[0] = 0.5 * duj[0][i][j] * Minmod_Flux(duj[0][i][jp1], duj[0][i][j]);
            deltl[1] = 0.5 * duj[1][i][jm1] * Minmod_Flux(duj[1][i][j], duj[1][i][jm1]);
            deltr[1] = 0.5 * duj[1][i][j] * Minmod_Flux(duj[1][i][jp1], duj[1][i][j]);
            deltl[2] = 0.5 * duj[2][i][jm1] * Minmod_Flux(duj[2][i][j], duj[2][i][jm1]);
            deltr[2] = 0.5 * duj[2][i][j] * Minmod_Flux(duj[2][i][jp1], duj[2][i][j]);
            deltl[3] = 0.5 * duj[3][i][jm1] * Minmod_Flux(duj[3][i][j], duj[3][i][jm1]);
            deltr[3] = 0.5 * duj[3][i][j] * Minmod_Flux(duj[3][i][jp1], duj[3][i][j]);

            /*deltl[0] = 0.5*duj[0][i][jm1]*VanAlbada_Limiter(duj[0][i][j], duj[0][i][jm1]);
            deltr[0] = 0.5*duj[0][i][j]*VanAlbada_Limiter(duj[0][i][jp1], duj[0][i][j]);
            deltl[1] = 0.5*duj[1][i][jm1]*VanAlbada_Limiter(duj[1][i][j], duj[1][i][jm1]);
            deltr[1] = 0.5*duj[1][i][j]*VanAlbada_Limiter(duj[1][i][jp1], duj[1][i][j]);
            deltl[2] = 0.5*duj[2][i][jm1]*VanAlbada_Limiter(duj[2][i][j], duj[2][i][jm1]);
            deltr[2] = 0.5*duj[2][i][j]*VanAlbada_Limiter(duj[2][i][jp1], duj[2][i][j]);
            deltl[3] = 0.5*duj[3][i][jm1]*VanAlbada_Limiter(duj[3][i][j], duj[3][i][jm1]);
            deltr[3] = 0.5*duj[3][i][j]*VanAlbada_Limiter(duj[3][i][jp1], duj[3][i][j]);*/

            /*deltl[0] = 0.5 * Venki_Limiter(duj[0][i][j], duj[0][i][jm1], dy);
            deltr[0] = 0.5 * Venki_Limiter(duj[0][i][jp1], duj[0][i][j], dy);
            deltl[1] = 0.5 * Venki_Limiter(duj[1][i][j], duj[1][i][jm1], dy);
            deltr[1] = 0.5 * Venki_Limiter(duj[1][i][jp1], duj[1][i][j], dy);
            deltl[2] = 0.5 * Venki_Limiter(duj[2][i][j], duj[2][i][jm1], dy);
            deltr[2] = 0.5 * Venki_Limiter(duj[2][i][jp1], duj[2][i][j], dy);
            deltl[3] = 0.5 * Venki_Limiter(duj[3][i][j], duj[3][i][jm1], dy);
            deltr[3] = 0.5 * Venki_Limiter(duj[3][i][jp1], duj[3][i][j], dy);*/


            rhol = dv[0][i][jm1] + deltl[0];
            ul   = dv[1][i][jm1] + deltl[1];
            vl   = dv[2][i][jm1] + deltl[2];
            pl   = dv[3][i][jm1] + deltl[3];
            hl   = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            qsll = (ul * sj[0][i][j] + vl * sj[1][i][j]) * rhol;

            rhor = dv[0][i][j] - deltr[0];
            ur   = dv[1][i][j] - deltr[1];
            vr   = dv[2][i][j] - deltr[2];
            pr   = dv[3][i][j] - deltr[3];
            hr   = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            qslr = (ur * sj[0][i][j] + vr * sj[1][i][j]) * rhor;

            pavg =  0.5 * (pl + pr);
            fc[0] = 0.5 * (qsll + qslr);
            fc[1] = 0.5 * (qsll * ul + qslr * ur) + pavg * sj[0][i][j];
            fc[2] = 0.5 * (qsll * vl + qslr * vr) + pavg * sj[1][i][j];
            fc[3] = 0.5 * (qsll * hl + qslr * hr);

            rhs[0][i][j] = rhs[0][i][j] + fc[0];
            rhs[1][i][j] = rhs[1][i][j] + fc[1];
            rhs[2][i][j] = rhs[2][i][j] + fc[2];
            rhs[3][i][j] = rhs[3][i][j] + fc[3];

            rhs[0][i][j - 1] = rhs[0][i][j - 1] - fc[0];
            rhs[1][i][j - 1] = rhs[1][i][j - 1] - fc[1];
            rhs[2][i][j - 1] = rhs[2][i][j - 1] - fc[2];
            rhs[3][i][j - 1] = rhs[3][i][j - 1] - fc[3];
        }
    }

    /*Flux_Symmetry(1, 1, 2, 9, cv, dv, si, sj, rhs);
    Flux_Boundary(2, id1, 2, jb, cv, dv, si, sj, rhs);
    Flux_Boundary(3, jd1, 2, ib, cv, dv, si, sj, rhs);
    Flux_Boundary(4, 1, 2, jb, cv, dv, si, sj, rhs);
    Flux_ViscWall(1, 1, 10, ib, cv, dv, si, sj, rhs);*/

    /*for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            std::cout<<"the computed rhs="<<rhs[1][i][j]<<std::endl;
        }
    }*/

    //at boundaries

    delete[] fc;
    delete[] deltl;
    delete[] deltr;

}

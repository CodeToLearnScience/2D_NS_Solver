#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

template <typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny);

template <typename T>
T Minmod_Flux(T del_p, T del_m);

template <typename T>
T VanAlbada_Limiter(T del_plus, T del_minus);

template <typename T>
T Venki_Limiter(T del_plus, T del_minus, T dx);

void Diss_Roe2_Prim(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                    double ***&sj, double ***&diss)
{

    int im1, jm1, ip1, jp1;
    double nx, ny, dx, dy, ds, rhoinv, rhol, ul, vl, rhor, ur, vr, pl, pr, hl, hr, gam1, ggm1;
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
    for (int j = 2; j <= jb; j++)
    {
        for (int i = 2; i <= id1; i++)
        {
            im1 = i - 1;
            ip1 = i + 1;

            ds = sqrt(si[0][i][j] * si[0][i][j] + si[1][i][j] * si[1][i][j]);
            nx = si[0][i][j] / ds;
            ny = si[1][i][j] / ds;
            dx = 0.5 * (x[i + 1][j] - x[i - 1][j]);

            //area_avg = 0.5*(area[im1][j]+area[i][j]);
            //left and right states

            for (int k = 0; k < 4; k++)
            {
                // deltl[k] = 0.5 * dui[k][im1][j] * Minmod_Flux(dui[k][i][j], dui[k][im1][j]);
                // deltr[k] = 0.5 * dui[k][i][j] * Minmod_Flux(dui[k][ip1][j], dui[k][i][j]);

                deltl[k] = 0.5*dui[k][im1][j]*VanAlbada_Limiter(dui[k][i][j], dui[k][im1][j]);
                deltr[k] = 0.5*dui[k][i][j]*VanAlbada_Limiter(dui[k][ip1][j], dui[k][i][j]);

                /*deltl[k] = 0.5 * Venki_Limiter(dui[k][i][j], dui[k][im1][j], dx);
                deltr[k] = 0.5 * Venki_Limiter(dui[k][ip1][j], dui[k][i][j], dx);*/
            }

            rhoinv = 1.0 / dv[0][im1][j];
            rhol = dv[0][im1][j] + deltl[0];
            ul = dv[1][im1][j] + deltl[1];
            vl = dv[2][im1][j] + deltl[2];
            pl = dv[3][im1][j] + deltl[3];
            hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            al = sqrt(Gamma * pl / rhol);

            rhoinv = 1.0 / dv[0][i][j];
            rhor = dv[0][i][j] - deltr[0];
            ur = dv[1][i][j] - deltr[1];
            vr = dv[2][i][j] - deltr[2];
            pr = dv[3][i][j] - deltr[3];
            hr = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
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

            double w, rho_roe, u_roe, v_roe, a_roe, Vt, Vn, lambda1, lambda2;
            double r1_M, r1_x, r1_y, r1_E, r2_M, r2_x, r2_y, r2_E;
            double alpha1, alpha2, dVn, dVt;

            //roe avgs
            w = sqrt(rhor / rhol);
            rho_roe = w * rhol;
            u_roe = (ul + w * ur) / (1.0 + w);
            v_roe = (vl + w * vr) / (1.0 + w);
            a_roe = sqrt((al * al * sqrt(rhol) + ar * ar * sqrt(rhor)) / (sqrt(rhol) + sqrt(rhor)));

            Vn = u_roe * nx + v_roe * ny;
            Vt = -u_roe * ny + v_roe * nx;

            lambda1 = fabs(Vn);
            lambda2 = fabs(sqrt((Gamma - 1) / Gamma) * a_roe);

            r1_M = 0.0;
            r1_x = nx;
            r1_y = ny;
            r1_E = Vn - a_roe / sqrt(Gamma * (Gamma - 1));

            r2_M = 1.0;
            r2_x = nx;
            r2_y = ny;
            r2_E = Vn + a_roe / sqrt(Gamma * (Gamma - 1));

            //Wave strengths
            dVn = (ur - ul) * nx + (vr - vl) * ny;
            dVt = (ur - ul) * ny - (vr - vl) * nx;

            alpha1 = ((rho_roe * dVn) / 2.0) - sqrt(Gamma / (Gamma - 1.0)) * ((pr - pl) / (2.0 * a_roe));
            alpha2 = ((rho_roe * dVn) / 2.0) + sqrt(Gamma / (Gamma - 1.0)) * ((pr - pl) / (2.0 * a_roe));

            double term1, term2, term3;

            fd[0] = 0.5 * (lambda1 * (rhor - rhol) + lambda2 * 0.0);
            fd[1] = 0.5 *
                    (lambda1 * (rho_roe * (ur - ul) + u_roe * (rhor - rhol)) + lambda2 * (alpha1 * nx + alpha2 * nx));
            fd[2] = 0.5 *
                    (lambda1 * (rho_roe * (vr - vl) + v_roe * (rhor - rhol)) + lambda2 * (alpha1 * ny + alpha2 * ny));
            term1 = (pr - pl) / (Gamma - 1.0);
            term2 = 0.5 * (u_roe * u_roe + v_roe * v_roe) * (rhor - rhol);
            term3 = rho_roe * (u_roe * (ur - ul) + v_roe * (vr - vl));
            fd[3] = 0.5 * (lambda1 * (term1 + term2 + term3) + lambda2 * (alpha1 * r1_E + alpha2 * r2_E));

            //final dissipation terms
            for (int k = 0; k < 4; k++)
            {
                diss[k][i][j] = diss[k][i][j] - fd[k] * ds;
                diss[k][im1][j] = diss[k][im1][j] + fd[k] * ds;
            }

            /*std::cout<<"diss flux:"<<i<<"\t"<<j<<"\t"<<diss[0][i][j]<<"\t"<<diss[1][i][j]<<"\t"<<diss[2][i][j]<<"\t"
                     <<diss[3][i][j]<<std::endl;*/
        }
    }

    //compute in j-direction
    for (int i = 2; i <= ib; i++)
    {
        for (int j = 2; j <= jd1; j++)
        {
            jm1 = j - 1;
            jp1 = j + 1;

            ds = sqrt(sj[0][i][j] * sj[0][i][j] + sj[1][i][j] * sj[1][i][j]);
            nx = sj[0][i][j] / ds;
            ny = sj[1][i][j] / ds;
            dy = 0.5 * (y[i][j + 1] - y[i][j - 1]);

            // area_avg = 0.5*(area[i][jm1]+area[i][j]);

            //left and right states

            for (int k = 0; k < 4; k++)
            {
                // deltl[k] = 0.5 * duj[k][i][jm1] * Minmod_Flux(duj[k][i][j], duj[k][i][jm1]);
                // deltr[k] = 0.5 * duj[k][i][j] * Minmod_Flux(duj[k][i][jp1], duj[k][i][j]);

                deltl[k] = 0.5*duj[k][i][jm1]*VanAlbada_Limiter(duj[k][i][j], duj[k][i][jm1]);
                deltr[k] = 0.5*duj[k][i][j]*VanAlbada_Limiter(duj[k][i][jp1], duj[k][i][j]);

                /*deltl[k] = 0.5 * Venki_Limiter(duj[k][i][j], duj[k][i][jm1], dy);
                deltr[k] = 0.5 * Venki_Limiter(duj[k][i][jp1], duj[k][i][j], dy);*/
            }

            rhol = dv[0][i][jm1] + deltl[0];
            ul = dv[1][i][jm1] + deltl[1];
            vl = dv[2][i][jm1] + deltl[2];
            pl = dv[3][i][jm1] + deltl[3];
            hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            al = sqrt(Gamma * pl / rhol);

            rhor = dv[0][i][j] - deltr[0];
            ur = dv[1][i][j] - deltr[1];
            vr = dv[2][i][j] - deltr[2];
            pr = dv[3][i][j] - deltr[3];
            hr = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
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

            double w, rho_roe, u_roe, v_roe, a_roe, Vt, Vn, lambda1, lambda2;
            double r1_M, r1_x, r1_y, r1_E, r2_M, r2_x, r2_y, r2_E;
            double alpha1, alpha2, dVn, dVt;

            //roe avgs
            w = sqrt(rhor / rhol);
            rho_roe = w * rhol;
            u_roe = (ul + w * ur) / (1.0 + w);
            v_roe = (vl + w * vr) / (1.0 + w);
            a_roe = sqrt((al * al * sqrt(rhol) + ar * ar * sqrt(rhor)) / (sqrt(rhol) + sqrt(rhor)));

            Vn = u_roe * nx + v_roe * ny;
            Vt = -u_roe * ny + v_roe * nx;

            lambda1 = fabs(Vn);
            lambda2 = fabs(sqrt((Gamma - 1) / Gamma) * a_roe);

            r1_M = 0.0;
            r1_x = nx;
            r1_y = ny;
            r1_E = Vn - a_roe / sqrt(Gamma * (Gamma - 1));

            r2_M = 1.0;
            r2_x = nx;
            r2_y = ny;
            r2_E = Vn + a_roe / sqrt(Gamma * (Gamma - 1));

            //Wave strengths
            dVn = (ur - ul) * nx + (vr - vl) * ny;
            dVt = (ur - ul) * ny - (vr - vl) * nx;

            alpha1 = ((rho_roe * dVn) / 2.0) - sqrt(Gamma / (Gamma - 1.0)) * ((pr - pl) / (2.0 * a_roe));
            alpha2 = ((rho_roe * dVn) / 2.0) + sqrt(Gamma / (Gamma - 1.0)) * ((pr - pl) / (2.0 * a_roe));

            double term1, term2, term3;

            fd[0] = 0.5 * (lambda1 * (rhor - rhol) + lambda2 * 0.0);
            fd[1] = 0.5 *
                    (lambda1 * (rho_roe * (ur - ul) + u_roe * (rhor - rhol)) + lambda2 * (alpha1 * nx + alpha2 * nx));
            fd[2] = 0.5 *
                    (lambda1 * (rho_roe * (vr - vl) + v_roe * (rhor - rhol)) + lambda2 * (alpha1 * ny + alpha2 * ny));
            term1 = (pr - pl) / (Gamma - 1.0);
            term2 = 0.5 * (u_roe * u_roe + v_roe * v_roe) * (rhor - rhol);
            term3 = rho_roe * (u_roe * (ur - ul) + v_roe * (vr - vl));
            fd[3] = 0.5 * (lambda1 * (term1 + term2 + term3) + lambda2 * (alpha1 * r1_E + alpha2 * r2_E));

            //final dissipation terms
            for (int k = 0; k < 4; k++)
            {
                diss[k][i][j] = diss[k][i][j] - fd[k] * ds;
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

/*********************Naveens_ROE_TV Scheme************************************/

void Diss_Roe2_TV_Prim(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                       double ***&sj, double ***&diss)
{

    int im1, jm1, ip1, jp1;
    double nx, ny, dx, dy, ds, rhoinv, rhol, ul, vl, rhor, ur, vr, pl, pr, hl, hr, El, Er, gam1, ggm1;
    double ar, al, max_eig, p_avg;
    double *fd, *fr, *fl, *Ur, *Ul, *deltl, *deltr;
    //double area_avg;

    double w, rho_roe, u_roe, v_roe, h_roe, a_roe, Vt, Vn, lambda1, lambda2, lambda3, lambda2_a, lambda3_a;
    double r1_M, r1_x, r1_y, r1_E, r2_M, r2_x, r2_y, r2_E;
    double alpha1, alpha2, dVn, dVt;
    double beta = 0.0;
    double term1, term2;

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
    for (int j = 2; j <= jb; j++)
    {
        for (int i = 2; i <= id1; i++)
        {
            im1 = i - 1;
            ip1 = i + 1;

            ds = sqrt(si[0][i][j] * si[0][i][j] + si[1][i][j] * si[1][i][j]);
            nx = si[0][i][j] / ds;
            ny = si[1][i][j] / ds;
            dx = 0.5 * (x[i + 1][j] - x[i - 1][j]);

            //area_avg = 0.5*(area[im1][j]+area[i][j]);
            //left and right states

            for (int k = 0; k < 4; k++)
            {
                deltl[k] = 0.5 * dui[k][im1][j] * Minmod_Flux(dui[k][i][j], dui[k][im1][j]);
                deltr[k] = 0.5 * dui[k][i][j] * Minmod_Flux(dui[k][ip1][j], dui[k][i][j]);

                /*deltl[k] = 0.5*dui[k][im1][j]*VanAlbada_Limiter(dui[k][i][j], dui[k][im1][j]);
                deltr[k] = 0.5*dui[k][i][j]*VanAlbada_Limiter(dui[k][ip1][j], dui[k][i][j]);*/

                /*deltl[k] = 0.5 * Venki_Limiter(dui[k][i][j], dui[k][im1][j], dx);
                deltr[k] = 0.5 * Venki_Limiter(dui[k][ip1][j], dui[k][i][j], dx);*/
            }

            rhoinv = 1.0 / dv[0][im1][j];
            rhol = dv[0][im1][j] + deltl[0];
            ul = dv[1][im1][j] + deltl[1];
            vl = dv[2][im1][j] + deltl[2];
            pl = dv[3][im1][j] + deltl[3];
            hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            El = pl / (gam1) + 0.5 * rhol * (ul * ul + vl * vl);
            al = sqrt(Gamma * pl / rhol);

            rhoinv = 1.0 / dv[0][i][j];
            rhor = dv[0][i][j] - deltr[0];
            ur = dv[1][i][j] - deltr[1];
            vr = dv[2][i][j] - deltr[2];
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

            //roe avgs
            w = sqrt(rhor / rhol);
            rho_roe = w * rhol;
            u_roe = (ul + w * ur) / (1.0 + w);
            v_roe = (vl + w * vr) / (1.0 + w);
            //a_roe = sqrt((al * al * sqrt(rhol) + ar * ar * sqrt(rhor)) / (sqrt(rhol) + sqrt(rhor)));
            p_avg = 0.5 * (pl + pr);
            a_roe = sqrt(Gamma * p_avg / rho_roe);

            Vn = u_roe * nx + v_roe * ny;
            Vt = -u_roe * ny + v_roe * nx;
            beta = sqrt(Vn * Vn + 4.0 * a_roe * a_roe);

            lambda1 = fabs(Vn);
            lambda2 = 0.5 * fabs(Vn - beta);
            lambda2_a = 0.5 * (Vn - beta);
            lambda3 = 0.5 * fabs(Vn + beta);
            lambda3_a = 0.5 * (Vn + beta);

            r1_M = 0.0;
            r1_x = nx;
            r1_y = ny;
            r1_E = Vn + (1 / 0.4) * lambda2_a;

            r2_M = 0.0;
            r2_x = nx;
            r2_y = ny;
            r2_E = Vn + (1 / 0.4) * lambda3_a;

            //Wave strengths
            dVn = (ur - ul) * nx + (vr - vl) * ny;
            dVt = (ur - ul) * ny - (vr - vl) * nx;

            alpha1 = ((rho_roe * dVn) / 2.0) - ((pr - pl) / beta) + (1.0 / (2 * beta)) * (Vn * rho_roe * dVn);
            alpha2 = ((rho_roe * dVn) / 2.0) + ((pr - pl) / beta) - (1.0 / (2 * beta)) * (Vn * rho_roe * dVn);

            fd[0] = 0.5 * (lambda1 * (rhor - rhol));
            fd[1] = 0.5 * (lambda1 * (rho_roe * (ur - ul) + u_roe * (rhor - rhol)) +
                           (alpha1 * lambda2 + alpha2 * lambda3) * nx);
            fd[2] = 0.5 * (lambda1 * (rho_roe * (vr - vl) + v_roe * (rhor - rhol)) +
                           (alpha1 * lambda2 + alpha2 * lambda3) * ny);
            term1 = 0.5 * (u_roe * u_roe + v_roe * v_roe) * (rhor - rhol);
            term2 = rho_roe * (u_roe * (ur - ul) + v_roe * (vr - vl));
            fd[3] = 0.5 * (lambda1 * (term1 + term2) + lambda2 * alpha1 * r1_E + alpha2 * lambda3 * r2_E);

            //final dissipation terms
            for (int k = 0; k < 4; k++)
            {
                diss[k][i][j] = diss[k][i][j] - fd[k] * ds;
                diss[k][im1][j] = diss[k][im1][j] + fd[k] * ds;
            }
        }
    }

    //compute in j-direction
    for (int i = 2; i <= ib; i++)
    {
        for (int j = 2; j <= jd1; j++)
        {
            jm1 = j - 1;
            jp1 = j + 1;

            ds = sqrt(sj[0][i][j] * sj[0][i][j] + sj[1][i][j] * sj[1][i][j]);
            nx = sj[0][i][j] / ds;
            ny = sj[1][i][j] / ds;
            dy = 0.5 * (y[i][j + 1] - y[i][j - 1]);

            // area_avg = 0.5*(area[i][jm1]+area[i][j]);

            //left and right states
            for (int k = 0; k < 4; k++)
            {
                deltl[k] = 0.5 * duj[k][i][jm1] * Minmod_Flux(duj[k][i][j], duj[k][i][jm1]);
                deltr[k] = 0.5 * duj[k][i][j] * Minmod_Flux(duj[k][i][jp1], duj[k][i][j]);

                /*deltl[k] = 0.5*duj[k][i][jm1]*VanAlbada_Limiter(duj[k][i][j], duj[k][i][jm1]);
                deltr[k] = 0.5*duj[k][i][j]*VanAlbada_Limiter(duj[k][i][jp1], duj[k][i][j]);*/

                /*deltl[k] = 0.5 * Venki_Limiter(duj[k][i][j], duj[k][i][jm1], dy);
                deltr[k] = 0.5 * Venki_Limiter(duj[k][i][jp1], duj[k][i][j], dy);*/
            }

            rhol = dv[0][i][jm1] + deltl[0];
            ul = dv[1][i][jm1] + deltl[1];
            vl = dv[2][i][jm1] + deltl[2];
            pl = dv[3][i][jm1] + deltl[3];
            hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            El = pl / (gam1) + 0.5 * rhol * (ul * ul + vl * vl);
            al = sqrt(Gamma * pl / rhol);

            rhor = dv[0][i][j] - deltr[0];
            ur = dv[1][i][j] - deltr[1];
            vr = dv[2][i][j] - deltr[2];
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

            //roe avgs
            w = sqrt(rhor / rhol);
            rho_roe = w * rhol;
            u_roe = (ul + w * ur) / (1.0 + w);
            v_roe = (vl + w * vr) / (1.0 + w);
            //a_roe = sqrt((al * al * sqrt(rhol) + ar * ar * sqrt(rhor)) / (sqrt(rhol) + sqrt(rhor)));
            p_avg = 0.5 * (pl + pr);
            a_roe = sqrt(Gamma * p_avg / rho_roe);

            Vn = u_roe * nx + v_roe * ny;
            Vt = -u_roe * ny + v_roe * nx;
            beta = sqrt(Vn * Vn + 4.0 * a_roe * a_roe);

            lambda1 = fabs(Vn);
            lambda2 = 0.5 * fabs(Vn - beta);
            lambda2_a = 0.5 * (Vn - beta);
            lambda3 = 0.5 * fabs(Vn + beta);
            lambda3_a = 0.5 * (Vn + beta);

            r1_M = 0.0;
            r1_x = nx;
            r1_y = ny;
            r1_E = Vn + (1 / 0.4) * lambda2_a;

            r2_M = 0.0;
            r2_x = nx;
            r2_y = ny;
            r2_E = Vn + (1 / 0.4) * lambda3_a;

            //Wave strengths
            dVn = (ur - ul) * nx + (vr - vl) * ny;
            dVt = (ur - ul) * ny - (vr - vl) * nx;

            alpha1 = ((rho_roe * dVn) / 2.0) - ((pr - pl) / beta) + (1.0 / (2 * beta)) * (Vn * rho_roe * dVn);
            alpha2 = ((rho_roe * dVn) / 2.0) + ((pr - pl) / beta) - (1.0 / (2 * beta)) * (Vn * rho_roe * dVn);

            fd[0] = 0.5 * (lambda1 * (rhor - rhol));
            fd[1] = 0.5 * (lambda1 * (rho_roe * (ur - ul) + u_roe * (rhor - rhol)) +
                           (alpha1 * lambda2 + alpha2 * lambda3) * nx);
            fd[2] = 0.5 * (lambda1 * (rho_roe * (vr - vl) + v_roe * (rhor - rhol)) +
                           (alpha1 * lambda2 + alpha2 * lambda3) * ny);
            term1 = 0.5 * (u_roe * u_roe + v_roe * v_roe) * (rhor - rhol);
            term2 = rho_roe * (u_roe * (ur - ul) + v_roe * (vr - vl));
            fd[3] = 0.5 * (lambda1 * (term1 + term2) + lambda2 * alpha1 * r1_E + alpha2 * lambda3 * r2_E);

            //final dissipation terms
            for (int k = 0; k < 4; k++)
            {
                diss[k][i][j] = diss[k][i][j] - fd[k] * ds;
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

/*******************************End of function for Naveen's roe scheme*********************************/

void Diss_ZB2(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
              double ***&sj, double ***&diss)
{

    int im1, jm1, ip1, jp1;
    double nx, ny, dx, dy, ds, rhoinv, rhol, ul, vl, rhor, ur, vr, pl, pr, hl, hr, El, Er, gam1, ggm1;
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
    for (int j = 2; j <= jb; j++)
    {
        for (int i = 2; i <= id1; i++)
        {
            im1 = i - 1;
            ip1 = i + 1;

            ds = sqrt(si[0][i][j] * si[0][i][j] + si[1][i][j] * si[1][i][j]);
            nx = si[0][i][j] / ds;
            ny = si[1][i][j] / ds;
            dx = 0.5 * (x[i + 1][j] - x[i - 1][j]);

            //area_avg = 0.5*(area[im1][j]+area[i][j]);
            //left and right states

            /*deltl[0] = 0.5*dui[0][im1][j]*Minmod_Flux(dui[0][i][j], dui[0][im1][j]);
            deltr[0] = 0.5*dui[0][i][j]*Minmod_Flux(dui[0][ip1][j], dui[0][i][j]);
            deltl[1] = 0.5*dui[1][im1][j]*Minmod_Flux(dui[1][i][j], dui[1][im1][j]);
            deltr[1] = 0.5*dui[1][i][j]*Minmod_Flux(dui[1][ip1][j], dui[1][i][j]);
            deltl[2] = 0.5*dui[2][im1][j]*Minmod_Flux(dui[2][i][j], dui[2][im1][j]);
            deltr[2] = 0.5*dui[2][i][j]*Minmod_Flux(dui[2][ip1][j], dui[2][i][j]);
            deltl[3] = 0.5*dui[3][im1][j]*Minmod_Flux(dui[3][i][j], dui[3][im1][j]);
            deltr[3] = 0.5*dui[3][i][j]*Minmod_Flux(dui[3][ip1][j], dui[3][i][j]);*/

            /*deltl[0] = 0.5*dui[0][im1][j]*VanAlbada_Limiter(dui[0][i][j], dui[0][im1][j]);
            deltr[0] = 0.5*dui[0][i][j]*VanAlbada_Limiter(dui[0][ip1][j], dui[0][i][j]);
            deltl[1] = 0.5*dui[1][im1][j]*VanAlbada_Limiter(dui[1][i][j], dui[1][im1][j]);
            deltr[1] = 0.5*dui[1][i][j]*VanAlbada_Limiter(dui[1][ip1][j], dui[1][i][j]);
            deltl[2] = 0.5*dui[2][im1][j]*VanAlbada_Limiter(dui[2][i][j], dui[2][im1][j]);
            deltr[2] = 0.5*dui[2][i][j]*VanAlbada_Limiter(dui[2][ip1][j], dui[2][i][j]);
            deltl[3] = 0.5*dui[3][im1][j]*VanAlbada_Limiter(dui[3][i][j], dui[3][im1][j]);
            deltr[3] = 0.5*dui[3][i][j]*VanAlbada_Limiter(dui[3][ip1][j], dui[3][i][j]);*/

            deltl[0] = 0.5 * Venki_Limiter(dui[0][i][j], dui[0][im1][j], dx);
            deltr[0] = 0.5 * Venki_Limiter(dui[0][ip1][j], dui[0][i][j], dx);
            deltl[1] = 0.5 * Venki_Limiter(dui[1][i][j], dui[1][im1][j], dx);
            deltr[1] = 0.5 * Venki_Limiter(dui[1][ip1][j], dui[1][i][j], dx);
            deltl[2] = 0.5 * Venki_Limiter(dui[2][i][j], dui[2][im1][j], dx);
            deltr[2] = 0.5 * Venki_Limiter(dui[2][ip1][j], dui[2][i][j], dx);
            deltl[3] = 0.5 * Venki_Limiter(dui[3][i][j], dui[3][im1][j], dx);
            deltr[3] = 0.5 * Venki_Limiter(dui[3][ip1][j], dui[3][i][j], dx);

            rhoinv = 1.0 / cv[0][im1][j];
            rhol = cv[0][im1][j] + deltl[0];
            ul = cv[1][im1][j] * rhoinv + deltl[1];
            vl = cv[2][im1][j] * rhoinv + deltl[2];
            pl = dv[0][im1][j] + deltl[3];
            hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            El = pl / (gam1) + 0.5 * rhol * (ul * ul + vl * vl);
            al = sqrt(Gamma * pl / rhol);

            rhoinv = 1.0 / cv[0][i][j];
            rhor = cv[0][i][j] - deltr[0];
            ur = cv[1][i][j] * rhoinv - deltr[1];
            vr = cv[2][i][j] * rhoinv - deltr[2];
            pr = dv[0][i][j] - deltr[3];
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

            double w, rho_roe, u_roe, v_roe, a_roe, Vt, Vn, lambda1, lambda2;
            double r1_M, r1_x, r1_y, r1_E, r2_M, r2_x, r2_y, r2_E;
            double alpha1, alpha2, dVn, dVt;

            //roe avgs
            w = sqrt(rhor / rhol);
            rho_roe = w * rhol;
            u_roe = (ul + w * ur) / (1.0 + w);
            v_roe = (vl + w * vr) / (1.0 + w);
            a_roe = sqrt((al * al * sqrt(rhol) + ar * ar * sqrt(rhor)) / (sqrt(rhol) + sqrt(rhor)));

            Vn = u_roe * nx + v_roe * ny;
            Vt = -u_roe * ny + v_roe * nx;

            lambda1 = fabs(Vn);
            lambda2 = fabs(sqrt((Gamma - 1) / Gamma) * a_roe);

            r1_M = 0.0;
            r1_x = nx;
            r1_y = ny;
            r1_E = Vn - a_roe / sqrt(Gamma * (Gamma - 1));

            r2_M = 1.0;
            r2_x = nx;
            r2_y = ny;
            r2_E = Vn + a_roe / sqrt(Gamma * (Gamma - 1));

            //Wave strengths
            dVn = (ur - ul) * nx + (vr - vl) * ny;
            dVt = (ur - ul) * ny - (vr - vl) * nx;

            alpha1 = ((rho_roe * dVn) / 2.0) - sqrt(Gamma / (Gamma - 1.0)) * ((pr - pl) / (2.0 * a_roe));
            alpha2 = ((rho_roe * dVn) / 2.0) + sqrt(Gamma / (Gamma - 1.0)) * ((pr - pl) / (2.0 * a_roe));

            double term1, term2, term3;

            fd[0] = 0.5 * (lambda1 * (rhor - rhol) + lambda2 * 0.0);
            fd[1] = 0.5 *
                    (lambda1 * (rho_roe * (ur - ul) + u_roe * (rhor - rhol)) + lambda2 * (alpha1 * nx + alpha2 * nx));
            fd[2] = 0.5 *
                    (lambda1 * (rho_roe * (vr - vl) + v_roe * (rhor - rhol)) + lambda2 * (alpha1 * ny + alpha2 * ny));
            term1 = (pr - pl) / (Gamma - 1.0);
            term2 = 0.5 * (u_roe * u_roe + v_roe * v_roe) * (rhor - rhol);
            term3 = rho_roe * (u_roe * (ur - ul) + v_roe * (vr - vl));
            fd[3] = 0.5 * (lambda1 * (term1 + term2 + term3) + lambda2 * (alpha1 * r1_E + alpha2 * r2_E));

            //final dissipation terms
            diss[0][i][j] = diss[0][i][j] - fd[0] * ds;
            diss[1][i][j] = diss[1][i][j] - fd[1] * ds;
            diss[2][i][j] = diss[2][i][j] - fd[2] * ds;
            diss[3][i][j] = diss[3][i][j] - fd[3] * ds;

            diss[0][im1][j] = diss[0][im1][j] + fd[0] * ds;
            diss[1][im1][j] = diss[1][im1][j] + fd[1] * ds;
            diss[2][im1][j] = diss[2][im1][j] + fd[2] * ds;
            diss[3][im1][j] = diss[3][im1][j] + fd[3] * ds;
        }
    }

    //compute in j-direction
    for (int i = 2; i <= ib; i++)
    {
        for (int j = 2; j <= jd1; j++)
        {
            jm1 = j - 1;
            jp1 = j + 1;

            ds = sqrt(sj[0][i][j] * sj[0][i][j] + sj[1][i][j] * sj[1][i][j]);
            nx = sj[0][i][j] / ds;
            ny = sj[1][i][j] / ds;
            dy = 0.5 * (y[i][j + 1] - y[i][j - 1]);

            // area_avg = 0.5*(area[i][jm1]+area[i][j]);

            //left and right states

            /*deltl[0] = 0.5*duj[0][i][jm1]*Minmod_Flux(duj[0][i][j], duj[0][i][jm1]);
            deltr[0] = 0.5*duj[0][i][j]*Minmod_Flux(duj[0][i][jp1], duj[0][i][j]);
            deltl[1] = 0.5*duj[1][i][jm1]*Minmod_Flux(duj[1][i][j], duj[1][i][jm1]);
            deltr[1] = 0.5*duj[1][i][j]*Minmod_Flux(duj[1][i][jp1], duj[1][i][j]);
            deltl[2] = 0.5*duj[2][i][jm1]*Minmod_Flux(duj[2][i][j], duj[2][i][jm1]);
            deltr[2] = 0.5*duj[2][i][j]*Minmod_Flux(duj[2][i][jp1], duj[2][i][j]);
            deltl[3] = 0.5*duj[3][i][jm1]*Minmod_Flux(duj[3][i][j], duj[3][i][jm1]);
            deltr[3] = 0.5*duj[3][i][j]*Minmod_Flux(duj[3][i][jp1], duj[3][i][j]);*/

            /*deltl[0] = 0.5*duj[0][i][jm1]*VanAlbada_Limiter(duj[0][i][j], duj[0][i][jm1]);
            deltr[0] = 0.5*duj[0][i][j]*VanAlbada_Limiter(duj[0][i][jp1], duj[0][i][j]);
            deltl[1] = 0.5*duj[1][i][jm1]*VanAlbada_Limiter(duj[1][i][j], duj[1][i][jm1]);
            deltr[1] = 0.5*duj[1][i][j]*VanAlbada_Limiter(duj[1][i][jp1], duj[1][i][j]);
            deltl[2] = 0.5*duj[2][i][jm1]*VanAlbada_Limiter(duj[2][i][j], duj[2][i][jm1]);
            deltr[2] = 0.5*duj[2][i][j]*VanAlbada_Limiter(duj[2][i][jp1], duj[2][i][j]);
            deltl[3] = 0.5*duj[3][i][jm1]*VanAlbada_Limiter(duj[3][i][j], duj[3][i][jm1]);
            deltr[3] = 0.5*duj[3][i][j]*VanAlbada_Limiter(duj[3][i][jp1], duj[3][i][j]);*/

            deltl[0] = 0.5 * Venki_Limiter(duj[0][i][j], duj[0][i][jm1], dy);
            deltr[0] = 0.5 * Venki_Limiter(duj[0][i][jp1], duj[0][i][j], dy);
            deltl[1] = 0.5 * Venki_Limiter(duj[1][i][j], duj[1][i][jm1], dy);
            deltr[1] = 0.5 * Venki_Limiter(duj[1][i][jp1], duj[1][i][j], dy);
            deltl[2] = 0.5 * Venki_Limiter(duj[2][i][j], duj[2][i][jm1], dy);
            deltr[2] = 0.5 * Venki_Limiter(duj[2][i][jp1], duj[2][i][j], dy);
            deltl[3] = 0.5 * Venki_Limiter(duj[3][i][j], duj[3][i][jm1], dy);
            deltr[3] = 0.5 * Venki_Limiter(duj[3][i][jp1], duj[3][i][j], dy);

            rhoinv = 1.0 / cv[0][i][jm1];
            rhol = cv[0][i][jm1] + deltl[0];
            ul = cv[1][i][jm1] * rhoinv + deltl[1];
            vl = cv[2][i][jm1] * rhoinv + deltl[2];
            pl = dv[0][i][jm1] + deltl[3];
            hl = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            El = pl / (gam1) + 0.5 * rhol * (ul * ul + vl * vl);
            al = sqrt(Gamma * pl / rhol);

            rhoinv = 1.0 / cv[0][i][j];
            rhor = cv[0][i][j] - deltr[0];
            ur = cv[1][i][j] * rhoinv - deltr[1];
            vr = cv[2][i][j] * rhoinv - deltr[2];
            pr = dv[0][i][j] - deltr[3];
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

            double w, rho_roe, u_roe, v_roe, a_roe, Vt, Vn, lambda1, lambda2;
            double r1_M, r1_x, r1_y, r1_E, r2_M, r2_x, r2_y, r2_E;
            double alpha1, alpha2, dVn, dVt;

            //roe avgs
            w = sqrt(rhor / rhol);
            rho_roe = w * rhol;
            u_roe = (ul + w * ur) / (1.0 + w);
            v_roe = (vl + w * vr) / (1.0 + w);
            a_roe = sqrt((al * al * sqrt(rhol) + ar * ar * sqrt(rhor)) / (sqrt(rhol) + sqrt(rhor)));

            Vn = u_roe * nx + v_roe * ny;
            Vt = -u_roe * ny + v_roe * nx;

            lambda1 = fabs(Vn);
            lambda2 = fabs(sqrt((Gamma - 1) / Gamma) * a_roe);

            r1_M = 0.0;
            r1_x = nx;
            r1_y = ny;
            r1_E = Vn - a_roe / sqrt(Gamma * (Gamma - 1));

            r2_M = 1.0;
            r2_x = nx;
            r2_y = ny;
            r2_E = Vn + a_roe / sqrt(Gamma * (Gamma - 1));

            //Wave strengths
            dVn = (ur - ul) * nx + (vr - vl) * ny;
            dVt = (ur - ul) * ny - (vr - vl) * nx;

            alpha1 = ((rho_roe * dVn) / 2.0) - sqrt(Gamma / (Gamma - 1.0)) * ((pr - pl) / (2.0 * a_roe));
            alpha2 = ((rho_roe * dVn) / 2.0) + sqrt(Gamma / (Gamma - 1.0)) * ((pr - pl) / (2.0 * a_roe));

            double term1, term2, term3;

            fd[0] = 0.5 * (lambda1 * (rhor - rhol) + lambda2 * 0.0);
            fd[1] = 0.5 *
                    (lambda1 * (rho_roe * (ur - ul) + u_roe * (rhor - rhol)) + lambda2 * (alpha1 * nx + alpha2 * nx));
            fd[2] = 0.5 *
                    (lambda1 * (rho_roe * (vr - vl) + v_roe * (rhor - rhol)) + lambda2 * (alpha1 * ny + alpha2 * ny));
            term1 = (pr - pl) / (Gamma - 1.0);
            term2 = 0.5 * (u_roe * u_roe + v_roe * v_roe) * (rhor - rhol);
            term3 = rho_roe * (u_roe * (ur - ul) + v_roe * (vr - vl));
            fd[3] = 0.5 * (lambda1 * (term1 + term2 + term3) + lambda2 * (alpha1 * r1_E + alpha2 * r2_E));

            //final dissipation terms
            diss[0][i][j] = diss[0][i][j] - fd[0] * ds;
            diss[1][i][j] = diss[1][i][j] - fd[1] * ds;
            diss[2][i][j] = diss[2][i][j] - fd[2] * ds;
            diss[3][i][j] = diss[3][i][j] - fd[3] * ds;

            diss[0][i][jm1] = diss[0][i][jm1] + fd[0] * ds;
            diss[1][i][jm1] = diss[1][i][jm1] + fd[1] * ds;
            diss[2][i][jm1] = diss[2][i][jm1] + fd[2] * ds;
            diss[3][i][jm1] = diss[3][i][jm1] + fd[3] * ds;
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

template <typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny)
{
    double L1r, L2r, L3r, L1l, L2l, L3l;
    L1r = fabs(ur * nx + vr * ny + ar);
    L1l = fabs(ul * nx + vl * ny + al);
    L2r = fabs(ur * nx + vr * ny);
    L2l = fabs(ul * nx + vl * ny);
    L3r = fabs(ur * nx + vr * ny - ar);
    L3l = fabs(ul * nx + vl * ny - al);

    //return Max3(max(L1l,L1r), max(L2l,L2r), max(L3l,L3r));
    return Max2(Max3(L1r, L2r, L3r), Max3(L1l, L2l, L3l));
}

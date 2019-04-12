#include "global_declarations.h"
#include "basic_functions.h"

template<typename T>
T ECCS(T Fr, T Fl, T Ur, T Ul,T VnL,T VnR, T dP, T Rho_I, T P_I);


void Diss_ECCS(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                     double ***&sj, double ***&diss) {

    int im1, jm1, ip1, jp1;
    double nx, ny, ds, rhol, ul, vl, rhor, ur, vr, pl, pr, hl, hr, El, Er, gam1, ggm1, Vnl, Vnr;
    double Rho_I, P_I, dP;
    double *fd, *fr, *fl, *Ur, *Ul, *deltl, *deltr, *mdiss;

    //double area_avg;

    fd = new double[nconv];
    fr = new double[nconv];
    fl = new double[nconv];
    Ur = new double[nconv];
    Ul = new double[nconv];
    deltl = new double[nconv];
    deltr = new double[nconv];
    mdiss = new double[nconv];

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
            ul   = dv[1][im1][j] + deltl[1];
            vl   = dv[2][im1][j] + deltl[2];
            pl   = dv[3][im1][j] + deltl[3];
            hl   = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            El   = pl / (gam1) + 0.5 * rhol * (ul * ul + vl * vl);
            Vnl  = ul * nx + vl * ny;

            rhor = dv[0][i][j] - deltr[0];
            ur   = dv[1][i][j]  - deltr[1];
            vr   = dv[2][i][j]  - deltr[2];
            pr   = dv[3][i][j] - deltr[3];
            hr   = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            Er = pr / (gam1) + 0.5 * rhor * (ur * ur + vr * vr);
            Vnr  = ur * nx + vr * ny;

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
			
			dP = pr - pl;
			Rho_I = 0.5*(rhor + rhol);
			P_I = 0.5*(pr+pl);

            for (int k=0; k<4; k++){
                fd[k] = 0.5*ECCS(fr[k], fl[k], Ur[k], Ul[k], Vnl, Vnr, dP, Rho_I, P_I) * (Ur[k] - Ul[k]);
            }

            //final dissipation terms
            for(int k=0; k<4; k++){
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
            ul   = dv[1][i][jm1] + deltl[1];
            vl   = dv[2][i][jm1] + deltl[2];
            pl   = dv[3][i][jm1] + deltl[3];
            hl   = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            El   = pl / (gam1) + 0.5 * rhol * (ul * ul + vl * vl);
            Vnl  = ul * nx + vl * ny;

            rhor = dv[0][i][j] - deltr[0];
            ur   = dv[1][i][j] - deltr[1];
            vr   = dv[2][i][j] - deltr[2];
            pr   = dv[3][i][j] - deltr[3];
            hr   = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            Er = pr / (gam1) + 0.5 * rhor * (ur * ur + vr * vr);
            Vnr  = ur * nx + vr * ny;

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
			
			dP = pr - pl;
			Rho_I = 0.5*(rhor + rhol);
			P_I = 0.5*(pr+pl);

			for (int k=0; k<4; k++){
                fd[k] = 0.5*ECCS(fr[k], fl[k], Ur[k], Ul[k], Vnl, Vnr, dP, Rho_I, P_I) * (Ur[k] - Ul[k]);
			}

            //final dissipation terms
            for(int k=0; k<4; k++){
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
    delete[] mdiss;
}


template<typename T>
T ECCS(T Fr, T Fl, T Ur, T Ul, T Vn_L, T Vn_R,T dP, T Rho_I, T P_I ) {

    double S = 0.0;

    double epsilon = 1e-10, Max=0.0;
//***************************************************************************************************************
	if(dP<epsilon)
		dP=0.0;
	
	Max = std::max(fabs(Vn_L),fabs(Vn_R));
	
    if((fabs(Fr-Fl) < epsilon and fabs(Ur - Ul) < epsilon))
    {
      S = 0.5*(fabs(Vn_L) + fabs(Vn_R));
    }
	else
	{
		//S =  Max + Sign(dP)*sqrt(Gamma*P_I/Rho_I);
        S = 0.5*(fabs(Vn_L) + fabs(Vn_R))+Sign(dP)*sqrt(Gamma*P_I/Rho_I);

	}
    
    return S;
}


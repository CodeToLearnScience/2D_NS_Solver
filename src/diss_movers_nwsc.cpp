#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"


template<typename T>
T Entropy_Corr(T z, T d);
template<typename T>
T MoversNWSC(T Fr, T Fl, T Ur, T Ul, T L_min);

template<typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny);

template<typename T>
T MinEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny);


void Diss_MOVERS_NWSC(int ib, int id1, int jb, int jd1, double ***&cv, double ***&dv, double ***&si, double ***&sj,
                 double ***&diss) {

    int im1, jm1, ip1, jp1;
    double nx, ny, ds, rhol, ul, vl, rhor, ur, vr, pl, pr, hl, hr, El, Er, gam1, ggm1;
    double ar, al, max_eig, min_eig,dP,p_I,beta=0.21,Sensor,Alpha_P,u_I,Vnl,Vnr,epsilon = 1e-8,du,max2,a;
    double *fd, *fr, *fl, *Ur, *Ul, *deltl, *deltr;



    //std::cout<<"in this function"<<std::endl;
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
            //left and right states

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


            al = sqrt(Gamma*pl/rhol);
            ar = sqrt(Gamma*pr/rhor);

            min_eig = MinEigVal(ur,ul,vr,vl,ar,al,nx,ny);
            max_eig = MaxEigVal(ur,ul,vr,vl,ar,al,nx,ny);

            u_I = 0.5*(fabs(Vnr*ds) + fabs(Vnl*ds));
            dP = pr - pl;
            p_I = 0.5*(pr + pl);

            du = ur - ul;

            a= 1.0;// - beta;

            max2 = a*max_eig + (1.0-a)*u_I;

            if(fabs(dP) < epsilon)
                dP = 0.0;

            Sensor = fabs(dP/(2.0*p_I));


            for(int k=0; k<4; k++){
                Alpha_P = MoversNWSC(fr[k], fl[k], Ur[k], Ul[k], min_eig);
                if(dP <= 0.0){Entropy_Corr(u_I,max2);}
                fd[k] = 0.5 * (beta*Sensor*Alpha_P + u_I*(Ur[k] - Ul[k]));
            }


            //final dissipation terms
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


            u_I = 0.5*(fabs(Vnr*ds) + fabs(Vnl*ds));
			dP = pr - pl;
			p_I = 0.5*(pr + pl);


            al = sqrt(Gamma*pl/rhol);
            ar = sqrt(Gamma*pr/rhor);

            min_eig = MinEigVal(ur,ul,vr,vl,ar,al,nx,ny);
            max_eig = MaxEigVal(ur,ul,vr,vl,ar,al,nx,ny);

			du = ur - ul;

			a= 1.0;// - beta;

            max2 = a*max_eig + (1.0-a)*u_I;

			if(fabs(dP) < epsilon)
			    dP = 0.0;

			Sensor = fabs(dP/(2.0*p_I));


            for(int k=0; k<4; k++){
                Alpha_P = MoversNWSC(fr[k], fl[k], Ur[k], Ul[k], min_eig);
               if(dP <= 0.0){Entropy_Corr(u_I,max2);}
                fd[k] = 0.5 * (beta*Sensor*Alpha_P + u_I*(Ur[k] - Ul[k]));
            }
            //final dissipation terms
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

}


template <typename T>
T MoversNWSC(T Fr, T Fl, T Ur,T Ul, T L_min){

    double S=0.0,dU=0.0,dF=0.0,epsilon = 1e-8;

    dU = Ur - Ul;
    dF = Fr - Fl;

    if (fabs(dU) < epsilon and fabs(dF) < epsilon)
        S=L_min;
    else
        S=Sign(dU)*fabs(dF);

    return S;
}

template<typename T>
T Entropy_Corr(T z, T d) {

    double k = 1.0 , d1 = k*d;

    if (fabs(z) > d)
        return d;
    else
        return 0.5 * (z * z + d1 * d1) / d1;
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
    return Max2(Max3(L1r, L2r, L3r), Max3(L1l, L2l, L3l));
}

template<typename T>
T MinEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny) {
    double L1r, L2r, L3r, L1l, L2l, L3l;
    L1r = fabs(ur * nx + vr * ny + ar);
    L1l = fabs(ul * nx + vl * ny + al);
    L2r = fabs(ur * nx + vr * ny);
    L2l = fabs(ul * nx + vl * ny);
    L3r = fabs(ur * nx + vr * ny - ar);
    L3l = fabs(ul * nx + vl * ny - al);

    //return Max3(min(L1l,L1r), min(L2l,L2r), min(L3l,L3r));
    return Max2(Min3(L1r, L2r, L3r), Min3(L1l, L2l, L3l));
}




/****************************************************MOVERS with No Wave Speed Correction ENDS ************************/

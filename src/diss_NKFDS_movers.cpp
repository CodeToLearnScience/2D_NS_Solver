/*****************************************************SHRINATHs scheme*******************************************************/


#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

template<typename T>
T Movers(T Fr, T Fl, T Ur, T Ul, T L_max, T L_min);

template<typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny);

template<typename T>
T MinEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny);

template <typename T>
T Mahalanobis(T rhoL, T rhoR, T unL, T unR, T pL, T pR, T maxa, T emax);

template<typename T>
T rounding(T dec, T n);


void Diss_NKFDS_MOVERS(int ib, int id1, int jb, int jd1, double ***&cv, double ***&dv, double ***&si, double ***&sj,
                 double ***&diss) {

    int im1, jm1, ip1, jp1;
    double nx, ny, ds, rhol, ul, vl, rhor, ur, vr, pl, pr, hl, hr, El, Er, gam1, ggm1;
    double ar, al, max_eig, min_eig, fed, emax=0,temp;
    double *fd,  *fr, *fl, *Ur, *Ul, *deltl, *deltr;

    fd = new double[nconv];
    fr = new double[nconv];
    fl = new double[nconv];
    Ur = new double[nconv];
    Ul = new double[nconv];
    deltl = new double[nconv];
    deltr = new double[nconv];

    gam1 = Gamma - 1.0;
    ggm1 = Gamma / gam1;

    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= id1; i++) {

          temp = dv[3][i][j] / (gam1) + 0.5 * dv[0][i][j] * (dv[1][i][j] * dv[1][i][j] + dv[2][i][j] * dv[2][i][j]);
          if (temp > emax){
              emax = temp;
          }
        }
    }


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

           /* for(int k=0; k<4; k++){
                deltl[k] = 0.5*dui[k][im1][j]*Minmod_Flux(dui[k][i][j], dui[k][im1][j]);
                deltr[k] = 0.5*dui[k][i][j]*Minmod_Flux(dui[k][ip1][j], dui[k][i][j]);

                /*deltl[k] = 0.5*dui[k][im1][j]*VanAlbada_Limiter(dui[k][i][j], dui[k][im1][j]);
                deltr[k] = 0.5*dui[k][i][j]*VanAlbada_Limiter(dui[k][ip1][j], dui[k][i][j]);*/

                /*deltl[k] = 0.5 * Venki_Limiter(dui[k][i][j], dui[k][im1][j], dx);
                deltr[k] = 0.5 * Venki_Limiter(dui[k][ip1][j], dui[k][i][j], dx);
            }*/


            //left and right states

            rhol = dv[0][im1][j] ;
            ul   = dv[1][im1][j] ;
            vl   = dv[2][im1][j] ;
            pl   = dv[3][im1][j];
            al   = sqrt(Gamma*pl/rhol);
            hl   = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            El   = pl / (gam1) + 0.5 * rhol * (ul * ul + vl * vl);

            rhor = dv[0][i][j];
            ur   = dv[1][i][j];
            vr   = dv[2][i][j];
            pr   = dv[3][i][j];
            ar   = sqrt(Gamma*pr/rhor);
            hr   = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            Er = pr / (gam1) + 0.5 * rhor * (ur * ur + vr * vr);

            max_eig = MaxEigVal(ur, ul, vr, vl, ar, al, nx, ny);
            min_eig = MinEigVal(ur, ul, vr, vl, ar, al, nx, ny);

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

            fed = Mahalanobis(rhol, rhor, ul, ur, pl,pr,max_eig,emax);
            for(int k=0; k<4; k++){

                fd[k] = 0.5 * Movers(fr[k], fl[k], Ur[k], Ul[k], max_eig, min_eig) ;
                if (fed > 0){
                    fd[k] = fed;
                }
                fd[k] = fd[k]* (Ur[k] - Ul[k]);
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



            rhol = dv[0][i][jm1];
            ul   = dv[1][i][jm1];
            vl   = dv[2][i][jm1];
            pl   = dv[3][i][jm1];
            al   = sqrt(Gamma*pl/rhol);
            hl   = ggm1 * pl / rhol + 0.5 * (ul * ul + vl * vl);
            El   = pl / (gam1) + 0.5 * rhol * (ul * ul + vl * vl);


            rhor = dv[0][i][j];
            ur   = dv[1][i][j];
            vr   = dv[2][i][j];
            pr   = dv[3][i][j];
            ar   = sqrt(Gamma*pr/rhor);
            hr   = ggm1 * pr / rhor + 0.5 * (ur * ur + vr * vr);
            Er = pr / (gam1) + 0.5 * rhor * (ur * ur + vr * vr);

            max_eig = MaxEigVal(ur, ul, vr, vl, ar, al, nx, ny);
            min_eig = MinEigVal(ur, ul, vr, vl, ar, al, nx, ny);

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

            fed = Mahalanobis(rhol, rhor, ul, ur, pl,pr,max_eig,emax);
            for(int k=0; k<4; k++){

                fd[k] = 0.5 * Movers(fr[k], fl[k], Ur[k], Ul[k], max_eig, min_eig) ;
                if (fed > 0){
                    fd[k] = fed;
                }
                fd[k] = fd[k]* (Ur[k] - Ul[k]);
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


//template <typename T>
//T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){
//
//	double S;
//	const double epsilon=1e-4;
//
//	if (fabs(Fr-Fl)<epsilon) return 0.0;
//	else if (fabs(Ur-Ul)<epsilon) return L_min;
//	else if (fabs(Ur-Ul)>epsilon and fabs(Fr-Fl)>epsilon) S=fabs(((Fr-Fl)/(Ur-Ul)));
//	else S=L_min;
//
//	if (S<epsilon)	return 0;
//	else if ((S)>=L_max) return (L_max);
//	else if((S)<=L_min) return (L_min);
//	else return (S);
//	return S;
//}

/*template<typename T>
T Movers(T Fr, T Fl, T Ur, T Ul, T L_max, T L_min) {

    double S, S1;
    const double eps = 10e-10;

    if (fabs(Fr - Fl) <= eps && fabs(Ur - Ul) <= eps) S = L_min;
    else if (fabs(Fr - Fl) <= eps && fabs(Ur - Ul) > eps) S = fabs(((Fr - Fl) / (Ur - Ul)));
    else if (fabs(Fr - Fl) > eps && fabs(Ur - Ul) <= eps) S = L_min;

    if (fabs(Fr - Fl) > eps && fabs(Ur - Ul) > eps) S1 = fabs(((Fr - Fl) / (Ur - Ul)));

    if (fabs(S1) >= L_max) return L_max;
    else if (fabs(S1) <= L_min) return L_min;
    else return fabs(((Fr - Fl) / (Ur - Ul)));
}*/

template <typename T>
T Moversold(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){

	double S=0.0;

	if (Ur!=Ul) S=fabs(((Fr-Fl)/(Ur-Ul)));
	else S=L_min;

	if (fabs(S)>=L_max) return L_max;
	else if(fabs(S)<=L_min) return L_min;
	else return S;
}

template <typename T>
T Movers(T Fr, T Fl, T Ur,T Ul,T L_max, T L_min){
    double deltaF, deltaU,epsilon = 1e-16;
    double S=0.0;
    deltaF = rounding(Fr- Fl,3.0);
    deltaU = rounding(Ur -Ul,3.0);
    if (fabs(deltaF)<epsilon) return 0.0;
    else if (fabs(deltaU)<epsilon) return L_min;
    else if (fabs(deltaU)>epsilon and fabs(deltaF)>epsilon) S=fabs(deltaF/deltaU);
    else S=L_min;


    if (fabs(S)>=L_max) return L_max;
    else if(fabs(S)<=L_min) return L_min;
    else return S;
}


template<typename T>
T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny) {
    double L1r, L2r, L3r, L1l, L2l, L3l;
    L1r = fabs(ur * nx + vr * ny);
    L1l = fabs(ul * nx + vl * ny);
    L2r=fabs(ur*nx+vr*ny); 	   L2l=fabs(ul*nx+vl*ny);
    L3r=fabs(ur*nx+vr*ny)-ar;  L3l=fabs(ul*nx+vl*ny)-al;

    return std::max(Max3(L1l, L2l, L3l), Max3(L1r, L2r, L3r));

}


template<typename T>
T MinEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny) {
    double L1r, L2r, L3r, L1l, L2l, L3l;
    L1r = fabs(ur * nx + vr * ny);
    L1l = fabs(ul * nx + vl * ny);
    L2r=fabs(ur*nx+vr*ny); 	   L2l=fabs(ul*nx+vl*ny);
    L3r=fabs(ur*nx+vr*ny)-ar;  L3l=fabs(ul*nx+vl*ny)-al;

    return std::max(Min3(L1l, L2l, L3l), Min3(L1r, L2r, L3r));
}

template<typename T>
T Mahalanobis(T rhoL, T rhoR, T unL, T unR, T pL, T pR, T maxa, T emax) {
    double cv = cp/Gamma;
    double RR = cp-cv;
    double TL = pL/rhoL/RR;
    double TR = pR/rhoR/RR;
    double SL, SR, delS;
    double mnb, DD,ric;
    DD = (rhoR-rhoL)*log((rhoR/rhoL)*pow((TL/TR),2.5))+rhoL/(2*RR*TL)*((rhoR/rhoL)+ (TL/TR))*(unL-unR)*(unL-unR)+ 2.5*(rhoL*(TL-TR)/TR +rhoR*(TR-TL)/TL);
    SL = cv*log(pL/(rhoL*(Gamma-1))) - RR*log(rhoL);
    SR = cv*log(pR/(rhoR*(Gamma-1))) - RR*log(rhoR);
    delS = SL - SR;
    // if (fabs(delS/emax)<1) {delS = 0;}
    //ric = (0.5*(fabs(unL)+fabs(unR))); // maxa;   //
         ric =  maxa;   //

    /* if ric = maxa it is LLF. If ric = (0.5*(fabs(unL)+fabs(unR))) then it is RICCA type */
    if ((DD >0) and (fabs(delS)<=1*emax)) {mnb = 1.0*(-ric);}
    else mnb =0;
    return mnb;
}

template<typename T>
T rounding(T dec, T n){
    double sp;
    sp = round(pow(10,n)*dec);
    return sp/pow(10,n);

}

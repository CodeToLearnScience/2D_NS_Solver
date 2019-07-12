#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

template <typename T>
T Entropy_Corr(T z, T d);

//void Diss_Roe1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si, double ***&sj, double ***&diss){
//
//    int im1, jm1;
//    double nx, ny, ds, rhoinv, rhol, ul, vl, rhor, ur, vr, pl, pr, hl, hr, gam1, ggm1;
//    double *fd;
//    double rav, dd, dd1, uav, vav, hav, q2a, c2a, cav, uv, du, h1, h2, h3, h4, h5, delta, eabs1, eabs2, eabs4;
//   fd = new double [nconv];
//
//   gam1 = Gamma - 1.0;
//   ggm1 = Gamma/gam1;
//
//    //compute in i-direction
//    for(int j=2;j<=jb;j++){
//        for(int i=2;i<=id1;i++){
//
//           im1=i-1;
//
//            ds = sqrt(si[0][i][j]*si[0][i][j]+si[1][i][j]*si[1][i][j]);
//            nx = si[0][i][j]/ds;
//            ny = si[1][i][j]/ds;
//
//            //left and right states
//
//            rhoinv = 1.0/cv[0][im1][j];
//            rhol   = cv[0][im1][j];
//            ul     = cv[1][im1][j]*rhoinv;
//            vl     = cv[2][im1][j]*rhoinv;
//            pl     = dv[0][im1][j];
//            hl     = ggm1*pl/rhol + 0.5*(ul*ul+vl*vl);
//
//            rhoinv = 1.0/cv[0][i][j];
//            rhor   = cv[0][i][j];
//            ur     = cv[1][i][j]*rhoinv;
//            vr     = cv[2][i][j]*rhoinv;
//            pr     = dv[0][i][j];
//            hr     = ggm1*pr/rhor + 0.5*(ur*ur+vr*vr);
//
//            //roe avgs
//            rav   = sqrt(rhol*rhor);
//            dd    = rav/rhol;
//            dd1   = 1.0/(1.0+dd);
//            uav   = (ul+dd*ur)*dd1;
//            vav   = (vl+dd*vr)*dd1;
//            hav   = (hl+dd*hr)*dd1;
//            q2a   = 0.5*(uav*uav+vav*vav);
//            c2a   = gam1*(hav-q2a);
//            cav   = sqrt(c2a);
//            uv    = uav*nx + vav*ny;
//            du    = (ur-ul)*nx + (vr-vl)*ny;
//
//
//
//            //eigen-values
//            h1 = fabs(uv-cav);
//            h2 = fabs(uv);
//            h4 = fabs(uv+cav);
//            delta = 1.5*h4;
//
//            eabs1 = Entropy_Corr(h1, delta);
//            eabs2 = Entropy_Corr(h2, delta);
//            eabs4 = Entropy_Corr(h4, delta);
//
//            /*eabs1 = h1;
//            eabs2 = h2;
//            eabs4 = h4;*/
//
//            h1 = rav*cav*du;
//            h2 = eabs1*(pr-pl - h1)/(2.0*c2a);
//            h3 = eabs2*(rhor-rhol - (pr-pl)/c2a);
//            h4 = eabs2*rav;
//            h5 = eabs4*(pr-pl + h1)/(2.0*c2a);
//
//            fd[0] = 0.5*(h2 + h3 + h5);
//            fd[1] = 0.5*(h2*(uav-cav*nx) + h3*uav + h4*(ur-ul-du*nx) + h5*(uav+cav*nx));
//            fd[2] = 0.5*(h2*(vav-cav*ny) + h3*vav + h4*(vr-vl-du*ny) + h5*(vav+cav*ny));
//            fd[3] = 0.5*(h2*(hav-cav*uv) + h3*q2a + h4*(uav*(ur-ul)  + vav*(vr-vl)-uv*du) + h5*(hav+cav*uv));
//
//
//            //final dissipation terms
//
//            diss[0][i][j] = diss[0][i][j] - fd[0]*ds;
//            diss[1][i][j] = diss[1][i][j] - fd[1]*ds;
//            diss[2][i][j] = diss[2][i][j] - fd[2]*ds;
//            diss[3][i][j] = diss[3][i][j] - fd[3]*ds;
//
//            diss[0][im1][j] = diss[0][im1][j] + fd[0]*ds;
//            diss[1][im1][j] = diss[1][im1][j] + fd[1]*ds;
//            diss[2][im1][j] = diss[2][im1][j] + fd[2]*ds;
//            diss[3][im1][j] = diss[3][im1][j] + fd[3]*ds;
//        }
//    }
//
//       //compute in j-direction
//    for(int i=2;i<=ib;i++){
//        for(int j=2;j<=jd1;j++){
//
//            jm1=j-1;
//
//            ds = sqrt(sj[0][i][j]*sj[0][i][j]+sj[1][i][j]*sj[1][i][j]);
//            nx = sj[0][i][j]/ds;
//            ny = sj[1][i][j]/ds;
//
//            //left and right states
//
//            rhoinv = 1.0/cv[0][i][jm1];
//            rhol   = cv[0][i][jm1];
//            ul     = cv[1][i][jm1]*rhoinv;
//            vl     = cv[2][i][jm1]*rhoinv;
//            pl     = dv[0][i][jm1];
//            hl     = ggm1*pl/rhol + 0.5*(ul*ul+vl*vl);
//
//            rhoinv = 1.0/cv[0][i][j];
//            rhor   = cv[0][i][j];
//            ur     = cv[1][i][j]*rhoinv;
//            vr     = cv[2][i][j]*rhoinv;
//            pr     = dv[0][i][j];
//            hr     = ggm1*pr/rhor + 0.5*(ur*ur+vr*vr);
//
//            //roe avgs
//            rav   = sqrt(rhol*rhor);
//            dd    = rav/rhol;
//            dd1   = 1.0/(1.0+dd);
//            uav   = (ul+dd*ur)*dd1;
//            vav   = (vl+dd*vr)*dd1;
//            hav   = (hl+dd*hr)*dd1;
//            q2a   = 0.5*(uav*uav+vav*vav);
//            c2a   = gam1*(hav-q2a);
//            cav   = sqrt(c2a);
//            uv    = uav*nx + vav*ny;
//            du    = (ur-ul)*nx + (vr-vl)*ny;
//
//
//
//            //eigen-values
//            h1 = fabs(uv-cav);
//            h2 = fabs(uv);
//            h4 = fabs(uv+cav);
//            delta = 1.5*h4;
//
//            eabs1 = Entropy_Corr(h1, delta);
//            eabs2 = Entropy_Corr(h2, delta);
//            eabs4 = Entropy_Corr(h4, delta);
//
//            /*eabs1 = h1;
//            eabs2 = h2;
//            eabs4 = h4;*/
//
//            h1 = rav*cav*du;
//            h2 = eabs1*(pr-pl - h1)/(2.0*c2a);
//            h3 = eabs2*(rhor-rhol - (pr-pl)/c2a);
//            h4 = eabs2*rav;
//            h5 = eabs4*(pr-pl + h1)/(2.0*c2a);
//
//            fd[0] = 0.5*(h2 + h3 + h5);
//            fd[1] = 0.5*(h2*(uav-cav*nx) + h3*uav + h4*(ur-ul-du*nx) + h5*(uav+cav*nx));
//            fd[2] = 0.5*(h2*(vav-cav*ny) + h3*vav + h4*(vr-vl-du*ny) + h5*(vav+cav*ny));
//            fd[3] = 0.5*(h2*(hav-cav*uv) + h3*q2a + h4*(uav*(ur-ul)  + vav*(vr-vl)-uv*du) + h5*(hav+cav*uv));
//
//            //final dissipation terms
//
//            diss[0][i][j] = diss[0][i][j] - fd[0]*ds;
//            diss[1][i][j] = diss[1][i][j] - fd[1]*ds;
//            diss[2][i][j] = diss[2][i][j] - fd[2]*ds;
//            diss[3][i][j] = diss[3][i][j] - fd[3]*ds;
//
//            diss[0][i][jm1] = diss[0][i][jm1] + fd[0]*ds;
//            diss[1][i][jm1] = diss[1][i][jm1] + fd[1]*ds;
//            diss[2][i][jm1] = diss[2][i][jm1] + fd[2]*ds;
//            diss[3][i][jm1] = diss[3][i][jm1] + fd[3]*ds;
//        }
//    }
//
//    /*for(int j=2;j<=jb;j++){
//        for(int i=2;i<=ib;i++){
//            std::cout<<"the computed dissipation values="<<diss[3][i][j]<<std::endl;
//        }
//    }*/
//
//    delete [] fd;
//
//}
/*
void Diss_Roe1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si, double ***&sj, double ***&diss){
    
    int im1, jm1;
    double nx, ny, ds, rhoinv, rhol, ul, vl, rhor, ur, vr, pl, pr, hl, hr, gam1, ggm1;
    double *fd;
    double w, rho_roe, u_roe, v_roe, h_roe, a_roe, Vt, Vn, lambda1, lambda2, lambda3, lambda4;
    double r1_M, r1_x, r1_y, r1_E,r2_M, r2_x, r2_y, r2_E,r3_M, r3_x, r3_y, r3_E,r4_M, r4_x, r4_y, r4_E;
    double alpha1, alpha2, alpha3, alpha4, dVn, dVt, delta;
    
    fd = new double [nconv]; 
    
    gam1 = Gamma - 1.0;
    ggm1 = Gamma/gam1;
   
    //compute in i-direction
    for(int j=2;j<=jb;j++){
        for(int i=2;i<=id1;i++){
            
           im1=i-1;
            
            ds = sqrt(si[0][i][j]*si[0][i][j]+si[1][i][j]*si[1][i][j]);
            nx = si[0][i][j]/ds;
            ny = si[1][i][j]/ds;
            
            //left and right states
            
            rhoinv = 1.0/cv[0][im1][j];
            rhol   = cv[0][im1][j];
            ul     = cv[1][im1][j]*rhoinv;
            vl     = cv[2][im1][j]*rhoinv;
            pl     = dv[0][im1][j];
            hl     = ggm1*pl/rhol + 0.5*(ul*ul+vl*vl);
            
            rhoinv = 1.0/cv[0][i][j];
            rhor   = cv[0][i][j];
            ur     = cv[1][i][j]*rhoinv;
            vr     = cv[2][i][j]*rhoinv;
            pr     = dv[0][i][j];
            hr     = ggm1*pr/rhor + 0.5*(ur*ur+vr*vr);
            
            
            double w, rho_roe, u_roe, v_roe, h_roe, a_roe, Vt, Vn, lambda1, lambda2, lambda3, lambda4;
            double r1_M, r1_x, r1_y, r1_E,r2_M, r2_x, r2_y, r2_E,r3_M, r3_x, r3_y, r3_E,r4_M, r4_x, r4_y, r4_E;
            double alpha1, alpha2, alpha3, alpha4, dVn, dVt;	
            //roe avgs
            w       = sqrt(rhol*rhor);
            rho_roe = w*rhol;
            u_roe   = (ul+w*ur)/(1.0+w);
            v_roe   = (vl+w*vr)/(1.0+w);
            h_roe   = (hl+w*hr)/(1.0+w);
            a_roe   = sqrt((h_roe-0.5*(u_roe*u_roe+v_roe*v_roe))*gam1);
            
            Vn      = u_roe*nx+v_roe*ny;
            Vt      = -u_roe*ny+v_roe*nx;
            
            lambda1=fabs(Vn-a_roe);
            lambda2=fabs(Vn);
            lambda3=lambda2;
            lambda4=fabs(Vn+a_roe);
            
            delta = 0.5*lambda4;
            
            lambda1 = Entropy_Corr(lambda1, delta);
            lambda2 = Entropy_Corr(lambda2, delta);
            lambda3 = Entropy_Corr(lambda3, delta);
            lambda4 = Entropy_Corr(lambda4, delta);
            
            r1_M=1.0;
            r1_x=u_roe-a_roe*nx;
            r1_y=v_roe-a_roe*ny;
            r1_E=h_roe-Vn*a_roe;
	
            r2_M=1.0;
            r2_x=u_roe;
            r2_y=v_roe;
            r2_E=0.5*(u_roe*u_roe+v_roe*v_roe);
			
            r3_M=0.0;
            r3_x=-ny;	
            r3_y=nx;	
            r3_E=Vt;
				
            r4_M=1.0;
            r4_x=u_roe+a_roe*nx;
            r4_y=v_roe+a_roe*ny;
            r4_E=h_roe+Vn*a_roe;
            
            //Wave strengths	
            
            dVn = (ur-ul)*nx + (vr-vl)*ny;
            dVt = -(ur-ul)*ny + (vr-vl)*nx;
			
            alpha1=((pr-pl)-(rho_roe*a_roe*dVn))/(2.0*a_roe*a_roe);
            alpha2=-((pr-pl)-(rhor-rhol)*a_roe*a_roe)/(a_roe*a_roe);
            alpha3=rho_roe*dVt;
            alpha4=((pr-pl)+(rho_roe*a_roe*dVn))/(2.0*a_roe*a_roe);

            fd[0] = 0.5*(alpha1*lambda1*r1_M+alpha2*lambda2*r2_M+alpha3*lambda3*r3_M+alpha4*lambda4*r4_M);
            fd[1] = 0.5*(alpha1*lambda1*r1_x+alpha2*lambda2*r2_x+alpha3*lambda3*r3_x+alpha4*lambda4*r4_x);
            fd[2] = 0.5*(alpha1*lambda1*r1_y+alpha2*lambda2*r2_y+alpha3*lambda3*r3_y+alpha4*lambda4*r4_y);
            fd[3] = 0.5*(alpha1*lambda1*r1_E+alpha2*lambda2*r2_E+alpha3*lambda3*r3_E+alpha4*lambda4*r4_E);
            
            //final dissipation terms
            
            diss[0][i][j] = diss[0][i][j] - fd[0]*ds;
            diss[1][i][j] = diss[1][i][j] - fd[1]*ds;
            diss[2][i][j] = diss[2][i][j] - fd[2]*ds;
            diss[3][i][j] = diss[3][i][j] - fd[3]*ds;
            
            diss[0][im1][j] = diss[0][im1][j] + fd[0]*ds;
            diss[1][im1][j] = diss[1][im1][j] + fd[1]*ds;
            diss[2][im1][j] = diss[2][im1][j] + fd[2]*ds;                                 
            diss[3][im1][j] = diss[3][im1][j] + fd[3]*ds;                                 
        }
    }
    
       //compute in j-direction
    for(int i=2;i<=ib;i++){
        for(int j=2;j<=jd1;j++){
            
            jm1=j-1;
            
            ds = sqrt(sj[0][i][j]*sj[0][i][j]+sj[1][i][j]*sj[1][i][j]);
            nx = sj[0][i][j]/ds;
            ny = sj[1][i][j]/ds;
            
            //left and right states
            
            rhoinv = 1.0/cv[0][i][jm1];
            rhol   = cv[0][i][jm1];
            ul     = cv[1][i][jm1]*rhoinv;
            vl     = cv[2][i][jm1]*rhoinv;
            pl     = dv[0][i][jm1];
            hl     = ggm1*pl/rhol + 0.5*(ul*ul+vl*vl);
            
            rhoinv = 1.0/cv[0][i][j];
            rhor   = cv[0][i][j];
            ur     = cv[1][i][j]*rhoinv;
            vr     = cv[2][i][j]*rhoinv;
            pr     = dv[0][i][j];
            hr     = ggm1*pr/rhor + 0.5*(ur*ur+vr*vr);
                        
            w       = sqrt(rhol*rhor);
            rho_roe = w*rhol;
            u_roe   = (ul+w*ur)/(1.0+w);
            v_roe   = (vl+w*vr)/(1.0+w);
            h_roe   = (hl+w*hr)/(1.0+w);
            a_roe   = sqrt((h_roe-0.5*(u_roe*u_roe+v_roe*v_roe))*gam1);
            
            Vn      = u_roe*nx+v_roe*ny;
            Vt      = -u_roe*ny+v_roe*nx;
            
            lambda1=fabs(Vn-a_roe);
            lambda2=fabs(Vn);
            lambda3=lambda2;
            lambda4=fabs(Vn+a_roe);
            
            delta = 0.5*lambda4;
            
            lambda1 = Entropy_Corr(lambda1, delta);
            lambda2 = Entropy_Corr(lambda2, delta);
            lambda3 = Entropy_Corr(lambda3, delta);
            lambda4 = Entropy_Corr(lambda4, delta);
            
            r1_M=1.0;
            r1_x=u_roe-a_roe*nx;
            r1_y=v_roe-a_roe*ny;
            r1_E=h_roe-Vn*a_roe;
	
            r2_M=1.0;
            r2_x=u_roe;
            r2_y=v_roe;
            r2_E=0.5*(u_roe*u_roe+v_roe*v_roe);
			
            r3_M=0.0;
            r3_x=-ny;	
            r3_y=nx;	
            r3_E=Vt;
				
            r4_M=1.0;
            r4_x=u_roe+a_roe*nx;
            r4_y=v_roe+a_roe*ny;
            r4_E=h_roe+Vn*a_roe;
            
            //Wave strengths	
            
            dVn = (ur-ul)*nx + (vr-vl)*ny;
            dVt = -(ur-ul)*ny + (vr-vl)*nx;
			
            alpha1=((pr-pl)-(rho_roe*a_roe*dVn))/(2.0*a_roe*a_roe);
            alpha2=-((pr-pl)-(rhor-rhol)*a_roe*a_roe)/(a_roe*a_roe);
            alpha3=rho_roe*dVt;
            alpha4=((pr-pl)+(rho_roe*a_roe*dVn))/(2.0*a_roe*a_roe);

            fd[0] = 0.5*(alpha1*lambda1*r1_M+alpha2*lambda2*r2_M+alpha3*lambda3*r3_M+alpha4*lambda4*r4_M);
            fd[1] = 0.5*(alpha1*lambda1*r1_x+alpha2*lambda2*r2_x+alpha3*lambda3*r3_x+alpha4*lambda4*r4_x);
            fd[2] = 0.5*(alpha1*lambda1*r1_y+alpha2*lambda2*r2_y+alpha3*lambda3*r3_y+alpha4*lambda4*r4_y);
            fd[3] = 0.5*(alpha1*lambda1*r1_E+alpha2*lambda2*r2_E+alpha3*lambda3*r3_E+alpha4*lambda4*r4_E);
            
            //final dissipation terms
            
            diss[0][i][j] = diss[0][i][j] - fd[0]*ds;
            diss[1][i][j] = diss[1][i][j] - fd[1]*ds;
            diss[2][i][j] = diss[2][i][j] - fd[2]*ds;
            diss[3][i][j] = diss[3][i][j] - fd[3]*ds;
            
            diss[0][i][jm1] = diss[0][i][jm1] + fd[0]*ds;
            diss[1][i][jm1] = diss[1][i][jm1] + fd[1]*ds;
            diss[2][i][jm1] = diss[2][i][jm1] + fd[2]*ds;       
            diss[3][i][jm1] = diss[3][i][jm1] + fd[3]*ds;       
        }
    }*/

/*for(int j=2;j<=jb;j++){
    for(int i=2;i<=ib;i++){
        std::cout<<"the computed dissipation values="<<diss[3][i][j]<<std::endl;
    }
}*/

/*   delete [] fd;

}*/

/*******************************************************************************************/
/*void Diss_Roe1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si, double ***&sj, double ***&diss){
    
    int im1, jm1;
    double nx, ny, ds, rhoinv, rhol, ul, vl, rhor, ur, vr, pl, pr, hl, hr, gam1, ggm1;
    double *fd;
    double w, rho_roe, u_roe, v_roe, h_roe, a_roe, Vt, Vn, lambda1, lambda2, lambda3, lambda4;
    double r1_M, r1_x, r1_y, r1_E,r2_M, r2_x, r2_y, r2_E,r3_M, r3_x, r3_y, r3_E,r4_M, r4_x, r4_y, r4_E;
    double alpha1, alpha2, alpha3, alpha4, dVn, dVt, delta;
    
    fd = new double [nconv]; 
    
    gam1 = Gamma - 1.0;
    ggm1 = Gamma/gam1;
   
    //compute in i-direction
    for(int j=2;j<=jb;j++){
        for(int i=2;i<=id1;i++){
            
           im1=i-1;
            
            ds = sqrt(si[0][i][j]*si[0][i][j]+si[1][i][j]*si[1][i][j]);
            nx = si[0][i][j]/ds;
            ny = si[1][i][j]/ds;
            
            //left and right states
            
            rhoinv = 1.0/cv[0][im1][j];
            rhol   = cv[0][im1][j];
            ul     = cv[1][im1][j]*rhoinv;
            vl     = cv[2][im1][j]*rhoinv;
            pl     = dv[0][im1][j];
            hl     = ggm1*pl/rhol + 0.5*(ul*ul+vl*vl);
            
            rhoinv = 1.0/cv[0][i][j];
            rhor   = cv[0][i][j];
            ur     = cv[1][i][j]*rhoinv;
            vr     = cv[2][i][j]*rhoinv;
            pr     = dv[0][i][j];
            hr     = ggm1*pr/rhor + 0.5*(ur*ur+vr*vr);
            
            
            double w, rho_roe, u_roe, v_roe, h_roe, a_roe, Vt, Vn, lambda1, lambda2, lambda3, lambda4;
            double r1_M, r1_x, r1_y, r1_E,r2_M, r2_x, r2_y, r2_E,r3_M, r3_x, r3_y, r3_E,r4_M, r4_x, r4_y, r4_E;
            double alpha1, alpha2, alpha3, alpha4, dVn, dVt;	
            //roe avgs
            w       = sqrt(rhor/rhol);
            rho_roe = w*rhol;
            u_roe   = (ul+w*ur)/(1.0+w);
            v_roe   = (vl+w*vr)/(1.0+w);
            h_roe   = (hl+w*hr)/(1.0+w);
            a_roe   = sqrt((h_roe-0.5*(u_roe*u_roe+v_roe*v_roe))*gam1);
            
            Vn      = u_roe*nx+v_roe*ny;
            Vt      = -u_roe*ny+v_roe*nx;
            
            lambda1=fabs(Vn-a_roe);
            lambda2=fabs(Vn);
            lambda3=lambda2;
            lambda4=fabs(Vn+a_roe);
            
            delta = 0.5*lambda4;
            
            lambda1 = Entropy_Corr(lambda1, delta);
            lambda2 = Entropy_Corr(lambda2, delta);
            lambda3 = Entropy_Corr(lambda3, delta);
            lambda4 = Entropy_Corr(lambda4, delta);
            
            r1_M=1.0;
            r1_x=u_roe-a_roe*nx;
            r1_y=v_roe-a_roe*ny;
            r1_E=h_roe-Vn*a_roe;
	
            r2_M=1.0;
            r2_x=u_roe;
            r2_y=v_roe;
            r2_E=0.5*(u_roe*u_roe+v_roe*v_roe);
			
            r3_M=0.0;
            r3_x=ny;	
            r3_y=-nx;	
            //r3_E=Vt;
	    r3_E = u_roe*ny-v_roe*nx;		
				
            r4_M=1.0;
            r4_x=u_roe+a_roe*nx;
            r4_y=v_roe+a_roe*ny;
            r4_E=h_roe+Vn*a_roe;
            
            //Wave strengths	
            dVn = (ur-ul)*nx + (vr-vl)*ny;
            //dVt = -(ur-ul)*ny + (vr-vl)*nx;
	    dVt = (ur-ul)*ny - (vr-vl)*nx;
			
            alpha1=((pr-pl)-(rho_roe*a_roe*dVn))/(2.0*a_roe*a_roe);
            //alpha2=-((pr-pl)-(rhor-rhol)*a_roe*a_roe)/(a_roe*a_roe);
            alpha2=(rhor-rhol)-((pr-pl)/(a_roe*a_roe));
            alpha3=rho_roe*dVt;
            alpha4=((pr-pl)+(rho_roe*a_roe*dVn))/(2.0*a_roe*a_roe);

            fd[0] = 0.5*(alpha1*lambda1*r1_M + alpha2*lambda2*r2_M + alpha3*lambda3*r3_M + alpha4*lambda4*r4_M);
            fd[1] = 0.5*(alpha1*lambda1*r1_x + alpha2*lambda2*r2_x + alpha3*lambda3*r3_x + alpha4*lambda4*r4_x);
            fd[2] = 0.5*(alpha1*lambda1*r1_y + alpha2*lambda2*r2_y + alpha3*lambda3*r3_y + alpha4*lambda4*r4_y);
            fd[3] = 0.5*(alpha1*lambda1*r1_E + alpha2*lambda2*r2_E + alpha3*lambda3*r3_E + alpha4*lambda4*r4_E);
            
            //final dissipation terms
            
            diss[0][i][j] = diss[0][i][j] - fd[0]*ds;
            diss[1][i][j] = diss[1][i][j] - fd[1]*ds;
            diss[2][i][j] = diss[2][i][j] - fd[2]*ds;
            diss[3][i][j] = diss[3][i][j] - fd[3]*ds;
            
            diss[0][im1][j] = diss[0][im1][j] + fd[0]*ds;
            diss[1][im1][j] = diss[1][im1][j] + fd[1]*ds;
            diss[2][im1][j] = diss[2][im1][j] + fd[2]*ds;                                 
            diss[3][im1][j] = diss[3][im1][j] + fd[3]*ds;                                 
        }
    }
    
       //compute in j-direction
    for(int i=2;i<=ib;i++){
        for(int j=2;j<=jd1;j++){
            
            jm1=j-1;
            
            ds = sqrt(sj[0][i][j]*sj[0][i][j]+sj[1][i][j]*sj[1][i][j]);
            nx = sj[0][i][j]/ds;
            ny = sj[1][i][j]/ds;
            
            //left and right states
            
            rhoinv = 1.0/cv[0][i][jm1];
            rhol   = cv[0][i][jm1];
            ul     = cv[1][i][jm1]*rhoinv;
            vl     = cv[2][i][jm1]*rhoinv;
            pl     = dv[0][i][jm1];
            hl     = ggm1*pl/rhol + 0.5*(ul*ul+vl*vl);
            
            rhoinv = 1.0/cv[0][i][j];
            rhor   = cv[0][i][j];
            ur     = cv[1][i][j]*rhoinv;
            vr     = cv[2][i][j]*rhoinv;
            pr     = dv[0][i][j];
            hr     = ggm1*pr/rhor + 0.5*(ur*ur+vr*vr);
                        
            w       = sqrt(rhor/rhol);
            rho_roe = w*rhol;
            u_roe   = (ul+w*ur)/(1.0+w);
            v_roe   = (vl+w*vr)/(1.0+w);
            h_roe   = (hl+w*hr)/(1.0+w);
            a_roe   = sqrt((h_roe-0.5*(u_roe*u_roe+v_roe*v_roe))*gam1);
            
            Vn      = u_roe*nx+v_roe*ny;
            Vt      = -u_roe*ny+v_roe*nx;
            
            lambda1=fabs(Vn-a_roe);
            lambda2=fabs(Vn);
            lambda3=lambda2;
            lambda4=fabs(Vn+a_roe);
            
            delta = 0.5*lambda4;
            
            lambda1 = Entropy_Corr(lambda1, delta);
            lambda2 = Entropy_Corr(lambda2, delta);
            lambda3 = Entropy_Corr(lambda3, delta);
            lambda4 = Entropy_Corr(lambda4, delta);
            
            r1_M=1.0;
            r1_x=u_roe-a_roe*nx;
            r1_y=v_roe-a_roe*ny;
            r1_E=h_roe-Vn*a_roe;
	
            r2_M=1.0;
            r2_x=u_roe;
            r2_y=v_roe;
            r2_E=0.5*(u_roe*u_roe+v_roe*v_roe);
			
            r3_M=0.0;
            r3_x=ny;	
            r3_y=-nx;	
            //r3_E=Vt;
	    r3_E = u_roe*ny-v_roe*nx;
				
            r4_M=1.0;
            r4_x=u_roe+a_roe*nx;
            r4_y=v_roe+a_roe*ny;
            r4_E=h_roe+Vn*a_roe;
            
            //Wave strengths	
            dVn = (ur-ul)*nx + (vr-vl)*ny;
            //dVt = -(ur-ul)*ny + (vr-vl)*nx;
	    dVt = (ur-ul)*ny - (vr-vl)*nx;
			
            alpha1=((pr-pl)-(rho_roe*a_roe*dVn))/(2.0*a_roe*a_roe);
            //alpha2=-((pr-pl)-(rhor-rhol)*a_roe*a_roe)/(a_roe*a_roe);
            alpha2=(rhor-rhol)-((pr-pl)/(a_roe*a_roe));
            alpha3=rho_roe*dVt;
            alpha4=((pr-pl)+(rho_roe*a_roe*dVn))/(2.0*a_roe*a_roe);

            fd[0] = 0.5*(alpha1*lambda1*r1_M+alpha2*lambda2*r2_M+alpha3*lambda3*r3_M+alpha4*lambda4*r4_M);
            fd[1] = 0.5*(alpha1*lambda1*r1_x+alpha2*lambda2*r2_x+alpha3*lambda3*r3_x+alpha4*lambda4*r4_x);
            fd[2] = 0.5*(alpha1*lambda1*r1_y+alpha2*lambda2*r2_y+alpha3*lambda3*r3_y+alpha4*lambda4*r4_y);
            fd[3] = 0.5*(alpha1*lambda1*r1_E+alpha2*lambda2*r2_E+alpha3*lambda3*r3_E+alpha4*lambda4*r4_E);
            
            //final dissipation terms
            
            diss[0][i][j] = diss[0][i][j] - fd[0]*ds;
            diss[1][i][j] = diss[1][i][j] - fd[1]*ds;
            diss[2][i][j] = diss[2][i][j] - fd[2]*ds;
            diss[3][i][j] = diss[3][i][j] - fd[3]*ds;
            
            diss[0][i][jm1] = diss[0][i][jm1] + fd[0]*ds;
            diss[1][i][jm1] = diss[1][i][jm1] + fd[1]*ds;
            diss[2][i][jm1] = diss[2][i][jm1] + fd[2]*ds;       
            diss[3][i][jm1] = diss[3][i][jm1] + fd[3]*ds;       
        }
    }*/

/*for(int j=2;j<=jb;j++){
    for(int i=2;i<=ib;i++){
        std::cout<<"the computed dissipation values="<<diss[3][i][j]<<std::endl;
    }
}*/

/*    delete [] fd;

}*/

/*
template <typename T>
T Entropy_Corr(T z, T d){
	
    if(z>d) return z;
    else return 0.5*(z*z+d*d)/d;
}*/

/********************************************************************************Naveen's Upwing scheme*****************************************************************************/
void Diss_Roe1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
               double ***&sj, double ***&diss)
{

    int im1, jm1;
    double nx, ny, ds, rhoinv, rhol, ul, vl, rhor, ur, vr, pl, pr, al, ar, gam1, ggm1;
    double *fd;
    double w, rho_roe, u_roe, v_roe, h_roe, a_roe, Vt, Vn, lambda1, lambda2;
    double r1_M, r1_x, r1_y, r1_E, r2_M, r2_x, r2_y, r2_E;
    double alpha1, alpha2, dVn, dVt;

    fd = new double[nconv];

    //gam1 = Gamma - 1.0;
    //ggm1 = Gamma/gam1;

    //compute in i-direction
    for (int j = 2; j <= jb; j++)
    {
        for (int i = 2; i <= id1; i++)
        {

            im1 = i - 1;

            ds = sqrt(si[0][i][j] * si[0][i][j] + si[1][i][j] * si[1][i][j]);
            nx = si[0][i][j] / ds;
            ny = si[1][i][j] / ds;

            //left and right states

            rhoinv = 1.0 / dv[0][im1][j];
            rhol = dv[0][im1][j];
            ul = dv[1][im1][j];
            vl = dv[2][im1][j];
            pl = dv[3][im1][j];
            al = sqrt(Gamma * pl / rhol);

            rhoinv = 1.0 / dv[0][i][j];
            rhor = dv[0][i][j];
            ur = dv[1][i][j];
            vr = dv[2][i][j];
            pr = dv[3][i][j];
            ar = sqrt(Gamma * pr / rhor);

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
            fd[1] = 0.5 * (lambda1 * (rho_roe * (ur - ul) + u_roe * (rhor - rhol)) + lambda2 * (alpha1 * nx + alpha2 * nx));
            fd[2] = 0.5 * (lambda1 * (rho_roe * (vr - vl) + v_roe * (rhor - rhol)) + lambda2 * (alpha1 * ny + alpha2 * ny));
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

            ds = sqrt(sj[0][i][j] * sj[0][i][j] + sj[1][i][j] * sj[1][i][j]);
            nx = sj[0][i][j] / ds;
            ny = sj[1][i][j] / ds;

            //left and right states

            rhoinv = 1.0 / dv[0][i][jm1];
            rhol = dv[0][i][jm1];
            ul = dv[1][i][jm1];
            vl = dv[2][i][jm1];
            pl = dv[3][i][jm1];
            al = sqrt(Gamma * pl / rhol);

            rhoinv = 1.0 / dv[0][i][j];
            rhor = dv[0][i][j];
            ur = dv[1][i][j];
            vr = dv[2][i][j];
            pr = dv[3][i][j];
            ar = sqrt(Gamma * pr / rhor);

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
            std::cout<<"the computed dissipation values="<<diss[3][i][j]<<std::endl;
        }
    }*/

    delete[] fd;
}

/********************************************************************************Naveen's Upwing scheme*****************************************************************************/
void Diss_Roe1_TV(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                  double ***&sj, double ***&diss)
{

    int im1, jm1;
    double nx, ny, ds, rhoinv, rhol, ul, vl, rhor, ur, vr, pl, pr, al, ar, gam1, ggm1;
    double *fd, p_avg;
    double w, rho_roe, u_roe, v_roe, h_roe, a_roe, Vt, Vn, lambda1, lambda2, lambda3, lambda2_a, lambda3_a;
    double r1_M, r1_x, r1_y, r1_E, r2_M, r2_x, r2_y, r2_E;
    double alpha1, alpha2, dVn, dVt;
    double beta = 0.0;
    //double w, rho_roe, u_roe, v_roe, a_roe, Vt, Vn, lambda1, lambda2, lambda3, lambda2_a, lambda3_a;
    //double r1_M, r1_x, r1_y, r1_E,r2_M, r2_x, r2_y, r2_E;
    //double alpha1, alpha2, dVn, dVt;
    //double beta=0.0;

    fd = new double[nconv];

    //gam1 = Gamma - 1.0;
    //ggm1 = Gamma/gam1;

    //compute in i-direction
    for (int j = 2; j <= jb; j++)
    {
        for (int i = 2; i <= id1; i++)
        {

            im1 = i - 1;

            ds = sqrt(si[0][i][j] * si[0][i][j] + si[1][i][j] * si[1][i][j]);
            nx = si[0][i][j] / ds;
            ny = si[1][i][j] / ds;

            //left and right states

            rhoinv = 1.0 / dv[0][im1][j];
            rhol = dv[0][im1][j];
            ul = dv[1][im1][j];
            vl = dv[2][im1][j];
            pl = dv[3][im1][j];
            al = sqrt(Gamma * pl / rhol);

            rhoinv = 1.0 / dv[0][i][j];
            rhor = dv[0][i][j];
            ur = dv[1][i][j];
            vr = dv[2][i][j];
            pr = dv[3][i][j];
            ar = sqrt(Gamma * pr / rhor);

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

            double term1, term2;

            fd[0] = 0.5 * (lambda1 * (rhor - rhol));
            fd[1] = 0.5 * (lambda1 * (rho_roe * (ur - ul) + u_roe * (rhor - rhol)) +
                           (alpha1 * lambda2 + alpha2 * lambda3) * nx);
            fd[2] = 0.5 * (lambda1 * (rho_roe * (vr - vl) + v_roe * (rhor - rhol)) +
                           (alpha1 * lambda2 + alpha2 * lambda3) * ny);
            term1 = 0.5 * (u_roe * u_roe + v_roe * v_roe) * (rhor - rhol);
            term2 = rho_roe * (u_roe * (ur - ul) + v_roe * (vr - vl));
            fd[3] = 0.5 * (lambda1 * (term1 + term2) + lambda2 * alpha1 * r1_E + alpha2 * lambda3 * r2_E);

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

            double term1, term2;

            fd[0] = 0.5 * (lambda1 * (rhor - rhol));
            fd[1] = 0.5 * (lambda1 * (rho_roe * (ur - ul) + u_roe * (rhor - rhol)) +
                           (alpha1 * lambda2 + alpha2 * lambda3) * nx);
            fd[2] = 0.5 * (lambda1 * (rho_roe * (vr - vl) + v_roe * (rhor - rhol)) +
                           (alpha1 * lambda2 + alpha2 * lambda3) * ny);
            term1 = 0.5 * (u_roe * u_roe + v_roe * v_roe) * (rhor - rhol);
            term2 = rho_roe * (u_roe * (ur - ul) + v_roe * (vr - vl));
            fd[3] = 0.5 * (lambda1 * (term1 + term2) + lambda2 * alpha1 * r1_E + alpha2 * lambda3 * r2_E);

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
            std::cout<<"the computed dissipation values="<<diss[3][i][j]<<std::endl;
        }
    }*/

    delete[] fd;
}

template <typename T>
T Entropy_Corr(T z, T d)
{

    if (z > d)
        return z;
    else
        return 0.5 * (z * z + d * d) / d;
}

#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

void Time_Step_Euler(int ib, int id1, int id2, int jb, int jd1, int jd2, int nconv, int ndepv, double ***&cv,
                     double ***&dv, double **&area, double ***&si, double ***&sj, double **&tstep, double **&sri,
                     double **&srj, double ***&epsij) {

    double rh, u, v, sx, sy, nx, ny, ds, vc, cs, tsmin, f1, f2, fac, dtv, srvi, srvj, reval, reval1, cflrat, ex, ey;
    const double cfac=2.0;

    time_step == "global";

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            tstep[i][j] = 0.0;
        }
    }

    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= ib; i++) {
            rh = 1.0/dv[0][i][j];
            u  = dv[1][i][j];
            v  = dv[2][i][j];

            //i-direction
            sx = 0.5 * (si[0][i][j] + si[0][i + 1][j]);
            sy = 0.5 * (si[1][i][j] + si[1][i + 1][j]);

            ds = sqrt(sx * sx + sy * sy);

            nx = sx / ds;
            ny = sx / ds;
            vc = u * nx + v * ny;
            sri[i][j] = (fabs(vc) + dv[5][i][j])*ds;


            //j-direction
            sx = 0.5 * (sj[0][i][j] + sj[0][i][j + 1]);
            sy = 0.5 * (sj[1][i][j] + sj[1][i][j + 1]);

            ds = sqrt(sx * sx + sy * sy);

            nx = sx / ds;
            ny = sx / ds;
            vc = u * nx + v * ny;
            srj[i][j] = (fabs(vc) + dv[5][i][j])*ds;


            tstep[i][j] = cfl*area[i][j] / (sri[i][j] + srj[i][j]);

            if(tstep[i][j]==0 or tstep[i][j]<0.0){
                std::cout<<"values inside time_step are:"<<"srvi="<<tstep[i][j]<<std::endl;
                exit(0);
            }

            if (tstep[i][j] < 0 or tstep[i][j] != tstep[i][j]) {
                std::cout << "t became negative at i=" << i << "\t" << "and at j=" << j << std::endl;
                std::cout << rh << "\t" << u << "\t" << v << "\t" << "\t" << cs << "\t" << sri[i][j] << "\t"
                          << srj[i][j] << "\t" << area[i][j] << std::endl;
                exit(0);
            }
        }
    }

    //std::cout<<"values inside time_step are:"<<"srvi="<<srvj<<std::endl;
    /****values of sri and srj at dummy nodes****/
    for (int i = 2; i <= ib; i++) {
        sri[i][1] = sri[i][2];
        sri[i][0] = sri[i][2];
        sri[i][jd1] = sri[i][jb];
        sri[i][jd2] = sri[i][jb];
        srj[i][1] = srj[i][2];
        srj[i][0] = srj[i][2];
        srj[i][jd1] = srj[i][jb];
        srj[i][jd2] = srj[i][jb];
    }

    for (int j = 2; j <= jb; j++) {
        sri[1][j] = sri[2][j];
        sri[0][j] = sri[2][j];
        sri[id1][j] = sri[ib][j];
        sri[id2][j] = sri[ib][j];
        srj[1][j] = srj[2][j];
        srj[0][j] = srj[2][j];
        srj[id1][j] = srj[ib][j];
        srj[id2][j] = srj[ib][j];
    }

    /*BC_Periodic(ib, jb, id1, id2, jd1, jd2, bind, sbind, ebind, pbind, spbind, epbind, 1, pbvar);

    for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            sri[i][j]=pbvar[1][i][j];
        }
    }*/

    /*BC_Periodic(ib, jb, id1, id2, jd1, jd2, bind, sbind, ebind, pbind, spbind, epbind, 1, pbvar);

    for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            srj[i][j]=pbvar[1][i][j];
        }
    }*/

    /****for global time-step*/

    /*std::ofstream time_all("time_all.dat");

     for(int j=2;j<=jb;j++){
            for(int i=2;i<=ib;i++){

                tsmin = std::min(tsmin, tstep[i][j]);
        //time_all<<tstep[i][j]<<std::endl;

            }
        }*/

    if(epsirs>0.0){
        cflrat = sqrt(1.0+4.0*epsirs);
        for (int j = 2; j <= jb; j++) {
            for (int i = 2; i <= ib; i++) {
                reval = srj[i][j]/sri[i][j];
                reval1 = 1.0/reval;
                ex = 0.25*(pow((cflrat/(1.0+0.125*reval)),2) - 1.0);
                ey = 0.25*(pow((cflrat/(1.0+0.125*reval1)),2) - 1.0);
                ex = std::min(epsirs,ex);
                ey = std::min(epsirs,ey);
                epsij[0][i][j] = std::max(0.0, ex);
                epsij[1][i][j] = std::max(0.0, ey);
            }
        }
    }


    tsmin = 1e32;
    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= ib; i++) {
            // tsmin = std::min(tsmin, tstep[i][j]);
            if (tstep[i][j] < tsmin) tsmin = tstep[i][j];
        }
    }

    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= ib; i++) {
            tstep[i][j] = tsmin;
        }
    }

    /*for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            std::cout<<"the computed time step values are="<<tstep[i][j]<<std::endl;
        }
    }*/
}


/**********************************************************************************************************************/
////////////////////////////////////////Time Step Compuatation for Non-Dimensional Equations////////////////////////////
/**********************************************************************************************************************/
void Time_Step_NS(int ib, int id1, int id2, int jb, int jd1, int jd2, int nconv, int ndepv, double ***&cv,
                     double ***&dv, double **&area, double ***&si, double ***&sj, double **&tstep, double **&sri,
                     double **&srj, double ***&epsij) {

    double rh, u, v, sx, sy, nx, ny, ds, vc, cs, tsmin, f1, f2, fac, dtv, srvi, srvj, reval, reval1, cflrat, ex, ey;
    const double cfac=2.0;

    time_step == "global";

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            tstep[i][j] = 0.0;
        }
    }

    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= ib; i++) {
            rh = 1.0/dv[0][i][j];
            u  = dv[1][i][j];
            v  = dv[2][i][j];

            //i-direction
            sx = 0.5 * (si[0][i][j] + si[0][i + 1][j]);
            sy = 0.5 * (si[1][i][j] + si[1][i + 1][j]);

            ds = sqrt(sx * sx + sy * sy);

            nx = sx / ds;
            ny = sx / ds;
            vc = u * nx + v * ny;
            sri[i][j] = (fabs(vc) + dv[5][i][j])*ds;

            //f1 = (4.0/3.0);
            f1 = (4.0/3.0)*rh;
            //f2 = Gamma/Pr;
            f2 = Gamma*rh;
            fac = std::max(f1, f2);
            //dtv = fac*(dv[6][i][j])*rh;
            dtv = fac*(dv[6][i][j]/Pr);
            srvi = (cfac*dtv*ds*ds)/area[i][j];

            //j-direction
            sx = 0.5 * (sj[0][i][j] + sj[0][i][j + 1]);
            sy = 0.5 * (sj[1][i][j] + sj[1][i][j + 1]);

            ds = sqrt(sx * sx + sy * sy);

            nx = sx / ds;
            ny = sx / ds;
            vc = u * nx + v * ny;
            srj[i][j] = (fabs(vc) + dv[5][i][j])*ds;

            //f1 = (4.0/3.0);
            f1 = (4.0/3.0)*rh;
            //f2 = Gamma/Pr;
            f2 = Gamma*rh;
            fac = std::max(f1, f2);
            //dtv = fac*(dv[6][i][j])*rh;
            dtv = fac*(dv[6][i][j]/Pr);
            srvj = (cfac*dtv*ds*ds)/area[i][j];

            tstep[i][j] = cfl*area[i][j] / (sri[i][j] + srj[i][j]+ (srvi + srvj));


            if(tstep[i][j]==0 or tstep[i][j]<0.0){
                std::cout<<"values inside time_step are:"<<"srvi="<<tstep[i][j]<<std::endl;
                exit(0);
            }

            if (tstep[i][j] < 0 or tstep[i][j] != tstep[i][j]) {
                std::cout << "t became negative at i=" << i << "\t" << "and at j=" << j << std::endl;
                std::cout << rh << "\t" << u << "\t" << v << "\t" << "\t" << cs << "\t" << sri[i][j] << "\t"
                          << srj[i][j] << "\t" << area[i][j] << std::endl;
                exit(0);
            }
        }
    }

    //std::cout<<"values inside time_step are:"<<"srvi="<<srvj<<std::endl;
    /****values of sri and srj at dummy nodes****/
    for (int i = 2; i <= ib; i++) {
        sri[i][1] = sri[i][2];
        sri[i][0] = sri[i][2];
        sri[i][jd1] = sri[i][jb];
        sri[i][jd2] = sri[i][jb];
        srj[i][1] = srj[i][2];
        srj[i][0] = srj[i][2];
        srj[i][jd1] = srj[i][jb];
        srj[i][jd2] = srj[i][jb];
    }

    for (int j = 2; j <= jb; j++) {
        sri[1][j] = sri[2][j];
        sri[0][j] = sri[2][j];
        sri[id1][j] = sri[ib][j];
        sri[id2][j] = sri[ib][j];
        srj[1][j] = srj[2][j];
        srj[0][j] = srj[2][j];
        srj[id1][j] = srj[ib][j];
        srj[id2][j] = srj[ib][j];
    }

    /*BC_Periodic(ib, jb, id1, id2, jd1, jd2, bind, sbind, ebind, pbind, spbind, epbind, 1, pbvar);

    for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            sri[i][j]=pbvar[1][i][j];
        }
    }*/

    /*BC_Periodic(ib, jb, id1, id2, jd1, jd2, bind, sbind, ebind, pbind, spbind, epbind, 1, pbvar);

    for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            srj[i][j]=pbvar[1][i][j];
        }
    }*/

    /****for global time-step*/

    /*std::ofstream time_all("time_all.dat");

     for(int j=2;j<=jb;j++){
            for(int i=2;i<=ib;i++){

                tsmin = std::min(tsmin, tstep[i][j]);
        //time_all<<tstep[i][j]<<std::endl;

            }
        }*/

    if(epsirs>0.0){
        cflrat = sqrt(1.0+4.0*epsirs);
        for (int j = 2; j <= jb; j++) {
            for (int i = 2; i <= ib; i++) {
                reval = srj[i][j]/sri[i][j];
                reval1 = 1.0/reval;
                ex = 0.25*(pow((cflrat/(1.0+0.125*reval)),2) - 1.0);
                ey = 0.25*(pow((cflrat/(1.0+0.125*reval1)),2) - 1.0);
                ex = std::min(epsirs,ex);
                ey = std::min(epsirs,ey);
                epsij[0][i][j] = std::max(0.0, ex);
                epsij[1][i][j] = std::max(0.0, ey);
            }
        }
    }


    tsmin = 1e32;
    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= ib; i++) {
            // tsmin = std::min(tsmin, tstep[i][j]);
            if (tstep[i][j] < tsmin) tsmin = tstep[i][j];
        }
    }

    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= ib; i++) {
            tstep[i][j] = tsmin;
        }
    }

    /*for(int j=2;j<=jb;j++){
        for(int i=2;i<=ib;i++){
            std::cout<<"the computed time step values are="<<tstep[i][j]<<std::endl;
        }
    }*/
}

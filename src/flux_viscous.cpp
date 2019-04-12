//
// Created by Maruthi on 16-02-2018.
//

#include "global_declarations.h"
#include "basic_functions.h"
void
Flux_Viscous(int nconv, int ndvar, double ***&cv, double ***&dv,
               double **&u, double **&v, double **&area, double ***&si, double ***&sj, double ***&gradfi, double ***&gradfj) {

    double two_by_3=2.0/3.0;
    double uav, vav, mav, kav, tauxx, tauyy, tauxy, phix, phiy, *fv;
    int im1, jm1;
    fv = new double[3];
    /**i-direction**/
    for (int j = 2; j <= jb; ++j) {
        for (int i = 2; i <= id1; ++i) {

            im1 = i-1;
            uav = 0.5*(dv[1][im1][j]+dv[1][i][j]);
            vav = 0.5*(dv[2][im1][j]+dv[2][i][j]);
            mav = 0.5*(dv[6][im1][j]+dv[6][i][j]);
            kav = 0.5*(dv[7][im1][j]+dv[7][i][j]);
            tauxx = two_by_3*mav*(2.0*gradfi[0][i][j]-gradfi[3][i][j]);
            tauyy = two_by_3*mav*(2.0*gradfi[3][i][j]-gradfi[0][i][j]);
            tauxy = mav*(gradfi[1][i][j]+gradfi[2][i][j]);
            phix = uav*tauxx + vav*tauxy + kav*gradfi[4][i][j];
            phiy = uav*tauxy + vav*tauyy + kav*gradfi[5][i][j];
            fv[0] = si[0][i][j]*tauxx + si[1][i][j]*tauxy;
            fv[1] = si[0][i][j]*tauxy + si[1][i][j]*tauyy;
            fv[2] = si[0][i][j]*phix + si[1][i][j]*phiy;

            /**dissipation term**/
            diss[1][i][j] = diss[1][i][j] + fv[0];
            diss[2][i][j] = diss[2][i][j] + fv[1];
            diss[3][i][j] = diss[3][i][j] + fv[2];

            diss[1][im1][j] = diss[1][im1][j] - fv[0];
            diss[2][im1][j] = diss[2][im1][j] - fv[1];
            diss[3][im1][j] = diss[3][im1][j] - fv[2];

            //std::cout<<" GradFaceI="<<i<<"\t"<<j<<"\t"<<diss[1][i][j]<<"\t"<<diss[2][i][j]<<"\t"<<diss[3][i][j]<<std::endl;
            //std::cout<<" valus of mav/mu and k are="<<i<<"\t"<<j<<"\t"<<mav<<"\t"<<kav<<std::endl;
        }
    }

/**j-direction**/
    for (int i = 2; i <= ib; ++i) {
        for (int j = 2; j <= jd1; ++j) {
            jm1 = j-1;
            uav = 0.5*(dv[1][i][jm1]+dv[1][i][j]);
            vav = 0.5*(dv[2][i][jm1]+dv[2][i][j]);
            mav = 0.5*(dv[6][i][jm1]+dv[6][i][j]);
            kav = 0.5*(dv[7][i][jm1]+dv[7][i][j]);
            tauxx = two_by_3*mav*(2.0*gradfj[0][i][j]-gradfj[3][i][j]);
            tauyy = two_by_3*mav*(2.0*gradfj[3][i][j]-gradfj[0][i][j]);
            tauxy = mav*(gradfj[1][i][j]+gradfj[2][i][j]);
            phix = uav*tauxx + vav*tauxy + kav*gradfj[4][i][j];
            phiy = uav*tauxy + vav*tauyy + kav*gradfj[5][i][j];
            fv[0] = sj[0][i][j]*tauxx + sj[1][i][j]*tauxy;
            fv[1] = sj[0][i][j]*tauxy + sj[1][i][j]*tauyy;
            fv[2] = sj[0][i][j]*phix  + sj[1][i][j]*phiy;

            /**dissipation term**/
            diss[1][i][j] = diss[1][i][j] + fv[0];
            diss[2][i][j] = diss[2][i][j] + fv[1];
            diss[3][i][j] = diss[3][i][j] + fv[2];

            diss[1][i][jm1] = diss[1][i][jm1] - fv[0];
            diss[2][i][jm1] = diss[2][i][jm1] - fv[1];
            diss[3][i][jm1] = diss[3][i][jm1] - fv[2];

            //std::cout<<" GradFaceI="<<i<<"\t"<<j<<"\t"<<diss[1][i][j]<<"\t"<<diss[2][i][j]<<"\t"<<diss[3][i][j]<<std::endl;
            //std::cout<<" valus of mav/mu and k are="<<i<<"\t"<<j<<"\t"<<mav<<"\t"<<kav<<std::endl;
        }
    }
}
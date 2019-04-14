#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

void
Gradient_FaceI(int nconv, int ndvar, double ***&cv, double ***&dv,
               double **&u, double **&v, double **&area, double ***&si, double ***&sj, double ***&gradfi) {

    //computes gradients of u, v, Temperature


    double sx, sy, uav, vav, tav, rvol;
    double *fgx, *fgy;

    fgx = new double[3];
    fgy = new double[3];

    /****right face of the auxiliary control volume****/
    for (int j = 2; j <= jd1; j++) {
        for (int i = 2; i <= ib; i++) {

            sx = 0.5 * (si[0][i][j] + si[0][i + 1][j]);
            sy = 0.5 * (si[1][i][j] + si[1][i + 1][j]);
            fgx[0] = dv[1][i][j] * sx;
            fgx[1] = dv[2][i][j] * sx;
            fgx[2] = dv[4][i][j] * sx;
            fgy[0] = dv[1][i][j] * sy;
            fgy[1] = dv[2][i][j] * sy;
            fgy[2] = dv[4][i][j] * sy;

            gradfi[0][i][j] = gradfi[0][i][j] - fgx[0];
            gradfi[1][i][j] = gradfi[1][i][j] - fgy[0];
            gradfi[2][i][j] = gradfi[2][i][j] - fgx[1];
            gradfi[3][i][j] = gradfi[3][i][j] - fgy[1];
            gradfi[4][i][j] = gradfi[4][i][j] - fgx[2];
            gradfi[5][i][j] = gradfi[5][i][j] - fgy[2];

            gradfi[0][i + 1][j] = gradfi[0][i + 1][j] + fgx[0];
            gradfi[1][i + 1][j] = gradfi[1][i + 1][j] + fgy[0];
            gradfi[2][i + 1][j] = gradfi[2][i + 1][j] + fgx[1];
            gradfi[3][i + 1][j] = gradfi[3][i + 1][j] + fgy[1];
            gradfi[4][i + 1][j] = gradfi[4][i + 1][j] + fgx[2];
            gradfi[5][i + 1][j] = gradfi[5][i + 1][j] + fgy[2];

            /****bottom face of auxiliary cv****/
            sx = 0.5 * (sj[0][i][j] + sj[0][i - 1][j]);
            sy = 0.5 * (sj[1][i][j] + sj[1][i - 1][j]);
            uav = 0.25 * (dv[1][i][j] + dv[1][i-1][j] + dv[1][i][j-1] + dv[1][i-1][j-1]);
            vav = 0.25 * (dv[2][i][j] + dv[2][i-1][j] + dv[2][i][j-1] + dv[2][i-1][j-1]);
            tav = 0.25 * (dv[4][i][j] + dv[4][i-1][j] + dv[4][i][j-1] + dv[4][i-1][j-1]);

            fgx[0] = uav * sx;
            fgx[1] = vav * sx;
            fgx[2] = tav * sx;

            fgy[0] = uav * sy;
            fgy[1] = vav * sy;
            fgy[2] = tav * sy;

            gradfi[0][i][j] = gradfi[0][i][j] + fgx[0];
            gradfi[1][i][j] = gradfi[1][i][j] + fgy[0];
            gradfi[2][i][j] = gradfi[2][i][j] + fgx[1];
            gradfi[3][i][j] = gradfi[3][i][j] + fgy[1];
            gradfi[4][i][j] = gradfi[4][i][j] + fgx[2];
            gradfi[5][i][j] = gradfi[5][i][j] + fgy[2];

            gradfi[0][i][j - 1] = gradfi[0][i][j - 1] - fgx[0];
            gradfi[1][i][j - 1] = gradfi[1][i][j - 1] - fgy[0];
            gradfi[2][i][j - 1] = gradfi[2][i][j - 1] - fgx[1];
            gradfi[3][i][j - 1] = gradfi[3][i][j - 1] - fgy[1];
            gradfi[4][i][j - 1] = gradfi[4][i][j - 1] - fgx[2];
            gradfi[5][i][j - 1] = gradfi[5][i][j - 1] - fgy[2];
        }
            //at boundary i=2
            sx = si[0][2][j];
            sy = si[1][2][j];

            fgx[0] = 0.5 * (dv[1][1][j] + dv[1][2][j]) * sx;
            fgx[1] = 0.5 * (dv[2][1][j] + dv[2][2][j]) * sx;
            fgx[2] = 0.5 * (dv[4][1][j] + dv[4][2][j]) * sx;
            fgy[0] = 0.5 * (dv[1][1][j] + dv[1][2][j]) * sy;
            fgy[1] = 0.5 * (dv[2][1][j] + dv[2][2][j]) * sy;
            fgy[2] = 0.5 * (dv[4][1][j] + dv[4][2][j]) * sy;

            gradfi[0][2][j] = gradfi[0][2][j] + fgx[0];
            gradfi[1][2][j] = gradfi[1][2][j] + fgy[0];
            gradfi[2][2][j] = gradfi[2][2][j] + fgx[1];
            gradfi[3][2][j] = gradfi[3][2][j] + fgy[1];
            gradfi[4][2][j] = gradfi[4][2][j] + fgx[2];
            gradfi[5][2][j] = gradfi[5][2][j] + fgy[2];

            //at boundary i=id1
            sx = si[0][id1][j];
            sy = si[1][id1][j];

            fgx[0] = 0.5 * (dv[1][ib][j] + dv[1][id1][j]) * sx;
            fgx[1] = 0.5 * (dv[2][ib][j] + dv[2][id1][j]) * sx;
            fgx[2] = 0.5 * (dv[4][ib][j] + dv[4][id1][j]) * sx;
            fgy[0] = 0.5 * (dv[1][ib][j] + dv[1][id1][j]) * sy;
            fgy[1] = 0.5 * (dv[2][ib][j] + dv[2][id1][j]) * sy;
            fgy[2] = 0.5 * (dv[4][ib][j] + dv[4][id1][j]) * sy;

            gradfi[0][id1][j] = gradfi[0][id1][j] - fgx[0];
            gradfi[1][id1][j] = gradfi[1][id1][j] - fgy[0];
            gradfi[2][id1][j] = gradfi[2][id1][j] - fgx[1];
            gradfi[3][id1][j] = gradfi[3][id1][j] - fgy[1];
            gradfi[4][id1][j] = gradfi[4][id1][j] - fgx[2];
            gradfi[5][id1][j] = gradfi[5][id1][j] - fgy[2];

            sx = 0.5 * sj[0][ib][j];
            sy = 0.5 * sj[1][ib][j];
            uav = 0.25 * (dv[1][ib][j] + dv[1][id1][j] + dv[1][ib][j-1] + dv[1][id1][j-1]);
            vav = 0.25 * (dv[2][ib][j] + dv[2][id1][j] + dv[2][ib][j-1] + dv[2][id1][j-1]);
            tav = 0.25 * (dv[4][ib][j] + dv[4][id1][j] + dv[4][ib][j-1] + dv[4][id1][j-1]);

            fgx[0] = uav * sx;
            fgx[1] = vav * sx;
            fgx[2] = tav * sx;
            fgy[0] = uav * sy;
            fgy[1] = vav * sy;
            fgy[2] = tav * sy;

            gradfi[0][id1][j] = gradfi[0][id1][j] + fgx[0];
            gradfi[1][id1][j] = gradfi[1][id1][j] + fgy[0];
            gradfi[2][id1][j] = gradfi[2][id1][j] + fgx[1];
            gradfi[3][id1][j] = gradfi[3][id1][j] + fgy[1];
            gradfi[4][id1][j] = gradfi[4][id1][j] + fgx[2];
            gradfi[5][id1][j] = gradfi[5][id1][j] + fgy[2];

            gradfi[0][id1][j - 1] = gradfi[0][id1][j - 1] - fgx[0];
            gradfi[1][id1][j - 1] = gradfi[1][id1][j - 1] - fgy[0];
            gradfi[2][id1][j - 1] = gradfi[2][id1][j - 1] - fgx[1];
            gradfi[3][id1][j - 1] = gradfi[3][id1][j - 1] - fgy[1];
            gradfi[4][id1][j - 1] = gradfi[4][id1][j - 1] - fgx[2];
            gradfi[5][id1][j - 1] = gradfi[5][id1][j - 1] - fgy[2];

        /*std::cout<<"Test values are at j="<<j<<"\t"<<gradfi[0][2][j]<<"\t"<<gradfi[1][2][j]<<"\t"
                 <<gradfi[2][2][j]<<"\t"<<gradfi[3][2][j]<<"\t"<<gradfi[4][2][j]<<"\t"<<gradfi[5][2][j]<<std::endl;*/
    }
    /**divide with control volume**/
    for (int j = 2; j <= jb; j++) {
        for (int i = 3; i <= ib; i++) {
            rvol = 2.0 / (area[i][j] + area[i - 1][j]);
            gradfi[0][i][j] = gradfi[0][i][j] * rvol;
            gradfi[1][i][j] = gradfi[1][i][j] * rvol;
            gradfi[2][i][j] = gradfi[2][i][j] * rvol;
            gradfi[3][i][j] = gradfi[3][i][j] * rvol;
            gradfi[4][i][j] = gradfi[4][i][j] * rvol;
            gradfi[5][i][j] = gradfi[5][i][j] * rvol;

            /*std::cout<<" GradFaceI="<<i<<"\t"<<j<<"\t"<<gradfi[0][i][j]<<"\t"<<gradfi[1][i][j]
                     <<"\t"<<gradfi[2][i][j]<<"\t"<<gradfi[3][i][j]<<"\t"<<gradfi[4][i][j]<<"\t"<<gradfi[5][i][j]<<std::endl;*/
        }
        rvol = 2.0 / (area[2][j]);
        gradfi[0][2][j] = gradfi[0][2][j] * rvol;
        gradfi[1][2][j] = gradfi[1][2][j] * rvol;
        gradfi[2][2][j] = gradfi[2][2][j] * rvol;
        gradfi[3][2][j] = gradfi[3][2][j] * rvol;
        gradfi[4][2][j] = gradfi[4][2][j] * rvol;
        gradfi[5][2][j] = gradfi[5][2][j] * rvol;

        rvol = 2.0 / (area[ib][j]);
        gradfi[0][id1][j] = gradfi[0][id1][j] * rvol;
        gradfi[1][id1][j] = gradfi[1][id1][j] * rvol;
        gradfi[2][id1][j] = gradfi[2][id1][j] * rvol;
        gradfi[3][id1][j] = gradfi[3][id1][j] * rvol;
        gradfi[4][id1][j] = gradfi[4][id1][j] * rvol;
        gradfi[5][id1][j] = gradfi[5][id1][j] * rvol;

        /*std::cout<<"Test values are at j="<<j<<"\t"<<gradfi[0][2][j]<<"\t"<<gradfi[1][2][j]<<"\t"
                 <<gradfi[2][2][j]<<"\t"<<gradfi[3][2][j]<<"\t"<<gradfi[4][2][j]<<"\t"<<gradfi[5][2][j]<<std::endl;*/
    }
     /*for (int j = 2; j <= jb; j++) {
         for (int i = 2; i <= ib; i++) {

             std::cout<<"the gradient values are at GradFaceI="<<i<<"\t"<<j<<"\t"<<gradfi[0][i][j]<<"\t"<<gradfi[1][i][j]
                      <<"\t"<<gradfi[2][i][j]<<"\t"<<gradfi[3][i][j]<<"\t"<<gradfi[4][i][j]<<"\t"<<gradfi[5][i][j]<<std::endl;
         }
     }*/

    delete fgx;
    delete fgy;
}

void
Gradient_FaceJ(int nconv, int ndvar, double ***&cv, double ***&dv,
               double **&u, double **&v, double **&area, double ***&si, double ***&sj, double ***&gradfj) {

    //computes gradients of u, v, Temperature

    double sx, sy, uav, vav, tav, rvol;
    double *fgx, *fgy;

    fgx = new double[3];
    fgy = new double[3];

    /****bottom face of the auxiliary control volume****/
    for (int i = 2; i <= id1; i++) {
        for (int j = 2; j <= jb; j++) {
            sx = 0.5 * (sj[0][i][j] + sj[0][i][j + 1]);
            sy = 0.5 * (sj[1][i][j] + sj[1][i][j + 1]);
            fgx[0] = dv[1][i][j] * sx;
            fgx[1] = dv[2][i][j] * sx;
            fgx[2] = dv[4][i][j] * sx;
            fgy[0] = dv[1][i][j] * sy;
            fgy[1] = dv[2][i][j] * sy;
            fgy[2] = dv[4][i][j] * sy;

            gradfj[0][i][j] = gradfj[0][i][j] - fgx[0];
            gradfj[1][i][j] = gradfj[1][i][j] - fgy[0];
            gradfj[2][i][j] = gradfj[2][i][j] - fgx[1];
            gradfj[3][i][j] = gradfj[3][i][j] - fgy[1];
            gradfj[4][i][j] = gradfj[4][i][j] - fgx[2];
            gradfj[5][i][j] = gradfj[5][i][j] - fgy[2];

            gradfj[0][i][j + 1] = gradfj[0][i][j + 1] + fgx[0];
            gradfj[1][i][j + 1] = gradfj[1][i][j + 1] + fgy[0];
            gradfj[2][i][j + 1] = gradfj[2][i][j + 1] + fgx[1];
            gradfj[3][i][j + 1] = gradfj[3][i][j + 1] + fgy[1];
            gradfj[4][i][j + 1] = gradfj[4][i][j + 1] + fgx[2];
            gradfj[5][i][j + 1] = gradfj[5][i][j + 1] + fgy[2];

            /****left face of auxiliary cv****/
            sx = 0.5 * (si[0][i][j] + si[0][i][j - 1]);
            sy = 0.5 * (si[1][i][j] + si[1][i][j - 1]);
            uav = 0.25 * (dv[1][i][j] + dv[1][i][j-1] + dv[1][i-1][j-1] + dv[1][i-1][j]);
            vav = 0.25 * (dv[2][i][j] + dv[2][i][j-1] + dv[2][i-1][j-1] + dv[2][i-1][j]);
            tav = 0.25 * (dv[4][i][j] + dv[4][i][j-1] + dv[4][i-1][j-1] + dv[4][i-1][j]);

            fgx[0] = uav * sx;
            fgx[1] = vav * sx;
            fgx[2] = tav * sx;

            fgy[0] = uav * sy;
            fgy[1] = vav * sy;
            fgy[2] = tav * sy;

            gradfj[0][i][j] = gradfj[0][i][j] + fgx[0];
            gradfj[1][i][j] = gradfj[1][i][j] + fgy[0];
            gradfj[2][i][j] = gradfj[2][i][j] + fgx[1];
            gradfj[3][i][j] = gradfj[3][i][j] + fgy[1];
            gradfj[4][i][j] = gradfj[4][i][j] + fgx[2];
            gradfj[5][i][j] = gradfj[5][i][j] + fgy[2];

            gradfj[0][i - 1][j] = gradfj[0][i - 1][j] - fgx[0];
            gradfj[1][i - 1][j] = gradfj[1][i - 1][j] - fgy[0];
            gradfj[2][i - 1][j] = gradfj[2][i - 1][j] - fgx[1];
            gradfj[3][i - 1][j] = gradfj[3][i - 1][j] - fgy[1];
            gradfj[4][i - 1][j] = gradfj[4][i - 1][j] - fgx[2];
            gradfj[5][i - 1][j] = gradfj[5][i - 1][j] - fgy[2];
        }
            //at boundary j=2
            sx = sj[0][i][2];
            sy = sj[1][i][2];

            fgx[0] = 0.5 * (dv[1][i][1] + dv[1][i][2]) * sx;
            fgx[1] = 0.5 * (dv[2][i][1] + dv[2][i][2]) * sx;
            fgx[2] = 0.5 * (dv[4][i][1] + dv[4][i][2]) * sx;
            fgy[0] = 0.5 * (dv[1][i][1] + dv[1][i][2]) * sy;
            fgy[1] = 0.5 * (dv[2][i][1] + dv[2][i][2]) * sy;
            fgy[2] = 0.5 * (dv[4][i][1] + dv[4][i][2]) * sy;

            gradfj[0][i][2] = gradfj[0][i][2] + fgx[0];
            gradfj[1][i][2] = gradfj[1][i][2] + fgy[0];
            gradfj[2][i][2] = gradfj[2][i][2] + fgx[1];
            gradfj[3][i][2] = gradfj[3][i][2] + fgy[1];
            gradfj[4][i][2] = gradfj[4][i][2] + fgx[2];
            gradfj[5][i][2] = gradfj[5][i][2] + fgy[2];

            //at boundary j=jd1
            sx = sj[0][i][jd1];
            sy = sj[1][i][jd1];

            fgx[0] = 0.5 * (dv[1][i][jb] + dv[1][i][jd1]) * sx;
            fgx[1] = 0.5 * (dv[2][i][jb] + dv[2][i][jd1]) * sx;
            fgx[2] = 0.5 * (dv[4][i][jb] + dv[4][i][jd1]) * sx;
            fgy[0] = 0.5 * (dv[1][i][jb] + dv[1][i][jd1]) * sy;
            fgy[1] = 0.5 * (dv[2][i][jb] + dv[2][i][jd1]) * sy;
            fgy[2] = 0.5 * (dv[4][i][jb] + dv[4][i][jd1]) * sy;

            gradfj[0][i][jd1] = gradfj[0][i][jd1] - fgx[0];
            gradfj[1][i][jd1] = gradfj[1][i][jd1] - fgy[0];
            gradfj[2][i][jd1] = gradfj[2][i][jd1] - fgx[1];
            gradfj[3][i][jd1] = gradfj[3][i][jd1] - fgy[1];
            gradfj[4][i][jd1] = gradfj[4][i][jd1] - fgx[2];
            gradfj[5][i][jd1] = gradfj[5][i][jd1] - fgy[2];

            sx = 0.5 * si[0][i][jb];
            sy = 0.5 * si[1][i][jb];
            uav = 0.25 * (dv[1][i][jb] + dv[1][i][jd1] + dv[1][i-1][jb] + dv[1][i-1][jd1]);
            vav = 0.25 * (dv[2][i][jb] + dv[2][i][jd1] + dv[2][i-1][jb] + dv[2][i-1][jd1]);
            tav = 0.25 * (dv[4][i][jb] + dv[4][i][jd1] + dv[4][i-1][jb] + dv[4][i-1][jd1]);

            fgx[0] = uav * sx;
            fgx[1] = vav * sx;
            fgx[2] = tav * sx;
            fgy[0] = uav * sy;
            fgy[1] = vav * sy;
            fgy[2] = tav * sy;

            gradfj[0][i][jd1] = gradfj[0][i][jd1] + fgx[0];
            gradfj[1][i][jd1] = gradfj[1][i][jd1] + fgy[0];
            gradfj[2][i][jd1] = gradfj[2][i][jd1] + fgx[1];
            gradfj[3][i][jd1] = gradfj[3][i][jd1] + fgy[1];
            gradfj[4][i][jd1] = gradfj[4][i][jd1] + fgx[2];
            gradfj[5][i][jd1] = gradfj[5][i][jd1] + fgy[2];

            gradfj[0][i - 1][jd1] = gradfj[0][i - 1][jd1] - fgx[0];
            gradfj[1][i - 1][jd1] = gradfj[1][i - 1][jd1] - fgy[0];
            gradfj[2][i - 1][jd1] = gradfj[2][i - 1][jd1] - fgx[1];
            gradfj[3][i - 1][jd1] = gradfj[3][i - 1][jd1] - fgy[1];
            gradfj[4][i - 1][jd1] = gradfj[4][i - 1][jd1] - fgx[2];
            gradfj[5][i - 1][jd1] = gradfj[5][i - 1][jd1] - fgy[2];
        }

    for (int i = 2; i <= ib; i++) {
        for (int j = 3; j <= jb; j++) {
            rvol = 2.0 / (area[i][j] + area[i][j - 1]);
            gradfj[0][i][j] = gradfj[0][i][j] * rvol;
            gradfj[1][i][j] = gradfj[1][i][j] * rvol;
            gradfj[2][i][j] = gradfj[2][i][j] * rvol;
            gradfj[3][i][j] = gradfj[3][i][j] * rvol;
            gradfj[4][i][j] = gradfj[4][i][j] * rvol;
            gradfj[5][i][j] = gradfj[5][i][j] * rvol;

            /*std::cout << "    GradFaceJ=" << i << "\t" << j << "\t" << gradfj[0][i][j] << "\t" << gradfj[1][i][j]
                      << "\t" << gradfj[2][i][j] << "\t" << gradfj[3][i][j] << "\t" << gradfj[4][i][j] << "\t"
                      << gradfj[5][i][j] << std::endl;*/
        }

        rvol = 2.0 / (area[i][2]);
        gradfj[0][i][2] = gradfj[0][i][2] * rvol;
        gradfj[1][i][2] = gradfj[1][i][2] * rvol;
        gradfj[2][i][2] = gradfj[2][i][2] * rvol;
        gradfj[3][i][2] = gradfj[3][i][2] * rvol;
        gradfj[4][i][2] = gradfj[4][i][2] * rvol;
        gradfj[5][i][2] = gradfj[5][i][2] * rvol;

        rvol = 2.0 / (area[i][jb]);
        gradfj[0][i][jd1] = gradfj[0][i][jd1] * rvol;
        gradfj[1][i][jd1] = gradfj[1][i][jd1] * rvol;
        gradfj[2][i][jd1] = gradfj[2][i][jd1] * rvol;
        gradfj[3][i][jd1] = gradfj[3][i][jd1] * rvol;
        gradfj[4][i][jd1] = gradfj[4][i][jd1] * rvol;
        gradfj[5][i][jd1] = gradfj[5][i][jd1] * rvol;
    }
    /*for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= ib; i++) {

            std::cout<<"the gradient values are at GradFaceJ="<<i<<"\t"<<j<<"\t"<<gradfj[0][i][j]<<"\t"<<gradfj[1][i][j]
                     <<"\t"<<gradfj[2][i][j]<<"\t"<<gradfj[3][i][j]<<"\t"<<gradfj[4][i][j]<<"\t"<<gradfj[5][i][j]<<std::endl;
        }
    }*/
    delete fgx;
    delete fgy;
}
#include "global_declarations.h"
#include "basic_functions.h"

void Initialize_Gradients(int id2, int jd2, int nconv, double ***&cv, double **&u, double **&v, double ***&gradfi,
                          double ***&gradfj){

/*    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            u[i][j] = cv[1][i][j] / cv[0][i][j];
            v[i][j] = cv[2][i][j] / cv[0][i][j];

            //std::cout<<"the u and v values in Initialize Grad="<<i<<"\t"<<j<<"\t"<<u[i][j]<<"\t"<<v[i][j]<<std::endl;
        }
        //std::cout<<"Test values are at j="<<j<<"\t"<<"j="<<j<<"\t"<<u[2][j]<<"\t"<<v[2][j]<<std::endl;
    }*/

    /*for (int j = 0; j <= jd1; j++) {
        for (int i = 0; i <= id1; i++) {
            std::cout<<"the u and v values in Initialize Grad="<<i<<"\t"<<j<<"\t"<<u[i][j]<<"\t"<<v[i][j]<<std::endl;
        }
    }*/

    //Initialize gradients to zero values: used later for viscous grad cal
    for (int n = 0; n < 6; n++) {
        for (int j = 0; j <= jd2; j++) {
            for (int i = 0; i <= id2; i++) {
                gradfi[n][i][j] = 0;
                gradfj[n][i][j] = 0;
            }
        }
    }
}

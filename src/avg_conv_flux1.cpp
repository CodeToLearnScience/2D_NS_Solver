#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

void Avg_Conv_Flux1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double **&area,
                    double ***&si, double ***&sj, double ***&diss, double ***&rhs) {

    double *fc;
    double rhoa=0.0, rhoua=0.0, rhova=0.0, rhoea=0.0, pa=0.0, vcont=0.0;

    fc = new double[nconv];

    for(int i=0; i<nconv; i++){
        fc[i] = 0.0;
    }

    for(int k=0; k<nconv; k++){
        for (int j = 2; j <= jb; j++) {
            for (int i = 2; i <= ib; i++) {

                rhs[k][i][j] = -diss[k][i][j];
            }
        }
    }

    for (int j = 2; j <= jb; j++) {
        //flux in i direction except at the boundaries
        for (int i = 2; i <= id1; i++) {

            rhoa =  0.5 * (dv[0][i-1][j] + dv[0][i][j]);
            rhoua = 0.5 * (dv[1][i-1][j] * dv[0][i-1][j] + dv[1][i][j] * dv[0][i][j]);
            rhova = 0.5 * (dv[2][i-1][j] * dv[0][i-1][j] + dv[2][i][j] * dv[0][i][j]);
            rhoea = 0.5 * (dv[3][i-1][j]/(Gamma-1.0)+0.5*dv[0][i-1][j]*(dv[1][i-1][j]*dv[1][i-1][j]+dv[2][i-1][j] *
                    dv[2][i-1][j]) + dv[3][i][j]/(Gamma-1.0)+0.5*dv[0][i][j]*(dv[1][i][j]*dv[1][i][j]+dv[2][i][j]*
                    dv[2][i][j]));
            pa =    0.5 * (dv[3][i-1][j] + dv[3][i][j]);

            if(pa<0){
                std::cout<<"the avg val of p in i-dir in avg_con_flux1 is="<<pa<<std::endl;
                exit(0);
            }

            vcont = (rhoua * si[0][i][j] + rhova * si[1][i][j]) / rhoa;

            fc[0] = vcont * rhoa;
            fc[1] = vcont * rhoua + pa * si[0][i][j];
            fc[2] = vcont * rhova + pa * si[1][i][j];
            fc[3] = vcont * (rhoea + pa);

            for(int k=0; k<nconv; k++){
                //rhs[k][i][j] = rhs[k][i][j] + fc[k];
                //rhs[k][i - 1][j] = rhs[k][i - 1][j] - fc[k];
                rhs[k][i][j] = rhs[k][i][j] + iAvgFlux[k][i][j];
                rhs[k][i - 1][j] = rhs[k][i - 1][j] - iAvgFlux[k][i][j];
            }

            //}
            /*std::cout<<"avg flux:"<<i<<"\t"<<j<<"\t"<<rhs[0][i][j]<<"\t"<<rhs[1][i][j]<<"\t"<<rhs[2][i][j]<<"\t"
                     <<rhs[3][i][j]<<std::endl;*/
        }
    }
        for (int j = 2; j <= jd1; j++) {
        //if (j > 1) {
            //flux in j-direction except at boundaries
            for (int i = 2; i <= ib; i++) {

                rhoa =  0.5 * (dv[0][i][j-1] + dv[0][i][j]);
                rhoua = 0.5 * (dv[1][i][j-1]*dv[0][i][j-1] + dv[1][i][j]*dv[0][i][j]);
                rhova = 0.5 * (dv[2][i][j-1]*dv[0][i][j-1] + dv[2][i][j]*dv[0][i][j]);
                rhoea = 0.5 * (dv[3][i][j-1]/(Gamma-1.0)+0.5*dv[0][i][j-1]*(dv[1][i][j-1]*dv[1][i][j-1]+dv[2][i][j-1] *
                        dv[2][i][j-1]) + dv[3][i][j]/(Gamma-1.0)+0.5*dv[0][i][j]*(dv[1][i][j]*dv[1][i][j]+dv[2][i][j] *
                        dv[2][i][j]));
                pa =    0.5 * (dv[3][i][j-1] + dv[3][i][j]);
                if(pa<0){
                    std::cout<<"the avg val of p in i-dir in avg_con_flux1 is="<<pa<<std::endl;
                    exit(0);
                }
                vcont = (rhoua * sj[0][i][j] + rhova * sj[1][i][j]) / rhoa;

                fc[0] = vcont * rhoa;
                fc[1] = vcont * rhoua + pa * sj[0][i][j];
                fc[2] = vcont * rhova + pa * sj[1][i][j];
                fc[3] = vcont * (rhoea + pa);

                for(int k=0; k<nconv; k++){
                    //rhs[k][i][j] = rhs[k][i][j] + fc[k];
                    //rhs[k][i][j - 1] = rhs[k][i][j - 1] - fc[k];

                    rhs[k][i][j] = rhs[k][i][j] + jAvgFlux[k][i][j];
                    rhs[k][i][j - 1] = rhs[k][i][j - 1] - jAvgFlux[k][i][j];
                }
            }
        }

    /*Flux_Symmetry(1, 1, 2, 9, cv, dv, si, sj, rhs);
    Flux_Boundary(2, id1, 2, jb, cv, dv, si, sj, rhs);
    Flux_Boundary(3, jd1, 2, ib, cv, dv, si, sj, rhs);
    Flux_Boundary(4, 1, 2, jb, cv, dv, si, sj, rhs);
    Flux_ViscWall(1, 1, 10, ib, cv, dv, si, sj, rhs);*/
    delete[] fc;
}


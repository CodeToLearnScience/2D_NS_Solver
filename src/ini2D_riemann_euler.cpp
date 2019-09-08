//
// Created by maruthinh on 27/4/19.
//

#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"
void Init2D_RiemannEuler(int id2, int jd2, int interface_flag, int interface_ind, double rho_l, double u_l, double v_l,
        double p_l, double rho_r, double u_r, double v_r, double p_r, double ***&cv, double ***&dv) {

    std::cout<<"the values to the left and right of the interface in 2D Riemann proble are: rho_l="<<rho_l<<"\t"
    <<"u_l="<<u_l<<"\t"<<"v_l="<<v_l<<"\t"<<"p_l="<<p_l<<"\t"<<"rho_r="<<rho_r<<"\t"<<"u_r="<<u_r<<"\t"<<"v_r="<<v_r<<"\t"
    <<"p_r="<<p_r<<std::endl;

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            cv[0][i][j] = 0.0;
            cv[1][i][j] = 0.0;
            cv[2][i][j] = 0.0;
            cv[3][i][j] = 0.0;

            dv[0][i][j] = 0.0; //rho
            dv[1][i][j] = 0.0; //u
            dv[2][i][j] = 0.0; //v
            dv[3][i][j] = 0.0; //p
            dv[4][i][j] = 0.0; //T
            dv[5][i][j] = 0.0; //a
            dv[6][i][j] = 0.0; //Mu
            dv[7][i][j] = 0.0; //K
        }
    }

    if(interface_flag==0){
        std::cout<<"diaphragm is located along x-direction"<<std::endl;

        for (int j = 0; j <= jd2; j++) {
            for (int i = 0; i <= interface_ind; i++) {
                dv[0][i][j] = rho_l; //rho
                dv[1][i][j] = u_l; //u
                dv[2][i][j] = v_l; //v
                dv[3][i][j] = p_l; //p
                dv[4][i][j] = p_l/(rho_l*Rgas); //T
                dv[5][i][j] = sqrt(Gamma*p_l/rho_l); //a

                /*std::cout <<"I.Cs: Dep Var"<<"\t"<<i<<"\t"<<j<<"\t"<<dv[0][i][j]<<"\t"<<dv[1][i][j]<<"\t"<< dv[2][i][j]
                          <<"\t"<<dv[3][i][j]<<"\t"<<dv[4][i][j]<<"\t"<<dv[5][i][j]<<"\t"<<dv[6][i][j]<<"\t"<<dv[7][i][j]
                          << std::endl;*/
            }
        }

        for (int j = 0; j <= jd2; j++) {
            for (int i = interface_ind; i <= id2; i++) {
                dv[0][i][j] = rho_r; //rho
                dv[1][i][j] = u_r; //u
                dv[2][i][j] = v_r; //v
                dv[3][i][j] = p_r; //p
                dv[4][i][j] = p_r/(rho_r*Rgas); //T
                dv[5][i][j] = sqrt(Gamma*p_r/rho_r); //a
                /*std::cout <<"I.Cs: Dep Var"<<"\t"<<i<<"\t"<<j<<"\t"<<dv[0][i][j]<<"\t"<<dv[1][i][j]<<"\t"<< dv[2][i][j]
                          <<"\t"<<dv[3][i][j]<<"\t"<<dv[4][i][j]<<"\t"<<dv[5][i][j]<<"\t"<<dv[6][i][j]<<"\t"<<dv[7][i][j]
                          << std::endl;*/
            }
        }
    }
    else if (interface_flag==1){
        std::cout<<"diaphragm is located along y-direction"<<std::endl;

        for (int j = 0; j <= interface_ind; j++) {
            for (int i = 0; i <= id2; i++) {
                dv[0][i][j] = rho_l; //rho
                dv[1][i][j] = u_l; //u
                dv[2][i][j] = v_l; //v
                dv[3][i][j] = p_l; //p
                dv[4][i][j] = p_l/(rho_l*Rgas); //T
                dv[5][i][j] = sqrt(Gamma*p_l/rho_l); //a

                /*std::cout <<"I.Cs: Dep Var"<<"\t"<<i<<"\t"<<j<<"\t"<<dv[0][i][j]<<"\t"<<dv[1][i][j]<<"\t"<< dv[2][i][j]
                          <<"\t"<<dv[3][i][j]<<"\t"<<dv[4][i][j]<<"\t"<<dv[5][i][j]<<"\t"<<dv[6][i][j]<<"\t"<<dv[7][i][j]
                          << std::endl;*/
            }
        }

        for (int j = interface_ind+1; j <= jd2; j++) {
            for (int i = 0; i <= id2; i++) {
                dv[0][i][j] = rho_r; //rho
                dv[1][i][j] = u_r; //u
                dv[2][i][j] = v_r; //v
                dv[3][i][j] = p_r; //p
                dv[4][i][j] = p_r/(rho_r*Rgas); //T
                dv[5][i][j] = sqrt(Gamma*p_r/rho_r); //a
                /*std::cout <<"I.Cs: Dep Var"<<"\t"<<i<<"\t"<<j<<"\t"<<dv[0][i][j]<<"\t"<<dv[1][i][j]<<"\t"<< dv[2][i][j]
                          <<"\t"<<dv[3][i][j]<<"\t"<<dv[4][i][j]<<"\t"<<dv[5][i][j]<<"\t"<<dv[6][i][j]<<"\t"<<dv[7][i][j]
                          << std::endl;*/
            }
        }
    }

    else{
        std::cout<<"wrong ind for the interface location. It shouble be 0 for x and 1 for y"<<std::endl;
        exit(0);
    }
        for (int j = 0; j <= jd2; j++) {
            for (int i = 0; i <= id2; i++) {
                cv[0][i][j] = dv[0][i][j];
                cv[1][i][j] = dv[0][i][j] * dv[1][i][j];
                cv[2][i][j] = dv[0][i][j] * dv[2][i][j];
                cv[3][i][j] = dv[3][i][j]/(Gamma-1.0)+0.5*dv[0][i][j]*(dv[1][i][j]*dv[1][i][j]+dv[2][i][j]*dv[2][i][j]);

                /*std::cout<<"I.Cs: Con Var" << "\t" << i<<"\t"<<j<<"\t"<<cv[0][i][j] << "\t" << cv[1][i][j]<< "\t"
                         << cv[2][i][j]<< "\t" << cv[3][i][j] << std::endl;*/
            }
        }
}

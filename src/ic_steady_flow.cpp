#include "global_declarations.h"
#include "basic_functions.h"


void Ic_Steady_Flow(int id2, int jd2, double rho_inf, double u_inf, double v_inf, double p_inf, double ***&cv, double ***&dv) {


    Re_inf=100000, Lref=0.8, Machinf=0.15, rhoinf = 1.0, tinf=1.0, qinf = 1.0, alpha=0.0;
    Tref = 288.16;

    Rgas = 1.0/(Gamma*Machinf*Machinf);
    pinf = rhoinf*tinf*Rgas;
    Cp=1.0/((Gamma-1.0)*Machinf*Machinf);
    refvisc = Lref/Re_inf;
    alpha = alpha/rad;
    uinf = qinf*cos(alpha);
    vinf = qinf*sin(alpha);
    rhoref = rhoinf;
    velref = qinf;

    /*Cp = 1005.0;
    Rgas = ((Gamma-1.0)/Gamma)*Cp;
    rhoinf = pinf/(Rgas*tinf);
    qinf = Machinf*sqrt((Gamma-1.0)*Cp*tinf);
    alpha = alpha/rad;
    uinf = qinf*cos(alpha);
    vinf = qinf*sin(alpha);
    refrho = rhoinf;
    refvel = qinf;
    refvisc = rhoinf*qinf*Lref/Re;
    //refvisc = 1.0/Re;*/

    std::cout<<"the values are at Ic_cond:"<<rhoinf<<"\t"<<Rgas<<"\t"<<qinf<<"\t"<<uinf<<"\t"<<pinf<<"\t"<<refvisc<<
             std::endl;
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

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            dv[0][i][j] = rhoinf; //rho
            dv[1][i][j] = uinf; //u
            dv[2][i][j] = vinf; //v
            dv[3][i][j] = pinf; //p
            dv[4][i][j] = tinf; //T
            dv[5][i][j] = sqrt(Gamma*pinf/rhoinf); //a
            //dv[6][i][j] = ((C0+C1)/(C1+dv[4][i][j]))*pow((dv[4][i][j]/C0),1.5)*refvisc; //Mu
            dv[6][i][j] = ((1.0+C1/C0)/(C1/C0+dv[4][i][j]))*pow(dv[4][i][j],1.5); //Mu
            dv[7][i][j] = dv[6][i][j]*Cp/Pr; //K

            /*std::cout <<"I.Cs: Dep Var"<<"\t"<<i<<"\t"<<j<<"\t"<<dv[0][i][j]<<"\t"<<dv[1][i][j]<<"\t"<< dv[2][i][j]
                      <<"\t"<<dv[3][i][j]<<"\t"<<dv[4][i][j]<<"\t"<<dv[5][i][j]<<"\t"<<dv[6][i][j]<<"\t"<<dv[7][i][j]
                      << std::endl;*/
        }
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


/*
void Ic_Steady_Flow(int id2, int jd2, double rho_inf, double u_inf, double v_inf, double p_inf, double ***&cv, double ***&dv) {

    //Re=100000, Lref=1.0, Machinf=0.15, rhoinf = 1.13235, tinf=297.62, qinf = 69.1687, alpha=0.0;
    Re=100000, Lref=1.0, Machinf=0.15, tinf=293.0, pinf=101325.0, alpha=0.0;
    Tref = 288.16;


    Cp = 1005.0;
    Rgas = ((Gamma-1.0)/Gamma)*Cp;
    rhoinf = pinf/(Rgas*tinf);
    qinf = Machinf*sqrt((Gamma-1.0)*Cp*tinf);
    alpha = alpha/rad;
    uinf = qinf*cos(alpha);
    vinf = qinf*sin(alpha);
    refrho = rhoinf;
    refvel = qinf;
    refvisc = rhoinf*qinf*Lref/Re;
    //refvisc = 1.0/Re;

    std::cout<<"the values are at Ic_cond:"<<rhoinf<<"\t"<<Rgas<<"\t"<<qinf<<"\t"<<uinf<<"\t"<<pinf<<"\t"<<refvisc<<
             std::endl;
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

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            dv[0][i][j] = rhoinf; //rho
            dv[1][i][j] = uinf; //u
            dv[2][i][j] = vinf; //v
            dv[3][i][j] = pinf; //p
            dv[4][i][j] = tinf; //T
            dv[5][i][j] = sqrt(Gamma*pinf/rhoinf); //a
            dv[6][i][j] = ((C0+C1)/(C1+dv[4][i][j]))*pow((dv[4][i][j]/C0),1.5)*refvisc; //Mu
            dv[7][i][j] = dv[6][i][j]*Cp/Pr; //K

            std::cout <<"I.Cs: Dep Var"<<"\t"<<i<<"\t"<<j<<"\t"<<dv[0][i][j]<<"\t"<<dv[1][i][j]<<"\t"<< dv[2][i][j]
                      <<"\t"<<dv[3][i][j]<<"\t"<<dv[4][i][j]<<"\t"<<dv[5][i][j]<<"\t"<<dv[6][i][j]<<"\t"<<dv[7][i][j]
                      << std::endl;
        }
    }

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            cv[0][i][j] = dv[0][i][j];
            cv[1][i][j] = dv[0][i][j] * dv[1][i][j];
            cv[2][i][j] = dv[0][i][j] * dv[2][i][j];
            cv[3][i][j] = dv[3][i][j]/(Gamma-1.0)+0.5*dv[0][i][j]*(dv[1][i][j]*dv[1][i][j]+dv[2][i][j]*dv[2][i][j]);

            std::cout<<"I.Cs: Con Var" << "\t" << i<<"\t"<<j<<"\t"<<cv[0][i][j] << "\t" << cv[1][i][j]<< "\t"
                     << cv[2][i][j]<< "\t" << cv[3][i][j] << std::endl;
        }
    }
}*/


/*
void Ic_Steady_Flow_SWBLI(int id2, int jd2, double rho_inf, double u_inf, double v_inf, double p_inf, double ***&cv,
                          double ***&dv) {

    double Einf;
    //Re=1301233.166, Lref=0.3048, Machinf=0.2, pinf=100000, tinf=300.0, alpha=0.0;
    //Re=100000, Lref=0.8, Machinf=2.15, pinf=100000, tinf=288.15, alpha=0.0;
    Re=96000, Lref=0.8, Machinf=2.2, pinf=107000, tinf=293.0, alpha=0.0;
    //Re=100000, Lref=0.8, Machinf=2.2, pinf=1.0, tinf=1.0, alpha=0.0;
    //Re=96000, Lref=0.8, Machinf=2.2, pinf=123000, tinf=303.0, alpha=0.0;
    //Re=100000, Lref=0.8, Machinf=2.15, pinf=10821.09160855605, tinf=152.24733697064173, alpha=0.0;
    //Re=100000, Lref=0.8, Machinf=2.15, pinf=12439.198764975645, tinf=157.4434918160561, alpha=0.0;
    //Re=100000, Lref=0.8, Machinf=2.15, pinf=107000.0/1.1, tinf=157.0003669264406, alpha=0.0;

    Cp = 1005.0;
    Rgas = ((Gamma-1.0)/Gamma)*Cp;
    rhoinf = pinf/(Rgas*tinf);
    qinf = Machinf*sqrt((Gamma-1.0)*Cp*tinf);
    alpha = alpha/rad;
    uinf = qinf*cos(alpha);
    vinf = qinf*sin(alpha);
    refrho = rhoinf;
    refvel = qinf;
    refvisc = rhoinf*qinf*Lref/Re;
    //refvisc = 1.0/Re;
    std::cout<<"the values are at Ic_cond :"<<rhoinf<<"\t"<<Rgas<<"\t"<<qinf<<"\t"<<uinf<<"\t"<<vinf<<"\t"<<refvisc<<std::endl;
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

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            dv[0][i][j] = rhoinf; //rho
            dv[1][i][j] = uinf; //u
            dv[2][i][j] = vinf; //v
            dv[3][i][j] = pinf; //p
            dv[4][i][j] = tinf; //T
            dv[5][i][j] = sqrt(Gamma*pinf/rhoinf); //a
            dv[6][i][j] = ((C0+C1)/(C1+dv[4][i][j]))*pow((dv[4][i][j]/C0),1.5)*refvisc; //Mu
            dv[7][i][j] = dv[6][i][j]*Cp/Pr; //K

            std::cout <<"I.Cs: Dep Var"<<"\t"<<i<<"\t"<<j<<"\t"<<dv[0][i][j]<<"\t"<<dv[1][i][j]<<"\t"<< dv[2][i][j]
                      <<"\t"<<dv[3][i][j]<<"\t"<<dv[4][i][j]<<"\t"<<dv[5][i][j]<<"\t"<<dv[6][i][j]<<"\t"<<dv[7][i][j]
                      << std::endl;
        }
    }

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            cv[0][i][j] = dv[0][i][j];
            cv[1][i][j] = dv[0][i][j] * dv[1][i][j];
            cv[2][i][j] = dv[0][i][j] * dv[2][i][j];
            cv[3][i][j] = dv[3][i][j]/(Gamma-1.0)+0.5*dv[0][i][j]*(dv[1][i][j]*dv[1][i][j]+dv[2][i][j]*dv[2][i][j]);

            std::cout<<"I.Cs: Con Var" << "\t" << i<<"\t"<<j<<"\t"<<cv[0][i][j] << "\t" << cv[1][i][j]<< "\t"
                     << cv[2][i][j]<< "\t" << cv[3][i][j] << std::endl;
        }
    }
}*/


void Ic_Steady_Flow_SWBLI(int id2, int jd2, double rho_inf, double u_inf, double v_inf, double p_inf, double ***&cv, double ***&dv) {

    double Einf;

    Re_inf=100000, Lref=0.8, Machinf=2.15, rhoinf = 1.0, tinf=1.0, qinf = 1.0, alpha=0.0;
    //Re=96000, Lref=0.08, Machinf=2.2, rhoinf = 1.0, tinf=1.0, qinf = 1.0, alpha=0.0;
    Tref = 288.16;

    refvisc = (rhoinf*qinf*Lref)/Re_inf;
    Rgas    = 1.0/(Gamma*Machinf*Machinf);
    Cp      = 1.0/((Gamma-1.0)*Machinf*Machinf);
    pinf    = rhoinf*Rgas*tinf;
    alpha   = alpha/rad;
    uinf    = qinf*cos(alpha);
    vinf    = qinf*sin(alpha);
    rhoref  = rhoinf;
    velref  = qinf;
    
    std::cout<<"the values are at Ic_cond:"<<rhoinf<<"\t"<<Rgas<<"\t"<<qinf<<"\t"<<uinf<<"\t"<<vinf<<"\t"<<refvisc<<std::endl;
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

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            dv[0][i][j] = rhoinf; //rho
            dv[1][i][j] = uinf; //u
            dv[2][i][j] = vinf; //v
            dv[3][i][j] = pinf; //p
            dv[4][i][j] = tinf; //T
            dv[5][i][j] = sqrt(Gamma*pinf/rhoinf); //a
            //dv[6][i][j] = ((C0+C1)/(C1+dv[4][i][j]))*pow((dv[4][i][j]/C0),1.5)*refvisc; //Mu
            //dv[6][i][j] = ((1.0+C1/C0)/(C1/C0+dv[4][i][j]))*pow(dv[4][i][j],1.5); //Mu
            dv[6][i][j] = ((1.0+C1/C0)/(C1/C0+dv[4][i][j]))*pow(dv[4][i][j],1.5)/Re_inf; //Mu
            //dv[6][i][j] = ((1.0+C1/C0)/(C1/C0+dv[4][i][j]))*pow(dv[4][i][j],1.5); //Mu
            //dv[7][i][j] = dv[6][i][j]*Cp/Pr; //K
            dv[7][i][j] = dv[6][i][j]/(GammaMinus1*Machinf*Machinf*Pr);

            std::cout <<"I.Cs: Dep Var"<<"\t"<<i<<"\t"<<j<<"\t"<<dv[0][i][j]<<"\t"<<dv[1][i][j]<<"\t"<< dv[2][i][j]
                      <<"\t"<<dv[3][i][j]<<"\t"<<dv[4][i][j]<<"\t"<<dv[5][i][j]<<"\t"<<dv[6][i][j]<<"\t"<<dv[7][i][j]
                      << std::endl;
        }
    }

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            cv[0][i][j] = dv[0][i][j];
            cv[1][i][j] = dv[0][i][j] * dv[1][i][j];
            cv[2][i][j] = dv[0][i][j] * dv[2][i][j];
            cv[3][i][j] = dv[3][i][j]/(Gamma-1.0)+0.5*dv[0][i][j]*(dv[1][i][j]*dv[1][i][j]+dv[2][i][j]*dv[2][i][j]);

            //std::cout<<"I.Cs: Con Var" << "\t" << i<<"\t"<<j<<"\t"<<cv[0][i][j] << "\t" << cv[1][i][j]<< "\t"
              //       << cv[2][i][j]<< "\t" << cv[3][i][j] << std::endl;
        }
    }
}

/**********************************************************************************************************************/
/**********************************************************************************************************************/
/*void Ic_Steady_Flow_SWBLI(int id2, int jd2, double rho_inf, double u_inf, double v_inf, double p_inf, double ***&cv,
                          double ***&dv) {

    Re=100000, Lref=0.8, Machinf=2.15, rhoinf = 1.0, tinf=1.0, qinf = 1.0, alpha=0.0;
    Tref = 288.16;
    const double RefVisc=1.716E-5;

    pinf = 1.0/(Gamma*Machinf*Machinf);
    Rgas = 1.0/(Gamma*Machinf*Machinf);
    Cp=1.0/((Gamma-1.0)*Machinf*Machinf);
    refvisc = Lref/Re;
    alpha = alpha/rad;
    uinf = qinf*cos(alpha);
    vinf = qinf*sin(alpha);
    //refrho = rhoinf;
    //refvel = qinf;


    std::cout<<"the values are at Ic_cond:"<<rhoinf<<"\t"<<Rgas<<"\t"<<qinf<<"\t"<<uinf<<"\t"<<pinf<<"\t"<<refvisc
             <<std::endl;
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

     for (int j = 0; j <= jd2; j++) {
         for (int i = 0; i <= id2; i++) {
             dv[0][i][j] = rhoinf; //rho
             dv[1][i][j] = uinf; //u
             dv[2][i][j] = vinf; //v
             dv[3][i][j] = pinf; //p
             dv[4][i][j] = tinf; //T
             dv[5][i][j] = sqrt(Gamma*pinf/rhoinf); //a
             dv[6][i][j] = ((C0+C1)/(C1+dv[4][i][j]))*pow((dv[4][i][j]/C0),1.5)*refvisc; //Mu
             dv[7][i][j] = dv[6][i][j]*Cp/Pr; //K

             std::cout <<"I.Cs: Dep Var"<<"\t"<<i<<"\t"<<j<<"\t"<<dv[0][i][j]<<"\t"<<dv[1][i][j]<<"\t"<< dv[2][i][j]
                       <<"\t"<<dv[3][i][j]<<"\t"<<dv[4][i][j]<<"\t"<<dv[5][i][j]<<"\t"<<dv[6][i][j]<<"\t"<<dv[7][i][j]
                       << std::endl;
         }
     }

     for (int j = 0; j <= jd2; j++) {
         for (int i = 0; i <= id2; i++) {
             cv[0][i][j] = dv[0][i][j];
             cv[1][i][j] = dv[0][i][j] * 6.471499999999999dv[1][i][j];
             cv[2][i][j] = dv[0][i][j] * dv[2][i][j];
             cv[3][i][j] = dv[3][i][j]/(Gamma-1.0)+0.5*dv[0][i][j]*(dv[1][i][j]*dv[1][i][j]+dv[2][i][j]*dv[2][i][j]);

             std::cout<<"I.Cs: Con Var" << "\t" << i<<"\t"<<j<<"\t"<<cv[0][i][j] << "\t" << cv[1][i][j]<< "\t"
                      << cv[2][i][j]<< "\t" << cv[3][i][j] << std::endl;
                    }
                }
 }*/

/**********************************************************************************************************************/
/**********************************************************************************************************************/
            /*
            void Ic_Steady_Flow_SWBLI(int id2, int jd2, double rho_inf, double u_inf, double v_inf, double p_inf, double ***&cv) {

                double Einf;

                //Re=1301233.166, Lref=0.3048, Machinf=0.2, pinf=100000, tinf=300.0, alpha=0.0;
                Re=100000, Lref=0.08, Machinf=2.15, pinf=1.0, rhoinf = 1.0, tinf=1.0, alpha=0.0;
                Tref = 288.16;
                //Re=96000, Lref=0.8, Machinf=2.2, pinf=107000, tinf=293.0, alpha=0.0;
                //Re=96000, Lref=0.8, Machinf=2.2, tinf=1.0, alpha=0.0;
                //Re=100000, Lref=0.8, Machinf=2.15, pinf=10821.09160855605, tinf=152.24733697064173, alpha=0.0;
                //Re=100000, Lref=0.8, Machinf=2.15, pinf=12439.198764975645, tinf=157.4434918160561, alpha=0.0;
                //Re=100000, Lref=0.8, Machinf=2.15, pinf=107000.0/1.1, tinf=157.0003669264406, alpha=0.0;

                Rgas = 1.0;
                qinf = Machinf*sqrt(Gamma);
                Cp=Gamma/(Gamma-1.0);
                refvisc = (qinf*Lref)/Re;
                qinf = 1.0/(Gamma*Machinf*Machinf);
                alpha = alpha/rad;
                uinf = qinf*cos(alpha);
                vinf = qinf*sin(alpha);
                //refrho = rhoinf;
                //refvel = qinf;


                std::cout<<"the values are at Ic_cond:"<<rhoinf<<"\t"<<Rgas<<"\t"<<qinf<<"\t"<<uinf<<"\t"<<pinf<<"\t"<<refvisc<<std::endl;
                for (int j = 0; j <= jd2; j++) {
                    for (int i = 0; i <= id2; i++) {
                        cv[0][i][j] = 0.0;
                        cv[1][i][j] = 0.0;
                        cv[2][i][j] = 0.0;
                        cv[3][i][j] = 0.0;
                    }
                }
                for (int j = 0; j <= jd2; j++) {
                    for (int i = 0; i <= id2; i++) {
                        cv[0][i][j] = rhoinf;
                        cv[1][i][j] = rhoinf * uinf;
                        cv[2][i][j] = rhoinf * vinf;
                        cv[3][i][j] = pinf/(Gamma-1.0) + 0.5*rhoinf*qinf*qinf;
                    }
                }

                for (int j = 0; j <= jd2; j++) {
                    for (int i = 0; i <= id2; i++) {

                        std::cout << "I.Cs" << "\t" << i<<"\t"<<j<<"\t"<<cv[0][i][j] << "\t" << cv[1][i][j]<< "\t" << cv[2][i][j]<< "\t" << cv[3][i][j] << std::endl;
                    }
                }
            }*/


void Ic_HyperSonicFlow(int id2, int jd2, double rho_inf, double u_inf, double v_inf, double p_inf, double ***&cv, double ***&dv) {

    //Re=514997.7, Lref=2.0*0.0381, Machinf=8.03, tinf=111.56, pinf=985.01, alpha=0.0;
    //Re=6758500.0, Lref=0.0381, Machinf=8.03, tinf=111.56, pinf=985.01, alpha=0.0;
    //Re=1.835e5, Lref=0.0381, Machinf=8.03, tinf=124.94, pinf=985.01, alpha=0.0;
    Re_inf=257354.0, Lref=0.0381, Machinf=8.03, tinf=111.56, pinf=985.01, alpha=0.0;

    Tref = 288.16;


    //Cp = 1005.0;
    //Rgas = ((Gamma-1.0)/Gamma)*Cp;
    Rgas = 286.92;
    Cp = Rgas*GamaByGamaMin1;
    rhoinf = pinf/(Rgas*tinf);
    qinf = Machinf*sqrt((Gamma-1.0)*Cp*tinf);
    alpha = alpha/rad;
    uinf = qinf*cos(alpha);
    vinf = qinf*sin(alpha);
    rhoref = rhoinf;
    velref = qinf;
    refvisc = rhoinf*qinf*Lref/Re_inf;
    //refvisc = 1.0/Re;

    std::cout<<"the values are at Ic_cond:"<<rhoinf<<"\t"<<Rgas<<"\t"<<qinf<<"\t"<<uinf<<"\t"<<pinf<<"\t"<<refvisc<<
             std::endl;
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

    for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {
            dv[0][i][j] = rhoinf; //rho
            dv[1][i][j] = uinf; //u
            dv[2][i][j] = vinf; //v
            dv[3][i][j] = pinf; //p
            dv[4][i][j] = tinf; //T
            dv[5][i][j] = sqrt(Gamma*pinf/rhoinf); //a
            dv[6][i][j] = ((C0+C1)/(C1+dv[4][i][j]))*pow((dv[4][i][j]/C0),1.5)*refvisc; //Mu
            dv[7][i][j] = dv[6][i][j]*Cp/Pr; //K

            /*std::cout <<"I.Cs: Dep Var"<<"\t"<<i<<"\t"<<j<<"\t"<<dv[0][i][j]<<"\t"<<dv[1][i][j]<<"\t"<< dv[2][i][j]
                      <<"\t"<<dv[3][i][j]<<"\t"<<dv[4][i][j]<<"\t"<<dv[5][i][j]<<"\t"<<dv[6][i][j]<<"\t"<<dv[7][i][j]
                      << std::endl;*/
        }
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
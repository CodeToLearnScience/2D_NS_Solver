#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

void Solver(int ib, int id1, int id2, int jb, int jd1, int jd2, int nconv, int ndepv, double **&x, double **&y,
    double ***&cv, double ***&dv, double ***&dui, double ***&duj, double **&area, double ***&si, double ***&sj,
    double **&tstep, double **&sri, double **&srj, double ***&gradfi, double ***&tgradfj, double ***&cvold,
    double ***&diss, double ***&rhs, double ***&epsij) {
    double pressure;
    double adtv=0;

    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= ib; i++) {
            cvold[0][i][j] = cv[0][i][j];
            cvold[1][i][j] = cv[1][i][j];
            cvold[2][i][j] = cv[2][i][j];
            cvold[3][i][j] = cv[3][i][j];

            diss[0][i][j] = 0.0;
            diss[1][i][j] = 0.0;
            diss[2][i][j] = 0.0;
            diss[3][i][j] = 0.0;
        }
    }

    if(eqn_type==0){
        Time_Step_Euler(ib, id1, id2, jb, jd1, jd2, nconv, ndepv, cv, dv, area, si, sj, tstep, sri, srj, epsij);
    }
    else if(eqn_type==1){
        Time_Step_NS(ib, id1, id2, jb, jd1, jd2, nconv, ndepv, cv, dv, area, si, sj, tstep, sri, srj, epsij);
    }


    switch (time_accuracy){
        case 0:
            //std::cout<<"I am  here in explicit time stepping"<<std::endl;
            ComputeFluxes(x, y, cv, dv, dui, duj, area, si, sj, tstep, sri,  srj, gradfi,  gradfj, cvold, diss, rhs,
                          epsij);

            for (int j = 2; j <= jb; j++) {
                for (int i = 2; i <= ib; i++) {
                    adtv = tstep[i][j] / area[i][j];
                    rhs[0][i][j] = adtv * rhs[0][i][j];
                    rhs[1][i][j] = adtv * rhs[1][i][j];
                    rhs[2][i][j] = adtv * rhs[2][i][j];
                    rhs[3][i][j] = adtv * rhs[3][i][j];
                }
            }

            //Imp_Res_Smoo(dirs, epsij, rhs);

            for (int j = 2; j <= jb; j++) {
                for (int i = 2; i <= ib; i++) {
                    cv[0][i][j] = cvold[0][i][j] - rhs[0][i][j];
                    cv[1][i][j] = cvold[1][i][j] - rhs[1][i][j];
                    cv[2][i][j] = cvold[2][i][j] - rhs[2][i][j];
                    cv[3][i][j] = cvold[3][i][j] - rhs[3][i][j];

                    pressure = (cv[3][i][j]-0.5*(cv[1][i][j]*cv[1][i][j]+cv[2][i][j]*cv[2][i][j])/cv[0][i][j])*(Gamma-1.0);

                    if ((cv[0][i][j] != cv[0][i][j]) or (cv[1][i][j] != cv[1][i][j]) or (cv[2][i][j] != cv[2][i][j]) or
                               (cv[3][i][j] != cv[3][i][j]) or pressure<0.0) {
                        std::cout << " In solver_explicit: the values of conserved var at i, j=" << i << "\t" << j << "\t"
                                  << cv[0][i][j] << "\t" << cv[1][i][j] << "\t" << cv[2][i][j] << "\t" << cv[3][i][j]
                                  <<"\t"<<"pressure="<<pressure<< std::endl;
                        exit(0);
                    }
                }
            }
            Dependent_Variables(id2, jd2, cv, dv);
            Boundary_Conditions(ib, id1, id2, jb, jd1, jd2, cv, dv, x, y, si, sj, pres_input_flag, bc_flag, bound_ind,
                                bound_cell, strt_bound_seg, end_bound_seg, pres_rhouvp);
            break;
        case 1:
            //SSPRK-2 loop
            for (int k=0;k<2;k++){

                ComputeFluxes(x, y, cv, dv, dui, duj, area, si, sj, tstep, sri,  srj, gradfi,  gradfj, cvold, diss,
                              rhs, epsij);

                for (int j = 2; j <= jb; j++) {
                    for (int i = 2; i <= ib; i++) {
                        //std::cout<<"ssprk2="<<i<<"\t"<<j<<"\t"<<k<<"\t"<<ConstSSPRK2a[k]<<"\t"<<ConstSSPRK2b[k]<<std::endl;
                        adtv = tstep[i][j] / area[i][j];
                        rhs[0][i][j] = adtv * rhs[0][i][j];
                        rhs[1][i][j] = adtv * rhs[1][i][j];
                        rhs[2][i][j] = adtv * rhs[2][i][j];
                        rhs[3][i][j] = adtv * rhs[3][i][j];
                    }
                }

                //Imp_Res_Smoo(dirs, epsij, rhs);

                for (int j = 2; j <= jb; j++) {
                    for (int i = 2; i <= ib; i++) {
                        cv[0][i][j] = ConstSSPRK2a[k] * cvold[0][i][j] + ConstSSPRK2b[k] * (cv[0][i][j] - rhs[0][i][j]);
                        cv[1][i][j] = ConstSSPRK2a[k] * cvold[1][i][j] + ConstSSPRK2b[k] * (cv[1][i][j] - rhs[1][i][j]);
                        cv[2][i][j] = ConstSSPRK2a[k] * cvold[2][i][j] + ConstSSPRK2b[k] * (cv[2][i][j] - rhs[2][i][j]);
                        cv[3][i][j] = ConstSSPRK2a[k] * cvold[3][i][j] + ConstSSPRK2b[k] * (cv[3][i][j] - rhs[3][i][j]);

                        pressure = (cv[3][i][j]-0.5*(cv[1][i][j]*cv[1][i][j]+cv[2][i][j]*cv[2][i][j])/cv[0][i][j])
                                   *(Gamma-1.0);

                        if ((cv[0][i][j] != cv[0][i][j]) or (cv[1][i][j] != cv[1][i][j]) or (cv[2][i][j] != cv[2][i][j])
                            or (cv[3][i][j] != cv[3][i][j]) or pressure<0.0) {
                            std::cout << " In solver ssprk2: the values of conserved var at i, j=" << i << "\t" << j
                                      << "\t" << cv[0][i][j] << "\t" << cv[1][i][j] << "\t" << cv[2][i][j] << "\t"
                                      << cv[3][i][j] <<"\t"<<"pressure="<<pressure<< std::endl;
                            exit(0);
                        }
                    }
                }
                Dependent_Variables(id2, jd2, cv, dv);
                Boundary_Conditions(ib, id1, id2, jb, jd1, jd2, cv, dv, x, y, si, sj, pres_input_flag, bc_flag,
                                    bound_ind, bound_cell, strt_bound_seg, end_bound_seg, pres_rhouvp);
            }
            break;
        case 2:
            //SSPRK-3 loop
            for (int k=0;k<3;k++){

                ComputeFluxes(x, y, cv, dv, dui, duj, area, si, sj, tstep, sri,  srj, gradfi,  gradfj, cvold, diss,
                              rhs, epsij);

                for (int j = 2; j <= jb; j++) {
                    for (int i = 2; i <= ib; i++) {
                        adtv = tstep[i][j] / area[i][j];
                        rhs[0][i][j] = adtv * rhs[0][i][j];
                        rhs[1][i][j] = adtv * rhs[1][i][j];
                        rhs[2][i][j] = adtv * rhs[2][i][j];
                        rhs[3][i][j] = adtv * rhs[3][i][j];
                    }
                }

                //Imp_Res_Smoo(dirs, epsij, rhs);

                for (int j = 2; j <= jb; j++) {
                    for (int i = 2; i <= ib; i++) {
                        cv[0][i][j] = ConstRK1[k] * cvold[0][i][j] + ConstRK2[k] * (cv[0][i][j] - rhs[0][i][j]);
                        cv[1][i][j] = ConstRK1[k] * cvold[1][i][j] + ConstRK2[k] * (cv[1][i][j] - rhs[1][i][j]);
                        cv[2][i][j] = ConstRK1[k] * cvold[2][i][j] + ConstRK2[k] * (cv[2][i][j] - rhs[2][i][j]);
                        cv[3][i][j] = ConstRK1[k] * cvold[3][i][j] + ConstRK2[k] * (cv[3][i][j] - rhs[3][i][j]);

                        pressure = (cv[3][i][j]-0.5*(cv[1][i][j]*cv[1][i][j]+cv[2][i][j]*cv[2][i][j])/cv[0][i][j])*
                                                                                                            (Gamma-1.0);

                        if ((cv[0][i][j] != cv[0][i][j]) or (cv[1][i][j] != cv[1][i][j]) or (cv[2][i][j] != cv[2][i][j])
                            or (cv[3][i][j] != cv[3][i][j]) or pressure<0.0) {
                            std::cout << " In solver ssprk3: the values of conserved var at i, j=" << i << "\t" << j << "\t"
                                      << cv[0][i][j] << "\t" << cv[1][i][j] << "\t" <<cv[2][i][j]<< "\t" << cv[3][i][j]
                                      <<"\t"<<"pressure="<<pressure<< std::endl;
                            exit(0);
                        }
                    }
                }
                Dependent_Variables(id2, jd2, cv, dv);
                Boundary_Conditions(ib, id1, id2, jb, jd1, jd2, cv, dv, x, y, si, sj, pres_input_flag, bc_flag,
                                    bound_ind, bound_cell, strt_bound_seg, end_bound_seg, pres_rhouvp);
            }
            break;
        default:
            std::cout<<"wrong time_accuray input: chose time_accuracy="<<time_accuracy<<std::endl;
            break;
    }
}

void ComputeFluxes(double **&x, double **&y,
                   double ***&cv, double ***&dv, double ***&dui, double ***&duj, double **&area, double ***&si,
                   double ***&sj, double **&tstep, double **&sri, double **&srj, double ***&gradfi, double ***&tgradfj,
                   double ***&cvold, double ***&diss, double ***&rhs, double ***&epsij){


    if(eqn_type==1){
        //std::cout<<"I am  here in viscous flux calculations"<<std::endl;
        Initialize_Gradients(id2, jd2, nconv, cv, u, v, gradfi, gradfj);
        Gradient_FaceI(nconv, ndvar, cv, dv, u, v, area, si, sj, gradfi);
        Gradient_FaceJ(nconv, ndvar, cv, dv, u, v, area, si, sj, gradfj);
        Flux_Viscous(nconv, ndvar, cv, dv, u, v, area, si, sj, gradfi, gradfj);
    }
    else if (eqn_type!=1 and eqn_type!=0){
        std::cout<<"error in equation type"<<std::endl;
        exit(0);
    }

    FluxConvDefns(id2, jd2, id1, jd1, ib, jb, cv, dv, iAvgFlux, jAvgFlux, iFluxDiff, jFluxDiff, iConVarDiff,
                  jConVarDiff);

    if(space_accuracy==1){
        //std::cout<<"I am  here in second order space"<<std::endl;
        Prim_Variables_Differences(ib, id2, jb, jd2, cv, dv, dui, duj);
        //Conv_Variables_Differences(ib, id2, jb, jd2, cv, dv, dui, duj);
        if(inviscid_scheme==20){
            Diss_LLF2_Prim(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
        }
        else if(inviscid_scheme==21){
            std::cout<<"second order movers is yet to be implemented"<<std::endl;
            //Diss_MOVERS_H2_Prim(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
            //Diss_MOVERS_H2_Conv(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
        }
        else if(inviscid_scheme==22){
            Diss_MOVERS_H2_Prim(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
            //Diss_MOVERS_H2_Conv(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
        }
        else if(inviscid_scheme==23){
            Diss_MOVERS_LE2_Prim(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
            //Diss_MOVERS_LE2_Conv(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
        }
        else if(inviscid_scheme==24){
            //std::cout<<"I am  here applying this scheme"<<std::endl;
            Diss_Roe2_Prim(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
        }
        else if(inviscid_scheme==25){
            Diss_Roe2_TV_Prim(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
        }
        else if(inviscid_scheme==26){
            FluxConvDefns(id2, jd2, id1, jd1, ib, jb, cv, dv, iAvgFlux, jAvgFlux, iFluxDiff, jFluxDiff, iConVarDiff,
                          jConVarDiff);
            KFDS_2ndOrder(id2, jd2, id1, jd1, ib, jb, cv, dv, iAvgFlux, jAvgFlux, iFluxDiff, jFluxDiff, iConVarDiff,
                          jConVarDiff);
        }
        else if(inviscid_scheme==27){
            Diss_ECCS(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
        }
        else if(inviscid_scheme==28){
            Diss_MOVERS(ib, id1, jb, jd1, cv, dv, si, sj, diss);
        }
        else if(inviscid_scheme==29){
            Diss_MOVERS_NWSC(ib, id1, jb, jd1, cv, dv, si, sj, diss);
        }
    }

    else if(space_accuracy==0){
        //std::cout<<"first order spatial accuracy is being used"<<std::endl;
        if(inviscid_scheme==0){
            Diss_LLF1(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
        }
        else if(inviscid_scheme==1){
            Diss_MOVERS(ib, id1, jb, jd1, cv, dv, si, sj, diss);
        }
        else if(inviscid_scheme==2){
            Diss_MOVERS_H1(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
            //MOVERS_H1(ib, id1, id2, jb, jd1, jd2, cv, dv, dui, duj, area, si, sj, diss, rhs);
        }
        else if(inviscid_scheme==3){
            Diss_MOVERS_LE1(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
        }
        else if(inviscid_scheme==5){
            //Diss_Roe1(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
            Diss_Roe1_TV(ib, id1, id2, jb, jd1, jd2, cv, dv, si, sj, diss);
        }
        else if(inviscid_scheme==6){
            FluxConvDefns(id2, jd2, id1, jd1, ib, jb, cv, dv, iAvgFlux, jAvgFlux, iFluxDiff, jFluxDiff, iConVarDiff,
                          jConVarDiff);
            KFDS_1stOrder(id2, jd2, id1, jd1, ib, jb, cv, dv, iAvgFlux, jAvgFlux, iFluxDiff, jFluxDiff, iConVarDiff,
                          jConVarDiff);
        }
    }

    //Avg_Flux1(ib, id1, jb, jd1, cv, dv, si, sj, diss, rhs);
    Avg_Conv_Flux1(ib, id1, id2, jb, jd1, jd2, cv, dv, area, si, sj, diss, rhs);
    //Avg_Conv_Flux2(ib, id1, id2, jb, jd1, jd2, cv, dv, area, si, sj, diss, rhs);
}

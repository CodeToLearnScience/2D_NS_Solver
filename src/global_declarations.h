#ifndef _Header_H
#define _Header_H

#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <string>
#include <iosfwd>
#include <sstream>
#include <limits>

/*****constants used throught the program******/
const double Gamma = 1.4, pi = 4.0 * atan(1.0), rad = 180.0 / pi, Pr=0.72;
const double C0=288.15, C2=1.458e-6, C1=110.4;
const double GammaMinus1 = Gamma - 1.0, GamaByGamaMin1=Gamma/GammaMinus1;
//const double Rgas = 286.92;

const int ndvar = 8, nconv = 4;
extern int mxeqn;

/*******SSPRK coefficients************/
static const double ConstRK1[]={0.0, 3.0/4.0, 1.0/3.0}, ConstRK2[]={1.0, 1.0/4.0, 2.0/3.0};
static const double ConstSSPRK2a[]={0.0, 1.0/2.0}, ConstSSPRK2b[]={1.0, 1.0/2.0};
/****global integer variables related to grid and domain*****/
extern int Nx, Ny, ib, jb, id1, id2, jd1, jd2, imax, jmax;

/****integers related to boundary conditions****/
extern int bind, sbind, ebind, pbind, spbind, epbind, iresmax, jresmax;
extern int no_boun_seg, no_bounseg_prescinp;
extern int *pres_input_flag, *bc_flag, *bound_ind, *bound_cell, *strt_bound_seg, *end_bound_seg;
extern double pres_rho, pres_u, pres_v, pres_p;
extern double **pres_rhouvp;

/****other variables used throughout the code****/
extern double Machinf, Re_inf, Lref, rhoinf, uinf, vinf, pinf, tinf, alpha, pref, rhoref, velref, Tref, qinf, cp, cl,
        cd, cm, Twall_iso,
        minject, tinject, maxwchg, maxichg, ttinl, ptinl, pout, betainl, betaout, Rgas, Cp, refvisc, Xref, Yref,
        epsirs, restart_t;
extern double iGlobalMaxEigVal, iGlobalMinEigVal, jGlobalMaxEigVal, jGlobalMinEigVal;

/****for residue calculation****/
extern double Resn1, Resn2, Resn3, Resn4;

/****global strings declared here, which will be used throughout the code*****/
extern std::string flo_intr_extr, flux_type, grid_file, vort_corr, time_step, RestartFile, MeshFile, test_case,
        MeshTopFile;
extern int eqn_type, flow_type, space_accuracy, time_accuracy, dimen, inviscid_scheme, visc_method, scaling_factor, restart,
        MaxIter, restart_iter;


/****global pointers declared here, which will be used throughout the code*****/
extern double cfl;
extern double *Icond, *Prim_var_L, *Prim_var_M, *Prim_var_R, *limref;
extern double ***cv, ***cvold, ***dv, ***rhs, ***dui, ***duj, ***cui, ***cuj, ***si, ***sj, ***pbvar, ***gradfi,
        ***gradfj, ***diss, ***epsij;
extern double **x, **y, **xc, **yc, **area, **p, **sri, **srj, **tstep, **u, **v, **dirs;
extern double **L_state, **R_state;
extern  double ***iAvgFlux, ***jAvgFlux, ***iFluxDiff, ***jFluxDiff, ***iConVarDiff, ***jConVarDiff, ***irhs, ***jrhs;

extern std::ofstream Res;
/*********************forward declaration of functions called in the main program **********************/
void Declaration_1d_array();

void Declaration_2d_array();

void Declaration_3d_array();

void Read_grid(std::string &grid_file, int Nx, int id1, int id2, int Ny, int jd1, int jd2, double **&x, double **&y);

void Grid_Computations(int ib, int jb, int id1, int id2, int jd1, int jd2, double **&x, double **&y, double **&area,
                       double ***&si, double ***&sj);

void Initialize_Gradients(int id2, int jd2, int nconv, double ***&cv, double **&u, double **&v, double ***&gradfi,
                          double ***&gradfj);
void ReadSolnFileToRestart(const std::string& RestartFile, int Nx, int id1, int id2, int Ny, int jd1, int jd2,
                           double ***&dv, double ***&cv);
void Ic_Steady_Flow(int id2, int jd2, double rhoinf, double uinf, double vinf, double pinf, double ***&dv,
                    double ***&cv);
void Ic_Steady_Flow_SWBLI(int id2, int jd2, double rhoinf, double uinf, double vinf, double pinf, double ***&cv,
                          double ***&dv);

void
Ic_Unsteady_Flow(int id1, int id2, int jd1, int jd2, int index, int loc, double rhol, double rhor, double ul, double ur,
                 double vl, double vr, double pl, double pr, double **&x, double **&y, double ***&cv);

void Dependent_Variables(int id2, int jd2, double ***&cv, double ***&dv);

void Dependent_Variables_One(int i, int j, double ***&cv, double ***&dv);

void Forces(int bind, int bnode, int sbind, int ebind, int ndvar, double **&x, double **&y, double ***&dv,
            double ***&si, double ***&sj);

void Bc_Transmitive(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv);

void Bc_Prescribed_Inflow(int bind, int bnode, int sbind, int ebind, double rhoinf, double uinf, double vinf,
                          double pinf, double ***&cv, double ***&dv);

void Bc_Prescribed_Inflow2(int ib, int id1, int id2, int jb, int jd1, int jd2, double loc, double rhol, double rhor,
                           double ul, double ur, double vl, double vr, double pl, double pr, int bind, int sbind,
                           int ebind, double **&x, double **&y, double ***&cv);
void Ic_HyperSonicFlow(int id2, int jd2, double rho_inf, double u_inf, double v_inf, double p_inf, double ***&cv,
                       double ***&dv);
void BC_wall(int bind, int bnode, int sbind, int ebind, double ***&cv,
             double ***&dv);

void BC_Eulerwall(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv);

void BC_Eulerwall_Interior(int int_dum, int bind, int sbind, int ebind, double ***&cv, double ***&dv);

void Flux_Wall(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv, double ***&si, double ***&sj,
               double ***&rhs);
void BC_NS_wall(int bind, int bnode, int sbind, int ebind, double ***&cv,
                double ***&dv);

void Flux_Interior_Wall(int int_dum, int bind, int sbind, int ebind, double ***&dv, double ***&rhs);

void BC_Interior_wall(int int_dum, int bind, int sbind, int ebind, double ***&cv, double ***&dv);
void BC_FarField(int bind, int bnode, int sbind, int ebind, double ***&si, double ***&sj, double ***&cv,
                 double ***&dv, double rhoinf, double uinf, double vinf, double pinf);

void Bc_Symmetry(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv);
void Bc_Cut(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv);
void BC_NS_wallIsothermal(int bind, int bnode, int sbind, int ebind, double ***&cv,
                          double ***&dv, double twall);

void Boundary_Conditions(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double **&x,
                         double **&y, double ***&si, double ***&sj, int *&pres_input_flag, int *&bc_flag, int *&bound_ind,
                         int *&bound_cell, int *&strt_bound_seg, int *&end_bound_seg, double **&pres_rhouvp);

void Prim_Variables_Differences(int ib, int id2, int jb, int jd2, double ***cv, double ***dv, double ***&dui,
                                double ***&duj);

void Conv_Variables_Differences(int ib, int id2, int jb, int jd2, double ***cv, double ***dv, double ***&dui,
                                double ***&duj);

void Time_Step_Euler(int ib, int id1, int id2, int jb, int jd1, int jd2, int nconv, int ndepv, double ***&cv,
                  double ***&dv, double **&area, double ***&si, double ***&sj, double **&tstep, double **&sri,
                  double **&srj, double ***&epsij);

void Time_Step_NS(int ib, int id1, int id2, int jb, int jd1, int jd2, int nconv, int ndepv, double ***&cv,
                     double ***&dv, double **&area, double ***&si, double ***&sj, double **&tstep, double **&sri,
                     double **&srj, double ***&epsij);

void Avg_Conv_Flux1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double **&area,
                    double ***&si, double ***&sj, double ***&diss, double ***&rhs);

void Avg_Conv_Flux2(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double **&area,
                    double ***&si, double ***&sj, double ***&diss, double ***&rhs);

void Diss_LLF1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
               double ***&sj, double ***&diss);

void Diss_LLF2_Prim(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                    double ***&sj, double ***&diss);

void Diss_Roe1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
               double ***&sj, double ***&diss);

void Diss_Roe1_TV(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                  double ***&sj, double ***&diss);

void Diss_Roe2_Prim(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                    double ***&sj, double ***&diss);

void Diss_Roe2_TV_Prim(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                       double ***&sj, double ***&diss);

void Diss_MOVERS(int ib, int id1, int jb, int jd1, double ***&cv, double ***&dv, double ***&si, double ***&sj,
                 double ***&diss);

void Diss_MOVERS2_Prim(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                       double ***&sj, double ***&diss);

void Diss_MOVERS_H1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                    double ***&sj, double ***&diss);

void Diss_MOVERS_LE1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                     double ***&sj, double ***&diss);

void
Diss_MOVERS_H2_Prim(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                    double ***&sj, double ***&diss);

void
Diss_MOVERS_H2_Conv(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                    double ***&sj, double ***&diss);

void
Diss_MOVERS_LE2_Conv(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                     double ***&sj, double ***&diss);

void
Diss_MOVERS_LE2_Prim(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
                     double ***&sj, double ***&diss);

void MOVERS_H1(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&dui,
               double ***&duj, double **&area, double ***&si, double ***&sj, double ***&diss, double ***&rhs);

void
Flux_Boundary(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv, double ***&si, double ***&sj,
              double ***&rhs);

void
Gradient_FaceI(int nconv, int ndvar, double ***&cv, double ***&dv,
               double **&u, double **&v, double **&area, double ***&si, double ***&sj, double ***&gradfi);
void
Gradient_FaceJ(int nconv, int ndvar, double ***&cv, double ***&dv, double **&u, double **&v, double **&area,
               double ***&si, double ***&sj, double ***&gradfi);
void Flux_Viscous(int nconv, int ndvar, double ***&cv, double ***&dv, double **&u, double **&v, double **&area,
                  double ***&si, double ***&sj, double ***&gradfi, double ***&gradfj);
void ComputeFluxes(double **&x, double **&y,
                   double ***&cv, double ***&dv, double ***&dui, double ***&duj, double **&area, double ***&si,
                   double ***&sj, double **&tstep, double **&sri, double **&srj, double ***&gradfi, double ***&tgradfj,
                   double ***&cvold, double ***&diss, double ***&rhs, double ***&epsij);

void Imp_Res_Smoo(double **&dirs, double ***&epsij, double ***&rhs);

void Solver(int ib, int id1, int id2, int jb, int jd1, int jd2, int nconv, int ndepv, double **&x, double **&y,
            double ***&cv, double ***&dv, double ***&dui, double ***&duj, double **&area, double ***&si, double ***&sj,
            double **&tstep, double **&sri, double **&srj, double ***&gradfi, double ***&tgradfj, double ***&cvold,
            double ***&diss, double ***&rhs, double ***&epsij);

void Residue_Cal(int iter, double dt, double t, double ***&cv, double ***&cvold, double ***&dv, double ***&si, double ***&sj);

void Write_Solution(int id1, int jd1, int iter, double t, double **&x, double **&y, double ***&cv);
void Interior_Solution(int bind, int indx, int beg_seg, int end_seg, int id1, int jd1, int iter, double t, double **&x,
                       double **&y, double ***&cv, double ***&dv, double ***&gradfi, double ***&gradfj);
void Write_Surf_Solution(int bind, int sur_ind, int beg_seg, int end_seg, int id1, int jd1, int iter, double t, double **&x, double **&y, double ***&cv,
                         double ***&dv, double ***&gradfi, double ***&gradfj);
void WriteRestartFile(int ib, int jb, int iter, double t, double ***&dv);
void Avg_Flux1(int ib, int id1, int jb, int jd1, double ***&cv, double ***&dv, double ***&si, double ***&sj,
               double ***&diss, double ***&rhs);
void Flux_ViscWall(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv, double ***&si, double ***&sj,
                   double ***&rhs);
void Flux_EulerWall(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv, double ***&si, double ***&sj,
                   double ***&rhs);
void Flux_Symmetry(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv, double ***&si, double ***&sj,
                   double ***&rhs);
void FluxConvDefns(int id2, int jd2, int id1, int jd1, int ib, int jb, double ***&cv, double ***&dv,
                   double ***&iAvgFlux, double ***&jAvgFlux, double ***&iFluxDiff, double ***&jFluxDiff,
                   double ***&iConVarDiff, double ***&jConVarDiff);
void KFDS_1stOrder(int id2, int jd2, int id1, int jd1, int ib, int jb, double ***&cv, double ***&dv,
                   double ***&iAvgFlux, double ***&jAvgFlux, double ***&iFluxDiff, double ***&jFluxDiff,
                   double ***&iConVarDiff, double ***&jConVarDiff);
void KFDS_2ndOrder(int id2, int jd2, int id1, int jd1, int ib, int jb, double ***&cv, double ***&dv,
              double ***&iAvgFlux, double ***&jAvgFlux, double ***&iFluxDiff, double ***&jFluxDiff,
              double ***&iConVarDiff, double ***&jConVarDiff);
void ReadSolverInput(const std::string& SolverInputFile);
void InitFlowDimensional(int id2, int jd2, double Re_inf, double Machinf, double Lref,  double alpha, double pinf,
                         double ***&cv, double ***&dv);
void InitFlowNonDimensional(int id2, int jd2, double Re_inf, double Machinf, double Lref,  double alpha, double ***&cv,
                            double ***&dv);
void ReadGridTopology(std::string &grid_top_file, int *&pres_input_flag, int *&bc_flag, int *&bound_ind,
                      int *&bound_cell, int *&strt_bound_seg, int *&end_bound_seg, double **&pres_rhouvp);
void extractIntegerWords(std::string str);
void Diss_ECCS(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double ***&si,
               double ***&sj, double ***&diss);
void Diss_MOVERS_NWSC(int ib, int id1, int jb, int jd1, double ***&cv, double ***&dv, double ***&si, double ***&sj,
                      double ***&diss);

#endif  //#ifndef _Header_H

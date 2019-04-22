#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

/****global integer variables related to grid and domain*****/

int Nx, Ny, ib, jb, id1, id2, jd1, jd2, imax, jmax;

int maxeqn, restart_iter;

/****global strings declared here, which will be used throughout the code*****/
std::string MeshFile, MeshTopFile, flo_intr_extr, flux_type, grid_file, vort_corr, time_step,
        RestartFile, test_case;
int eqn_type, flow_type, time_accuracy, space_accuracy, dimen, inviscid_scheme, visc_method, scaling_factor, restart,
        disp_freq, outp_freq, MaxIter;
double cfl;

/****variables related to boundary conditons****/
int bind, sbind, ebind, pbind, spbind, epbind, idum1, idum2, jdum1, jdum2,
        iins1, iins2, jins1, jins2, iresmax, jresmax;
int no_boun_seg, no_bounseg_prescinp;
int *pres_input_flag=NULL, *bc_flag=NULL, *bound_ind=NULL, *bound_cell=NULL, *strt_bound_seg=NULL, *end_bound_seg=NULL;
double iGlobalMaxEigVal, iGlobalMinEigVal, jGlobalMaxEigVal, jGlobalMinEigVal;
double pres_rho, pres_u, pres_v, pres_p;
double **pres_rhouvp;

/****other variables used throughout the code****/
double Machinf, Re_inf, Lref, rhoinf, uinf, vinf, pinf, tinf, alpha, pref, rhoref, velref, Tref, qinf, cp, cl, cd, cm,
        Twall_iso, minject, tinject, maxwchg, maxichg, ttinl, ptinl, pout, betainl, betaout, Rgas=286.92, Cp, refvisc,
        Xref, Yref, epsirs=0.3, restart_t;

/****for residue calculation****/
double Resn1, Resn2, Resn3, Resn4;

/****declaration of the global pointer variables used in the code****/
double *Icond, *Prim_var_L, *Prim_var_R, *limref;
double ***cv, ***cvold, ***dv, ***rhs, ***dui, ***duj, ***cui, ***cuj, ***si, ***sj, ***pbvar, ***diss, ***epsij, ***gradfi, ***gradfj;
double **x, **y, **xc, **yc, **p, **area, **sri, **srj, **tstep, **u, **v, **dirs;
double ***iAvgFlux, ***jAvgFlux, ***iFluxDiff, ***jFluxDiff, ***iConVarDiff, ***jConVarDiff, ***irhs, ***jrhs;

std::ofstream Res;
/**************  to allocate 1d array***********/
void Declaration_1d_array() {
    Icond = new double[nconv];
    Prim_var_L = new double[nconv];
    Prim_var_R = new double[nconv], limref = new double[nconv];
}

/**************  to allocate 2d array***********/
void Declaration_2d_array() {
    Allocate_2D(x, imax, jmax);
    Allocate_2D(y, imax, jmax);
    Allocate_2D(xc, imax, jmax);
    Allocate_2D(yc, imax, jmax);
    Allocate_2D(area, imax, jmax);
    Allocate_2D(p, imax, jmax);
    Allocate_2D(sri, imax, jmax);
    Allocate_2D(srj, imax, jmax);
    Allocate_2D(tstep, imax, jmax);
    Allocate_2D(u, imax, jmax);
    Allocate_2D(v, imax, jmax);
    Allocate_2D(dirs, imax, jmax);
}

/**************to allocate 3d array***********/
void Declaration_3d_array(){
    Allocate_3D(cv, nconv, imax, jmax);
    Allocate_3D(cvold, nconv, imax, jmax);
    Allocate_3D(dv, ndvar, imax, jmax);
    Allocate_3D(rhs, nconv, imax, jmax);
    Allocate_3D(dui, nconv, imax, jmax);
    Allocate_3D(duj, nconv, imax, jmax);
    Allocate_3D(si, 2, imax, jmax);
    Allocate_3D(sj, 2, imax, jmax);
    Allocate_3D(pbvar, 10, imax, jmax);
    Allocate_3D(gradfi, 6, imax, jmax);
    Allocate_3D(gradfj, 6, imax, jmax);
    Allocate_3D(diss, nconv, imax, jmax);
    Allocate_3D(epsij, 2, imax, jmax);
    Allocate_3D(cui, nconv, imax, jmax);
    Allocate_3D(cuj, nconv, imax, jmax);
    Allocate_3D(iAvgFlux, nconv, imax, jmax);
    Allocate_3D(jAvgFlux, nconv, imax, jmax);
    Allocate_3D(iFluxDiff, nconv, imax, jmax);
    Allocate_3D(jFluxDiff, nconv, imax, jmax);
    Allocate_3D(iConVarDiff, nconv, imax, jmax);
    Allocate_3D(jConVarDiff, nconv, imax, jmax);
    Allocate_3D(irhs, nconv, imax, jmax);
    Allocate_3D(jrhs, nconv, imax, jmax);
}

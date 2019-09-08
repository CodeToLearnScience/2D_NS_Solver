#include "../inc/global_declarations.h"

void Boundary_Conditions(int ib, int id1, int id2, int jb, int jd1, int jd2, double ***&cv, double ***&dv, double **&x,
                         double **&y, double ***&si, double ***&sj, int *&pres_input_flag, int *&bc_flag, int *&bound_ind,
                         int *&bound_cell, int *&strt_bound_seg, int *&end_bound_seg, double **&pres_rhouvp) {

    /****fill dummy nodes at the corners****/

    for (int neqn = 0; neqn < ndvar; neqn++) {
        dv[neqn][1][1]     = 0.5 * (dv[neqn][1][2] + dv[neqn][2][1]);
        dv[neqn][id1][1]   = 0.5 * (dv[neqn][ib][1] + dv[neqn][id1][2]);
        dv[neqn][1][jd1]   = 0.5 * (dv[neqn][2][jd1] + dv[neqn][1][jb]);
        dv[neqn][id1][jd1] = 0.5 * (dv[neqn][id1][jb] + dv[neqn][ib][jd1]);
    }

    Dependent_Variables_One(1, 1, cv, dv);
    Dependent_Variables_One(id1, 1, cv, dv);
    Dependent_Variables_One(1, jd1, cv, dv);
    Dependent_Variables_One(id1, jd1, cv, dv);

    int no_of_pres_count=-1;
    for(int bc=0; bc<no_boun_seg; bc++) {
        if(pres_input_flag[bc]==0 and bc_flag[bc]==10){ //prescribed inflow b.c
            //std::cout<<"I am applying b.c 10 without 1"<<std::endl;
            Bc_Prescribed_Inflow(bound_ind[bc], bound_cell[bc], strt_bound_seg[bc], end_bound_seg[bc], rhoinf, uinf,
                                 vinf, pinf, cv, dv);
        }
        else if(pres_input_flag[bc]==1 and bc_flag[bc]==10){ //prescribed inflow b.c with values other than free stream
            no_of_pres_count++;
            //::cout<<"I am applying b.c 10 with 1"<<std::endl;
            Bc_Prescribed_Inflow(bound_ind[bc], bound_cell[bc], strt_bound_seg[bc], end_bound_seg[bc],
                                 pres_rhouvp[0][no_of_pres_count], pres_rhouvp[1][no_of_pres_count],
                                 pres_rhouvp[2][no_of_pres_count], pres_rhouvp[3][no_of_pres_count], cv, dv);
        }
        else if(bc_flag[bc]==20){ //transmissive boundary condition
            //std::cout<<"I am applying b.c 20"<<std::endl;
            Bc_Transmitive(bound_ind[bc], bound_cell[bc], strt_bound_seg[bc], end_bound_seg[bc], cv, dv);
        }
        else if(bc_flag[bc]==30){ //Euler wall boundary condition
            // BC_wall(bound_ind[bc], bound_cell[bc], strt_bound_seg[bc], end_bound_seg[bc], cv, dv);
            BC_Eulerwall(bound_ind[bc], bound_cell[bc], strt_bound_seg[bc], end_bound_seg[bc], cv, dv);
        }
        else if(bc_flag[bc]==40){ //navier-stokes wall boundary condition
            //std::cout<<"I am applying b.c 40"<<std::endl;
            BC_NS_wall(bound_ind[bc], bound_cell[bc], strt_bound_seg[bc], end_bound_seg[bc], cv, dv);
        }
        else if(bc_flag[bc]==50){ //navier-stokes isothermal-wall boundary condition
            BC_NS_wallIsothermal(bound_ind[bc], bound_cell[bc], strt_bound_seg[bc], end_bound_seg[bc], cv, dv,
                                 Twall_iso);
        }
        else if(bc_flag[bc]==60){ //far-field boundary conditions
            //std::cout<<" I am here in farfield"<<std::endl;
            BC_FarField(bound_ind[bc], bound_cell[bc], strt_bound_seg[bc], end_bound_seg[bc], si, sj, cv, dv, rhoinf,
                        uinf, vinf, pinf);
        }
        else if(bc_flag[bc]==70){ //Symmetry boundary condition
            //std::cout<<"I am applying b.c 70"<<std::endl;
            Bc_Symmetry(bound_ind[bc], bound_cell[bc], strt_bound_seg[bc], end_bound_seg[bc], cv, dv);
        }
        else if(bc_flag[bc]==80){ //Cut wall boundary condition
            Bc_Cut(bound_ind[bc], bound_cell[bc], strt_bound_seg[bc], end_bound_seg[bc], cv, dv);
        }
    }
}
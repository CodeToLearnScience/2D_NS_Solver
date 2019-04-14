#include "../inc/global_declarations.h"

int main(int argc, char** argv) {

    std::cout << std::setprecision(8);
    std::cout << std::ios::dec;
    std::cout << std::scientific;
    
    ReadSolverInput(argv[1]);
    //ReadSolverInput("../input/InputHypersonicFlowEuler.dat");
    //ReadSolverInput("../input/InputRampEuler.dat");
    //ReadSolverInput("../input/InputSWBLI_NS.dat");
    //ReadSolverInput("../input/InputSWBLI_NS_ND1.dat");
    //ReadSolverInput("../input/InputSWBLI_NS_Dimen.dat");
    //ReadSolverInput("../input/InputBlasius.dat");
    ib = Nx + 1, jb = Ny + 1, id1 = Nx + 2, id2 = Nx + 3, jd1 = Ny + 2, jd2 = Ny + 3, imax = Nx + 4, jmax = Ny + 4;

    ReadGridTopology(MeshTopFile, pres_input_flag, bc_flag, bound_ind, bound_cell, strt_bound_seg, end_bound_seg,
                     pres_rhouvp);
    Declaration_1d_array();
    Declaration_2d_array();
    Declaration_3d_array();

    Read_grid(MeshFile, Nx, id1, id2, Ny, jd1, jd2, x, y);
    Grid_Computations(ib, jb, id1, id2, jd1, jd2, x, y, area, si, sj);

    //exit(0);

    if(dimen==0){
        InitFlowNonDimensional(id2, jd2, Re_inf, Machinf, Lref,  alpha, cv, dv);
    }
    else if(dimen==1) {
        InitFlowDimensional(id2, jd2, Re_inf, Machinf, Lref,  alpha, pinf, cv, dv);
    }
    else{
        std::cout<<"wrong input for dimensionality of the eqns"<<std::endl;
        exit(0);
    }

    if (restart==0){
        Res.open("../output/Residue_"+test_case+std::to_string(Nx)+"_"+std::to_string(Ny)+".dat");
        Res.flags(std::ios::dec | std::ios::scientific);
        Res.precision(10);

        if (!Res) {
            std::cerr << "File couldn't be opened to write the residue" << std::endl;
            exit(1);
        }
    }
    else if (restart==1){
        ReadSolnFileToRestart(RestartFile, Nx, id1, id2, Ny, jd1, jd2, dv, cv);
        Res.open("../output/Residue_Iter_"+test_case+std::to_string(restart_iter)+"_"+std::to_string(Nx)+"_"
                 +std::to_string(Ny)+".dat", std::ios_base::app);
        Res.flags(std::ios::dec | std::ios::scientific);
        Res.precision(10);

        if (!Res) {
            std::cerr << "File couldn't be opened to write the solution" << std::endl;
            exit(1);
        }
    }

    Dependent_Variables(id2, jd2, cv, dv);
    Boundary_Conditions(ib, id1, id2, jb, jd1, jd2, cv, dv, x, y, si, sj, pres_input_flag, bc_flag, bound_ind,
                        bound_cell, strt_bound_seg, end_bound_seg, pres_rhouvp) ;

    int iter;
    double time;

    if(restart==0){
        iter=0;
        time = 0.0;
    }
    else if(restart==1){
        iter = restart_iter;
        time = restart_t;
    }
    if(flow_type==0){
        while(iter<MaxIter){
            time=time+tstep[2][2];
            iter=iter+1;

            Solver(ib, id1, id2, jb, jd1, jd2, nconv, ndvar, x, y, cv, dv, dui, duj, area, si, sj, tstep, sri,  srj,
                   gradfi,  gradfj, cvold, diss, rhs, epsij);
            Residue_Cal(iter, tstep[2][2], time, cv, cvold, dv, si, sj);
            if(iter%1000==0){
                Write_Solution(id1, jd1, iter, time, x, y, cv);
                //Write_Surf_Solution(4, 1, 2, jd1, id1, jd1, iter, time, x, y, cv, dv, gradfi, gradfj); //fpr cylinder
                Write_Surf_Solution(1, 1, 10, id1, id1, jd1, iter, time, x, y, cv, dv, gradfi, gradfj);
                Interior_Solution(1, 100, 2, jd1, id1, jd1, iter, time, x, y, cv, dv, gradfi, gradfj);
                //WriteRestartFile(ib, jb, iter, time, dv);
            }
        }
    }
    else if (flow_type==1){
        int iter = 0;
        double t = 0, tmax = 1.958e-3;
        Write_Solution(id1, jd1, iter, time, x, y, cv);
        Write_Surf_Solution(4, 1, 2, jd1, id1, jd1, iter, time, x, y, cv, dv, gradfi, gradfj);
        Interior_Solution(1, 69, 2, jd1, id1, jd1, iter, time, x, y, cv, dv, gradfi, gradfj);
        while (t < tmax) {

            t = t + tstep[2][2];
            iter = iter + 1;
            std::cout << "iter number" << iter << "\t" << "current time=" << tstep[2][2] << "\t" << "total time=" << t
                      << std::endl;

            Solver(ib, id1, id2, jb, jd1, jd2, nconv, ndvar, x, y, cv, dv, dui, duj, area, si, sj, tstep, sri,  srj,
                   gradfi,  gradfj, cvold, diss, rhs, epsij);
            if(iter%1000==0){
                Write_Solution(id1, jd1, iter, time, x, y, cv);
                Write_Surf_Solution(1, 1, 10, 10, id1, jd1, iter, time, x, y, cv, dv, gradfi, gradfj);
                Interior_Solution(1, 69, 2, jd1, id1, jd1, iter, time, x, y, cv, dv, gradfi, gradfj);
                WriteRestartFile(ib, jb, iter, time, dv);
            }
        }
    }
    else{
        std::cout<<"wrong input for the flow type: it should be 0 for steady and 1 for unsteady"<<std::endl;
    }

    Write_Solution(id1, jd1, iter, time, x, y, cv);
    //Write_Surf_Solution(4, 1, 2, jd1, id1, jd1, iter, time, x, y, cv, dv, gradfi, gradfj);
    Write_Surf_Solution(1, 1, 2, id1, id1, jd1, iter, time, x, y, cv, dv, gradfi, gradfj);
    Interior_Solution(1, 56, 2, jd1, id1, jd1, iter, time, x, y, cv, dv, gradfi, gradfj);
    delete gradfi, gradfj;
    return 0;
}

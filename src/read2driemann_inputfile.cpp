//
// Created by maruthinh on 27/4/19.
//

#include "../inc/global_declarations.h"

void Read2DEulerRiemannProbInput(const std::string& RiemannSolverInputFile){

    std::string TitleLine, StartDelim, EndDelim, junk;

    std::ifstream infile;
    infile.open(RiemannSolverInputFile);

    std::cout<<"string solverInput="<<RiemannSolverInputFile<<std::endl;

    if (infile.fail()) {
        std::cerr << "Solver Input File couldn't be opened to read" << std::endl;
        std::exit(1);
    }

    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>test_case;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"Test case being solved..."<<test_case<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>flow_type;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if(flow_type==0) std::cout<<"solving steady state problem"<<std::endl;
    else std::cout<<"solving unsteady state problem"<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>time_accuracy;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if(time_accuracy==0) std::cout<<"using explicit time stepping method"<<std::endl;
    else if(time_accuracy==1) std::cout<<"using ssprk1 time stepping method"<<std::endl;
    else std::cout<<"using ssprk2 time stepping method"<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>inviscid_scheme;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if(inviscid_scheme==0) std::cout<<"using llf method"<<std::endl;
    else if(time_accuracy==1) std::cout<<"using movers method"<<std::endl;
    else std::cout<<"using kfds method"<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>space_accuracy;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"spatial accuracy is="<<space_accuracy<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>cfl;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"cfl is="<<cfl<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>restart;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>MeshFile;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"MeshFile is="<<MeshFile<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>Nx>>Ny;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"Grid size, Nx and Ny="<<Nx<<"\t"<<Ny<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>MeshTopFile;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"Mesh topology File is="<<MeshTopFile<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>RestartFile;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if(restart==0) std::cout<<"Starting solution afresh"<<std::endl;
    else std::cout<<"Restart File is="<<RestartFile<<std::endl;
    //parameter to read and control display frequency of results
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>disp_freq;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"display frequency is=="<<disp_freq<<std::endl;
    //parameter to read and control output frequency of results
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>outp_freq;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"output frequency is=="<<outp_freq<<std::endl;
    //parameter to control maximum iterations results
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>MaxIter;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"Maximum number of iterations are="<<MaxIter<<std::endl;
    //to read the maximum time at which unsteady flow will be stopped
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>tot_time;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"Time at which unsteady flow will be stopped="<<tot_time<<std::endl;
/******************Reading flow parameters****************************************/
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>interface_flag;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"Dir along which diaphragm is placed ="<<interface_flag<<std::endl;

    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>interface_ind;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"location at which diaphragm is placed ="<<interface_ind<<std::endl;

    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>rho_l;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"density to the left of interface ="<<rho_l<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>rho_r;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"density to the right of interface ="<<rho_r<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>u_l;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"u-velocity to the left of interface ="<<u_l<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>u_r;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"u-velocity to the right of interface ="<<u_r<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>v_l;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"v-velocity to the left of interface ="<<v_l<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>v_r;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"v-velocity to the right of interface ="<<v_l<<std::endl;

    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>p_l;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"pressure to the left of interface ="<<p_l<<std::endl;

    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>p_r;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"pressure to the right of interface ="<<p_r<<std::endl;
    std::cout<<TitleLine<<std::endl;
}
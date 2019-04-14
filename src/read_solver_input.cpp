//
// Created by Maruthi NH on 13-07-2018.
//


#include "../inc/global_declarations.h"

void ReadSolverInput(const std::string& SolverInputFile){

    std::string TitleLine, StartDelim, EndDelim, junk;

    std::ifstream infile;
    infile.open(SolverInputFile);

    std::cout<<"string solverInput="<<SolverInputFile<<std::endl;

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
    infile>>eqn_type;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if(eqn_type==0) std::cout<<"solving Euler equations"<<std::endl;
    else std::cout<<"solving Navier-Stokes equations"<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>dimen;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if(dimen==0) std::cout<<"solving in Non-dimensional form"<<std::endl;
    else std::cout<<"solving in dimensional form"<<std::endl;
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
    infile>>visc_method;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if(visc_method==0) std::cout<<"using Green-Gauss method"<<std::endl;
    else if(visc_method==1) std::cout<<"using Least Square method"<<std::endl;
    else std::cout<<"wrong input"<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>cfl;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"cfs is="<<cfl<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>scaling_factor;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if(scaling_factor==1) std::cout<<"Grid is not scaled"<<std::endl;
    else std::cout<<"grid is scaled with factor="<<scaling_factor<<std::endl;
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
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>MaxIter;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"Maximum number of iterations are="<<MaxIter<<std::endl;
/******************Reading flow parameters*******************************/
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>Re_inf;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"Flow Reynolds number is="<<Re_inf<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>Machinf;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"Flow Mach number is="<<Machinf<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>Lref;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"Reference length is="<<Lref<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>alpha;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"Inlet angle of attack is="<<alpha<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>pinf;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"Inlet pressure is="<<pinf<<std::endl;
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>tinf;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"Inlet temperature is="<<tinf<<std::endl;
    /******************Reading flow parameters for Non-dimensional simulations*******************************/
    if(dimen==0){
        getline(infile, TitleLine);
        std::cout<<TitleLine<<std::endl;
        getline(infile, TitleLine);
        std::cout<<TitleLine<<std::endl;
        infile>>pref;
        infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout<<"ref pressure is="<<pref<<std::endl;
        getline(infile, TitleLine);
        std::cout<<TitleLine<<std::endl;
        infile>>Tref;
        infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout<<"ref temperature is="<<Tref<<std::endl;
        getline(infile, TitleLine);
        std::cout<<TitleLine<<std::endl;
        infile>>rhoref;
        infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout<<"ref density is="<<rhoref<<std::endl;
        getline(infile, TitleLine);
        std::cout<<TitleLine<<std::endl;
        infile>>velref;
        infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout<<"ref velocity is="<<velref<<std::endl;
        getline(infile, TitleLine);
        std::cout<<TitleLine<<std::endl;
    }
    else{
        getline(infile, TitleLine);
        getline(infile, TitleLine);
        getline(infile, TitleLine);
        getline(infile, TitleLine);
        getline(infile, TitleLine);
        getline(infile, TitleLine);
        getline(infile, TitleLine);
        getline(infile, TitleLine);
        getline(infile, TitleLine);
        getline(infile, TitleLine);
        std::cout<<TitleLine<<std::endl;
    }
}

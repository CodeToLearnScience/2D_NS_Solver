//
// Created by maruthinh on 15/5/18.
//

#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

void ReadSolnFileToRestart(const std::string& RestartFile, int Nx, int id1, int id2, int Ny, int jd1, int jd2,
                           double ***&dv, double ***&cv) {
    int k;
    std::string TitleLine, StartDelim, EndDelim, junk;

    std::ifstream infile;

    infile.open("../output/"+RestartFile);

    if (infile.fail()) {
        std::cerr << "File couldn't be opened to read Restart file" << std::endl;
        std::exit(1);
    }
    getline(infile, TitleLine);
    std::cout<<"read string is="<<TitleLine<<std::endl;
    StartDelim = "Iter=";
    EndDelim = ",";
    restart_iter = std::stoi(get_str_between_two_str(TitleLine, StartDelim, EndDelim));

    StartDelim = "time=";
    EndDelim = "\"";
    restart_t = std::stod(get_str_between_two_str(TitleLine, StartDelim, EndDelim));

    std::cout<<"restart iter = "<<restart_iter<<"\t"<<"restart time = "<<restart_t<< std::endl;

    std::cout << "---------------------------------------------------------" << std::endl
              << "Reading solution to restart from the file" << std::endl
              << "---------------------------------------------------------" << std::endl;
    getline(infile, junk);
    std::cout<<junk<<std::endl;
    getline(infile, junk);
    std::cout<<junk<<std::endl;
    int count=0;
    //while (!infile.eof()) {
        for (int j = 2; j <= jd1; j++) {
            for (int i = 2; i <= id1; i++) {
                count++;
                infile >> junk >> junk >> dv[0][i][j] >> dv[1][i][j] >> dv[2][i][j] >> dv[3][i][j] >> junk >> junk;
                std::cout << "read values are=" << count << "\t" << dv[0][i][j] << "\t" << dv[1][i][j] << "\t"
                          << dv[2][i][j] << "\t" <<
                          dv[3][i][j] << std::endl;
                dv[4][i][j] = dv[3][i][j]/(dv[0][i][j]*Rgas);
                dv[5][i][j] = sqrt(Gamma*dv[3][i][j]/dv[0][i][j]); //a

                if(dimen==0){
                    dv[6][i][j] = ((1.0+C1/C0)/(C1/C0+dv[4][i][j]))*pow(dv[4][i][j],1.5)/Re_inf; //Mu
                    dv[7][i][j] = dv[6][i][j]/(GammaMinus1*Machinf*Machinf*Pr); //K
                }
                else if(dimen==1){
                    dv[6][i][j] = ((C0+C1)/(C1+dv[4][i][j]))*pow((dv[4][i][j]/C0),1.5)*refvisc; //Mu
                    dv[7][i][j] = dv[6][i][j]*Cp/Pr; //K
                }
            }
        }
    //}

    infile.close();

    //grid points for dummy cells

    for (k=0;k<=ndvar;k++){
        for (int i = 2; i <= id1; i++) {
            dv[k][i][0]   = dv[k][i][2];
            dv[k][i][1]   = dv[k][i][2];
            //dv[k][i][jd1] = dv[k][i][jb];
            dv[k][i][jd2] = dv[k][i][jd1];
        }

        for (int j = 0; j <= jd2; j++) {
            dv[k][0][j]   = dv[k][2][j];
            dv[k][1][j]   = dv[k][2][j];
            //dv[k][id1][j] = dv[k][ib][j];
            dv[k][id2][j] = dv[k][id1][j];
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

    /*for (int j = 0; j <= jd2; j++) {
        for (int i = 0; i <= id2; i++) {

            //std::cout << k << "\t" << x[i][j] << "\t" << y[i][j] << std::endl;
            std::cout<<"i="<< i<<"\t"<<"j="<<j<<"\t"<< x[i][j] << "\t" << y[i][j] << std::endl;
        }
    }*/
}




#include "global_declarations.h"

void
Ic_Unsteady_Flow(int id1, int id2, int jd1, int jd2, int index, int loc, double rhol, double rhor, double ul, double ur,
                 double vl, double vr, double pl, double pr, double **&x, double **&y, double ***&cv) {

    double El, Er;

    El = pl / (Gamma - 1.0) + 0.5 * rhol * (ul * ul + vl * vl);
    Er = pr / (Gamma - 1.0) + 0.5 * rhor * (ur * ur + vr * vr);

    if (index == 1) {
        for (int j = 2; j <= jd1; j++) {
            for (int i = 2; i <= loc; i++) {

                cv[0][i][j] = rhol;
                cv[1][i][j] = rhol * ul;
                cv[2][i][j] = rhol * vl;
                cv[3][i][j] = El;
            }
        }

        for (int j = 2; j <= jd1; j++) {
            for (int i = loc + 1; i <= id1; i++) {

                cv[0][i][j] = rhor;
                cv[1][i][j] = rhor * ur;
                cv[2][i][j] = rhor * vr;
                cv[3][i][j] = Er;
            }
        }
    } else if (index == 2) {

        for (int j = 2; j <= loc; j++) {
            for (int i = 2; i <= id1; i++) {

                cv[0][i][j] = rhol;
                cv[1][i][j] = rhol * ul;
                cv[2][i][j] = rhol * vl;
                cv[3][i][j] = El;
            }
        }

        for (int j = loc + 1; j <= jd1; j++) {
            for (int i = 2; i <= id1; i++) {

                cv[0][i][j] = rhor;
                cv[1][i][j] = rhor * ur;
                cv[2][i][j] = rhor * vr;
                cv[3][i][j] = Er;
            }
        }
    }


    for (int i = 2; i <= id1; i++) {

        cv[0][i][0] = cv[0][i][2];
        cv[1][i][0] = cv[1][i][2];
        cv[2][i][0] = cv[2][i][2];
        cv[3][i][0] = cv[3][i][2];

        cv[0][i][1] = cv[0][i][2];
        cv[1][i][1] = cv[1][i][2];
        cv[2][i][1] = cv[2][i][2];
        cv[3][i][1] = cv[3][i][2];

        cv[0][i][jd2] = cv[0][i][jd1];
        cv[1][i][jd2] = cv[1][i][jd1];
        cv[2][i][jd2] = cv[2][i][jd1];
        cv[3][i][jd2] = cv[3][i][jd1];

    }

    for (int j = 0; j <= jd2; j++) {

        cv[0][0][j] = cv[0][2][j];
        cv[1][0][j] = cv[1][2][j];
        cv[2][0][j] = cv[2][2][j];
        cv[3][0][j] = cv[3][2][j];

        cv[0][1][j] = cv[0][2][j];
        cv[1][1][j] = cv[1][2][j];
        cv[2][1][j] = cv[2][2][j];
        cv[3][1][j] = cv[3][2][j];

        cv[0][id2][j] = cv[0][id1][j];
        cv[1][id2][j] = cv[1][id1][j];
        cv[2][id2][j] = cv[2][id1][j];
        cv[3][id2][j] = cv[3][id1][j];
    }

    for (int j = 0; j <= jd2; j++) {
        //for(int i=0; i<=id2; i++){

        //std::cout<<"the initialized conserved variables are="<<"\t"<<cv[0][i][j]<<"\t"<<cv[1][i][j]<<"\t"<<cv[2][i][j]<<"\t"<<cv[3][i][j]<<std::endl;
        std::cout << "the initialized conserved variables are=" << "\t" << cv[0][23][j] << "\t" << cv[1][22][j] << "\t"
                  << cv[2][22][j] << "\t" << cv[3][22][j] << std::endl;

        //}
    }

}






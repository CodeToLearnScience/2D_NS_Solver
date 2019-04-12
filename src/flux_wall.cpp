#include "global_declarations.h"
#include "basic_functions.h"

void Flux_EulerWall(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv, double ***&si, double ***&sj,
               double ***&rhs) {

    double sx, sy, pa;
    int dum1, ins1, ins2;

    if (bind == 1 or bind == 4) {

        dum1 = bnode;
        ins1 = bnode + 1;
        ins2 = bnode + 2;
    } else if (bind == 2 or bind == 3) {

        dum1 = bnode;
        ins1 = bnode - 1;
        ins2 = bnode - 2;

    }else {
        std::cout << "Flux at wall: please check the boundary index, 1=bottom, 2=right, 3=top, 4=left and bind is =" <<bind<<std::endl;
        exit(0);
    }

    if (bind == 1 or bind == 3) {
        for (int i = sbind; i <= ebind; i++) {
            if (bind == 1) {
                sx = sj[0][i][ins1];
                sy = sj[1][i][ins1];
            } else {
                sx = -sj[0][i][dum1];
                sy = -sj[1][i][dum1];
            }

            pa = 0.5 * (3.0 * dv[3][i][ins1] - dv[3][i][ins2]);

            rhs[1][i][ins1] = rhs[1][i][ins1] + pa * sx;
            rhs[2][i][ins1] = rhs[2][i][ins1] + pa * sy;
        }
    } else if (bind == 2 or bind == 4) {
        for (int j = sbind; j <= ebind; j++) {
            if (bind == 4) {
                sx = si[0][ins1][j];
                sy = si[1][ins1][j];
            } else {
                sx = -si[0][dum1][j];
                sy = -si[1][dum1][j];
            }

            pa = 0.5 * (3.0 * dv[3][ins1][j] - dv[3][ins2][j]);

            rhs[1][ins1][j] = rhs[1][ins1][j] + pa * sx;
            rhs[2][ins1][j] = rhs[2][ins1][j] + pa * sy;
        }

        if (pa != pa) {
            std::cout << "In fun Flux_wall: Value of p bacame negative at the wall i and j are" << bind << "\t" << bnode << "\t" << pa
                      << std::endl;
            exit(0);
        }
    }
}

void Flux_ViscWall(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv, double ***&si, double ***&sj,
                    double ***&rhs) {

    double sx, sy, pa;
    int dum1, ins1, ins2;

    if (bind == 1 or bind == 4) {

        dum1 = bnode;
        ins1 = bnode + 1;
        ins2 = bnode + 2;
    } else if (bind == 2 or bind == 3) {

        dum1 = bnode;
        ins1 = bnode - 1;
        ins2 = bnode - 2;

    }else {
        std::cout << "Flux at wall: please check the boundary index, 1=bottom, 2=right, 3=top, 4=left and bind is =" <<bind<<std::endl;
        exit(0);
    }

    if (bind == 1 or bind == 3) {
        for (int i = sbind; i <= ebind; i++) {
            if (bind == 1) {
                sx = sj[0][i][ins1];
                sy = sj[1][i][ins1];
            } else {
                sx = -sj[0][i][dum1];
                sy = -sj[1][i][dum1];
            }

            pa = dv[3][i][ins1];

            rhs[1][i][ins1] = rhs[1][i][ins1] + pa * sx;
            rhs[2][i][ins1] = rhs[2][i][ins1] + pa * sy;
        }
    } else if (bind == 2 or bind == 4) {
        for (int j = sbind; j <= ebind; j++) {
            if (bind == 4) {
                sx = si[0][ins1][j];
                sy = si[1][ins1][j];
            } else {
                sx = -si[0][dum1][j];
                sy = -si[1][dum1][j];
            }

            pa = dv[3][ins1][j];

            rhs[1][ins1][j] = rhs[1][ins1][j] + pa * sx;
            rhs[2][ins1][j] = rhs[2][ins1][j] + pa * sy;
        }

        if (pa != pa) {
            std::cout << "In fun Flux_wall: Value of p bacame negative at the wall i and j are" << bind << "\t" << bnode << "\t" << pa
                      << std::endl;
            exit(0);
        }
    }
}

void Flux_Symmetry(int bind, int bnode, int sbind, int ebind, double ***&cv, double ***&dv, double ***&si, double ***&sj,
                    double ***&rhs) {

    double sx, sy, pa;
    int dum1, ins1, ins2;

    if (bind == 1 or bind == 4) {

        dum1 = bnode;
        ins1 = bnode + 1;
        ins2 = bnode + 2;
    } else if (bind == 2 or bind == 3) {

        dum1 = bnode;
        ins1 = bnode - 1;
        ins2 = bnode - 2;

    }else {
        std::cout << "Flux at wall: please check the boundary index, 1=bottom, 2=right, 3=top, 4=left and bind is =" <<bind<<std::endl;
        exit(0);
    }

    if (bind == 1 or bind == 3) {
        for (int i = sbind; i <= ebind; i++) {
            if (bind == 1) {
                sx = sj[0][i][ins1];
                sy = sj[1][i][ins1];
            } else {
                sx = -sj[0][i][dum1];
                sy = -sj[1][i][dum1];
            }

            pa = 0.5 * (dv[3][i][ins1] + dv[3][i][dum1]);

            rhs[1][i][ins1] = rhs[1][i][ins1] + pa * sx;
            rhs[2][i][ins1] = rhs[2][i][ins1] + pa * sy;
        }
    } else if (bind == 2 or bind == 4) {
        for (int j = sbind; j <= ebind; j++) {
            if (bind == 4) {
                sx = si[0][ins1][j];
                sy = si[1][ins1][j];
            } else {
                sx = -si[0][dum1][j];
                sy = -si[1][dum1][j];
            }

            pa = 0.5 * (dv[3][ins1][j] + dv[3][dum1][j]);

            rhs[1][ins1][j] = rhs[1][ins1][j] + pa * sx;
            rhs[2][ins1][j] = rhs[2][ins1][j] + pa * sy;
        }

        if (pa != pa) {
            std::cout << "In fun Flux_wall: Value of p bacame negative at the wall i and j are" << bind << "\t" << bnode << "\t" << pa
                      << std::endl;
            exit(0);
        }
    }
}
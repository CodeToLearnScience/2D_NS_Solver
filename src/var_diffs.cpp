#include "global_declarations.h"
#include "basic_functions.h"


void Prim_Variables_Differences(int ib, int id2, int jb, int jd2, double ***cv, double ***dv, double ***&dui,
                                double ***&duj) {

    for (int j = 2; j <= jb; j++) {
        for (int i = 1; i <= id2; i++) {

            dui[0][i][j] = dv[0][i][j] - dv[0][i - 1][j];
            dui[1][i][j] = dv[1][i][j] - dv[1][i - 1][j];
            dui[2][i][j] = dv[2][i][j] - dv[2][i - 1][j];
            dui[3][i][j] = dv[3][i][j] - dv[3][i - 1][j];
        }

        dui[0][0][j] = dui[0][1][j];
        dui[1][0][j] = dui[1][1][j];
        dui[2][0][j] = dui[2][1][j];
        dui[3][0][j] = dui[3][1][j];
    }

    for (int i = 2; i <= ib; i++) {
        for (int j = 1; j <= jd2; j++) {

            duj[0][i][j] = dv[0][i][j] - dv[0][i][j - 1];
            duj[1][i][j] = dv[1][i][j] - dv[1][i][j - 1];
            duj[2][i][j] = dv[2][i][j] - dv[2][i][j - 1];
            duj[3][i][j] = dv[3][i][j] - dv[3][i][j - 1];
        }

        duj[0][i][0] = duj[0][i][1];
        duj[1][i][0] = duj[1][i][1];
        duj[2][i][0] = duj[2][i][1];
        duj[3][i][0] = duj[3][i][1];
    }
}

void Conv_Variables_Differences(int ib, int id2, int jb, int jd2, double ***cv, double ***dv, double ***&dui,
                                double ***&duj) {

    for (int j = 2; j <= jb; j++) {
        for (int i = 1; i <= id2; i++) {

            cui[0][i][j] = cv[0][i][j] - cv[0][i - 1][j];
            cui[1][i][j] = cv[1][i][j] - cv[1][i - 1][j];
            cui[2][i][j] = cv[2][i][j] - cv[2][i - 1][j];
            cui[3][i][j] = cv[3][i][j] - cv[3][i - 1][j];
        }

        cui[0][0][j] = cui[0][1][j];
        cui[1][0][j] = cui[1][1][j];
        cui[2][0][j] = cui[2][1][j];
        cui[3][0][j] = cui[3][1][j];
    }

    for (int i = 2; i <= ib; i++) {
        for (int j = 1; j <= jd2; j++) {

            cuj[0][i][j] = cv[0][i][j] - cv[0][i][j - 1];
            cuj[1][i][j] = cv[1][i][j] - cv[1][i][j - 1];
            cuj[2][i][j] = cv[2][i][j] - cv[2][i][j - 1];
            cuj[3][i][j] = cv[3][i][j] - cv[3][i][j - 1];
        }

        cuj[0][i][0] = cuj[0][i][1];
        cuj[1][i][0] = cuj[1][i][1];
        cuj[2][i][0] = cuj[2][i][1];
        cuj[3][i][0] = cuj[3][i][1];
    }
}

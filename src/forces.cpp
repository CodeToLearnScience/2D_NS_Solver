#include "global_declarations.h"
#include "basic_functions.h"

void Forces(int bind, int bnode, int sbind, int ebind, int ndvar, double **&x, double **&y, double ***&dv,
            double ***&si, double ***&sj) {

    double cx, cy, sx, sy, pwall, xa, ya, dcy, dcx;
    double xref, yref, cref;

    cx = 0;
    cy = 0;
    cm = 0;
    xref = 0.25;
    yref = 0.0;
    cref = 1.0;
    int dum1, dum2, ins1, ins2;
    if (bind == 1 or bind == 4) {
        dum1 = bnode;
        dum2 = bnode - 1;
        ins1 = bnode + 1;
        ins2 = bnode + 2;
    } else if (bind == 2 or bind == 3) {
        dum1 = bnode;
        dum2 = bnode + 1;
        ins1 = bnode - 1;
        ins2 = bnode - 2;
    } else {
        std::cout << "Inflo b.c., please check the boundary index, 1=bottom, 2=right, 3=top, 4=left" << std::endl;
        exit(0);
    }

    if (bind == 1 or bind == 3) {

        //then sbind, ebind will be i

        for (int i = sbind; i <= ebind; i++) {

            if (bind == 1) {

                sx = sj[0][i][ins1];
                sy = sj[1][i][ins1];
            } else {
                sx = -sj[0][i][dum1];
                sy = -sj[1][i][dum1];
            }

            pwall = 0.5 * (dv[3][i][dum1] + dv[3][i][ins1]);
            cp = 2.0 * (pwall - pinf) / (rhoinf * qinf * qinf);
            xa = (0.5 * (x[i][ins1] + x[i + 1][ins1]) - xref) / cref;
            ya = (0.5 * (y[i][ins1] + y[i + 1][ins1]) - yref) / cref;
            dcy = sy * cp;
            dcx = sx * cp;
            cy = cy + dcy;
            cx = cx + dcx;
            cm = cm + dcx * ya - dcy * xa;

        }
    }
    if (bind == 2 or bind == 4  ) {

        //then sbind, ebind will be j

        for (int j = sbind; j <= ebind; j++) {

            if (bind == 4) {

                sx = si[0][ins1][j];
                sy = si[1][ins1][j];
            } else {
                sx = -si[0][dum1][j];
                sy = -si[1][dum1][j];
            }

            pwall = 0.5 * (dv[3][dum1][j] + dv[3][ins1][j]);
            cp = 2.0 * (pwall - pinf) / (rhoinf * qinf * qinf);
            xa = (0.5 * (x[ins1][j] + x[ins1][j + 1]) - xref) / cref;
            ya = (0.5 * (y[ins1][j] + y[ins1][j + 1]) - yref) / cref;
            dcy = sy * cp;
            dcx = sx * cp;
            cy = cy + dcy;
            cx = cx + dcx;
            cm = cm + dcx * ya - dcy * xa;
        }
    }

    cl = cy * cos(alpha) - cx * sin(alpha);
    cd = cy * sin(alpha) - cx * cos(alpha);
    //std::cout << "the values of cl and cd in forces" << "\t" << cl << "\t" << cd << "\t" << std::endl;
}


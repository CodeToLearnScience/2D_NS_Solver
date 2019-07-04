#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"

void Write_SolutionVtk(int id1, int jd1, int iter, double t, double **&x, double **&y, double ***&cv)
{

    double rho, u, v, p, T, a, Mach, mav, kav, tauxx, tauyy, tauxy, two_by_3, phix, phiy, sxx, syx, dsx, nxx, nyx, sxy,
        syy, dsy, nxy, nyy;

    std::string dataType = "double";

    std::string oup;
    oup = "../output/";
    std::ofstream other(oup + test_case + "_IsoCont" + std::to_string(Nx) + "_" + std::to_string(Ny) + "Iter" + std::to_string(iter) +
                        ".vtk");
    other.flags(std::ios::dec | std::ios::scientific);
    other.precision(16);

    if (!other)
    {
        std::cerr << "File couldn't be opened to write the solution" << std::endl;
        exit(1);
    }

    int Nz = 1;
    other << "# vtk DataFile Version 3.0" << std::endl;
    other << "Field at"<< " "<< "T="<< ToStringWithPrecision(t, 16)<< std::endl;
    other << "ASCII" << std::endl;
    other << "DATASET STRUCTURED_GRID" << std::endl;
    other << "DIMENSIONS"
          << " " << Nx + 1 << " " << Ny + 1 << " " << 1 << std::endl;
    other << "POINTS"
          << " " << (Nx + 1) * (Ny + 1) * Nz << " "
          << "double" << std::endl;

    for (int j = 2; j <= jd1; j++)
    {
        for (int i = 2; i <= id1; i++)
        {
            other << x[i][j] << "\t" << y[i][j] << "\t" << 1 << std::endl;
        }
    }
    other << "POINT_DATA"
          << " " << (Nx + 1) * (Ny + 1) * Nz << std::endl;
    other << "SCALARS density"
          << " "
          << dataType << std::endl;
    other << "LOOKUP_TABLE default" << std::endl;

    for (int j = 2; j <= jd1; j++)
    {
        for (int i = 2; i <= id1; i++)
        {
            rho = 0.25 * (dv[0][i][j] + dv[0][i - 1][j] + dv[0][i - 1][j - 1] + dv[0][i][j - 1]);
            u = 0.25 * (dv[1][i][j] + dv[1][i - 1][j] + dv[1][i - 1][j - 1] + dv[1][i][j - 1]);
            v = 0.25 * (dv[2][i][j] + dv[2][i - 1][j] + dv[2][i - 1][j - 1] + dv[2][i][j - 1]);
            p = 0.25 * (dv[3][i][j] + dv[3][i - 1][j] + dv[3][i - 1][j - 1] + dv[3][i][j - 1]);
            T = 0.25 * (dv[4][i][j] + dv[4][i - 1][j] + dv[4][i - 1][j - 1] + dv[4][i][j - 1]);
            mav = 0.25 * (dv[6][i][j] + dv[6][i - 1][j] + dv[6][i - 1][j - 1] + dv[6][i][j - 1]);
            kav = 0.25 * (dv[7][i][j] + dv[7][i - 1][j] + dv[7][i - 1][j - 1] + dv[7][i][j - 1]);

            Mach = sqrt(u * u + v * v) / a;

            other << rho << std::endl;
        }
    }

    other << "SCALARS pressure"
          << " "
          << dataType << std::endl;
    other << "LOOKUP_TABLE default" << std::endl;

    for (int j = 2; j <= jd1; j++)
    {
        for (int i = 2; i <= id1; i++)
        {
            rho = 0.25 * (dv[0][i][j] + dv[0][i - 1][j] + dv[0][i - 1][j - 1] + dv[0][i][j - 1]);
            u = 0.25 * (dv[1][i][j] + dv[1][i - 1][j] + dv[1][i - 1][j - 1] + dv[1][i][j - 1]);
            v = 0.25 * (dv[2][i][j] + dv[2][i - 1][j] + dv[2][i - 1][j - 1] + dv[2][i][j - 1]);
            p = 0.25 * (dv[3][i][j] + dv[3][i - 1][j] + dv[3][i - 1][j - 1] + dv[3][i][j - 1]);
            T = 0.25 * (dv[4][i][j] + dv[4][i - 1][j] + dv[4][i - 1][j - 1] + dv[4][i][j - 1]);
            mav = 0.25 * (dv[6][i][j] + dv[6][i - 1][j] + dv[6][i - 1][j - 1] + dv[6][i][j - 1]);
            kav = 0.25 * (dv[7][i][j] + dv[7][i - 1][j] + dv[7][i - 1][j - 1] + dv[7][i][j - 1]);

            Mach = sqrt(u * u + v * v) / a;

            other << p << std::endl;
        }
    }

    other << "VECTORS velocity"
          << " " << dataType << std::endl;

    for (int j = 2; j <= jd1; j++)
    {
        for (int i = 2; i <= id1; i++)
        {
            rho = 0.25 * (dv[0][i][j] + dv[0][i - 1][j] + dv[0][i - 1][j - 1] + dv[0][i][j - 1]);
            u = 0.25 * (dv[1][i][j] + dv[1][i - 1][j] + dv[1][i - 1][j - 1] + dv[1][i][j - 1]);
            v = 0.25 * (dv[2][i][j] + dv[2][i - 1][j] + dv[2][i - 1][j - 1] + dv[2][i][j - 1]);
            p = 0.25 * (dv[3][i][j] + dv[3][i - 1][j] + dv[3][i - 1][j - 1] + dv[3][i][j - 1]);
            T = 0.25 * (dv[4][i][j] + dv[4][i - 1][j] + dv[4][i - 1][j - 1] + dv[4][i][j - 1]);
            mav = 0.25 * (dv[6][i][j] + dv[6][i - 1][j] + dv[6][i - 1][j - 1] + dv[6][i][j - 1]);
            kav = 0.25 * (dv[7][i][j] + dv[7][i - 1][j] + dv[7][i - 1][j - 1] + dv[7][i][j - 1]);

            Mach = sqrt(u * u + v * v) / a;

            other << u << " " << v << " " << 0.0 << std::endl;
        }
    }

    other.close();
}

#include "global_declarations.h"
#include "basic_functions.h"

void Write_Solution(int id1, int jd1, int iter, double t, double **&x, double **&y, double ***&cv) {

    double rho, u, v, p, T, a, Mach, mav, kav, tauxx, tauyy, tauxy, two_by_3, phix, phiy, sxx, syx, dsx, nxx, nyx, sxy,
            syy, dsy, nxy, nyy;

    std::string oup;
    oup = "../output_kfds/";
    std::ofstream other(oup+test_case+"_IsoCont"+std::to_string(Nx)+"_"+std::to_string(Ny)+"Iter"+std::to_string(iter)+
                        ".dat");
    other.flags(std::ios::dec | std::ios::scientific);
    other.precision(16);

    if (!other) {
        std::cerr << "File couldn't be opened to write the solution" << std::endl;
        exit(1);
    }

    /*other << "TITLE = \"Flow\"" << std::endl << "VARIABLES = \"xc\", \"yc\", \"rho\", \"u\", \"v\", \"p\", \"T\", \"Mach\", \"mav\""
            ", \"kav\", \"nxx\", \"nyx\", \"nxy\", \"nyy\",\"tauxx\", \"tauyy\", \"tauxy\", \"phix\", \"phiy\"" << std::endl;*/

    other << "TITLE = "<<"\"Iter="<<iter<<", "<<"time="<<ToStringWithPrecision(t, 16)<<"\""<< std::endl << "VARIABLES = "
            "\"xc\", \"yc\", \"rho\", \"u\", \"v\", \"p\", \"T\", ""\"Mach\"" << std::endl;
    other << "Zone"<<" "<<"T="<<"\""<<ToStringWithPrecision(t, 16)<<"\""<<","<<" "<<"I=" << Nx + 1 <<","<<" "<<"J="
          << Ny + 1 << std::endl;

    for (int j = 2; j <= jd1; j++) {
        for (int i = 2; i <= id1; i++) {
           rho  = 0.25*(dv[0][i][j]+dv[0][i-1][j]+dv[0][i-1][j-1]+dv[0][i][j-1]);
            u   = 0.25*(dv[1][i][j]+dv[1][i-1][j]+dv[1][i-1][j-1]+dv[1][i][j-1]);
            v   = 0.25*(dv[2][i][j]+dv[2][i-1][j]+dv[2][i-1][j-1]+dv[2][i][j-1]);
            p   = 0.25*(dv[3][i][j]+dv[3][i-1][j]+dv[3][i-1][j-1]+dv[3][i][j-1]);
            T   = 0.25*(dv[4][i][j]+dv[4][i-1][j]+dv[4][i-1][j-1]+dv[4][i][j-1]);
            a   = 0.25*(dv[5][i][j]+dv[5][i-1][j]+dv[5][i-1][j-1]+dv[5][i][j-1]);
            //mav = 0.25*(dv[6][i][j]+dv[6][i-1][j]+dv[6][i-1][j-1]+dv[6][i][j-1]);
            //kav = 0.25*(dv[7][i][j]+dv[7][i-1][j]+dv[7][i-1][j-1]+dv[7][i][j-1]);

            Mach = sqrt(u*u+v*v)/a;

            /*sxx = si[0][i][j];
            syx = si[1][i][j];
            dsx = sqrt(sxx*sxx+syx*syx);
            nxx = sxx/dsx;
            nyx = syx/dsx;
            sxy = sj[0][i][j];
            syy = sj[1][i][j];
            dsy = sqrt(sxy*sxy+syy*syy);
            nxy = sxy/dsy;
            nyy = syy/dsy;
            tauxx = two_by_3*(mav/Re)*(2.0*gradfj[0][i][j]-gradfj[3][i][j]);
            tauyy = two_by_3*(mav/Re)*(2.0*gradfj[3][i][j]-gradfj[0][i][j]);
            tauxy = (mav/Re)*(gradfj[1][i][j]+gradfj[2][i][j]);
            phix = u*tauxx + v*tauxy + kav*gradfj[4][i][j];
            phiy = u*tauxy + v*tauyy + kav*gradfj[5][i][j];*/


            /*rho = cv[0][i][j];
            u = (cv[1][i][j] / cv[0][i][j]);
            v = (cv[2][i][j] / cv[0][i][j]);
            p = (dv[0][i][j]);
            a = (dv[2][i][j]);
            Mach = sqrt(u * u + v * v) / a;*/

            /*other << x[i][j]<<"\t"<<y[i][j]<<"\t"<<rho<<"\t"<<u<<"\t"<<v<<"\t"<<p<<"\t"<<T<<"\t"<<Mach<<"\t"<<mav<<"\t"<<kav
                  <<"\t"<<nxx<<"\t"<<nyx<<"\t"<<nxy<<"\t"<<nyy<<"\t"<<tauxx<<"\t"<<tauyy<<"\t"<<tauxy<<"\t"<<phix<<"\t"
                  << phiy<<std::endl;*/

            other << x[i][j]<<"\t"<<y[i][j]<<"\t"<<rho<<"\t"<<u<<"\t"<<v<<"\t"<<p<<"\t"<<T<<"\t"<<Mach<<std::endl;

        }
    }

    other.close();
}

void WriteRestartFile(int ib, int jb, int iter, double t, double ***&dv) {

    std::string oup;
    oup = "../output_kfds/";
    //std::ofstream other(oup+"RestartFile"+"Iter"+std::to_string(iter)+"_"+"t"+std::to_string(t)+".dat");
    std::ofstream other(oup+test_case+"RestartFile"+std::to_string(Nx)+"_"+std::to_string(Ny)+".dat");
    other.flags(std::ios::dec | std::ios::scientific);
    other.precision(16);

    if (!other) {
        std::cerr << "File couldn't be opened to write the solution" << std::endl;
        exit(1);
    }
    other << "iter"<<"\t"<<iter<<std::endl;
    other << "t"<<"\t"<<t<<std::endl;

    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= ib; i++) {

            other << dv[0][i][j] <<"\t"<<dv[1][i][j]<<"\t"<<dv[2][i][j]<<"\t"<<dv[3][i][j]<<std::endl;

        }
    }
    other.close();
}

void Write_Surf_Solution(int bind, int sur_ind, int beg_seg, int end_seg, int id1, int jd1, int iter, double t,
                         double **&x, double **&y, double ***&cv, double ***&dv, double ***&gradfi, double ***&gradfj) {

    double rho, u, v, p, a, Mach, T, Cf, mav, kav, *gradavg, sx, sy, ds, nx, ny, gradnx, gradny, gradnn, dvdnx,
            dvdny, sgn, dvdna;
    int ip1, jp1, id, jd;
    int counter;
    gradavg = new double[4];

    std::string oup;
    oup = "../output_kfds/";
    std::ofstream sur(oup + test_case+"_Sur_"+std::to_string(Nx)+"_"+std::to_string(Ny)+"Iter"+std::to_string(iter)
                      +".dat");
    sur.flags(std::ios::dec | std::ios::scientific);
    sur.precision(16);

    if (!sur) {
        std::cerr << "File couldn't be opened to write the solution" << std::endl;
        exit(1);
    }

    sur << "TITLE = "<<"\"Iter="<<iter<<", "<<"time="<<ToStringWithPrecision(t, 16)<<"\""<< std::endl
        << "VARIABLES = counter, xc, yc, rho, u, v, p, T, a, Mach, visc, k, Cf" << std::endl;

    if(bind==1){
        sur << "Zone"<<" "<<"T="<<"\""<<ToStringWithPrecision(t, 16)<<"\""<<","<<" "<<"I=" << ib <<","<<"\t"
            << " F = POINT" << std::endl;
        counter=0;
        for (int i = beg_seg; i <= end_seg; i++) {
            int j = sur_ind;
            counter++;

            if(i==id1) id=ib;
            else id = i;
            /*rho = 0.25 * (dv[0][i][j] + dv[0][i + 1][j] + dv[0][i + 1][j + 1] + dv[0][i][j + 1]);
            u   = 0.25 * (dv[1][i][j] + dv[1][i + 1][j] + dv[1][i + 1][j + 1] + dv[1][i][j + 1]);
            v   = 0.25 * (dv[2][i][j] + dv[2][i + 1][j] + dv[2][i + 1][j + 1] + dv[2][i][j + 1]);
            p   = 0.25 * (dv[3][i][j] + dv[3][i + 1][j] + dv[3][i + 1][j + 1] + dv[3][i][j + 1]);
            T   = 0.25 * (dv[4][i][j] + dv[4][i + 1][j] + dv[4][i + 1][j + 1] + dv[4][i][j + 1]);
            a   = 0.25 * (dv[5][i][j] + dv[5][i + 1][j] + dv[5][i + 1][j + 1] + dv[5][i][j + 1]);
            mav = 0.25 * (dv[6][i][j] + dv[6][i + 1][j] + dv[6][i + 1][j + 1] + dv[6][i][j + 1]);
            kav = 0.25 * (dv[7][i][j] + dv[7][i + 1][j] + dv[7][i + 1][j + 1] + dv[7][i][j + 1]);*/

            rho = 0.25 * (dv[0][id][j] + dv[0][id][j + 1] + dv[0][id - 1][j + 1] + dv[0][id - 1][j]);
            u   = 0.25 * (dv[1][id][j] + dv[1][id][j + 1] + dv[1][id - 1][j + 1] + dv[1][id - 1][j]);
            v   = 0.25 * (dv[2][id][j] + dv[2][id][j + 1] + dv[2][id - 1][j + 1] + dv[2][id - 1][j]);
            p   = 0.25 * (dv[3][id][j] + dv[3][id][j + 1] + dv[3][id - 1][j + 1] + dv[3][id - 1][j]);
            T   = 0.25 * (dv[4][id][j] + dv[4][id][j + 1] + dv[4][id - 1][j + 1] + dv[4][id - 1][j]);
            a   = 0.25 * (dv[5][id][j] + dv[5][id][j + 1] + dv[5][id - 1][j + 1] + dv[5][id - 1][j]);
            mav = 0.25 * (dv[6][id][j] + dv[6][id][j + 1] + dv[6][id - 1][j + 1] + dv[6][id - 1][j]);
            kav = 0.25 * (dv[7][id][j] + dv[7][id][j + 1] + dv[7][id - 1][j + 1] + dv[7][id - 1][j]);

            /*//to compute Cf
            if(i==ib) ip1=ib;
            else ip1 = i + 1;

            gradavg[0] = 0.5 * (gradfj[0][ip1][j] + gradfj[0][ip1][j + 1]);
            gradavg[1] = 0.5 * (gradfj[1][ip1][j] + gradfj[1][ip1][j + 1]);
            gradavg[2] = 0.5 * (gradfj[2][ip1][j] + gradfj[2][ip1][j + 1]);
            gradavg[3] = 0.5 * (gradfj[3][ip1][j] + gradfj[3][ip1][j + 1]);*/

            //to compute Cf
            gradavg[0] = 0.5 * (gradfj[0][id][j] + gradfj[0][id][j + 1]);
            gradavg[1] = 0.5 * (gradfj[1][id][j] + gradfj[1][id][j + 1]);
            gradavg[2] = 0.5 * (gradfj[2][id][j] + gradfj[2][id][j + 1]);
            gradavg[3] = 0.5 * (gradfj[3][id][j] + gradfj[3][id][j + 1]);

            /*sx = -0.5 * (sj[0][ip1][j] + sj[0][ip1][j + 1]);
            sy = -0.5 * (sj[1][ip1][j] + sj[1][ip1][j + 1]);
            ds = sqrt(sx * sx + sy * sy);
            nx = sx / ds;
            ny = sy / ds;
            
            gradnx = gradavg[0] * nx + gradavg[1] * ny;
            gradny = gradavg[2] * nx + gradavg[3] * ny;
            gradnn = gradnx * nx + gradny * ny;
            dvdnx  = gradnx - gradnn * nx;
            dvdny  = gradny - gradnn * ny;
            sgn    = Sign(gradnx);
            dvdna  = sqrt(dvdnx * dvdnx + dvdny * dvdny);
            Cf     = 2.0 * sgn * mav * dvdna / (rhoinf * qinf * qinf);*/

            sx = -0.5 * (sj[0][id][j] + sj[0][id][j + 1]);
            sy = -0.5 * (sj[1][id][j] + sj[1][id][j + 1]);
            ds = sqrt(sx * sx + sy * sy);
            nx = sx / ds;
            ny = sy / ds;

            gradnx = gradavg[0] * nx + gradavg[1] * ny;
            gradny = gradavg[2] * nx + gradavg[3] * ny;
            gradnn = gradnx * nx + gradny * ny;
            dvdnx  = gradnx - gradnn * nx;
            dvdny  = gradny - gradnn * ny;
            sgn    = Sign(gradnx);
            dvdna  = sqrt(dvdnx * dvdnx + dvdny * dvdny);
            Cf     = 2.0 * sgn * mav * dvdna / (rhoinf * qinf * qinf);

            sur << counter << "\t" << x[i][j] << "\t" << y[i][j] << "\t" << rho << "\t" << u << "\t" << v << "\t"
                << p << "\t" << T << "\t" << a << "\t" << Mach << "\t"
                << mav << "\t" << kav << "\t" << Cf << "\t" << std::endl;
        }
    }
    else if(bind==2){
        counter=0;
        sur << "Zone T = \"2\"," <<"I = " << jb <<","<<"\t"<< " F = POINT" << std::endl;

        for (int j = beg_seg; j <= end_seg; j++) {
            int i = sur_ind;
            counter++;

            if(j==jd1) jd = jb;
            else jd = j;

            rho = 0.25 * (dv[0][i][jd] + dv[0][i-1][jd] + dv[0][i - 1][jd - 1] + dv[0][i][jd - 1]);
            u   = 0.25 * (dv[1][i][jd] + dv[1][i-1][jd] + dv[1][i - 1][jd - 1] + dv[1][i][jd - 1]);
            v   = 0.25 * (dv[2][i][jd] + dv[2][i-1][jd] + dv[2][i - 1][jd - 1] + dv[2][i][jd - 1]);
            p   = 0.25 * (dv[3][i][jd] + dv[3][i-1][jd] + dv[3][i - 1][jd - 1] + dv[3][i][jd - 1]);
            T   = 0.25 * (dv[4][i][jd] + dv[4][i-1][jd] + dv[4][i - 1][jd - 1] + dv[4][i][jd - 1]);
            a   = 0.25 * (dv[5][i][jd] + dv[5][i-1][jd] + dv[5][i - 1][jd - 1] + dv[5][i][jd - 1]);
            mav = 0.25 * (dv[6][i][jd] + dv[6][i-1][jd] + dv[6][i - 1][jd - 1] + dv[6][i][jd - 1]);
            kav = 0.25 * (dv[7][i][jd] + dv[7][i-1][jd] + dv[7][i - 1][jd - 1] + dv[7][i][jd - 1]);

            gradavg[0] = 0.5 * (gradfj[0][i][jd] + gradfj[0][i - 1][jd]);
            gradavg[1] = 0.5 * (gradfj[1][i][jd] + gradfj[1][i - 1][jd]);
            gradavg[2] = 0.5 * (gradfj[2][i][jd] + gradfj[2][i - 1][jd]);
            gradavg[3] = 0.5 * (gradfj[3][i][jd] + gradfj[3][i - 1][jd]);

            sx = -0.5 * (si[0][i][jd] + si[0][i - 1][jd]);
            sy = -0.5 * (si[1][i][jd] + si[1][i - 1][jd]);
            ds = sqrt(sx * sx + sy * sy);
            nx = sx / ds;
            ny = sy / ds;
            gradnx = gradavg[0] * nx + gradavg[1] * ny;
            gradny = gradavg[2] * nx + gradavg[3] * ny;
            gradnn = gradnx * nx + gradny * ny;
            dvdnx = gradnx - gradnn * nx;
            dvdny = gradny - gradnn * ny;
            sgn = Sign(gradnx);
            dvdna = sqrt(dvdnx * dvdnx + dvdny * dvdny);
            Cf = 2.0 * sgn * mav * dvdna / (rhoinf * qinf * qinf);

            sur << counter << "\t" << x[i][j] << "\t" << y[i][j] << "\t" << rho << "\t" << u << "\t" << v << "\t"
                << p << "\t" << T << "\t" << a << "\t" << Mach << "\t"
                << mav << "\t" << kav << "\t" << Cf << "\t" << std::endl;
        }
    }

    else if(bind==3){
        counter=0;
        sur << "Zone T = \"3\"," <<"I = " << ib <<","<<"\t"<< " F = POINT" << std::endl;

        for (int i = beg_seg; i <= end_seg; i++) {
            int j = sur_ind;

            //to compute Cf
            if(i==id1) id = ib;
            else id = i;

            counter++;
            rho = 0.25 * (dv[0][id][j] + dv[0][id - 1][j] + dv[0][id - 1][j - 1] + dv[0][id][j - 1]);
            u   = 0.25 * (dv[1][id][j] + dv[1][id - 1][j] + dv[1][id - 1][j - 1] + dv[1][id][j - 1]);
            v   = 0.25 * (dv[2][id][j] + dv[2][id - 1][j] + dv[2][id - 1][j - 1] + dv[2][id][j - 1]);
            p   = 0.25 * (dv[3][id][j] + dv[3][id - 1][j] + dv[3][id - 1][j - 1] + dv[3][id][j - 1]);
            T   = 0.25 * (dv[4][id][j] + dv[4][id - 1][j] + dv[4][id - 1][j - 1] + dv[4][id][j - 1]);
            a   = 0.25 * (dv[5][id][j] + dv[5][id - 1][j] + dv[5][id - 1][j - 1] + dv[5][id][j - 1]);
            mav = 0.25 * (dv[6][id][j] + dv[6][id - 1][j] + dv[6][id - 1][j - 1] + dv[6][id][j - 1]);
            kav = 0.25 * (dv[7][id][j] + dv[7][id - 1][j] + dv[7][id - 1][j - 1] + dv[7][id][j - 1]);

            gradavg[0] = 0.5 * (gradfj[0][id][j] + gradfj[0][id][j - 1]);
            gradavg[1] = 0.5 * (gradfj[1][id][j] + gradfj[1][id][j - 1]);
            gradavg[2] = 0.5 * (gradfj[2][id][j] + gradfj[2][id][j - 1]);
            gradavg[3] = 0.5 * (gradfj[3][id][j] + gradfj[3][id][j - 1]);

            sx = -0.5 * (sj[0][id][j] + sj[0][id][j - 1]);
            sy = -0.5 * (sj[1][id][j] + sj[1][id][j - 1]);
            ds = sqrt(sx * sx + sy * sy);
            nx = sx / ds;
            ny = sy / ds;
            gradnx = gradavg[0] * nx + gradavg[1] * ny;
            gradny = gradavg[2] * nx + gradavg[3] * ny;
            gradnn = gradnx * nx + gradny * ny;
            dvdnx = gradnx - gradnn * nx;
            dvdny = gradny - gradnn * ny;
            sgn = Sign(gradnx);
            dvdna = sqrt(dvdnx * dvdnx + dvdny * dvdny);
            Cf = 2.0 * sgn * mav * dvdna / (rhoinf * qinf * qinf);

            sur << counter << "\t" << x[i][j] << "\t" << y[i][j] << "\t" << rho << "\t" << u << "\t" << v << "\t"
                << p << "\t" << T << "\t" << a << "\t" << Mach << "\t"
                << mav << "\t" << kav << "\t" << Cf << "\t" << std::endl;
        }
    }

    else if(bind==4){
        counter=0;
        sur << "Zone T = \"4\"," <<"I = " << jb <<","<<"\t"<< " F = POINT" << std::endl;

        for (int j = beg_seg; j <= end_seg; j++) {
            int i = sur_ind;
            counter++;
            if(j==jd1) jd = jb;
            else jd = j;

            rho = 0.25 * (dv[0][i][jd] + dv[0][i ][jd - 1] + dv[0][i + 1][jd - 1] + dv[0][i + 1][jd]);
            u   = 0.25 * (dv[1][i][jd] + dv[1][i ][jd - 1] + dv[1][i + 1][jd - 1] + dv[1][i + 1][jd]);
            v   = 0.25 * (dv[2][i][jd] + dv[2][i ][jd - 1] + dv[2][i + 1][jd - 1] + dv[2][i + 1][jd]);
            p   = 0.25 * (dv[3][i][jd] + dv[3][i ][jd - 1] + dv[3][i + 1][jd - 1] + dv[3][i + 1][jd]);
            T   = 0.25 * (dv[4][i][jd] + dv[4][i ][jd - 1] + dv[4][i + 1][jd - 1] + dv[4][i + 1][jd]);
            a   = 0.25 * (dv[5][i][jd] + dv[5][i ][jd - 1] + dv[5][i + 1][jd - 1] + dv[5][i + 1][jd]);
            mav = 0.25 * (dv[6][i][jd] + dv[6][i ][jd - 1] + dv[6][i + 1][jd - 1] + dv[6][i + 1][jd]);
            kav = 0.25 * (dv[7][i][jd] + dv[7][i ][jd - 1] + dv[7][i + 1][jd - 1] + dv[7][i + 1][jd]);


            gradavg[0] = 0.5 * (gradfj[0][i][jd] + gradfj[0][i + 1][jd]);
            gradavg[1] = 0.5 * (gradfj[1][i][jd] + gradfj[1][i + 1][jd]);
            gradavg[2] = 0.5 * (gradfj[2][i][jd] + gradfj[2][i + 1][jd]);
            gradavg[3] = 0.5 * (gradfj[3][i][jd] + gradfj[3][i + 1][jd]);

            sx = -0.5 * (si[0][i][jd] + si[0][i + 1][jd]);
            sy = -0.5 * (si[1][i][jd] + si[1][i + 1][jd]);
            ds = sqrt(sx * sx + sy * sy);
            nx = sx / ds;
            ny = sy / ds;
            gradnx = gradavg[0] * nx + gradavg[1] * ny;
            gradny = gradavg[2] * nx + gradavg[3] * ny;
            gradnn = gradnx * nx + gradny * ny;
            dvdnx = gradnx - gradnn * nx;
            dvdny = gradny - gradnn * ny;
            sgn = Sign(gradnx);
            dvdna = sqrt(dvdnx * dvdnx + dvdny * dvdny);
            Cf = 2.0 * sgn * mav * dvdna / (rhoinf * qinf * qinf);

            sur << counter << "\t" << x[i][j] << "\t" << y[i][j] << "\t" << rho << "\t" << u << "\t" << v << "\t"
                << p << "\t" << T << "\t" << a << "\t" << Mach << "\t" << mav << "\t" << kav << "\t" << Cf << "\t"
                << std::endl;
        }
    }
    else{
        std::cout<<"wrong input fot the bind: Can't write surface values"<<std::endl;
        exit(0);
    }
    sur.close();
    delete gradavg;
}

void Interior_Solution(int bind, int indx, int beg_seg, int end_seg, int id1, int jd1, int iter, double t, double **&x,
                       double **&y, double ***&cv, double ***&dv, double ***&gradfi, double ***&gradfj){

    double rho, u, v, p, a, Mach, T, mav, kav, *gradavg;
    int counter;
    gradavg = new double[4];

    std::string oup;
    oup = "../output_kfds/";
    std::ofstream sur(oup + test_case+"_Inter_"+std::to_string(Nx)+"_"+std::to_string(Ny)+"Iter"+std::to_string(iter)
                      +".dat");
    sur.flags(std::ios::dec | std::ios::scientific);
    sur.precision(16);

    if (!sur) {
        std::cerr << "File couldn't be opened to write the interior solution" << std::endl;
        exit(1);
    }

    sur << "TITLE = "<<"\"Iter="<<iter<<", "<<"time="<<ToStringWithPrecision(t, 16)<<"\""<< std::endl
        << "VARIABLES = counter, x, y, rho, u, v, p, T, a, Mach, visc, k" << std::endl;

    if(bind==1){
        sur << "Zone"<<" "<<"T="<<"\""<<ToStringWithPrecision(t, 16)<<"\""<<","<<" "<<"I=" << ib <<","<<"\t"
            << " F = POINT" << std::endl;
        counter=0;
        for (int j = beg_seg; j <= end_seg; j++) {
            counter++;
            rho = 0.25 * (dv[0][indx][j] + dv[0][indx][j+1] + dv[0][indx-1][j + 1] + dv[0][indx][j - 1]);
            u   = 0.25 * (dv[1][indx][j] + dv[1][indx][j+1] + dv[1][indx-1][j + 1] + dv[1][indx][j - 1]);
            v   = 0.25 * (dv[2][indx][j] + dv[2][indx][j+1] + dv[2][indx-1][j + 1] + dv[2][indx][j - 1]);
            p   = 0.25 * (dv[3][indx][j] + dv[3][indx][j+1] + dv[3][indx-1][j + 1] + dv[3][indx][j - 1]);
            T   = 0.25 * (dv[4][indx][j] + dv[4][indx][j+1] + dv[4][indx-1][j + 1] + dv[4][indx][j - 1]);
            a   = 0.25 * (dv[5][indx][j] + dv[5][indx][j+1] + dv[5][indx-1][j + 1] + dv[5][indx][j - 1]);
            mav = 0.25 * (dv[6][indx][j] + dv[6][indx][j+1] + dv[6][indx-1][j + 1] + dv[6][indx][j - 1]);
            kav = 0.25 * (dv[7][indx][j] + dv[7][indx][j+1] + dv[7][indx-1][j + 1] + dv[7][indx][j - 1]);

            Mach = sqrt(u*u+v*v)/a;

             sur << counter << "\t" << x[indx][j] << "\t" << y[indx][j] << "\t" << rho << "\t" << u << "\t" << v << "\t"
                << p << "\t" << T << "\t" << a << "\t" << Mach << "\t"<< mav << "\t" << kav << std::endl;
        }
    }

    else if(bind==2){
        counter=0;
        sur << "Zone T = \"2\"," <<"I = " << jb <<","<<"\t"<< " F = POINT" << std::endl;

        for (int i = beg_seg; i <= beg_seg; i--) {
            counter++;
            rho = 0.25 * (dv[0][i][indx] + dv[0][i-1][indx] + dv[0][i - 1][indx-1] + dv[0][i][indx-1]);
            u   = 0.25 * (dv[1][i][indx] + dv[1][i-1][indx] + dv[1][i - 1][indx-1] + dv[1][i][indx-1]);
            v   = 0.25 * (dv[2][i][indx] + dv[2][i-1][indx] + dv[2][i - 1][indx-1] + dv[2][i][indx-1]);
            p   = 0.25 * (dv[3][i][indx] + dv[3][i-1][indx] + dv[3][i - 1][indx-1] + dv[3][i][indx-1]);
            T   = 0.25 * (dv[4][i][indx] + dv[4][i-1][indx] + dv[4][i - 1][indx-1] + dv[4][i][indx-1]);
            a   = 0.25 * (dv[5][i][indx] + dv[5][i-1][indx] + dv[5][i - 1][indx-1] + dv[5][i][indx-1]);
            mav = 0.25 * (dv[6][i][indx] + dv[6][i-1][indx] + dv[6][i - 1][indx-1] + dv[6][i][indx-1]);
            kav = 0.25 * (dv[7][i][indx] + dv[7][i-1][indx] + dv[7][i - 1][indx-1] + dv[7][i][indx-1]);

            Mach = sqrt(u*u+v*v)/a;

            sur << counter << "\t" << x[i][indx] << "\t" << y[i][indx] << "\t" << rho << "\t" << u << "\t" << v << "\t"
                << p << "\t" << T << "\t" << a << "\t" << Mach << "\t"<< mav << "\t" << kav << std::endl;
        }
    }

    else if(bind==3){
        counter=0;
        sur << "Zone T = \"3\"," <<"I = " << ib <<","<<"\t"<< " F = POINT" << std::endl;

        for (int j = beg_seg; j <= beg_seg; j--) {
            counter++;
            rho = 0.25 * (dv[0][indx][j] + dv[0][indx-1][j] + dv[0][indx-1][j - 1] + dv[0][indx][j - 1]);
            u   = 0.25 * (dv[1][indx][j] + dv[1][indx-1][j] + dv[1][indx-1][j - 1] + dv[1][indx][j - 1]);
            v   = 0.25 * (dv[2][indx][j] + dv[2][indx-1][j] + dv[2][indx-1][j - 1] + dv[2][indx][j - 1]);
            p   = 0.25 * (dv[3][indx][j] + dv[3][indx-1][j] + dv[3][indx-1][j - 1] + dv[3][indx][j - 1]);
            T   = 0.25 * (dv[4][indx][j] + dv[4][indx-1][j] + dv[4][indx-1][j - 1] + dv[4][indx][j - 1]);
            a   = 0.25 * (dv[5][indx][j] + dv[5][indx-1][j] + dv[5][indx-1][j - 1] + dv[5][indx][j - 1]);
            mav = 0.25 * (dv[6][indx][j] + dv[6][indx-1][j] + dv[6][indx-1][j - 1] + dv[6][indx][j - 1]);
            kav = 0.25 * (dv[7][indx][j] + dv[7][indx-1][j] + dv[7][indx-1][j - 1] + dv[7][indx][j - 1]);

            Mach = sqrt(u*u+v*v)/a;

            sur << counter << "\t" << x[indx][j] << "\t" << y[indx][j] << "\t" << rho << "\t" << u << "\t" << v << "\t"
                << p << "\t" << T << "\t" << a << "\t" << Mach << "\t"<< mav << "\t" << kav << std::endl;
        }
    }

    else if(bind==4){
        counter=0;
        sur << "Zone T = \"4\"," <<"I = " << jb <<","<<"\t"<< " F = POINT" << std::endl;

        for (int i = beg_seg; i <= end_seg; i++) {
            counter++;
            rho = 0.25 * (dv[0][i][indx] + dv[0][i][indx - 1] + dv[0][i + 1][indx - 1] + dv[0][i + 1][indx]);
            u   = 0.25 * (dv[1][i][indx] + dv[1][i][indx - 1] + dv[1][i + 1][indx - 1] + dv[1][i + 1][indx]);
            v   = 0.25 * (dv[2][i][indx] + dv[2][i][indx - 1] + dv[2][i + 1][indx - 1] + dv[2][i + 1][indx]);
            p   = 0.25 * (dv[3][i][indx] + dv[3][i][indx - 1] + dv[3][i + 1][indx - 1] + dv[3][i + 1][indx]);
            T   = 0.25 * (dv[4][i][indx] + dv[4][i][indx - 1] + dv[4][i + 1][indx - 1] + dv[4][i + 1][indx]);
            a   = 0.25 * (dv[5][i][indx] + dv[5][i][indx - 1] + dv[5][i + 1][indx - 1] + dv[5][i + 1][indx]);
            mav = 0.25 * (dv[6][i][indx] + dv[6][i][indx - 1] + dv[6][i + 1][indx - 1] + dv[6][i + 1][indx]);
            kav = 0.25 * (dv[7][i][indx] + dv[7][i][indx - 1] + dv[7][i + 1][indx - 1] + dv[7][i + 1][indx]);

            Mach = sqrt(u*u+v*v)/a;

            sur << counter << "\t" << x[i][indx] << "\t" << y[i][indx] << "\t" << rho << "\t" << u << "\t" << v << "\t"
                << p << "\t" << T << "\t" << a << "\t" << Mach << "\t"<< mav << "\t" << kav << std::endl;
        }
    }

    else{
        std::cout<<"wrong input fot the bind: Can't write surface values"<<std::endl;
        exit(0);
    }

    sur.close();
    delete gradavg;
}

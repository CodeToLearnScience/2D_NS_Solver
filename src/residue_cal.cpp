#include "../inc/global_declarations.h"

void Residue_Cal(int iter, double dt, double t, double ***&cv, double ***&cvold, double ***&dv, double ***&si, double ***&sj) {

    Resn1 = 0.0, Resn2 = 0.0, Resn3 = 0.0, Resn4 = 0.0;
    double dcv1=0.0, dcv1max=0.0, dcv11=0.0;


    for (int j = 2; j <= jb; j++) {
        for (int i = 2; i <= ib; i++) {

            dcv1 = cv[0][i][j] - cvold[0][i][j];
            Resn1 = Resn1 + dcv1*dcv1;
            if(fabs(dcv1)>=dcv1max){
                dcv1max = fabs(dcv1);
                iresmax = i;
                jresmax = j;
            }
            Resn2 = Resn2 + pow((cv[1][i][j] - cvold[1][i][j]), 2);
            Resn3 = Resn3 + pow((cv[2][i][j] - cvold[2][i][j]), 2);
            Resn4 = Resn4 + pow((cv[3][i][j] - cvold[3][i][j]), 2);
            //std::cout<<iter<<"\t"<<cv[0][i][j]<<"\t"<<cvold[0][i][j]<<std::endl;
            //std::cout<<iter<<"\t"<<dcv1<<"\t"<<dcv1max<<"\t"<<iresmax<<"\t"<<jresmax<<"\t"<<cl<<"\t"<<cd<<"\t"<<cm<<std::endl;
        }
    }

    if(iter==1){
        dcv11 = sqrt(Resn1)+1e-32;
        Resn1 = 1.0;
    }
    else{
        Resn1 = sqrt(Resn1)/dcv11;
    }

    Forces(4, 1, 2, jb, ndvar, x, y, dv, si, sj);

    //std::cout<<iter<<"\t"<<Resn1<<"\t"<<dcv1max<<"\t"<<iresmax<<"\t"<<jresmax<<"\t"<<cl<<"\t"<<cd<<"\t"<<cm<<std::endl;
    //Res<<iter<<"\t"<<Resn1<<"\t"<<dcv1max<<"\t"<<iresmax<<"\t"<<jresmax<<"\t"<<cl<<"\t"<<cd<<"\t"<<cm<<std::endl;
    //Res.close();
    if(iter%500==0) {
        std::cout << iter << "\t" << dt << "\t" << t << "\t" << dcv1max << "\t" << iresmax << "\t" << jresmax << "\t"
                  << cl << "\t" << cd << "\t" << cm << std::endl;
    }
    Res<<iter<<"\t"<<dt<<"\t"<<t<<"\t"<<dcv1max<<"\t"<<iresmax<<"\t"<<jresmax<<"\t"<<cl<<"\t"<<cd<<"\t"<<cm<<std::endl;
}

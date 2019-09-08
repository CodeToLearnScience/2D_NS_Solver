//
// Created by Maruthi on 14-07-2018.
//
#include "../inc/global_declarations.h"
#include "../inc/basic_functions.h"
void ReadGridTopology(std::string &grid_top_file, int *&pres_input_flag, int *&bc_flag, int *&bound_ind,
                        int *&bound_cell, int *&strt_bound_seg, int *&end_bound_seg, double **&pres_rhouvp) {
    int k, tot_seg_loop;
    std::string TitleLine;
    std::ifstream infile;

    infile.open("../input/" + MeshTopFile);
    std::cout<<"read top file="<<"../input/" + MeshTopFile<<std::endl;
    if(!infile){
        std::cerr << "File couldn't be opened to read the grid topology" << std::endl;
        std::exit(1);
    }

    std::cout << "---------------------------------------------------------" << std::endl
              << "Reading grid topology file" << std::endl
              << "---------------------------------------------------------" << std::endl;


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
    
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>no_boun_seg;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;
    infile>>no_bounseg_prescinp;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout<<"no of boun seg are="<<no_boun_seg<<"\t"<<"no of bound seg with prescr input="<<no_bounseg_prescinp
             <<std::endl;
    
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;

    tot_seg_loop=no_boun_seg+no_bounseg_prescinp;
    pres_input_flag=new int[no_boun_seg], bc_flag=new int[no_boun_seg], bound_ind=new int[no_boun_seg],
    bound_cell=new int[no_boun_seg], strt_bound_seg=new int[no_boun_seg], end_bound_seg=new int[no_boun_seg];

    Allocate_2D(pres_rhouvp, 4, no_bounseg_prescinp);
    int pres_val_coun=-1;
    for(int j=0; j<no_boun_seg; j++){
        infile>>pres_input_flag[j]>>bc_flag[j]>>bound_ind[j]>>bound_cell[j]>>strt_bound_seg[j]>>end_bound_seg[j];
        infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        if(pres_input_flag[j]==1){
            pres_val_coun++;
            std::cout<<"I am here and j="<<pres_val_coun<<std::endl;
            //infile>>pres_rho>>pres_u>>pres_v>>pres_p;
            infile>>pres_rhouvp[0][pres_val_coun]>>pres_rhouvp[1][pres_val_coun]>>pres_rhouvp[2][pres_val_coun]>>
                  pres_rhouvp[3][pres_val_coun];
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        else if(bc_flag[j]==50){
            infile>>Twall_iso;
            std::cout<<"Twall_iso="<<Twall_iso<<std::endl;
        }
    }
    for (int j = 0; j<no_boun_seg; j++) {
        std::cout<<pres_input_flag[j]<<" "<<bc_flag[j]<<" "<<bound_ind[j]<<" "<<bound_cell[j]<<" "
                 <<strt_bound_seg[j]<<" "<<end_bound_seg[j]<<std::endl;
    }

    for (int k=0; k<4; k++){
        for (int j=0; j<no_bounseg_prescinp; j++){
            std::cout<<"rho, u, v, p="<<pres_rhouvp[k][0]<<"\t";
        }
        std::cout<<std::endl;
    }

    
    getline(infile, TitleLine);
    std::cout<<TitleLine<<std::endl;

    infile.close();
}


void extractIntegerWords(std::string str){
    std::stringstream ss;

    /* Storing the whole string into string stream */
    ss << str;

    /* Running loop till the end of the stream */
    std::string temp;
    int found;
    while (!ss.eof()) {

        /* extracting word by word from stream */
        ss >> temp;

        /* Checking the given word is integer or not */
        if (std::stringstream(temp) >> found)
            std::cout << found << " ";

        /* To save from space at the end of string */
        temp = "";
    }
}

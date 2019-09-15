#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>


using namespace std;

enum stretch_flag{str_of_line, end_of_line};

int Nx, Ny;
double **x, **y;
double dx, dy, xll, xul, yll, yul, len, D, A, num, den;

template <typename T>
T **Allocate_2D(T ** &m, int t1, int t2) {
	m=new T* [t1];
        for (int i=0; i<t1; ++i) {
                m[i]=new T [t2];
                for (int j=0; j<t2; ++j)
                        m[i][j]=0.0;
        }
        return m;
}

template <typename T>
T *linePartition(T *&t_line_to_part, const T t_start_of_line, const T t_end_of_line, const int t_part_size){

    T part_len = (t_end_of_line - t_start_of_line)/t_part_size;

    for(int i=0;i<=t_part_size; i++){
        t_line_to_part[i] = t_start_of_line + i*part_len;
    }

    return t_line_to_part;
}

void WriteToFile(const int Nx, const int Ny, double **&x, double **&y){

	//ofstream GridFFS("SWBLIGrid"+std::to_string(Nx)+std::to_string(Ny)+".dat");
  	ofstream GridFFS("Blasius"+std::to_string(Nx)+std::to_string(Ny)+".dat");

	ofstream GridFFSTecPlot("SWBLIGridTecPlot"+std::to_string(Nx)+std::to_string(Ny)+".dat");
	//GridFFS12040.flags(ios::dec | ios::scientific);
	//GridFFS12040.precision(5);

	if(!GridFFS){
		cerr<<"File couldn't be opened to write the solution"<<endl;
		exit(1);
	}
	
	GridFFSTecPlot << "TITLE = Flow" << endl << "VARIABLES = xc, yc, rho" << endl;
	GridFFSTecPlot << "Zone T = Omega I = " << Nx+1 << " J = " << Ny+1 << endl ;
	int count=0;

	for(int j=0;j<=Ny;j++){
		for(int i=0;i<=Nx;i++){
			count=count+1;
			GridFFS<<count<<"\t"<<x[i][j]<<"\t"<<y[i][j]<<"\t"<<endl;
			GridFFSTecPlot<<x[i][j]<<"\t"<<y[i][j]<<"\t"<<0.0<<endl;

		}
	}

	GridFFS.close();
	GridFFSTecPlot.close();

}

void Write_SolutionVtk(const int Nx, const int Ny, double **&x, double **&y)
{

    std::string dataType = "double";

     std::ofstream other("Blasius" + std::to_string(Nx) + "_" + std::to_string(Ny)+".vtk");
    other.flags(std::ios::dec | std::ios::scientific);
    other.precision(16);

    if (!other)
    {
        std::cerr << "File couldn't be opened to write the solution" << std::endl;
        exit(1);
    }

    int Nz = 1;
    other << "# vtk DataFile Version 3.0" << std::endl;
    other << "Blasius grid to visualize"<< std::endl;
    other << "ASCII" << std::endl;
    other << "DATASET STRUCTURED_GRID" << std::endl;
    other << "DIMENSIONS"
          << " " << Nx + 1<< " " << Ny +1 << " " << 1 << std::endl;
    other << "POINTS"
          << " " << (Nx+1) * (Ny+1) * Nz << " "
          << "double" << std::endl;

    for (int j = 0; j <= Ny; j++)
    {
        for (int i = 0; i <= Nx; i++)
        {
            other << x[i][j] << "\t" << y[i][j] << "\t" << 1 << std::endl;
        }
    }
}

template <typename T>
T *strechedGrid(T *&vec, T vec_max, T vec_min, T vec_size){

    T *eta, *A, B;
    const float beta = 1.01;

    B = (beta + 1.0)/(beta - 1.0);
    eta = new T [vec_size];
    A = new T [vec_size];
    
    for(int j=0;j<=vec_size;j++){

             A[j] = (beta + (1 - vec[j]/vec_max))/(beta - (1 - vec[j]/vec_max));
             std::cout<<"val of A="<<A[j]<<std::endl;
    
    }

	    for(int j=0;j<=vec_size;j++){
            
             eta[j] = 1.0 - (log(A[j])/log(B));

    }
    
	    for(int j=0;j<=vec_size;j++){
                //  vec[j] = yMax*(((beta+1.0)- (beta-1.0)*pow(B, 1-eta[j]))/(1.0+pow(B, 1-eta[j])));
               vec[j] = (((beta+1.0)- (beta-1.0)*pow(B, 1-vec[j]))/(1.0+pow(B, 1-vec[j])));

                std::cout<<"val of vec="<<vec[j]<<std::endl;

        }
   
    delete eta, A;
    return vec;
}

// Function to reverse arr[] from start to end*/
template <typename T>
void reverseArray(T *&arr, int start, int end) 
{ 
    if (start >= end) 
    return; 
      
    double temp = arr[start];  
    arr[start] = arr[end]; 
    arr[end] = temp; 
      
    // Recursive Function calling 
    reverseArray(arr, start + 1, end - 1);  
}  

// Function to reverse arr[] from start to end*/
template <typename T>
void rvereseArray(T *&arr, int start, int end) 
{ 
    while (start < end) 
    { 
        double temp = arr[start];  
        arr[start] = arr[end]; 
        arr[end] = temp; 
        start++; 
        end--; 
    }  
}      

template <typename T>
T *strechedGridUsingTanh(T *&vec, T vec_max, T vec_min, const int vec_size, const double delta, stretch_flag flag){

    T *eta;

    std::cout<<"vec_size="<<vec_size<<std::endl;
    
    eta = new T [vec_size+1];
    
    
    for(int j=0;j<=vec_size;j++){

            if(flag == str_of_line){
                eta[j] = 1.0 + (tanh(delta*(vec[j]/vec_max - 1.0))/ tanh(delta));//(((beta+1.0)- (beta-1.0)*pow(B, 1-vec[j]))/(1.0+pow(B, 1-vec[j])));
            }
            else if(flag == end_of_line) 
            {
                eta[j] =  (tanh(delta*(vec[j]/vec_max - 1.0))/ tanh(delta));//(((beta+1.0)- (beta-1.0)*pow(B, 1-vec[j]))/(1.0+pow(B, 1-vec[j])));

            }
            
            // std::cout<<"val of eta="<<eta[j]<<std::endl;
    
    }

                // reverseArray(eta, 0, vec_size);
            if(flag == end_of_line) 
            {
                        rvereseArray(eta, 0, vec_size);

            }

	    for(int j=0;j<=vec_size;j++){
            if(flag == str_of_line){
                vec[j] = vec_max*eta[j];
            }
            else if(flag==end_of_line){
              vec[j] = -vec_max*eta[j];



              std::cout<<"val of vec="<<vec[j]<<std::endl;

            }

            // std::cout<<"val of vec="<<vec[j]<<std::endl;
    }
    delete[] eta;
    return vec;
    
}

    

int main(){
	
	
	// dx=(xul-xll)/(Nx); dy=(yul-yll)/(Ny);
	// len = (xul-xll);	
	double *yStretch, *xStretch1, *xStretch2;
    // yStretch = new double [Ny+1];

    // double xll_1 = 0.0;
    // double xul_1 = 1.0;

    // yll = 0.0;
    // yul = 2.0;

    int Nx1 = 40; Ny = 80;
    int Nx2 = 20;

    xStretch1 = new double [Nx1+1];
    xStretch2 = new double [Nx2+1];
    yStretch = new double [Ny+1];

    Allocate_2D(x, Nx1+Nx2+1, Ny+1); Allocate_2D(y, Nx1+Nx2+1, Ny+1);
    std::cout<<"---------------first stretch in x---------------"<<std::endl;
    linePartition(xStretch1, double(0.0), double(1.0), Nx1);
    strechedGridUsingTanh(xStretch1, 1.0, 0.0, Nx1, 1.1, end_of_line);

    linePartition(yStretch, double (0.0), double(2.0), Ny);
    strechedGridUsingTanh(yStretch, 2.0, 1.0, Ny, 1.5, str_of_line);

    std::cout<<"---------------second stretch in x---------------"<<std::endl;
        linePartition(xStretch2, double(0.0), double(0.5), Nx2);
        strechedGridUsingTanh(xStretch2, 0.5, 0.0, Nx2, 1.1, str_of_line);

    for (int i=0;i<=Nx1;i++){
	    for(int j=0;j<=Ny;j++){
		        x[i][j]=xStretch1[i];
        		y[i][j]=yStretch[j];
        
        }
    }


    for (int i=0;i<=Nx2;i++){
        xStretch2[i] = xStretch2[i] + 1;
        // std::cout<<"xStretch2 val="<<xStretch2[i]<<std::endl;
    }
    for (int i=Nx1+1;i<=(Nx1+Nx2);i++){
	    for(int j=0;j<=Ny;j++){
		        x[i][j]=xStretch2[i-(Nx1)];
        		y[i][j]=yStretch[j];
        
        }
    }


	WriteToFile(Nx1+Nx2, Ny, x, y);
    Write_SolutionVtk(Nx1+Nx2, Ny, x, y);

    delete [] xStretch1;
    delete [] xStretch2;
    delete [] yStretch;
	return 0;
}

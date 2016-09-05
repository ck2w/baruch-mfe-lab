#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <cholesky.h>
#include <load.h>
#include <Eigen/Dense>

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " 
                  << argv[0] 
                  << " <matrix.csv> -n " 
                  << std::endl;
        return 1;
    }

    int p=6;
    if (argc == 3) {
        p = atoi(&argv[2][1]);
    }
                
    int N = countRows(argv[1]);
    assert(N>0);
    
    Eigen::MatrixXd mat = readCSV(argv[1], N, N);
    Eigen::MatrixXd U = cholesky(mat);

    for (int i=0; i<N; i++) {
        std::cout << " | ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << mat(i,j) 
                      << " ";
        }
        std::cout << " | = | ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << U(j,i) 
                      << " ";
        }
        std::cout << " | ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << U(i,j) 
                      << " ";
        }
        std::cout << " | " << std::endl;
    }

    return 0;
}


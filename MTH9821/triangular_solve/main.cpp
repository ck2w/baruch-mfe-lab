#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <triangular_solve.h>
#include <load.h>
#include <Eigen/Dense>

int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: " 
                  << argv[0] 
                  << " -l/u <matrix.csv> <vector.csv> -n " 
                  << std::endl;
        return 1;
    }

    int p=6;
    if (argc == 5) {
        p = atoi(&argv[4][1]);
    }

    int N = countRows(argv[2]);
    assert(N>0);
    
    Eigen::MatrixXd mat = readCSV(argv[2], N, N);
    Eigen::VectorXd vec = readCSV(argv[3], N);

    Eigen::VectorXd x;
    std::string mode = argv[1];
    if (mode == "-l") {
        mat.triangularView<Eigen::StrictlyUpper>().setZero();
        x = forward_subst(mat, vec);
    }
    else if (mode == "-u") {
        mat.triangularView<Eigen::StrictlyLower>().setZero();
        x = backward_subst(mat, vec);
    }

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
        std::cout  << "| " 
                   << std::right 
                   << std::setw(p+3) 
                   << x(i) 
                   << " = " 
                   << vec(i) 
                   << std::endl;
    }

    return 0;
}


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <load.h>
#include <linear_solve.h>
#include <Eigen/Dense>

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: " 
                  << argv[0] 
                  << " <matrix.csv> <vector.csv> -n " 
                  << std::endl;
        return 1;
    }

    int p=6;
    if (argc == 4) {
        p = atoi(&argv[3][1]);
    }

    int N = countRows(argv[1]);
    assert(N>0);
    
    Eigen::MatrixXd mat = readCSV(argv[1], N, N);
    Eigen::VectorXd vec = readCSV(argv[2], N);

    Eigen::VectorXd x;
    std::string mode = argv[1];
    x = spd_solve(mat, vec);

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


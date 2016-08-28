#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <lu.h>
#include <load.h>
#include <Eigen/Dense>

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " 
                  << argv[0] 
                  << " <matrix.csv> -p -n " 
                  << std::endl;
        return 1;
    }

    bool pivoting = false;
    int p=6;
    if (argc > 2) {
        if (std::string(argv[2]) == "-p") {
            pivoting = true;
        }
        else {
            p = atoi(&argv[3][1]);
        }

        if (argc > 3) {
            if ( std::string(argv[3]) == "-p" ) {
                pivoting = true;
            }
            else {
                p = atoi(&argv[3][1]);
            }
        }
    }
                
    int N = countRows(argv[1]);
    assert(N>0);
    
    Eigen::MatrixXd mat = readCSV(argv[1], N, N);

    Eigen::MatrixXd L;
    Eigen::MatrixXd U;
    Eigen::MatrixXd P;

    if (pivoting) {
        std::tuple<Eigen::VectorXi, Eigen::MatrixXd, Eigen::MatrixXd> res 
            = lu_row_pivoting(mat);
        Eigen::VectorXi p = std::get<0>(res);
        L = std::get<1>(res);
        U = std::get<2>(res);
        P = Eigen::MatrixXd::Zero(N,N);
        for (int i=0; i<N; i++) {
            P(i,p(i)-1) = 1;
        }
    }
    else {
        std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> res = lu_no_pivoting(mat);
        L = std::get<0>(res);
        U = std::get<1>(res);
    }

    for (int i=0; i<N; i++) {
        if (pivoting) {
            std::cout << " | ";
            for (int j=0; j<N; j++) {
                std::cout << std::right 
                          << std::setw(p+3) 
                          << std::fixed 
                          << std::setprecision(p) 
                          << P(i,j)
                          << " ";
            }
        }
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
                      << L(i,j) 
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


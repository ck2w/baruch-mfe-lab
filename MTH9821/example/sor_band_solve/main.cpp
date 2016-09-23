#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <Eigen/Dense>
#include <band_mult.h>
#include <band_conversion.h>
#include <jacobi.h>
#include <sor.h>

int main(int argc, char* argv[])
{
    // Relaxation parameters
    double omega_start = 0.95;
    double omega_delta = 0.25;
    int    omega_count = 2;
    // Tolerance factor
    double tol=1e-6;
    // Output precision
    int p=9;
    // Matrix dimension
    int N=8;
    // Band width
    int band=3;
    // Matrix specification
    Eigen::ArrayXXd A = Eigen::ArrayXXd::Zero(N,2*band+1);

    //diagonal
    for (int i=0; i<N; i++) {A(i,band) = 8.0;}
    //lower diagonal
    for (int i=1; i<N; i++) { 
        //if (i>=2) A(i,band-2) = 3;
        //if (i>=3) A(i,band-3) = 1;
        A(i,band-2) = 2.0;
        A(i,band-3) = 1.0;
    }
    //upper diagonal
    for (int i=0; i<N-1; i++) 
    {
        //if (i<=5) A(i,band+2) = -2;
        //if (i<=4) A(i,band+3) = -1;
        A(i,band+2) = -1.0;
        A(i,band+3) = -2.0;
    }

    Eigen::MatrixXd b = Eigen::VectorXd::Zero(N);
    for (int i=0; i<N; i++) {
        b(i) = (2.0*i-3.0)/(2.0*i*i+1.0);
    }
    
    // Jacobi iterative solve
    std::tuple<Eigen::VectorXd, int> res = jacobi(A, band, b, tol);
    Eigen::VectorXd x = std::get<0>(res);
    int ic = std::get<1>(res);
    Eigen::VectorXd r = b - band_mult(A,band,band,x);

    std::cout << std::fixed
              << std::setprecision(p)
              << "Jacobi: ,"
              << ", Iterations =, " << ic 
              << ", residual =, " << r.norm() 
              << std::endl;

    std::cout << "x = " << std::endl;
    for (int i=0; i<N; i++) {
        std::cout << std::fixed << std::setprecision(p) << x(i) << std::endl;
    }
    
    // Gauss-Seidel iterative solve
    res = gs(A, band, b, tol);
    x = std::get<0>(res);
    ic = std::get<1>(res);
    r = b - band_mult(A,band,band,x);

    std::cout << std::fixed
              << std::setprecision(p)
              << "Gauss-Seidel: ,"
              << ", Iterations =, " << ic 
              << ", residual =, " << r.norm() 
              << std::endl;

    std::cout << "x = " << std::endl;
    for (int i=0; i<N; i++) {
        std::cout << std::fixed << std::setprecision(p) << x(i) << std::endl;
    }

    // SOR iterative solve for all omega values
    double omega = omega_start;
    for (int k=0; k<omega_count; k++) {
        std::tuple<Eigen::VectorXd, int> res = sor(omega, A, band, b, tol);
        Eigen::VectorXd x = std::get<0>(res);
        int ic = std::get<1>(res);
        Eigen::VectorXd r = b - band_mult(A,band,band,x);

        std::cout << std::fixed
                  << std::setprecision(p)
                  << "Omega = ," << omega
                  << ", Iterations =, " << ic 
                  << ", residual =, " << r.norm() 
                  << std::endl;

        std::cout << "x = " << std::endl;
        for (int i=0; i<N; i++) {
            std::cout << std::fixed << std::setprecision(p) << x(i) << std::endl;
        }
        
        omega += omega_delta;
    }

    return 0;
}


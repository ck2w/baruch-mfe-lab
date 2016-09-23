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
    double omega_delta = 0.26;
    int    omega_count = 2;
    // Tolerance factor
    double tol=1e-6;
    // Output precision
    int p=9;
    // Matrix dimension
    int N=8;
    // Matrix specification
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N,N);

    //diagonal
    for (int i=0; i<N; i++) {A(i,i) = 9.0;}
    //lower diagonal
    for (int i=0; i<N; i++) { 
        if (i<=5) A(i,i+2) = -2.0;
        if (i>=2) A(i,i-2) = 4.0;
        if (i<=4) A(i,i+3) = -2.0;
        if (i>=3) A(i,i-3) = -1.0; 
    }

    Eigen::MatrixXd b = Eigen::VectorXd::Zero(N);
    for (int i=0; i<N; i++) {
        b(i) = (4.0*i-3.0)/(2.0*i*i+1.0);
    }
    
    // Jacobi iterative solve
    std::tuple<Eigen::VectorXd, int> res = jacobi(A, b, tol);
    Eigen::VectorXd x = std::get<0>(res);
    int ic = std::get<1>(res);
    Eigen::VectorXd r = b - A*x;
    
    std::cout << "x solution = " << std::endl;
    for (int i=0; i<N; i++) {
        std::cout << std::fixed << std::setprecision(p) << x(i) << std::endl;
    }

    std::cout << std::fixed
              << std::setprecision(p)
              << "Jacobi: ,"
              << ", Iterations =, " << ic 
              << ", residual =, " << r.norm() 
              << std::endl;
    
    // Gauss-Seidel iterative solve
    res = gs(A, b, tol);
    x = std::get<0>(res);
    ic = std::get<1>(res);
    r = b - A*x; 
    
    std::cout << "x solution = " << std::endl;
    for (int i=0; i<N; i++) {
        std::cout << std::fixed << std::setprecision(p) << x(i) << std::endl;
    }

    std::cout << std::fixed
              << std::setprecision(p)
              << "Gauss-Seidel: ,"
              << ", Iterations =, " << ic 
              << ", residual =, " << r.norm() 
              << std::endl;

    // SOR iterative solve for all omega values
    double omega = omega_start;
    for (int k=0; k<omega_count; k++) {
        std::tuple<Eigen::VectorXd, int> res = sor(omega, A, b, tol);
        Eigen::VectorXd x = std::get<0>(res);
        int ic = std::get<1>(res);
        Eigen::VectorXd r = b - A*x; 
        
        std::cout << "x solution = " << std::endl;
        for (int i=0; i<N; i++) {
            std::cout << std::fixed << std::setprecision(p) << x(i) << std::endl;
        }

        std::cout << std::fixed
                  << std::setprecision(p)
                  << "Omega = ," << omega
                  << ", Iterations =, " << ic 
                  << ", residual =, " << r.norm() 
                  << std::endl;
        
        omega += omega_delta;
    }

    return 0;
}


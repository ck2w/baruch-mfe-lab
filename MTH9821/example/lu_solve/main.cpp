#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <cmath>
#include <lu.h>
#include <triangular_solve.h>
#include <Eigen/Dense>

int main(int argc, char* argv[])
{
    // LU without pivoting or with row-pivoting
    bool pivoting = true;
    // Output precision
    int p=9;
    // Matrix dimension
    int N=9;
    // Matrix specification
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N,N); 
    for (int i=0; i<N; i++) {
        A(i,i) = 4.0;
        if (i>=1) A(i,i-1) = -1.0;
        if (i>=2) A(i,i-2) = 3.0;
        if (i<=6) A(i,i+2) = -2.0;
    }
    // Vector specification
    Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
    for (int i=0; i<N; i++) {
        b(i) = std::sqrt(0.5*i*i-i+6.0);
    }

    // Implementation details
    Eigen::VectorXd Pb = b;
    Eigen::MatrixXd L;
    Eigen::MatrixXd U;
    Eigen::MatrixXd P;

    if (pivoting) {
        std::tuple<Eigen::VectorXi, Eigen::MatrixXd, Eigen::MatrixXd> res 
            = lu_row_pivoting(A);
        Eigen::VectorXi p = std::get<0>(res);
        L = std::get<1>(res);
        U = std::get<2>(res);
        P = Eigen::MatrixXd::Zero(N,N);
        for (int i=0; i<N; i++) {
            P(i,p(i)-1) = 1;
        }

        Pb = P*b;
    }
    else {
        std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> res = lu_no_pivoting(A);
        L = std::get<0>(res);
        U = std::get<1>(res);
    }

    // Perform forward substitution and backward substitution
    Eigen::VectorXd y = forward_subst(L, Pb);
    Eigen::VectorXd x = backward_subst(U, y);

    // Compute residual error
    double R = (A*x-b).norm();
    std::cout << "Residual error =, "
              << std::setprecision(p)
              << R
              << std::endl;

    // Compute L-inverse, U-inverse, and A-inverse
    Eigen::MatrixXd L_inv = Eigen::MatrixXd::Zero(N,N);
    Eigen::MatrixXd U_inv = Eigen::MatrixXd::Zero(N,N);
    Eigen::MatrixXd A_inv = Eigen::MatrixXd::Zero(N,N);
    // L-inverse
    for (int i=0; i<N; i++) {
        Eigen::VectorXd ei = Eigen::VectorXd::Zero(N);
        ei(i) = 1;
        Eigen::VectorXd x = forward_subst(L, ei);
        L_inv.col(i) = x;
    }
    // U-inverse
    for (int i=0; i<N; i++) {
        Eigen::VectorXd ei = Eigen::VectorXd::Zero(N);
        ei(i) = 1;
        Eigen::VectorXd x = backward_subst(U, ei);
        U_inv.col(i) = x;
    }
    // A-inverse
    for (int i=0; i<N; i++) {
        Eigen::VectorXd ei = Eigen::VectorXd::Zero(N);
        ei(i) = 1;
        if (pivoting) {ei = P*ei;}
        Eigen::VectorXd y = forward_subst(L, ei);
        Eigen::VectorXd x = backward_subst(U, y);
        A_inv.col(i) = x;
    }
    
    // Compute residual error by matrix inverse
    Eigen::VectorXd v = A_inv*b;
    R = (b-A*v).norm();
    std::cout << "Residual error by inverse =, "
              << std::setprecision(p)
              << R
              << std::endl;

    std::cout << "LU - decomposition" << std::endl;
    for (int i=0; i<N; i++) {
        if (pivoting) {
            std::cout << " |, ";
            for (int j=0; j<N; j++) {
                std::cout << std::right 
                          << std::setw(p+3) 
                          << std::fixed 
                          << std::setprecision(p) 
                          << P(i,j)
                          << ", ";
            }
        }
        std::cout << " |, ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << A(i,j) 
                      << ", ";
        }
        std::cout << " |, =, |, ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << L(i,j) 
                      << ", ";
        }
        std::cout << " |, ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << U(i,j) 
                      << ", ";
        }
        std::cout << " |, " << std::endl;
    }

    std::cout << "Linear solve - x:" << std::endl;
    for (int i=0; i<N; i++) {
        std::cout << std::right 
                  << std::setw(p+3) 
                  << std::fixed 
                  << std::setprecision(p) 
                  << x(i)
                  << std::endl; 
    }
    
    std::cout << "Inverse relations:" << std::endl;
    for (int i=0; i<N; i++) {
        std::cout << " |, ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << A_inv(i,j) 
                      << ", ";
        }
        if (pivoting) {
            std::cout << " |, ";
            for (int j=0; j<N; j++) {
                std::cout << std::right 
                          << std::setw(p+3) 
                          << std::fixed 
                          << std::setprecision(p) 
                          << P.transpose()(i,j)
                          << ", ";
            }
        }
        std::cout << " |, =, |, ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << U_inv(i,j) 
                      << ", ";
        }
        std::cout << " |, ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << L_inv(i,j) 
                      << ", ";
        }
        std::cout << " |, " << std::endl;
    }
    
    std::cout << "Linear solve by inversion- x:" << std::endl;
    for (int i=0; i<N; i++) {
        std::cout << std::right 
                  << std::setw(p+3) 
                  << std::fixed 
                  << std::setprecision(p) 
                  << v(i)
                  << std::endl; 
    }
    
    // Verify accuracy of matrix inverse
    std::cout << "L * L-inverse residual error =, "
              << std::scientific
              << std::setprecision(p)
              << (L*L_inv-Eigen::MatrixXd::Identity(N,N)).norm() 
              << std::endl;
    std::cout << "U * U-inverse residual error =, "
              << std::scientific
              << std::setprecision(p)
              << (U*U_inv-Eigen::MatrixXd::Identity(N,N)).norm() 
              << std::endl;
    std::cout << "A * A-inverse residual error =, "
              << std::scientific
              << std::setprecision(p)
              << (A*A_inv-Eigen::MatrixXd::Identity(N,N)).norm() 
              << std::endl;

    return 0;
}


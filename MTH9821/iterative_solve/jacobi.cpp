#include <jacobi.h>
#include <linear_iterate.h>
#include <Eigen/Dense>
#include <cassert>
#include <iostream>

std::tuple<Eigen::VectorXd, int> jacobi(const Eigen::MatrixXd & A, 
                                        const Eigen::VectorXd & b,
                                        double tol)
{
    int n = A.rows();
    assert(n == A.cols());

    Eigen::VectorXd d = A.diagonal();
    Eigen::MatrixXd N = A;
    N.diagonal().setZero();
   
    // estimate of the solution 
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(n); 
    return linear_iterate_diagonal(x0,d,N,b,tol);
}

std::tuple<Eigen::VectorXd, int> jacobi(const Eigen::ArrayXXd & A, int m, 
                                        const Eigen::VectorXd & b, 
                                        double tol)
{
    int n = A.rows();
    assert(2*m+1 == A.cols());
    Eigen::ArrayXd d = A.col(m);
    Eigen::ArrayXXd N = A;
    N.col(m).setZero();
    
    // estimate of the solution 
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(n); 
    return linear_iterate_diagonal_banded(x0,d,N,b,m,tol);
}

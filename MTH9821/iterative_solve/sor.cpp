#include <sor.h>
#include <linear_iterate.h>
#include <Eigen/Dense>
#include <tuple>

std::tuple<Eigen::VectorXd, int> sor(double omega,
                                     const Eigen::MatrixXd & A, 
                                     const Eigen::VectorXd & b,
                                     double tol)
{
    assert(omega>0);
    assert(omega<2);
    int n = A.rows();
    assert(n == A.cols());
    double w = 1.0/omega;

    Eigen::VectorXd d = A.diagonal();

    Eigen::MatrixXd M = A; 
    M.triangularView<Eigen::StrictlyUpper>().setZero();
    M += ((w-1)*d).asDiagonal();

    Eigen::MatrixXd N = A;
    N.triangularView<Eigen::StrictlyLower>().setZero();
    N -= (w*d).asDiagonal();

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(n);
    return linear_iterate_triangular(x0,M,N,b,tol);
}

std::tuple<Eigen::VectorXd, int> sor(double omega,
                                     const Eigen::ArrayXXd & A, int m, 
                                     const Eigen::VectorXd & b, 
                                     double tol)
{
    assert(omega>0);
    assert(omega<2);
    int nrow = A.rows();
    int ncol = A.cols();
    assert(2*m+1 == ncol);
    double w = 1.0/omega;

    Eigen::ArrayXd d = A.col(m);

    Eigen::ArrayXXd M = A.block(0,0,nrow,m+1);
    M.col(m) += (w-1)*d;

    Eigen::ArrayXXd N = A.block(0,m,nrow,ncol-m);
    N.col(0) -= w*d;
    
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(nrow);
    return linear_iterate_triangular_banded(x0,M,N,b,m,tol);
}

std::tuple<Eigen::VectorXd, int> gs(const Eigen::MatrixXd & A, 
                                    const Eigen::VectorXd & b,
                                    double tol)
{
    return sor(1,A,b,tol);
}

std::tuple<Eigen::VectorXd, int> gs(const Eigen::ArrayXXd & A, int m,
                                    const Eigen::VectorXd & b,
                                    double tol)
{
    return sor(1,A,m,b,tol);
}

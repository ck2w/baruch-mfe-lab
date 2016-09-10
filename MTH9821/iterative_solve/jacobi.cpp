#include <jacobi.h>
#include <linear_iterate.h>
#include <Eigen/Dense>
#include <cassert>

std::tuple<Eigen::VectorXd, int> jacobi(const Eigen::MatrixXd & A, 
                                        const Eigen::VectorXd & b,
                                        double tol)
{
    int n = A.rows();
    assert(n == A.cols());
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(n);

    Eigen::VectorXd d = A.diagonal();
    Eigen::MatrixXd N = A;
    N.diagonal().setZero();
    
    return linear_iterate(x0,d,N,b,tol);
}

#include <sor.h>
#include <linear_iterate.h>
#include <Eigen/Dense>
#include <tuple>

std::tuple<Eigen::VectorXd, int> sor(double omega,
                                     const Eigen::MatrixXd & A, 
                                     const Eigen::VectorXd & b,
                                     double tol)
{
    int n = A.rows();
    assert(n == A.cols());
    double w = 1.0/omega;
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(n);

    Eigen::VectorXd D = A.diagonal();

    Eigen::MatrixXd M = A; 
    M.triangularView<Eigen::StrictlyUpper>().setZero();
    M += ((w-1)*D).asDiagonal();

    Eigen::MatrixXd N = A;
    N.triangularView<Eigen::StrictlyLower>().setZero();
    N -= (w*D).asDiagonal();

    return linear_iterate(x0,M,N,b,tol);
}

std::tuple<Eigen::VectorXd, int> gs(const Eigen::MatrixXd & A, 
                                    const Eigen::VectorXd & b,
                                    double tol)
{
    return sor(1,A,b,tol);
}

//std::tuple<Eigen::VectorXd, int> sor_banded(double omega,
//                                            const Eigen::MatrixXd & A, int m, 
//                                            const Eigen::VectorXd & b, 
//                                            double tol)
//{
//
//}

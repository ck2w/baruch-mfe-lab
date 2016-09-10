#include <linear_iterate.h>
#include <Eigen/Dense>
#include <tuple>
#include <cassert>
#include <triangular_solve.h>

std::tuple<Eigen::VectorXd, int> linear_iterate(const Eigen::VectorXd & x0,
                                                const Eigen::MatrixXd & M, 
                                                const Eigen::MatrixXd & N, 
                                                const Eigen::VectorXd & b,
                                                double tol)
{
    Eigen::VectorXd c = forward_subst(M,b);
    Eigen::VectorXd x = x0;
    Eigen::VectorXd x_old = x;
    Eigen::VectorXd r = b-(M+N)*x0; 
    double stop_resid = tol*r.norm();

    int n = 0;
    while (r.norm() > stop_resid && n < MAX_ITERATION) {
        x = c-forward_subst(M, N*x_old);
        r = b-(M+N)*x;
        x_old = x;
        n++;
    }

    return std::make_tuple(x,n);
}

std::tuple<Eigen::VectorXd, int> linear_iterate(const Eigen::VectorXd & x0,
                                                const Eigen::VectorXd & d,
                                                const Eigen::MatrixXd & N, 
                                                const Eigen::VectorXd & b,
                                                double tol)
{
    Eigen::MatrixXd D = d.asDiagonal();
    Eigen::VectorXd c = forward_subst_banded(D,0,b);
    Eigen::VectorXd x = x0;
    Eigen::VectorXd x_old = x;
    Eigen::VectorXd r = b-(D+N)*x0; 
    double stop_resid = tol*r.norm();

    int n = 0;
    while (r.norm() > stop_resid && n < MAX_ITERATION) {
        x = c-forward_subst_banded(D,0,N*x_old);
        r = b-(D+N)*x;
        x_old = x;
        n++;
    }

    return std::make_tuple(x,n);
}

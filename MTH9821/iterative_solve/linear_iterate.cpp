#include <linear_iterate.h>
#include <Eigen/Dense>
#include <tuple>
#include <cassert>
#include <triangular_solve.h>
#include <band_mult.h>
#include <band_add.h>

std::tuple<Eigen::VectorXd, int> linear_iterate_triangular
                                               (const Eigen::VectorXd & x0,
                                                const Eigen::MatrixXd & M, 
                                                const Eigen::MatrixXd & N, 
                                                const Eigen::VectorXd & b,
                                                double tol)
{
    Eigen::VectorXd c = forward_subst(M,b);
    Eigen::VectorXd x = x0, x_old = x0;
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

std::tuple<Eigen::VectorXd, int> linear_iterate_triangular_banded
                                               (const Eigen::VectorXd & x0,
                                                const Eigen::ArrayXXd & M, 
                                                const Eigen::ArrayXXd & N, 
                                                const Eigen::VectorXd & b,
                                                int m, double tol)
{
    Eigen::ArrayXXd MN = band_add(M,m,0,N,0,m);
    Eigen::VectorXd c = forward_subst_banded(M,b);
    Eigen::VectorXd x = x0, x_old = x0;
    Eigen::VectorXd r = b-band_mult(MN,m,m,x0);
    double stop_resid = tol*r.norm();

    int n = 0;
    while (r.norm() > stop_resid && n < MAX_ITERATION) {
        x = c-forward_subst_banded(M, band_mult(N,0,m,x_old));
        r = b-band_mult(MN,m,m,x);
        x_old = x;
        n++;
    }

    return std::make_tuple(x,n);
}


std::tuple<Eigen::VectorXd, int> linear_iterate_diagonal
                                               (const Eigen::VectorXd & x0,
                                                const Eigen::VectorXd & d,
                                                const Eigen::MatrixXd & N, 
                                                const Eigen::VectorXd & b,
                                                double tol)
{
    Eigen::MatrixXd D = d.asDiagonal();
    Eigen::VectorXd c = forward_subst_banded(D,0,b);
    Eigen::VectorXd x = x0, x_old = x0;
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

std::tuple<Eigen::VectorXd, int> linear_iterate_diagonal_banded
                                               (const Eigen::VectorXd & x0,
                                                const Eigen::ArrayXd  & d,
                                                const Eigen::ArrayXXd & N, 
                                                const Eigen::VectorXd & b,
                                                int m, double tol)
{
    Eigen::ArrayXXd MN = N;
    MN.col(m) += d;
    Eigen::VectorXd c = b.array()/d;
    Eigen::VectorXd x = x0, x_old = x0;
    Eigen::VectorXd r = b-band_mult(MN,m,m,x0);
    double stop_resid = tol*r.norm();

    int n = 0;
    while (r.norm() > stop_resid && n < MAX_ITERATION) {
        x = c.array()-band_mult(N,m,m,x_old).array()/d;
        r = b-band_mult(MN,m,m,x);
        x_old = x;
        n++;
    }

    return std::make_tuple(x,n);
}

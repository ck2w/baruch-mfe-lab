#include <linear_iterate.h>
#include <Eigen/Dense>
#include <tuple>
#include <cassert>
#include <triangular_solve.h>
#include <band_mult.h>
#include <band_add.h>
#include <iostream>

std::tuple<Eigen::VectorXd, int> linear_iterate_triangular
                                               (const Eigen::VectorXd & x0,
                                                const Eigen::MatrixXd & M, 
                                                const Eigen::MatrixXd & N, 
                                                const Eigen::VectorXd & b,
                                                double tol,
                                                const Eigen::VectorXd & p)
{
    // IsAmerican if overriding vector for early exercise is given.
    bool isAmerican = (p.size()>0);
    Eigen::VectorXd x = x0, x_old = x0;
    // If isAmerican, use residual error to stop.
    Eigen::VectorXd r = b-(M+N)*x0;
    double stop_resid = tol*r.norm();
    // If not isAmerican, use diff of consecutive approximation to stop.
    Eigen::VectorXd diff = x0;  
    // This is prescribed by homework, might subject to change.
    double stop_diff = tol;

    int n = 0;
    //while ( ( ( !isAmerican && r.norm() > stop_resid ) || 
    //          ( isAmerican && diff.norm() > stop_diff ) ) &&
    while ( diff.norm() > stop_diff && n < MAX_ITERATION ) {
        x = forward_subst(M, b-N*x_old, p);
        r = b-(M+N)*x;
        diff = x-x_old;
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
                                                int m, double tol,
                                                const Eigen::VectorXd & p)
{
    // IsAmerican if overriding vector for early exercise is given.
    bool isAmerican = (p.size()>0);
    Eigen::ArrayXXd MN = band_add(M,m,0,N,0,m);
    Eigen::VectorXd x = x0, x_old = x0;
    // If isAmerican, use residual error to stop.
    Eigen::VectorXd r = b-band_mult(MN,m,m,x0);
    double stop_resid = tol*r.norm();
    // If not isAmerican, use diff of consecutive approximation to stop.
    Eigen::VectorXd diff = x0;  
    // This is prescribed by homework, might subject to change.
    double stop_diff = tol;

    int n = 0;
    while ( ( ( !isAmerican && r.norm() > stop_resid ) ||
              ( isAmerican && diff.norm() > stop_diff ) ) &&
            n < MAX_ITERATION ) {
        x = forward_subst_banded(M, b-band_mult(N,0,m,x_old), p);
        r = b-band_mult(MN,m,m,x);
        diff = x-x_old;
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

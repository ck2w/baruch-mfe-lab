#include <linear_solve.h>
#include <triangular_solve.h>
#include <cholesky.h>
#include <Eigen/Dense>
#include <cassert>

Eigen::VectorXd spd_solve(const Eigen::MatrixXd & A, 
                          const Eigen::VectorXd & b)
{
    int n = b.size();
    assert(A.rows() == n);
    assert(A.cols() == n);

    Eigen::MatrixXd U = cholesky(A);
    Eigen::VectorXd y = forward_subst(U.transpose(), b);
    Eigen::VectorXd x = backward_subst(U, y);
    
    return x;
}

Eigen::VectorXd banded_spd_solve(const Eigen::MatrixXd & A, int m,
                                 const Eigen::VectorXd & b)
{
    int n = b.size();
    assert(A.rows() == n);
    assert(A.cols() == n);

    Eigen::MatrixXd U = cholesky_banded(A, m);
    Eigen::VectorXd y = forward_subst_banded(U.transpose(), m, b);
    Eigen::VectorXd x = backward_subst_banded(U, m, y);
    
    return x;
}

Eigen::VectorXd tridiagonal_spd_solve(const Eigen::MatrixXd & A, 
                                      const Eigen::VectorXd & b)
{
    int n = b.size();
    assert(A.rows() == n);
    assert(A.cols() == n);

    Eigen::MatrixXd U = cholesky_tridiagonal(A);
    Eigen::VectorXd y = forward_subst_tridiagonal(U.transpose(), b);
    Eigen::VectorXd x = backward_subst_tridiagonal(U, y);
    
    return x;
}

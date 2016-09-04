#include <triangular_solve.h>
#include <Eigen/Dense>
#include <cassert>
#include <algorithm>

Eigen::VectorXd forward_subst_banded(const Eigen::MatrixXd & L, int m,
                                     const Eigen::VectorXd & b)
{
    int n = b.size();
    assert(L.rows() == n);
    assert(L.cols() == n);

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    x(0) = b(0)/L(0,0);
    for (int i=1; i<n; i++) {
        double sum = 0;
        int boundary = std::max(0,i-m);
        for (int j=boundary; j<i; j++) {
            sum += L(i,j)*x(j);
        }
        x(i) = (b(i)-sum)/L(i,i);
    }
    
    return x;
}

Eigen::VectorXd forward_subst(const Eigen::MatrixXd & L, 
                              const Eigen::VectorXd & b)
{
    int n = b.size();
    return forward_subst_banded(L, n, b);
}

Eigen::VectorXd forward_subst_tridiagonal(const Eigen::MatrixXd & L, 
                                          const Eigen::VectorXd & b)
{
    return forward_subst_banded(L, 1, b);
}

Eigen::VectorXd backward_subst_banded(const Eigen::MatrixXd & U, int m,
                                      const Eigen::VectorXd & b)
{
    int n = b.size();
    assert(U.rows() == n);
    assert(U.cols() == n);

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    x(n-1) = b(n-1)/U(n-1,n-1);
    for (int i=n-2; i>=0; i--) {
        double sum = 0;
        int boundary = std::min(n,i+m+1);
        for (int j=i+1; j<boundary; j++) {
            sum += U(i,j)*x(j);
        }
        x(i) = (b(i)-sum)/U(i,i);
    }

    return x;
}

Eigen::VectorXd backward_subst(const Eigen::MatrixXd & U, 
                               const Eigen::VectorXd & b)
{
    int n = b.size();
    return backward_subst_banded(U, n, b);
}

Eigen::VectorXd backward_subst_tridiagonal(const Eigen::MatrixXd & U, 
                                           const Eigen::VectorXd & b)
{
    return backward_subst_banded(U, 1, b);
}

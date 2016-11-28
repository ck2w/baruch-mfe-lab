#include <triangular_solve.h>
#include <Eigen/Dense>
#include <cassert>
#include <algorithm>
#include <iostream>

Eigen::VectorXd forward_subst_banded(const Eigen::MatrixXd & L, int m,
                                     const Eigen::VectorXd & b,
                                     const Eigen::VectorXd & p)
{
    int n = b.size();
    assert(L.rows() == n);
    assert(L.cols() == n);
    assert(p.size() == 0 || p.size() == n);

    bool isAmerican = (p.size()>0);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    for (int i=0; i<n; i++) {
        double sum = 0;
        int boundary = std::max(0,i-m);
        for (int j=boundary; j<i; j++) {
            sum += L(i,j)*x(j);
        }
        x(i) = (b(i)-sum)/L(i,i);
        // overriding with early exercise premium
        if (isAmerican && p(i)>x(i)) { x(i) = p(i); }
    }
    
    return x;
}

Eigen::VectorXd forward_subst(const Eigen::MatrixXd & L, 
                              const Eigen::VectorXd & b,
                              const Eigen::VectorXd & p)
{
    int n = b.size();
    return forward_subst_banded(L, n-1, b, p);
}

Eigen::VectorXd forward_subst_tridiagonal(const Eigen::MatrixXd & L, 
                                          const Eigen::VectorXd & b,
                                          const Eigen::VectorXd & p)
{
    return forward_subst_banded(L, 1, b, p);
}

// band storage implementation
Eigen::VectorXd forward_subst_banded(const Eigen::ArrayXXd & L,
                                     const Eigen::VectorXd & b,
                                     const Eigen::VectorXd & p)
{
    int n = b.size();
    assert(L.rows() == n);
    assert(p.size() == 0 || p.size() == n);
    int m = L.cols()-1;

    bool isAmerican = (p.size()>0);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    for (int i=0; i<n; i++) {
        double sum = 0;
        int boundary = std::max(0,i-m);
        for (int j=boundary;j<i;j++) {
            sum += L(i,m-i+j)*x(j);
        }
        x(i) = (b(i)-sum)/L(i,m);
        // overriding with early exercise premium
        if (isAmerican && p(i)>x(i)) { x(i) = p(i); }
    }

    return x;
}

Eigen::VectorXd backward_subst_banded(const Eigen::MatrixXd & U, int m,
                                      const Eigen::VectorXd & b,
                                      const Eigen::VectorXd & p)
{
    int n = b.size();
    assert(U.rows() == n);
    assert(U.cols() == n);
    assert(p.size() == 0 || p.size() == n);

    bool isAmerican = (p.size()>0);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    for (int i=n-1; i>=0; i--) {
        double sum = 0;
        int boundary = std::min(n,i+m+1);
        for (int j=i+1; j<boundary; j++) {
            sum += U(i,j)*x(j);
        }
        x(i) = (b(i)-sum)/U(i,i);
        // overriding with early exercise premium
        if (isAmerican && p(i)>x(i)) { x(i) = p(i); }
    }

    return x;
}

Eigen::VectorXd backward_subst(const Eigen::MatrixXd & U, 
                               const Eigen::VectorXd & b,
                               const Eigen::VectorXd & p)
{
    int n = b.size();
    return backward_subst_banded(U, n-1, b, p);
}

Eigen::VectorXd backward_subst_tridiagonal(const Eigen::MatrixXd & U, 
                                           const Eigen::VectorXd & b,
                                           const Eigen::VectorXd & p)
{
    return backward_subst_banded(U, 1, b, p);
}

// band storage implementation
Eigen::VectorXd backward_subst_banded(const Eigen::ArrayXXd & U,
                                      const Eigen::VectorXd & b,
                                      const Eigen::VectorXd & p)
{
    int n = b.size();
    assert(U.rows() == n);
    assert(p.size() == 0 || p.size() == n);
    int m = U.cols()-1;

    bool isAmerican = (p.size()>0);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    for (int i=n-1; i>=0; i--) {
        double sum = 0;
        int boundary = std::min(n,i+m+1);
        for (int j=i+1; j<boundary; j++) {
            sum += U(i,j-i)*x(j);
        }
        x(i) = (b(i)-sum)/U(i,0);
        // overriding with early exercise premium
        if (isAmerican && p(i)>x(i)) { x(i) = p(i); }
    }

    return x;
}

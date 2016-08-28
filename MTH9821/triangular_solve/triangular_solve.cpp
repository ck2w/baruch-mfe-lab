#include <triangular_solve.h>
#include <Eigen/Dense>
#include <cassert>

Eigen::VectorXd forward_subst(const Eigen::MatrixXd & L, 
                              const Eigen::VectorXd & b)
{
    int n = b.size();
    assert(L.rows() == n);
    assert(L.cols() == n);

    Eigen::VectorXd x(n);
    x(0) = b(0)/L(0,0);
    for (int i=1; i<n; i++) {
        double sum = 0;
        for (int j=0; j<i; j++) {
            sum += L(i,j)*x(j);
        }
        x(i) = (b(i)-sum)/L(i,i);
    }
    
    return x;
}

Eigen::VectorXd backward_subst(const Eigen::MatrixXd & U, 
                               const Eigen::VectorXd & b)
{
    int n = b.size();
    assert(U.rows() == n);
    assert(U.cols() == n);

    Eigen::VectorXd x(n);
    x(n-1) = b(n-1)/U(n-1,n-1);
    for (int i=n-2; i>=0; i--) {
        double sum = 0;
        for (int j=i+1; j<n; j++) {
            sum += U(i,j)*x(j);
        }
        x(i) = (b(i)-sum)/U(i,i);
    }

    return x;
}


#include <cholesky.h>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <cassert>

Eigen::MatrixXd cholesky_banded(const Eigen::MatrixXd & A, int m)
{
    Eigen::MatrixXd Acopy = A;
    int n = Acopy.rows();
    assert(n == Acopy.cols());
    assert(Acopy == Acopy.transpose());

    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(n,n);

    for (int k=0; k<n-1; k++) {
        assert(Acopy(k,k)>0);

        U(k,k) = std::sqrt(Acopy(k,k));

        int boundary = std::min(k+m+1,n);
        for (int i=k+1; i<boundary; i++) {
            U(k,i) = Acopy(k,i)/U(k,k);
        }

        for (int i=k+1; i<boundary; i++) {
            for (int j=i; j<boundary; j++) {
                Acopy(i,j) -= U(k,i)*U(k,j);
            }
        }
    }

    U(n-1,n-1) = std::sqrt(Acopy(n-1,n-1));

    return U; 
}

Eigen::MatrixXd cholesky(const Eigen::MatrixXd & A)
{
    int n = A.rows();
    return cholesky_banded(A, n); 
}

Eigen::MatrixXd cholesky_tridiagonal(const Eigen::MatrixXd & A)
{
    return cholesky_banded(A, 1);
}


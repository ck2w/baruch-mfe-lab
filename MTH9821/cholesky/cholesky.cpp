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

Eigen::MatrixXd cholesky_banded(const Eigen::ArrayXXd & A, int m)
{
    Eigen::ArrayXXd Acopy = A;
    int n = Acopy.rows();
    assert(m+1 == Acopy.cols());

    // upper triangular matrix
    // with lower band width = 0
    // with upper band width = m
    Eigen::ArrayXXd U = Eigen::ArrayXXd::Zero(n,m+1);

    for (int k=0; k<n-1; k++) {
        assert(Acopy(k,0)>0);

        U(k,0) = std::sqrt(Acopy(k,0));

        int boundary = std::min(k+m+1,n);
        for (int i=k+1; i<boundary; i++) {
            U(k,i-k) = Acopy(k,i-k)/U(k,0);
        }

        for (int i=k+1; i<boundary; i++) {
            for (int j=i; j<boundary; j++) {
                Acopy(i,j-i) -= U(k,i-k)*U(k,j-k);
            }
        }
    }

    U(n-1,0) = std::sqrt(Acopy(n-1,0));

    return U; 
}

Eigen::MatrixXd cholesky(const Eigen::MatrixXd & A)
{
    int n = A.rows();
    return cholesky_banded(A, n-1); 
}

Eigen::MatrixXd cholesky_tridiagonal(const Eigen::MatrixXd & A)
{
    return cholesky_banded(A, 1);
}


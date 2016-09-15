#include <band_conversion.h>
#include <Eigen/Dense>
#include <cassert>
#include <algorithm>

Eigen::MatrixXd dense_from_band(const Eigen::ArrayXXd & a, int m1, int m2)
{
    int n = a.rows();
    assert( a.cols() == m1+m2+1 );

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n,n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            if (i-j>m1 || j-i>m2) {
                continue;
            }
            A(i,j) = a(i,m1-i+j);
        }
    }

    return A;
}

Eigen::ArrayXXd band_from_dense(const Eigen::MatrixXd & A, int m1, int m2)
{
    int n = A.rows();
    assert(n == A.cols());

    Eigen::ArrayXXd a = Eigen::ArrayXXd::Zero(n,m1+m2+1);
    a.col(m1) = A.diagonal();

    for (int i=1; i<=m1; i++) {
        a.col(m1-i).block(i,0,n-i,1) = A.diagonal(-i);
    }

    for (int i=1; i<=m2; i++) {
        a.col(m1+i).block(0,0,n-i,1) = A.diagonal(i);
    }

    return a;
}

Eigen::ArrayXXd band_transpose(const Eigen::ArrayXXd & a, int m1, int m2)
{
    int n = a.rows();
    int m = m1+m2+1;
    assert( a.cols() == m );

    Eigen::ArrayXXd b = Eigen::ArrayXXd::Zero(n,m);

    b.col(m2) = a.col(m1);
    for (int i=1; i<=m1; i++) {
        b.col(m2+i).block(0,0,n-i,1) = a.col(m1-i).block(i,0,n-i,1);
    }
    for (int i=1; i<=m2; i++) {
        b.col(m2-i).block(i,0,n-i,1) = a.col(m1+i).block(0,0,n-i,1);
    }

    return b;
}

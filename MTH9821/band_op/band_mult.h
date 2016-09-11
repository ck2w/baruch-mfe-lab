#ifndef BAND_OP_H
#define BAND_OP_H 

#include <Eigen/Dense>

// The array a[0,...,n-1][0,...,m1+m2] stores
// a matrix A as follows: 
//
// The diagonal elements are in a[0,...n-1][m1]. 
//
// Subdiagonal elements are in a[j,...,n-1][0,...,m1-1] 
// with j>0 appropriate to the number of elements on each subdiagonal;
//
// Superdiagonal elements are in a[0,...,j][m1+1,...,m1+m2]
// with j<n-1 appropriate to the number of elements on each superdiagonal.
//
// m1: width of lower band
// m2: width of upper band
//
// x: a vector of size n
// y: matrix-vector multipliation y=Ax
Eigen::VectorXd band_mult(const Eigen::ArrayXXd & a, int m1, int m2,
                          const Eigen::VectorXd & x);

#endif /* BAND_OP_H */


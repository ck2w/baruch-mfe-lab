#ifndef BAND_CONVERSION_H
#define BAND_CONVERSION_H 

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
Eigen::MatrixXd dense_from_band(const Eigen::ArrayXXd & a, int m1, int m2);
Eigen::ArrayXXd band_from_dense(const Eigen::MatrixXd & A, int m1, int m2);

#endif /* BAND_CONVERSION_H */


#ifndef CHOLESKY_H
#define CHOLESKY_H 

#include <Eigen/Dense>
#include <tuple>

// A: symmetric positive definite matrix
// U: Cholesky decomposition of A=U'U
Eigen::MatrixXd cholesky(const Eigen::MatrixXd &A);

// A: symmetric positive definite matrix with band m
// U: Cholesky decomposition of A=U'U
Eigen::MatrixXd cholesky_banded(const Eigen::MatrixXd &A, int m);

// A: tridiagonal symmetric positive definite matrix
// U: Cholesky decomposition of A=U'U
Eigen::MatrixXd cholesky_tridiagonal(const Eigen::MatrixXd &A);

#endif /* CHOLESKY_H */

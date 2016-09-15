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

// A: symmetric positive definite matrix with band m
// U: Cholesky decomposition of A=U'U
// Note: only the upper triangular part of A and U are stored
Eigen::MatrixXd cholesky_banded(const Eigen::ArrayXXd &A, int m);

// A: tridiagonal symmetric positive definite matrix
// U: Cholesky decomposition of A=U'U
Eigen::MatrixXd cholesky_tridiagonal(const Eigen::MatrixXd &A);

#endif /* CHOLESKY_H */

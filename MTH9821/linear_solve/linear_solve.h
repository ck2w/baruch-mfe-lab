#ifndef LINEAR_SOLVE_H
#define LINEAR_SOLVE_H 

#include <Eigen/Dense>

// A: non-singular symmetric positive definite matrix
// m: band width of A if A is banded
// b: column vector
// x: solution to Ax=b
Eigen::VectorXd spd_solve(const Eigen::MatrixXd & A, 
                          const Eigen::VectorXd & b);

Eigen::VectorXd banded_spd_solve(const Eigen::MatrixXd & A, int m, 
                                 const Eigen::VectorXd & b);

// use upper-triangular band storage for A
Eigen::VectorXd banded_spd_solve(const Eigen::ArrayXXd & A, int m, 
                                 const Eigen::VectorXd & b);

Eigen::VectorXd tridiagonal_spd_solve(const Eigen::MatrixXd & A,
                                      const Eigen::VectorXd & b);

#endif /* LINEAR_SOLVE_H */

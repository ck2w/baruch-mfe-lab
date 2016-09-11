#ifndef JACOBI_H
#define JACOBI_H 

#include <Eigen/Dense>
#include <tuple>

// A: non-singular symmetric positive definite matrix
// m: band width of A if A is banded
// b: column vector
// x: solution to Ax=b
// tol: tolerance factor
std::tuple<Eigen::VectorXd, int> jacobi(const Eigen::MatrixXd & A, 
                                        const Eigen::VectorXd & b,
                                        double tol);

std::tuple<Eigen::VectorXd, int> jacobi(const Eigen::ArrayXXd & A, int m, 
                                        const Eigen::VectorXd & b, 
                                        double tol);

#endif /* JACOBI_H */

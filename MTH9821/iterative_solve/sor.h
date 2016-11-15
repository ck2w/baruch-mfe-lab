#ifndef SOR_H
#define SOR_H 

#include <Eigen/Dense>
#include <tuple>

// omega: relaxation parameter
// A: non-singular symmetric positive definite matrix
// m: band width of A if A is banded
// b: column vector
// x: solution to Ax=b
// tol: tolerance factor
std::tuple<Eigen::VectorXd, int> sor(double omega,
                                     const Eigen::MatrixXd & A, 
                                     const Eigen::VectorXd & b,
                                     double tol,
                                     const Eigen::VectorXd & x0
                                          = Eigen::VectorXd::Zero(0));

std::tuple<Eigen::VectorXd, int> sor(double omega,
                                     const Eigen::ArrayXXd & A, int m, 
                                     const Eigen::VectorXd & b, 
                                     double tol,
                                     const Eigen::VectorXd & x0
                                          = Eigen::VectorXd::Zero(0));

std::tuple<Eigen::VectorXd, int> gs(const Eigen::MatrixXd & A, 
                                    const Eigen::VectorXd & b,
                                    double tol,
                                    const Eigen::VectorXd & x0
                                          = Eigen::VectorXd::Zero(0));

std::tuple<Eigen::VectorXd, int> gs(const Eigen::ArrayXXd & A, int m,
                                    const Eigen::VectorXd & b,
                                    double tol,
                                    const Eigen::VectorXd & x0
                                          = Eigen::VectorXd::Zero(0));

#endif /* SOR_H */

#ifndef TRIANGULAR_SOLVE_H
#define TRIANGULAR_SOLVE_H 

#include <Eigen/Dense>

// L: non-singular lower triangular matrix
// b: column vector
// x: solution to Lx=b
Eigen::VectorXd forward_subst(const Eigen::MatrixXd & L, 
                              const Eigen::VectorXd & b);

// L: non-singular upper triangular matrix
// b: column vector
// x: solution to Ux=b
Eigen::VectorXd backward_subst(const Eigen::MatrixXd & U, 
                               const Eigen::VectorXd & b);

#endif /* TRIANGULAR_SOLVE_H */

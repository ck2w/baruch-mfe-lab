#ifndef TRIANGULAR_SOLVE_H
#define TRIANGULAR_SOLVE_H 

#include <Eigen/Dense>

// L: non-singular lower triangular matrix
// m: band width of L if L is banded
// b: column vector
// x: solution to Lx=b
Eigen::VectorXd forward_subst(const Eigen::MatrixXd & L, 
                              const Eigen::VectorXd & b);

Eigen::VectorXd forward_subst_banded(const Eigen::MatrixXd & L, int m,
                                     const Eigen::VectorXd & b);

Eigen::VectorXd forward_subst_tridiagonal(const Eigen::MatrixXd & L, 
                                          const Eigen::VectorXd & b);

// banded lower-triangular matrix, storage scheme defined in band_mult
Eigen::VectorXd forward_subst_banded(const Eigen::ArrayXXd & L,
                                     const Eigen::VectorXd & b);

// U: non-singular upper triangular matrix
// m: band width of U if U is banded
// b: column vector
// x: solution to Ux=b
Eigen::VectorXd backward_subst(const Eigen::MatrixXd & U, 
                               const Eigen::VectorXd & b);

Eigen::VectorXd backward_subst_banded(const Eigen::MatrixXd & U, int m,
                                      const Eigen::VectorXd & b);

Eigen::VectorXd backward_subst_tridiagonal(const Eigen::MatrixXd & U, 
                                           const Eigen::VectorXd & b);

// banded upper-triangular matrix, storage scheme defined in band_mult
Eigen::VectorXd backward_subst_banded(const Eigen::ArrayXXd & U,
                                      const Eigen::VectorXd & b);

#endif /* TRIANGULAR_SOLVE_H */

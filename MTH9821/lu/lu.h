#ifndef LU_H
#define LU_H 

#include <Eigen/Dense>
#include <tuple>

// A: non-singular matrix
// (L,U): LU decomposition of A
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> 
    lu_no_pivoting(const Eigen::MatrixXd &A);

// A: non-singular matrix
// (L,U): LU decomposition of A
std::tuple<Eigen::VectorXi, Eigen::MatrixXd, Eigen::MatrixXd> 
    lu_row_pivoting(const Eigen::MatrixXd &A);

#endif /* LU_H */

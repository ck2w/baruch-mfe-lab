#ifndef LAPLACE2D_H
#define LAPLACE2D_H

#include <Eigen/Dense>
#include <tuple>

// Successive over relaxation solution of 2-dimensional Laplace
// equation with Chebyshev acceleration.
//
// The arrays a, b, c, d, e, and f are input as the coefficients
// of th equation, each dimensioned to the grid size imax-jmax.
//
// a(i,j)u(i+1,j) + b(i,j)u(i-1,j) + c(i,j)u(i,j+1) + d(i,j)u(i,j-1)
// + e(i,j)u(i,j) = f(i,j)
//
// The array u is input as the initial guess to the solution,
// usually set to zero in the interior and boundary conditions 
// on the boudanry. Returns with the final value.
//
// The number omega is input as the relaxation parameter, typically
// determined from the spectral radius of the Jacobi iteration.
//
// The bool chebyshev determines whether Chebyshev acceleration
// is used.
//
// Returns (# of iterations, absolute residual, relative residual).

std::tuple<int, double, double> laplace2d(const Eigen::ArrayXXd & a,
                                          const Eigen::ArrayXXd & b,
                                          const Eigen::ArrayXXd & c,
                                          const Eigen::ArrayXXd & d,
                                          const Eigen::ArrayXXd & e,
                                          const Eigen::ArrayXXd & f,
                                          Eigen::ArrayXXd * u,
                                          double omega=1,
                                          double tol=1e-6,
                                          bool chebyshev=false);

#endif /* LAPLACE2D_H */

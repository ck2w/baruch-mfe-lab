#include <linear_iterate.h>
#include <sor.h>
#include <jacobi.h>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <tuple>
#include <iostream>

class IterativeSolveTest : public ::testing::Test
{
    protected:

        Eigen::MatrixXd get_diagonally_dominant_matrx(int n)
        {
            // create a stricly diagonally dominant matrix
            Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
            A.diagonal().setZero();
            Eigen::VectorXd a = A.cwiseAbs().rowwise().sum();
            A += a.asDiagonal();
            Eigen::VectorXd d = Eigen::VectorXd::Random(4).array().abs();
            A += d.asDiagonal();

            return A;
        }

        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(IterativeSolveTest, TriangularLinearIterate)
{
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(4);
    Eigen::VectorXd b = Eigen::VectorXd::Random(4);
    Eigen::MatrixXd A = get_diagonally_dominant_matrx(4);
            
    // Gauss-Seidel decomposition
    Eigen::MatrixXd L = A;
    L.triangularView<Eigen::StrictlyUpper>().setZero();
    Eigen::MatrixXd U = A; 
    U.triangularView<Eigen::Lower>().setZero();

    double tol = 1e-6;
    // linear iterate
    std::tuple<Eigen::VectorXd, int> res = linear_iterate(x0,L,U,b,tol);
    Eigen::VectorXd x = std::get<0>(res);
    int n = std::get<1>(res);

    EXPECT_LT(n, 20);
    EXPECT_NEAR((A*x-b).norm(), 0, tol*(A*x0-b).norm());

    // Gauss-Seidel iteration
    res = gs(A,b,tol);
    Eigen::VectorXd y = std::get<0>(res);
    int m = std::get<1>(res);
    
    EXPECT_LT(m, 20);
    EXPECT_NEAR((A*y-b).norm(), 0, tol*b.norm());
    EXPECT_NEAR((y-x).norm(), 0, tol*b.norm());

    // Gauss-Seidel is a special case of triangular iteration
    EXPECT_EQ(m,n);
}

TEST_F(IterativeSolveTest, DiagonalLinearIterate)
{
    Eigen::VectorXd x0 = Eigen::VectorXd::Random(4);
    Eigen::VectorXd b = Eigen::VectorXd::Random(4);
    Eigen::MatrixXd A = get_diagonally_dominant_matrx(4);

    // Jacobi decomposition
    Eigen::VectorXd dd = A.diagonal();
    Eigen::MatrixXd N = A;
    N.diagonal().setZero();

    double tol = 1e-6;
    // linear iterate
    std::tuple<Eigen::VectorXd, int> res = linear_iterate(x0,dd,N,b,tol);
    Eigen::VectorXd x = std::get<0>(res);
    int n = std::get<1>(res);

    EXPECT_LT(n, 20);
    EXPECT_NEAR((A*x-b).norm(), 0, tol*(A*x0-b).norm());
}

TEST_F(IterativeSolveTest, JacobiSolve)
{
    Eigen::VectorXd b = Eigen::VectorXd::Random(4);
    Eigen::MatrixXd A = get_diagonally_dominant_matrx(4);

    double tol = 1e-6;
    std::tuple<Eigen::VectorXd, int> res = jacobi(A,b,tol);
    Eigen::VectorXd x = std::get<0>(res);
    int n = std::get<1>(res);

    EXPECT_LT(n, 1000000);
    EXPECT_NEAR((A*x-b).norm(), 0, tol*b.norm());
}

TEST_F(IterativeSolveTest, SorSolve)
{
    Eigen::VectorXd b = Eigen::VectorXd::Random(4);
    Eigen::MatrixXd A = get_diagonally_dominant_matrx(4);

    double omega = 1.5;
    double tol = 1e-6;
    std::tuple<Eigen::VectorXd, int> res = sor(omega,A,b,tol);
    Eigen::VectorXd x = std::get<0>(res);
    int n = std::get<1>(res);

    EXPECT_LT(n, 1000000);
    EXPECT_NEAR((A*x-b).norm(), 0, tol*b.norm());
}


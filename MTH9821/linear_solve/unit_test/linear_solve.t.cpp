#include <linear_solve.h>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <iostream>

class LinearSolveTest : public ::testing::Test
{
    protected:

        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(LinearSolveTest, SpdLinearSolveVerification)
{
    Eigen::MatrixXd A(4,4);
    A <<  9, -3,  6, -3,
         -3,  5, -4,  7,
          6, -4, 21,  3,
         -3,  7,  3, 15;

    Eigen::VectorXd b = Eigen::VectorXd::Random(4);
    Eigen::VectorXd x = spd_solve(A, b);
    double tol = 1e-13;
    EXPECT_NEAR((A*x-b).norm(), 0, tol);
}

TEST_F(LinearSolveTest, BandedSpdLinearSolveVerification)
{
    Eigen::MatrixXd A(4,4);
    A <<  9, -3,  6,  0,
         -3,  5, -4,  7,
          6, -4, 21,  3,
          0,  7,  3, 15;

    Eigen::VectorXd b = Eigen::VectorXd::Random(4);
    Eigen::VectorXd x = banded_spd_solve(A, 2, b);

    double tol = 1e-13;
    EXPECT_NEAR((A*x-b).norm(), 0, tol);
}

TEST_F(LinearSolveTest, TridiagonalSpdLinearSolveVerification)
{
    Eigen::MatrixXd A(4,4);
    A <<  9, -3,  0,  0,
         -3,  5, -4,  0,
          0, -4, 21,  3,
          0,  0,  3, 15;

    Eigen::VectorXd b = Eigen::VectorXd::Random(4);
    Eigen::VectorXd x = tridiagonal_spd_solve(A, b);

    double tol = 1e-16;
    EXPECT_NEAR((A*x-b).norm(), 0, tol);
}


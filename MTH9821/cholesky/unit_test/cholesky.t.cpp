#include <cholesky.h>
#include <gtest/gtest.h>
#include <Eigen/Dense>

class CholeskyTest : public ::testing::Test
{
    protected:

        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(CholeskyTest, CholeskyDecomposition4X4Matrix)
{
    Eigen::MatrixXd A(4,4);
    A <<  9, -3,  6, -3,
         -3,  5, -4,  7,
          6, -4, 21,  3,
         -3,  7,  3, 15;

    Eigen::MatrixXd U = cholesky(A);

    double tol = 1e-16;
    EXPECT_NEAR((U.transpose()*U-A).norm(), 0, tol);
}

TEST_F(CholeskyTest, BandedCholeskyDecomposition4X4Matrix)
{
    Eigen::MatrixXd A(4,4);
    A <<  9, -3,  6,  0,
         -3,  5, -4,  7,
          6, -4, 21,  3,
          0,  7,  3, 15;

    Eigen::MatrixXd U = cholesky_banded(A,2);

    double tol = 1e-16;
    EXPECT_NEAR((U.transpose()*U-A).norm(), 0, tol);
}

TEST_F(CholeskyTest, TridiagonalCholeskyDecomposition4X4Matrix)
{
    Eigen::MatrixXd A(4,4);
    A <<  9, -3,  0,  0,
         -3,  5, -4,  0,
          0, -4, 21,  3,
          0,  0,  3, 15;

    Eigen::MatrixXd U = cholesky_tridiagonal(A);

    double tol = 1e-16;
    EXPECT_NEAR((U.transpose()*U-A).norm(), 0, tol);
}

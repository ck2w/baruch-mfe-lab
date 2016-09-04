#include <triangular_solve.h>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <iostream>

class TriangularSolveTest : public ::testing::Test
{
    protected:

        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(TriangularSolveTest, ForwardSubstMatchEigenComputedResults)
{
    Eigen::MatrixXd L = Eigen::MatrixXd::Random(2,2);
    Eigen::VectorXd b = Eigen::VectorXd::Random(2);
    L.triangularView<Eigen::StrictlyUpper>().setZero();
    Eigen::VectorXd x = forward_subst(L, b);
    Eigen::VectorXd y = L.triangularView<Eigen::Lower>().solve(b);
    double tol = 1e-16;
    EXPECT_NEAR((x-y).norm(), 0, tol);
}

TEST_F(TriangularSolveTest, BackwardSubstMatchEigenComputedResults)
{
    Eigen::MatrixXd U = Eigen::MatrixXd::Random(2,2);
    Eigen::VectorXd b = Eigen::VectorXd::Random(2);
    U.triangularView<Eigen::StrictlyLower>().setZero();
    Eigen::VectorXd x = backward_subst(U, b);
    Eigen::VectorXd y = U.triangularView<Eigen::Upper>().solve(b);
    double tol = 1e-16;
    EXPECT_NEAR((x-y).norm(), 0, tol);
}

TEST_F(TriangularSolveTest, ForwardSubstBandedMatchEigenComputedResults)
{
    Eigen::MatrixXd L = Eigen::MatrixXd::Random(6,6);
    Eigen::VectorXd b = Eigen::VectorXd::Random(6);
    L.triangularView<Eigen::StrictlyUpper>().setZero();
    int band = 3;
    for (int j=0; j<6; j++) {
        for (int i=j; i<6; i++) {
            if (i-j>band) {
                L(i,j) = 0;
            }
        }
    }

    Eigen::VectorXd x = forward_subst_banded(L, band, b);
    Eigen::VectorXd y = L.triangularView<Eigen::Lower>().solve(b);
    double tol = 1e-13;
    EXPECT_NEAR((x-y).norm(), 0, tol);
    EXPECT_NEAR((L*x-b).norm(), 0, tol);
}

TEST_F(TriangularSolveTest, BackwardSubstBandedMatchEigenComputedResults)
{
    Eigen::MatrixXd U = Eigen::MatrixXd::Random(6,6);
    Eigen::VectorXd b = Eigen::VectorXd::Random(6);
    U.triangularView<Eigen::StrictlyLower>().setZero();
    int band = 3;
    for (int i=0; i<6; i++) {
        for (int j=i; j<6; j++) {
            if (j-i>band) {
                U(i,j) = 0;
            }
        }
    }
    Eigen::VectorXd x = backward_subst_banded(U, band, b);
    Eigen::VectorXd y = U.triangularView<Eigen::Upper>().solve(b);
    double tol = 1e-13;
    EXPECT_NEAR((x-y).norm(), 0, tol);
    EXPECT_NEAR((U*x-b).norm(), 0, tol);
}

TEST_F(TriangularSolveTest, ForwardSubstTridiagonalMatchEigenComputedResults)
{
    Eigen::MatrixXd L = Eigen::MatrixXd::Random(6,6);
    Eigen::VectorXd b = Eigen::VectorXd::Random(6);
    L.triangularView<Eigen::StrictlyUpper>().setZero();
    int band = 1;
    for (int j=0; j<6; j++) {
        for (int i=j; i<6; i++) {
            if (i-j>band) {
                L(i,j) = 0;
            }
        }
    }

    Eigen::VectorXd x = forward_subst_tridiagonal(L, b);
    Eigen::VectorXd y = L.triangularView<Eigen::Lower>().solve(b);
    double tol = 1e-15;
    EXPECT_NEAR((x-y).norm(), 0, tol);
    EXPECT_NEAR((L*x-b).norm(), 0, tol);
}

TEST_F(TriangularSolveTest, BackwardSubstTridiagonalMatchEigenComputedResults)
{
    Eigen::MatrixXd U = Eigen::MatrixXd::Random(6,6);
    Eigen::VectorXd b = Eigen::VectorXd::Random(6);
    U.triangularView<Eigen::StrictlyLower>().setZero();
    int band = 1;
    for (int i=0; i<6; i++) {
        for (int j=i; j<6; j++) {
            if (j-i>band) {
                U(i,j) = 0;
            }
        }
    }

    Eigen::VectorXd x = backward_subst_tridiagonal(U, b);
    Eigen::VectorXd y = U.triangularView<Eigen::Upper>().solve(b);
    double tol = 1e-15;
    EXPECT_NEAR((x-y).norm(), 0, tol);
    EXPECT_NEAR((U*x-b).norm(), 0, tol);
}

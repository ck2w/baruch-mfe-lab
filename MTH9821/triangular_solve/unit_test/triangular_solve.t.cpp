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


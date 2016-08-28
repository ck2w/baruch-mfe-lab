#include <lu.h>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <tuple>

class LuTest : public ::testing::Test
{
    protected:

        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(LuTest, LuNoPivoting)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(4,4);
    std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> res = lu_no_pivoting(A);
    Eigen::MatrixXd L = std::get<0>(res);
    Eigen::MatrixXd U = std::get<1>(res);
    double tol = 1e-16;
    EXPECT_NEAR((L*U-A).norm(), 0, tol);
}

TEST_F(LuTest, LuRowPivoting)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(4,4);
    std::tuple<Eigen::VectorXi, Eigen::MatrixXd, Eigen::MatrixXd> res 
        = lu_row_pivoting(A);
    Eigen::VectorXi p = std::get<0>(res);
    Eigen::MatrixXd L = std::get<1>(res);
    Eigen::MatrixXd U = std::get<2>(res);

    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(4,4);
    for (int i=0; i<4; i++) {
        P(i,p(i)-1) = 1;
    }

    double tol = 2e-16;
    EXPECT_NEAR((L*U-P*A).norm(), 0, tol);
}


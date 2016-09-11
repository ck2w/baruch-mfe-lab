#include <band_mult.h>
#include <band_add.h>
#include <band_conversion.h>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <iostream>
#include <algorithm>

class BandOpTest : public ::testing::Test
{
    protected:
        
        void make_banded(Eigen::MatrixXd * A, int m1, int m2)
        {
            int n = A->rows();
            for (int i=0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    if ( i-j>m1 || j-i>m2 ) {
                        (*A)(i,j) = 0;
                    }
                }
            }
        }

        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(BandOpTest, CheckBandedMatrixAddition)
{
    int n=6;
    int m1=3;
    int m2=2;
    int n1=1;
    int n2=4;
    int l1=std::max(m1,n1);
    int l2=std::max(m2,n2);

    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(n,n);

    make_banded(&A,m1,m2);
    make_banded(&B,n1,n2);
    
    Eigen::ArrayXXd a = band_from_dense(A,m1,m2);
    Eigen::ArrayXXd b = band_from_dense(B,n1,n2);

    Eigen::ArrayXXd c = band_add(a,m1,m2,b,n1,n2);

    Eigen::MatrixXd C = dense_from_band(c,l1,l2);
    double tol = 1e-16;
    EXPECT_NEAR(C.norm(),(A+B).norm(), tol);
}

TEST_F(BandOpTest, CheckResultsWithDenseMult)
{
    int n=5;
    int m1=3;
    int m2=2;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
    Eigen::VectorXd b = Eigen::VectorXd::Random(n);

    // create band matrix
    make_banded(&A,m1,m2);

    Eigen::ArrayXXd a = band_from_dense(A,m1,m2);

    Eigen::VectorXd x = A*b;
    Eigen::VectorXd y = band_mult(a,m1,m2,b);
    
    double tol = 3e-16;
    EXPECT_NEAR((x-y).norm(), 0, tol);
}

TEST_F(BandOpTest, CheckBandMatrixConversion)
{
    int n=5;
    int m1=3;
    int m2=2;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
    Eigen::VectorXd b = Eigen::VectorXd::Random(n);

    // create band matrix
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            if (i-j>m1 || j-i>m2) {
                A(i,j) = 0;
            }
        }
    }

    Eigen::ArrayXXd a = band_from_dense(A,m1,m2);
    Eigen::MatrixXd Acopy = dense_from_band(a,m1,m2);

    double tol = 1e-16;
    EXPECT_NEAR((A-Acopy).norm(),0,tol);
}

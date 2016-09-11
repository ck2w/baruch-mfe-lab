#include <linear_iterate.h>
#include <sor.h>
#include <jacobi.h>
#include <band_mult.h>
#include <band_conversion.h>
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
            Eigen::VectorXd d = Eigen::VectorXd::Random(n).array().abs();
            A += d.asDiagonal();

            return A;
        }

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
    std::tuple<Eigen::VectorXd, int> res 
        = linear_iterate_triangular(x0,L,U,b,tol);
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

TEST_F(IterativeSolveTest, TriangularBandedLinearIterate)
{
    int n=6;
    int m=3;
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd b = Eigen::VectorXd::Random(n);
    Eigen::MatrixXd A = get_diagonally_dominant_matrx(n);
    make_banded(&A,m,m);
    
    // Gauss-Seidel decomposition, dense version
    Eigen::MatrixXd L = A;
    L.triangularView<Eigen::StrictlyUpper>().setZero();
    Eigen::MatrixXd U = A; 
    U.triangularView<Eigen::Lower>().setZero();

    double tol = 1e-6;
    // linear iterate
    std::tuple<Eigen::VectorXd, int> res1 
        = linear_iterate_triangular(x0,L,U,b,tol);
    Eigen::VectorXd x1 = std::get<0>(res1);
    int step1 = std::get<1>(res1);

    EXPECT_LT(step1, 20);
    EXPECT_NEAR((A*x1-b).norm(), 0, tol*(A*x0-b).norm());

    // banded version
    Eigen::ArrayXXd M = band_from_dense(L,m,0);
    Eigen::ArrayXXd N = band_from_dense(U,0,m);

    std::tuple<Eigen::VectorXd, int> res2 
        = linear_iterate_triangular_banded(x0,M,N,b,m,tol);
    Eigen::VectorXd x2 = std::get<0>(res2);
    int step2 = std::get<1>(res2);
    
    EXPECT_EQ(step1, step2);
    EXPECT_NEAR((A*x2-b).norm(), 0, tol*(A*x0-b).norm());
    EXPECT_NEAR((x2-x1).norm(), 0, 1e-16);
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
    std::tuple<Eigen::VectorXd, int> res 
        = linear_iterate_diagonal(x0,dd,N,b,tol);
    Eigen::VectorXd x = std::get<0>(res);
    int n = std::get<1>(res);

    EXPECT_LT(n, 20);
    EXPECT_NEAR((A*x-b).norm(), 0, tol*(A*x0-b).norm());
}

TEST_F(IterativeSolveTest, DiagonalBandedLinearIterate)
{
    int n=6;
    int m=2;
    Eigen::VectorXd x0 = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = Eigen::VectorXd::Random(n);
    Eigen::MatrixXd A = get_diagonally_dominant_matrx(n);
    make_banded(&A,m,m);

    // Jacobi decomposition
    Eigen::VectorXd dd = A.diagonal();
    Eigen::MatrixXd N = A;
    N.diagonal().setZero();

    double tol = 1e-6;
    // linear iterate
    std::tuple<Eigen::VectorXd, int> res1 
        = linear_iterate_diagonal(x0,dd,N,b,tol);
    Eigen::VectorXd x1 = std::get<0>(res1);
    int step1 = std::get<1>(res1);

    EXPECT_LT(n, 20);
    EXPECT_NEAR((A*x1-b).norm(), 0, tol*(A*x0-b).norm());
   
    // banded version
    Eigen::ArrayXd d = dd.array();
    Eigen::ArrayXXd B = band_from_dense(N,m,m);
    
    std::tuple<Eigen::VectorXd, int> res2 
        = linear_iterate_diagonal_banded(x0,d,B,b,m,tol);
    Eigen::VectorXd x2 = std::get<0>(res2);
    int step2 = std::get<1>(res2);
    
    EXPECT_EQ(step1, step2);
    EXPECT_NEAR((A*x2-b).norm(), 0, tol*(A*x0-b).norm());
    EXPECT_NEAR((x2-x1).norm(), 0, 1e-16);
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

TEST_F(IterativeSolveTest, JacobiBandedSolve)
{
    int n=6;
    int m=2;
    Eigen::VectorXd b = Eigen::VectorXd::Random(n);
    Eigen::MatrixXd A = get_diagonally_dominant_matrx(n);
    make_banded(&A,m,m);
    Eigen::ArrayXXd B = band_from_dense(A,m,m);

    double tol = 1e-6;
    std::tuple<Eigen::VectorXd, int> res = jacobi(B,m,b,tol);
    Eigen::VectorXd x = std::get<0>(res);
    int step = std::get<1>(res);

    EXPECT_LT(step, 1000000);
    EXPECT_NEAR((A*x-b).norm(), 0, tol*b.norm());
    EXPECT_NEAR((band_mult(B,m,m,x)-b).norm(), 0, tol*b.norm());
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

TEST_F(IterativeSolveTest, SorBandedSolve)
{
    int n=6;
    int m=2;
    Eigen::VectorXd b = Eigen::VectorXd::Random(n);
    Eigen::MatrixXd A = get_diagonally_dominant_matrx(n);
    make_banded(&A,m,m);
    Eigen::ArrayXXd B = band_from_dense(A,m,m);

    double omega = 1.5;
    double tol = 1e-6;
    std::tuple<Eigen::VectorXd, int> res = sor(omega,B,m,b,tol);
    Eigen::VectorXd x = std::get<0>(res);
    int step= std::get<1>(res);

    EXPECT_LT(step, 1000000);
    EXPECT_NEAR((A*x-b).norm(), 0, tol*b.norm());
}


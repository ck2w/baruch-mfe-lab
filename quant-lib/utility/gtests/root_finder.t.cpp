#include <root_finder.h>
#include <gtest/gtest.h>
#include <cmath>

double function(double x)
{
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x2*x2;

    return (x4-5*x2+4-1.0/(1+std::exp(x3)));
}

double derivative(double x)
{
    double x2 = x*x;
    double x3 = x2*x;
    double ex3 = std::exp(x3);
    double y = (1+ex3)*(1+ex3);
    
    return (4*x3-10*x+3*x2*ex3/y);
}

class RootFinderTest : public ::testing::Test
{
    protected:
        
        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(RootFinderTest, NewtonShouldReturnInitialGuessIfDerivativeSetToNull)
{
    double tol=1e-9;
    double eps=1e-6;

    double (*null_derivative)(double);
    null_derivative = 0;
    Root_finder rf(&function, null_derivative);
    
    double initial_guess = -3;
    double x1 = rf.newton_find_root(initial_guess, tol, eps, false);
    EXPECT_EQ(initial_guess, x1);
}

TEST_F(RootFinderTest, NewtonShouldReturnInitialGuessIfDerivativeNotProvided)
{
    double tol=1e-9;
    double eps=1e-6;

    Root_finder rf(&function);
    
    double initial_guess = -3;
    double x1 = rf.newton_find_root(initial_guess, tol, eps, false);
    EXPECT_EQ(initial_guess, x1);
}

TEST_F(RootFinderTest, NewtonShouldFindCorrectRootForMoreComplicatedEquation)
{
    double tol=1e-9;
    double eps=1e-6;
    Root_finder rf(&function, &derivative);
    
    double x1 = rf.newton_find_root(-3, tol, eps, false);
    EXPECT_NEAR( -2.074304, x1, eps);
    EXPECT_NEAR( 0.0, function(x1), tol);
    
    double x2 = rf.newton_find_root(-0.5, tol, eps, false);
    EXPECT_NEAR( -0.889642, x2, eps);
    EXPECT_NEAR( 0.0, function(x2), tol);

    double x3 = rf.newton_find_root(0.5, tol, eps, false);
    EXPECT_NEAR( 0.950748, x3, eps);
    EXPECT_NEAR( 0.0, function(x3), tol);
    
    double x4 = rf.newton_find_root(3, tol, eps, false);
    EXPECT_NEAR( 2.000028, x4, eps);
    EXPECT_NEAR( 0.0, function(x4), tol);
}

TEST_F(RootFinderTest, SecantShouldFindCorrectRootForMoreComplicatedEquation)
{
    double tol=1e-9;
    double eps=1e-6;
    Root_finder rf(&function);
    
    double x1_second = -3;
    double x1_first = x1_second - 0.01;
    double x1 = rf.secant_find_root(x1_first, x1_second, tol, eps, false);
    EXPECT_NEAR( -2.074304, x1, eps);
    EXPECT_NEAR( 0.0, function(x1), tol);

    double x2_second = -0.5;
    double x2_first = x2_second - 0.01; 
    double x2 = rf.secant_find_root(x2_first, x2_second, tol, eps, false);
    EXPECT_NEAR( -0.889642, x2, eps);
    EXPECT_NEAR( 0.0, function(x2), tol);

    double x3_second = 0.5;
    double x3_first = x3_second - 0.01;
    double x3 = rf.secant_find_root(x3_first, x3_second, tol, eps, false);
    EXPECT_NEAR( 0.950748, x3, eps);
    EXPECT_NEAR( 0.0, function(x3), tol);
   
    double x4_second = 3;
    double x4_first = x4_second - 0.01; 
    double x4 = rf.secant_find_root(x4_first, x4_second, tol, eps, false);
    EXPECT_NEAR( 2.000028, x4, eps);
    EXPECT_NEAR( 0.0, function(x4), tol);
}

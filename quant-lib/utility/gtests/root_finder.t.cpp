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

TEST_F(RootFinderTest, ShouldFindCorrectRootForMoreComplicatedEquation)
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

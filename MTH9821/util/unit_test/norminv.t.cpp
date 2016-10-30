#include <norminv.h>
#include <gtest/gtest.h>
#include <iostream>
#include <cmath>
#include <vector>

class InvCdfTest : public ::testing::Test
{
    protected:
        
        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(InvCdfTest, CheckBeasleySpringerMoroAlgorithm)
{
    double tol = 1e-6;
    EXPECT_NEAR(norminv(0.5), 0.0, tol);
    EXPECT_NEAR(norminv(0.4), -0.2533471, tol);
    EXPECT_NEAR(norminv(0.6), 0.2533471, tol);
    EXPECT_NEAR(norminv(0.3), -0.5244005, tol);
    EXPECT_NEAR(norminv(0.7), 0.5244005, tol);
}


#include <european_option.h>
#include <constant_parameters.h>
#include <gtest/gtest.h>
#include <iostream>
#include <iomanip>

class BlackScholesTest : public ::testing::Test
{
    protected:

        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(BlackScholesTest, BlackScholesFormulaForCallShouldReturnCorrectValue)
{
    European_call c(40.0, 0.5);
    Constant_parameters vol(0.3);
    Constant_parameters div(0.03);
    Constant_parameters rate(0.05);
    double value= c.Black_Scholes_evaluate( 42.0,
                                            vol,
                                            div,
                                            rate); 

    EXPECT_NEAR(4.705325, value, 1e-6);
}

TEST_F(BlackScholesTest, BlackScholesFormulaForPutShouldReturnCorrectValue)
{
    European_put p(40.0, 0.5);
    Constant_parameters vol(0.3);
    Constant_parameters div(0.03);
    Constant_parameters rate(0.05);
    double value= p.Black_Scholes_evaluate( 42.0,
                                            vol,
                                            div,
                                            rate); 

    EXPECT_NEAR(2.343020, value, 1e-6);
}

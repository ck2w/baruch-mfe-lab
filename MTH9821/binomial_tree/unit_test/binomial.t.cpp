#include <payoff.h>
#include <option_value.h>
#include <binomial.h>
#include <black_scholes.h>
#include <gtest/gtest.h>
#include <iostream>
#include <tuple>

class BinomialTest : public ::testing::Test
{
    protected:

        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(BinomialTest, EuropeanCallTestCompareWithKnownResult)
{
    double T = 0.5;
    double K = 800;
    double S = 810;
    double r = 0.05;
    double q = 0.02;
    double vol = 0.2;
    CallPayoff call_payoff(K);
    BinomialTree b1(call_payoff, T, S, r, q, vol);
    double v = b1.evaluate(2).price;
    double tol = 1e-2;
    EXPECT_NEAR(v, 53.39, tol);
}

TEST_F(BinomialTest, AmericanCallTestCompareWithKnownResult)
{
    double T = 0.25;
    double K = 0.60;
    double S = 0.61;
    double r = 0.05;
    double q = 0.07;
    double vol = 0.12;
    CallPayoff call_payoff(K);
    BinomialTree b1(call_payoff, T, S, r, q, vol);
    double v_AmericanCall = b1.evaluate(3, true).price;
    double tol = 1e-3;
    EXPECT_NEAR(v_AmericanCall, 0.019, tol);
}

TEST_F(BinomialTest, AmericanCallEqualEuropeanCallWithZeroDividend)
{
    double T = 0.25;
    double K = 0.60;
    double S = 0.61;
    double r = 0.05;
    double q = 0.0;
    double vol = 0.12;
    CallPayoff call_payoff(K);
    BinomialTree b1(call_payoff, T, S, r, q, vol);
    double v_AmericanCall = b1.evaluate(3, true).price;
    double v_EuropeanCall = b1.evaluate(3, false).price;
    double tol = 1e-16;
    EXPECT_NEAR(v_AmericanCall, v_EuropeanCall, tol);
}

TEST_F(BinomialTest, AmericanPutTestCompareWithKnownResult)
{
    double T = 0.75;
    double K = 30;
    double S = 31;
    double r = 0.05;
    double q = 0.05;
    double vol = 0.3;
    PutPayoff put_payoff(K);
    BinomialTree b1(put_payoff, T, S, r, q, vol);
    double v_AmericanPut = b1.evaluate(3, true).price;
    double tol = 1e-2;
    EXPECT_NEAR(v_AmericanPut, 2.84, tol);
}

TEST_F(BinomialTest, BlackScholesTest)
{
    double T = 1;
    double K = 40;
    double S = 41;
    double r = 0.03;
    double q = 0.00;
    double vol = 0.3;

    std::tuple<OptionValue,OptionValue> res = BlackScholes(T,K,S,r,q,vol);
    double BsCall = std::get<0>(res).price;
    double BsPut = std::get<1>(res).price;
    double tol = 1e-2;
    EXPECT_NEAR(BsCall, 5.93, tol);
    EXPECT_NEAR(BsPut, 3.75, tol);
}

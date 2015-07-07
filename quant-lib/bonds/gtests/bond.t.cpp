#include <bond.h>
#include <time_dependent_parameters.h>
#include <gtest/gtest.h>
#include <cmath>

double constant_zero_rate(double x)
{
    return 0.05;
}

double logarithmic_zero_rate(double x)
{
    return (0.0525+std::log(1+2*x)/200);
}

double an_instantaneous_rate(double x)
{
    return (0.0525+0.01/(1+std::exp(-x*x)));
}

// the following instantaneous rate corresponds 
// to the logarithmic zero rate
double another_instantaneous_rate(double x)
{
    return (0.0525+std::log(1+2*x)/200+0.01*x/(1+2*x));
}

class BondTest : public ::testing::Test
{
    protected:
    
        const static double eps=1e-9;
        
        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(BondTest, ShouldHaveNoCouponsAfterConstructionWithParValueOnly)
{
    Bond b(1,100); 
    EXPECT_EQ( 0, b.get_number_of_coupons() );
}

TEST_F(BondTest, ShouldReturnCorrectCashFlowsAfterConstructionWithCashFlows)
{
    Bond b(1,100);
    b.add_coupon(0.25, 5);
    b.add_coupon(0.50, 5);
    b.add_coupon(0.75, 5);
    b.add_coupon(1.00, 5);

    EXPECT_EQ( 4, b.get_number_of_coupons() );
    b.print_cashflow(std::cout, 6);
}

TEST_F(BondTest, 
       ShouldReturnCorrectCashFlowsAfterConstructionWithCouponRateAndFrequency)
{
    Bond b1(1,100, 0.05, 2);
    b1.print_cashflow(std::cout, 6);
    
    Bond b2(1.25,100, 0.05, 4);
    b2.print_cashflow(std::cout, 6);
}

TEST_F(BondTest, 
       ShouldCalculateCorrectBondPriceGivenZeroRateCurveWithoutCoupons)
{
    Bond b(1, 100);
    Time_dependent_parameters zero_rates(&constant_zero_rate);
    EXPECT_NEAR( 100*std::exp(-0.05), 
                 b.evaluate_from_zero_rate(zero_rates), 
                 eps );
}

TEST_F(BondTest,
       ShouldCalculateCorrectSemiAnnualBondPriceGivenZeroRateCurve)
{
    double expiry = 20.0/12; // twenty months in units of years
    double par_value = 100.0;
    double coupon_rate = 0.06;
    int coupon_frequency = 2;
    Bond b(expiry, par_value, coupon_rate, coupon_frequency);

    double tol = 1e-6;
    Time_dependent_parameters zero_rates(&logarithmic_zero_rate);
    // The correct result can be found in Primer, Chapter 2.
    EXPECT_NEAR(101.888216, b.evaluate_from_zero_rate(zero_rates), tol);
}

TEST_F(BondTest,
       ShouldCalculateCorrectSemiAnnualBondPriceGivenInstantaneousRate)
{
    double expiry = 20.0/12; // twenty months in units of years
    double par_value = 100.0;
    double coupon_rate = 0.06;
    int coupon_frequency = 2;
    Bond b(expiry, par_value, coupon_rate, coupon_frequency);

    double tol = 1e-6;
    Time_dependent_parameters instantaneous_rates(&an_instantaneous_rate);
    // The correct result can be found in Primer, Chapter 2.
    EXPECT_NEAR(101.954564, 
                b.evaluate_from_instantaneous_rate(instantaneous_rates), tol);
}

TEST_F(BondTest,
       ZeroRateAndInstantaneousRateOfTheSameRateShouldGiveSameResult)
{
    double expiry = 20.0/12; // twenty months in units of years
    double par_value = 100.0;
    double coupon_rate = 0.06;
    int coupon_frequency = 2;

    Bond b(expiry, par_value, coupon_rate, coupon_frequency);

    Time_dependent_parameters zero_rates(&logarithmic_zero_rate);
    Time_dependent_parameters instantaneous_rates(&another_instantaneous_rate);

    double tol=1e-8;
    EXPECT_NEAR(b.evaluate_from_zero_rate(zero_rates),
                b.evaluate_from_instantaneous_rate(another_instantaneous_rate),
                tol);
}

TEST_F(BondTest,
       ShouldCalculateCorrectBondPriceDurationConvexityFromYield)
{
    double expiry = 20.0/12; // twenty months in units of years
    double par_value = 100.0;
    double coupon_rate = 0.06;
    int coupon_frequency = 2;
    double yield = 0.065;
    
    double bond_price = 0.0;
    double bond_duration = 0.0;
    double bond_convexity = 0.0;

    Bond b(expiry, par_value, coupon_rate, coupon_frequency);
    bond_price = b.evaluate_from_yield(yield, &bond_duration, &bond_convexity);

    double tol=1e-6;
    EXPECT_NEAR(101.046193, bond_price,     tol);
    EXPECT_NEAR( 1.5804216, bond_duration,  tol);
    EXPECT_NEAR( 2.5916859, bond_convexity, tol);
}

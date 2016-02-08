#include <time_dependent_parameters.h>
#include <gtest/gtest.h>
#include <cmath>

double constant_function(double x)
{
    return 10;
}

double exponential_function(double x)
{
    return std::exp(x);
}

class TimeDependentParametersTest : public ::testing::Test
{
    protected:

        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(TimeDependentParametersTest,
       IntegralMethodsShouldReturnCorrectValuesForExponentialFunction)
{
    Time_dependent_parameters p(&exponential_function);

    double eps = 1e-9;
    double time1 = 1.0;
    double time2 = 2.0;
    EXPECT_NEAR(std::exp(2)-std::exp(1), p.integral(time1, time2), eps);
    EXPECT_NEAR(std::exp(1)-std::exp(2), p.integral(time2, time1), eps);
}

TEST_F(TimeDependentParametersTest, 
       TimeProductMethodsShouldReturnCorrectValuesForTimeIndependentCase)
{
    Time_dependent_parameters p(&constant_function);

    double eps = 1e-9;
    double time1 = 2.5;
    double time2 = 5.0;
    EXPECT_NEAR( 25, p.time_product_diff(time1, time2), eps);
    EXPECT_NEAR(-25, p.time_product_diff(time2, time1), eps);
}

TEST_F(TimeDependentParametersTest, 
       TimeProductMethodsShouldReturnCorrectValuesForExponentialFunction)
{
    Time_dependent_parameters p(&exponential_function);

    double eps = 1e-9;
    double time1 = 0.0;
    double time2 = 1.0;
    EXPECT_NEAR( std::exp(1), p.time_product_diff(time1, time2), eps);
}

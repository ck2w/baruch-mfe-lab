#include <cumulative_normal.h>
#include <gtest/gtest.h>
#include <cmath>

class CumulativeNormalTest : public ::testing::Test
{
    protected:
        
        static const double eps = 7.5e-7;
        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(CumulativeNormalTest, StaticMemberFunctionShouldReturnExpectedValues)
{
    EXPECT_NEAR(0.5, Cumulative_normal::approximate_evaluate(0.0), eps);
}

TEST_F(CumulativeNormalTest, 
       ApproximateEvaluationShouldReturnApproximationsOfExactValues)
{
    EXPECT_NEAR(Cumulative_normal::exact_evaluate(0.0),
                Cumulative_normal::approximate_evaluate(0.0), eps);
    EXPECT_NEAR(Cumulative_normal::exact_evaluate(0.1),
                Cumulative_normal::approximate_evaluate(0.1), eps);
    EXPECT_NEAR(Cumulative_normal::exact_evaluate(0.5),
                Cumulative_normal::approximate_evaluate(0.5), eps);
    EXPECT_NEAR(Cumulative_normal::exact_evaluate(1.0),
                Cumulative_normal::approximate_evaluate(1.0), eps);
    EXPECT_NEAR(Cumulative_normal::exact_evaluate(-0.1),
                Cumulative_normal::approximate_evaluate(-0.1), eps);
    EXPECT_NEAR(Cumulative_normal::exact_evaluate(-0.5),
                Cumulative_normal::approximate_evaluate(-0.5), eps);
    EXPECT_NEAR(Cumulative_normal::exact_evaluate(-1.0),
                Cumulative_normal::approximate_evaluate(-1.0), eps);
}

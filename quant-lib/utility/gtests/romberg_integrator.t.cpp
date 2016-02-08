#include <romberg_integrator.h>
#include <gtest/gtest.h>
#include <cmath>

double step_function(double x)
{
    if ( x >= 0 )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

double quadratic_function(double x)
{
    return x*x;
}

double cubic_function(double x)
{
    return x*x*x;
}

double quartic_function(double x)
{
    double square = x*x;
    return square*square;
}

double normal_function(double x)
{
    return std::exp(-x*x);
}

class RombergIntegratorTest : public ::testing::Test
{
    protected:
        
        virtual void SetUp() {}
        virtual void TearDown() {}
};

TEST_F(RombergIntegratorTest, DefaultConstructor)
{
    Romberg_integrator r;

    EXPECT_EQ( 0.0, r.lower_end() );
    EXPECT_EQ( 0.0, r.upper_end() );

    EXPECT_EQ( 0, r.bisection_order() );
    EXPECT_EQ( 1, r.expansion_order() );
    EXPECT_EQ( 1, r.interval_number() );

    EXPECT_EQ( 0.0, r.interval_width() );
}

TEST_F(RombergIntegratorTest, ConstructorWithIntegrationRangeAndIntegrand)
{
    Romberg_integrator r( &step_function, -1, 1 );
    
    EXPECT_EQ( -1, r.lower_end() );
    EXPECT_EQ(  1, r.upper_end() );

    EXPECT_EQ( 1, r.integrand(0.5) );
    EXPECT_EQ( 1, r.integrand(2.0) );
    EXPECT_EQ( 1, r.integrand(0) );
    EXPECT_EQ( 0, r.integrand(-1) );
    
    EXPECT_EQ( 0, r.bisection_order() );
    EXPECT_EQ( 1, r.expansion_order() );
    EXPECT_EQ( 1, r.interval_number() );
    EXPECT_EQ( 2.0, r.interval_width() );
}

TEST_F(RombergIntegratorTest, ShouldNotBeAbleToRefineAfterDefaultConstruction)
{
    Romberg_integrator r1;
    EXPECT_FALSE(r1.refine_bisection());
    
    Romberg_integrator r2;
    EXPECT_FALSE(r2.refine_expansion());
}

TEST_F(RombergIntegratorTest, 
       ShouldReturnCorrectIntervalNumberAndWidthAfterRefineBisection)
{
    Romberg_integrator r( &quadratic_function, 0, 1);
    EXPECT_EQ( 1, r.interval_number() );
    EXPECT_EQ( 1, r.interval_width() );

    r.refine_bisection();
    EXPECT_EQ( 2, r.interval_number() );
    EXPECT_EQ( 0.5, r.interval_width() );
    
    r.refine_bisection();
    EXPECT_EQ( 4, r.interval_number() );
    EXPECT_EQ( 0.25, r.interval_width() );
    
    r.refine_bisection();
    EXPECT_EQ( 8, r.interval_number() );
    EXPECT_EQ( 0.125, r.interval_width() );
}

TEST_F(RombergIntegratorTest, 
       ShouldNotBeAbleToRefineExpansionExceedingBisection)
{
    Romberg_integrator r( &quadratic_function, 0, 1);
    EXPECT_FALSE(r.refine_expansion());
    EXPECT_TRUE( r.refine_bisection());
    EXPECT_TRUE( r.refine_expansion());
    EXPECT_FALSE(r.refine_expansion());
}

TEST_F(RombergIntegratorTest, 
       ShouldReturnCorrectResultsWhenRefineOnQuadraticFunction)
{
    Romberg_integrator r( &quadratic_function, 0, 1);
    EXPECT_EQ( 0, r.bisection_order() );
    EXPECT_EQ( 1, r.expansion_order() );

    r.refine_bisection();
    EXPECT_EQ( 1, r.bisection_order() );
    EXPECT_EQ( 1, r.expansion_order() );

    r.refine_expansion();
    EXPECT_EQ( 1, r.bisection_order() );
    EXPECT_EQ( 2, r.expansion_order() );
    
    r.refine_bisection();
    EXPECT_EQ( 2, r.bisection_order() );
    EXPECT_EQ( 2, r.expansion_order() );

    r.refine_expansion();
    EXPECT_EQ( 2, r.bisection_order() );
    EXPECT_EQ( 3, r.expansion_order() );

    double eps = 1e-6;
    EXPECT_NEAR(0.250000, r.get_value(0,0), eps);
    EXPECT_NEAR(0.312500, r.get_value(1,0), eps);
    EXPECT_NEAR(0.328125, r.get_value(2,0), eps);
    EXPECT_NEAR(0.500000, r.get_value(0,1), eps);
    EXPECT_NEAR(0.375000, r.get_value(1,1), eps);
    EXPECT_NEAR(0.343750, r.get_value(2,1), eps);
    EXPECT_NEAR(0.333333, r.get_value(1,2), eps);
    EXPECT_NEAR(0.333333, r.get_value(2,2), eps);
    EXPECT_NEAR(0.333333, r.get_value(2,3), eps);

    EXPECT_NEAR(0.328125, r.midpoint_value(),  eps);
    EXPECT_NEAR(0.343750, r.trapezoid_value(), eps);
    EXPECT_NEAR(0.333333, r.simpson_value(),   eps);
}

TEST_F(RombergIntegratorTest, 
       ShouldReturnCorrectResultsWhenRefineOnCubicFunction)
{
    Romberg_integrator r( &cubic_function, 0, 1);
    EXPECT_EQ( 0, r.bisection_order() );
    EXPECT_EQ( 1, r.expansion_order() );

    r.refine_bisection();
    EXPECT_EQ( 1, r.bisection_order() );
    EXPECT_EQ( 1, r.expansion_order() );

    r.refine_bisection();
    EXPECT_EQ( 2, r.bisection_order() );
    EXPECT_EQ( 1, r.expansion_order() );

    r.refine_expansion();
    EXPECT_EQ( 2, r.bisection_order() );
    EXPECT_EQ( 2, r.expansion_order() );
    
    r.refine_expansion();
    EXPECT_EQ( 2, r.bisection_order() );
    EXPECT_EQ( 3, r.expansion_order() );

    double eps = 1e-6;
    EXPECT_NEAR(0.125000, r.get_value(0,0), eps);
    EXPECT_NEAR(0.218750, r.get_value(1,0), eps);
    EXPECT_NEAR(0.242188, r.get_value(2,0), eps);
    EXPECT_NEAR(0.500000, r.get_value(0,1), eps);
    EXPECT_NEAR(0.312500, r.get_value(1,1), eps);
    EXPECT_NEAR(0.265625, r.get_value(2,1), eps);
    EXPECT_NEAR(0.250000, r.get_value(1,2), eps);
    EXPECT_NEAR(0.250000, r.get_value(2,2), eps);
    EXPECT_NEAR(0.250000, r.get_value(2,3), eps);
    
    EXPECT_NEAR(0.242188, r.midpoint_value(),  eps);
    EXPECT_NEAR(0.265625, r.trapezoid_value(), eps);
    EXPECT_NEAR(0.250000, r.simpson_value(),   eps);
}

TEST_F(RombergIntegratorTest, 
       ShouldReturnCorrectResultsWhenRefineOnQuarticFunction)
{
    Romberg_integrator r( &quartic_function, 0, 1);

    r.refine_bisection();
    r.refine_bisection();
    r.refine_expansion();
    r.refine_expansion();

    double eps = 1e-6;
    EXPECT_NEAR(0.062500, r.get_value(0,0), eps);
    EXPECT_NEAR(0.160156, r.get_value(1,0), eps);
    EXPECT_NEAR(0.189697, r.get_value(2,0), eps);
    EXPECT_NEAR(0.500000, r.get_value(0,1), eps);
    EXPECT_NEAR(0.281250, r.get_value(1,1), eps);
    EXPECT_NEAR(0.220703, r.get_value(2,1), eps);
    EXPECT_NEAR(0.208333, r.get_value(1,2), eps);
    EXPECT_NEAR(0.200521, r.get_value(2,2), eps);
    EXPECT_NEAR(0.200000, r.get_value(2,3), eps);
    
    EXPECT_NEAR(0.189697, r.midpoint_value(),  eps);
    EXPECT_NEAR(0.220703, r.trapezoid_value(), eps);
    EXPECT_NEAR(0.200521, r.simpson_value(),   eps);
}

TEST_F(RombergIntegratorTest,
       DirectAndRecursiveMethodsShouldReturnTheSameResultForQuadratic)
{
    Romberg_integrator r( &quadratic_function, 0, 1);
    
    r.refine_bisection();
    r.refine_bisection();
    r.refine_expansion();
    r.refine_expansion();

    double eps = 1e-9;
    EXPECT_NEAR(r.midpoint_direct_evaluate(4),  r.midpoint_value(),  eps);
    EXPECT_NEAR(r.trapezoid_direct_evaluate(4), r.trapezoid_value(), eps);
    EXPECT_NEAR(r.simpson_direct_evaluate(2),   r.simpson_value(),   eps);
}

TEST_F(RombergIntegratorTest,
       DirectAndRecursiveMethodsShouldReturnTheSameResultForCubic)
{
    Romberg_integrator r( &cubic_function, 0, 1);
    
    r.refine_bisection();
    r.refine_bisection();
    r.refine_expansion();
    r.refine_expansion();

    double eps = 1e-9;
    EXPECT_NEAR(r.midpoint_direct_evaluate(4),  r.midpoint_value(),  eps);
    EXPECT_NEAR(r.trapezoid_direct_evaluate(4), r.trapezoid_value(), eps);
    EXPECT_NEAR(r.simpson_direct_evaluate(2),   r.simpson_value(),   eps);
}

TEST_F(RombergIntegratorTest,
       DirectAndRecursiveMethodsShouldReturnTheSameResultForQuartic)
{
    Romberg_integrator r( &quartic_function, 0, 1);
    
    r.refine_bisection();
    r.refine_bisection();
    r.refine_expansion();
    r.refine_expansion();

    double eps = 1e-9;
    EXPECT_NEAR(r.midpoint_direct_evaluate(4),  r.midpoint_value(),  eps);
    EXPECT_NEAR(r.trapezoid_direct_evaluate(4), r.trapezoid_value(), eps);
    EXPECT_NEAR(r.simpson_direct_evaluate(2),   r.simpson_value(),   eps);
}

TEST_F(RombergIntegratorTest, 
       ShouldReturnCorrectResultsWhenRefineOnNormalFunction)
{
    Romberg_integrator r( &normal_function, 0, 2);

    r.refine_bisection();  //2
    r.refine_bisection();  //4
    r.refine_bisection();  //8
    r.refine_bisection();  //16
    r.refine_bisection();  //32
    r.refine_bisection();  //64
    r.refine_bisection();  //128
    r.refine_bisection();  //256
    r.refine_bisection();  //512
    r.refine_expansion();
    r.refine_expansion();

    double eps = 1e-8;

    EXPECT_NEAR(0.88278895, r.get_value(2,0), eps);
    EXPECT_NEAR(0.88226870, r.get_value(3,0), eps);
    EXPECT_NEAR(0.88212887, r.get_value(4,0), eps);
    EXPECT_NEAR(0.88209330, r.get_value(5,0), eps);
    EXPECT_NEAR(0.88208437, r.get_value(6,0), eps);
    EXPECT_NEAR(0.88208214, r.get_value(7,0), eps);
    EXPECT_NEAR(0.88208158, r.get_value(8,0), eps);
    EXPECT_NEAR(0.88208144, r.get_value(9,0), eps);
    
    EXPECT_NEAR(0.88061863, r.get_value(2,1), eps);
    EXPECT_NEAR(0.88170379, r.get_value(3,1), eps);
    EXPECT_NEAR(0.88198624, r.get_value(4,1), eps);
    EXPECT_NEAR(0.88205756, r.get_value(5,1), eps);
    EXPECT_NEAR(0.88207543, r.get_value(6,1), eps);
    EXPECT_NEAR(0.88207990, r.get_value(7,1), eps);
    EXPECT_NEAR(0.88208102, r.get_value(8,1), eps);
    EXPECT_NEAR(0.88208130, r.get_value(9,1), eps);
    
    EXPECT_NEAR(0.88206551, r.get_value(3,2), eps);
    EXPECT_NEAR(0.88208040, r.get_value(4,2), eps);
    EXPECT_NEAR(0.88208133, r.get_value(5,2), eps);
    EXPECT_NEAR(0.88208139, r.get_value(6,2), eps);
    
    EXPECT_NEAR(r.midpoint_direct_evaluate(512),  r.midpoint_value(),  eps);
    EXPECT_NEAR(r.trapezoid_direct_evaluate(512), r.trapezoid_value(), eps);
    EXPECT_NEAR(r.simpson_direct_evaluate(256),   r.simpson_value(),   eps);
}

TEST_F(RombergIntegratorTest, 
       MidpointRecursiveShouldReturnResultAccurateUpToTolerance)
{
    double exact_value = 0.88208139076242;
    double tol = 1e-6;
    
    Romberg_integrator r( &normal_function, 0, 2);
    double value = r.midpoint_recursive_evaluate(tol);

    EXPECT_NEAR(exact_value, value, tol);
    EXPECT_EQ(8, r.bisection_order());
}

TEST_F(RombergIntegratorTest, 
       TrapezoidRecursiveShouldReturnResultAccurateUpToTolerance)
{
    double exact_value = 0.88208139076242;
    double tol = 1e-7;
    
    Romberg_integrator r( &normal_function, 0, 2);
    double value = r.trapezoid_recursive_evaluate(tol);

    EXPECT_NEAR(exact_value, value, tol);
    EXPECT_EQ(10, r.bisection_order());
}

TEST_F(RombergIntegratorTest, 
       SimpsonRecursiveShouldReturnResultAccurateUpToTolerance)
{
    double exact_value = 0.88208139076242;
    double tol = 1e-8;
    
    Romberg_integrator r( &normal_function, 0, 2);
    double value = r.simpson_recursive_evaluate(tol);

    EXPECT_NEAR(exact_value, value, tol);
    EXPECT_EQ(7, r.bisection_order());

    r.print(std::cout, 9);
}

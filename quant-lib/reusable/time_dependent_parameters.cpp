#include <time_dependent_parameters.h>
#include <romberg_integrator.h>

Time_dependent_parameters::Time_dependent_parameters( double (*func)(double) )
                                                    : d_func(func)
{}

Parameters* Time_dependent_parameters::clone() const
{
    return new Time_dependent_parameters(*this);
}

double Time_dependent_parameters::integral(double time1, double time2) const
{
    Romberg_integrator r(d_func, time1, time2);
    return r.simpson_recursive_evaluate(1e-9);
}

double Time_dependent_parameters::square_integral(double time1, 
                                                  double time2) const
{
    // to be implemented
    return 0.0;
}

double Time_dependent_parameters::time_product_diff(double time1, 
                                                    double time2) const
{
    return time2*(*d_func)(time2)-time1*(*d_func)(time1);
}

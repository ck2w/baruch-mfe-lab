#include <romberg_integrator.h>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

Romberg_integrator:: Romberg_integrator() : d_integrand(0),
                                            d_lower_end(0),
                                            d_upper_end(0),
                                            d_bisection_order(0),
                                            d_expansion_order(1),
                                            d_interval_number(1)
{}

Romberg_integrator:: Romberg_integrator( double (*func)(double), 
                                         double lower_end,
                                         double upper_end )
                                       : d_integrand(func),
                                         d_lower_end(lower_end),
                                         d_upper_end(upper_end),
                                         d_bisection_order(0),
                                         d_expansion_order(1),
                                         d_interval_number(1)
{
    // initialize the midpoint approximations
    double midpoint = (d_lower_end + d_upper_end)/2;
    double midpoint_value = (*d_integrand)(midpoint)*interval_width();
    std::vector<double> midpoint_values(1, midpoint_value);
    d_romberg_approximations.push_back( midpoint_values );

    // initialize the romberg approximations with trapezoidal rule
    double trapezoid_value = ( (*d_integrand)(d_lower_end) 
                             + (*d_integrand)(d_upper_end) )
                             * interval_width() * 0.5;
    std::vector<double> trapezoid_values(1, trapezoid_value);
    d_romberg_approximations.push_back( trapezoid_values );
}

double Romberg_integrator::lower_end() const
{
    return d_lower_end;
}

double Romberg_integrator::upper_end() const
{
    return d_upper_end;
}

double Romberg_integrator::integrand(double x) const
{
    return (*d_integrand)(x);
}

int Romberg_integrator::bisection_order() const
{
    return d_bisection_order;
}

int Romberg_integrator::expansion_order() const
{
    return d_expansion_order;
}

int Romberg_integrator::interval_number() const
{
    return d_interval_number;
}

double Romberg_integrator::interval_width() const
{
    return (d_upper_end-d_lower_end)/d_interval_number;
}

bool Romberg_integrator::refine_bisection()
{
    if ( d_romberg_approximations.size() < 2 )
    {
        return false;
    }

    // new trapezoidal value
    double curr_midpoint = d_romberg_approximations[MIDPOINT_INDEX].back();
    double curr_trapezoid = d_romberg_approximations[TRAPEZOID_INDEX].back();
    double next_trapezoid = (curr_midpoint+curr_trapezoid)*0.5;
    d_romberg_approximations[TRAPEZOID_INDEX].push_back(next_trapezoid);
  
    // new Romberg values
    int coefficient = 1;
    for (int i=SIMPSON_INDEX; i<=d_expansion_order; i++)
    {
        std::vector<double> prev_expansion = d_romberg_approximations[i-1];
        std::vector<double>::iterator it = prev_expansion.end();
        it--; double prev_value = (*it);
        it--; double prev_prev_value = (*it);
        coefficient *= 4;
        double curr_value 
            = (coefficient*prev_value-prev_prev_value)/(coefficient-1);
        d_romberg_approximations[i].push_back(curr_value);
    }

    d_bisection_order++;
    d_interval_number *= 2;

    // new midpoint value
    double midpoint_value = midpoint_direct_evaluate(d_interval_number);
    d_romberg_approximations.front().push_back( midpoint_value );

    return true;
}

bool Romberg_integrator::refine_expansion()
{
    if ( d_romberg_approximations.size() < 2 ||
         d_expansion_order - d_bisection_order > 0 )
    {
        return false;
    }

    std::vector<double> curr_expansion = d_romberg_approximations.back();
    std::vector<double> next_expansion;

    int coefficient = 1;
    for (int i=1; i<=d_expansion_order; i++) { coefficient *= 4; }
    
    for ( std::vector<double>::const_iterator it = curr_expansion.begin()+1;
          it != curr_expansion.end();
          ++it )
    {
        double prev_value = (*it);
        double prev_prev_value = (*(it-1));
        double curr_value 
            = (coefficient*prev_value-prev_prev_value)/(coefficient-1);
        next_expansion.push_back(curr_value);
    }

    d_romberg_approximations.push_back(next_expansion);
    d_expansion_order++;

    return true;
}

double Romberg_integrator::midpoint_value() const
{
    if ( d_romberg_approximations.size() > MIDPOINT_INDEX &&
         d_romberg_approximations[MIDPOINT_INDEX].size() > 0 )
    {
        return d_romberg_approximations[MIDPOINT_INDEX].back();
    }
    else
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

double Romberg_integrator::trapezoid_value() const
{
    if ( d_romberg_approximations.size() > TRAPEZOID_INDEX &&
         d_romberg_approximations[TRAPEZOID_INDEX].size() > 0 )
    {
        return d_romberg_approximations[TRAPEZOID_INDEX].back();
    }
    else
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

double Romberg_integrator::simpson_value() const
{
    if ( d_romberg_approximations.size() > SIMPSON_INDEX &&
         d_romberg_approximations[SIMPSON_INDEX].size() > 0 )
    {
        return d_romberg_approximations[SIMPSON_INDEX].back();
    }
    else
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

double Romberg_integrator::best_estimate() const
{
    if ( d_romberg_approximations.size() > 0 &&
         d_romberg_approximations.back().size() > 0 )
    {
        return d_romberg_approximations.back().back();
    }
    else
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

double Romberg_integrator::best_estimate(int expansion_order) const
{
    if ( d_romberg_approximations.size() > expansion_order &&
         d_romberg_approximations.back().size() > 0 )
    {
        return d_romberg_approximations[expansion_order].back();
    }
    else
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

double Romberg_integrator::get_value(int bisection_order, 
                                     int expansion_order) const
{
    int bisection_index 
        = get_bisection_index(bisection_order, expansion_order);
    if ( d_romberg_approximations.size() > expansion_order &&
         d_romberg_approximations[expansion_order].size() > bisection_index )
    {
        return d_romberg_approximations[expansion_order][bisection_index];
    }
    else
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

double Romberg_integrator::recursive_evaluate(int expansion_order, double tol)
{
    // refine bisection until the given expansion order is supported
    while ( d_bisection_order < expansion_order )
    {
        if ( !refine_bisection() ) 
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    // refine expansion until the given order is reached 
    while ( d_expansion_order < expansion_order )
    {
        if ( !refine_expansion() )
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    // refine bisection until the given tolerance is satisfied
    double value = best_estimate(expansion_order);
    double error = 0;
    do {
        refine_bisection();
        double next_value = best_estimate(expansion_order);
        error = std::fabs(next_value-value);
        value = next_value;
    } while ( error > tol );

    return value;
}

double Romberg_integrator::midpoint_recursive_evaluate(double tol)
{
    return recursive_evaluate(MIDPOINT_INDEX, tol);
}

double Romberg_integrator::trapezoid_recursive_evaluate(double tol)
{
    return recursive_evaluate(TRAPEZOID_INDEX, tol);
}

double Romberg_integrator::simpson_recursive_evaluate(double tol)
{
    return recursive_evaluate(SIMPSON_INDEX, tol);
}

void Romberg_integrator::print(std::ostream & os, int num_digits) const
{
    int spacing = num_digits+5;

    os << std::setw(spacing) << "Midpoint" 
       << std::setw(spacing) << "Trapezoid" 
       << std::setw(spacing) << "Simpson" 
       << std::setw(spacing) << "Rombergs" 
       << std::endl;

    for (int i=0; i<=d_bisection_order; i++)
    {
        for (int j=0; j<=d_expansion_order; j++)
        {
            double value = get_value(i,j);
            if ( !std::isnan(value) ) 
            { 
                os << std::setw(spacing) 
                   << std::setprecision(num_digits)
                   << std::fixed 
                   << value;
            }
        }

        os << std::endl;
    }
}

int Romberg_integrator::get_bisection_index(int bisection_order, 
                                            int expansion_order) const
{
    if ( expansion_order == 0 )
    {
        return bisection_order;
    }
    else
    {
        return (bisection_order-expansion_order+1);
    }
}

double Romberg_integrator::midpoint_direct_evaluate(int partition_number) const
{
    double h = (d_upper_end-d_lower_end)/partition_number;
    double midpoint_value = 0.0;
    for (int i=0; i<partition_number; i++)
    {
        double abscissa = d_lower_end+(i+0.5)*h;
        midpoint_value += (*d_integrand)(abscissa);
    }

    midpoint_value *= h;
    return midpoint_value;
}

double Romberg_integrator::trapezoid_direct_evaluate(int partition_number) const
{
    double h = (d_upper_end-d_lower_end)/partition_number;
    double trapezoid_value = ( (*d_integrand)(d_lower_end) 
                             + (*d_integrand)(d_upper_end) ) / 2; 

    for (int i=1; i<partition_number; i++)
    {
        double abscissa = d_lower_end + i*h;
        trapezoid_value += (*d_integrand)(abscissa);
    }

    trapezoid_value *= h;
    return trapezoid_value;
}

double Romberg_integrator::simpson_direct_evaluate(int partition_number) const
{
    // For Simpson's rule, partition_number = interval_number/2
    double h = (d_upper_end-d_lower_end)/partition_number;
    double simpson_value = ( (*d_integrand)(d_lower_end) 
                           + (*d_integrand)(d_upper_end) ) / 6; 

    for (int i=1; i<partition_number; i++)
    {
        double abscissa = d_lower_end + i*h;
        simpson_value += (*d_integrand)(abscissa)/3;
    }

    for (int i=0; i<partition_number; i++)
    {
        double abscissa = d_lower_end+(i+0.5)*h;
        simpson_value += 2*(*d_integrand)(abscissa)/3;
    }

    simpson_value *= h;
    return simpson_value;
}

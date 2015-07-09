#include <root_finder.h>
#include <cmath>
#include <iostream>
#include <iomanip>

static void print_one_iteration( int iteration,
                                 double (*func)(double), 
                                 double x_new, 
                                 double x_old )
{
    std::cout << iteration << std::setw(6)
              << std::setw(15) << std::setprecision(9) << std::fixed
              << x_new
              << std::setw(15) << std::setprecision(9) << std::fixed
              << (*func)(x_new) 
              << std::setw(15) << std::setprecision(9) << std::fixed
              << x_new-x_old
              << std::endl;
}

double Root_finder::newton_find_root(double initial_guess,
                                     double tol,
                                     double eps,
                                     bool print)
{
    double x = initial_guess;

    if ( d_function && d_derivative )
    {
        double x_old = x-1;
        int iteration = 0;
        while ( std::fabs((*d_function)(x)) > tol || std::fabs(x-x_old) > eps )
        {
            if (print)
            {
                print_one_iteration(iteration, d_function, x, x_old);
            }

            x_old = x;
            x = x_old - d_function(x_old)/d_derivative(x_old); 
            iteration++;
        }
    }

    return x;
}

double Root_finder::secant_find_root(double first_guess,
                                     double second_guess,
                                     double tol,
                                     double eps,
                                     bool print)
{
    double x = second_guess;

    if ( d_function )
    {
        double x_old = first_guess;
        int iteration = 0;
        while ( std::fabs((*d_function)(x)) > tol || std::fabs(x-x_old) > eps )
        {
            if (print)
            {
                print_one_iteration(iteration, d_function, x, x_old);
            }

            double x_oldest = x_old;
            x_old = x;
            x = x_old - d_function(x_old) * 
                (x_old-x_oldest)/(d_function(x_old)-d_function(x_oldest));
        }
    }

    return x;
}

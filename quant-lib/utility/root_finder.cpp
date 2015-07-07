#include <root_finder.h>
#include <cmath>
#include <iostream>
#include <iomanip>

double Root_finder::newton_find_root(double initial_guess,
                                     double tol,
                                     double eps,
                                     bool print)
{
    double x = initial_guess;
    double x_old = x-1;
    int iteration = 0;
    while ( std::fabs((*d_function)(x)) > tol || std::fabs(x-x_old) > eps )
    {
        if (print)
        {
            std::cout << iteration << std::setw(6)
                      << std::setw(15) << std::setprecision(9) << std::fixed
                      << x
                      << std::setw(15) << std::setprecision(9) << std::fixed
                      << (*d_function)(x) 
                      << std::setw(15) << std::setprecision(9) << std::fixed
                      << x-x_old
                      << std::endl;
        }

        x_old = x;
        x = x_old - d_function(x_old)/d_derivative(x_old); 
        iteration++;
    }

    return x;
}

double Root_finder::secant_find_root(double first_guess,
                                     double second_guess,
                                     double tol,
                                     double eps)
{
    // to be implemented
    return 0.0;
}

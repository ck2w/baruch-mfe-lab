#ifndef ROOT_FINDER_H
#define ROOT_FINDER_H

class Root_finder
{
    public:

        Root_finder( double (*function)(double), 
                     double (*derivative)(double)=0 )
                   : d_function(function), 
                     d_derivative(derivative)
        {}

        double newton_find_root(double initial_guess,
                                double tol,  // largest admissible |func(x)|
                                double eps,  // largest admissible |x'-x|
                                bool print=false);

        double secant_find_root(double first_guess,
                                double second_guess,
                                double tol,  // largest admissible |func(x)|
                                double eps,  // largest admissible |x'-x|
                                bool print);                        

    private:

        Root_finder();

        double (*d_function)(double);
        double (*d_derivative)(double);
};

#endif /* ROOT_FINDER_H */

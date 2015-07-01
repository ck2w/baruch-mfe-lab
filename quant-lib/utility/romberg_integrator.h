#ifndef ROMBERG_INTEGRATOR_H
#define ROMBERG_INTEGRATOR_H 
#include <vector>
#include <ostream>

class Romberg_integrator
{
    public:

        Romberg_integrator();

        Romberg_integrator( double (*func)(double), 
                            double lower_end,
                            double upper_end );

        double lower_end() const;
        double upper_end() const;
        double integrand(double x) const;

        int bisection_order() const;
        int expansion_order() const;

        int interval_number() const;
        double interval_width() const;
        
        bool refine_bisection();
        bool refine_expansion();

        double midpoint_value() const;
        double trapezoid_value() const;
        double simpson_value() const;
        double best_estimate() const;
        double best_estimate(int expansion_order) const;
        double get_value(int bisection_order, int expansion_order) const;

        double recursive_evaluate(int expansion_order, double tol);
        double midpoint_recursive_evaluate(double  tol);
        double trapezoid_recursive_evaluate(double tol);
        double simpson_recursive_evaluate(double   tol);
        
        void print(std::ostream & os, int num_digits) const;
        
        double midpoint_direct_evaluate(int  partition_number) const;
        double trapezoid_direct_evaluate(int partition_number) const;
        double simpson_direct_evaluate(int   partition_number) const;

    private:

        enum { MIDPOINT_INDEX  = 0 };
        enum { TRAPEZOID_INDEX = 1 };
        enum { SIMPSON_INDEX   = 2 };

        double (*d_integrand)(double);
        double d_lower_end;
        double d_upper_end;

        int d_bisection_order;
        int d_expansion_order;
        int d_interval_number;
        
        int get_bisection_index(int bisection_order, int expansion_order) const;

        /***********************************************************
         *                                                         *
         * The hierarchy of Romberg approximations                 *
         *                                                         *
         * A(0,0)=M(1)   A(1,0)=T(1)                               *
         * A(0,1)=M(2)   A(1,1)=T(2)  A(2,1)=S(1)                  *
         * A(0,2)=M(4)   A(1,2)=T(4)  A(2,2)=S(2)  A(3,2)          *
         * A(0,3)=M(8)   A(1,3)=T(8)  A(2,3)=S(4)  A(3,3)  A(4,3)  *
         *    .             .            .            .       .    *
         *    .             .            .            .       .    *
         *    .             .            .            .       .    *
         *                                                         *
         ***********************************************************/
        std::vector< std::vector<double> > d_romberg_approximations;
};

#endif /* ROMBERG_INTEGRATOR_H */

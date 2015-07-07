#ifndef BOND_H
#define BOND_H
#include <vector>
#include <ostream>
#include <constant_parameters.h>
#include <time_dependent_parameters.h>

typedef std::vector< std::pair<double, double> > Cashflow;

class Bond
{
    public:

        Bond(double expiry, double par_value);
        // coupon rate: a percentage of the par value per year
        // coupon frequency: the number of coupon payments per year
        // coupon period: the time interval between coupon payments
        Bond(double expiry, double par_value, 
             double coupon_rate, int coupon_frequency);
        Bond(double expiry, double par_value, 
             double coupon_rate, double coupon_period);

        double get_expiry() const;
        double get_par_value() const;

        int get_number_of_coupons() const;
        bool add_coupon(double time, double value);
        void sort_coupons();

        void print_cashflow(std::ostream & os, int num_digits) const;

        double evaluate_from_zero_rate
            ( const Time_dependent_parameters & zero_rate ) const;

        double evaluate_from_instantaneous_rate
            ( const Time_dependent_parameters & instantaneous_rate ) const;

        // evaluate bond price, duration, and convexity from yield
        double evaluate_from_yield( const Constant_parameters & yield,
                                    double * duration, 
                                    double * convexity ) const;
        
        // called by constructor to initialize coupons
        void init_coupon_list(double expiry, 
                              double coupon_payment, 
                              double coupon_period);

    private:

        Bond() {}

        double d_expiry;
        double d_par_value;
        Cashflow d_coupons;
};

#endif /* BOND_H */

#include <bond.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

Bond::Bond(double expiry, double par_value) : d_expiry(expiry), 
                                              d_par_value(par_value)
{}

void Bond::init_coupon_list(double expiry, 
                            double coupon_payment,
                            double coupon_period)
{
    for ( double residual = expiry; residual > 0; residual -= coupon_period)
    {
        add_coupon(residual, coupon_payment);
    }
}

Bond::Bond( double expiry, double par_value, 
            double coupon_rate, int coupon_frequency )
          : d_expiry(expiry), 
            d_par_value(par_value)
{
    double coupon_payment = par_value * coupon_rate / coupon_frequency;
    init_coupon_list(expiry, coupon_payment, 1.0/coupon_frequency);
}

Bond::Bond( double expiry, double par_value, 
            double coupon_rate, double coupon_period )
          : d_expiry(expiry), 
            d_par_value(par_value)
{
    double coupon_payment = par_value * coupon_rate * coupon_period;
    init_coupon_list(expiry, coupon_payment, coupon_period);
}

double Bond::get_expiry() const
{
    return d_expiry;
}

double Bond::get_par_value() const
{
    return d_par_value;
}

bool Bond::add_coupon(double time, double value)
{
    if ( time <= d_expiry )
    {
        bool need_coupon_sort = get_number_of_coupons() > 0 &&
                                time < d_coupons.back().first;
        d_coupons.push_back( std::make_pair(time, value) );
        if ( need_coupon_sort )
        {
            sort_coupons();
        }

        return true;
    }
    else
    {
        return false;
    }
}

int Bond::get_number_of_coupons() const
{
    return d_coupons.size();
}

void Bond::sort_coupons()
{
    std::sort(d_coupons.begin(), d_coupons.end());
}

void Bond::print_cashflow(std::ostream & os, int num_digits) const
{
    int spacing = num_digits+5;
    os << std::setw(spacing) << "Time" << std::setw(spacing) << "Cashflow"
       << std::endl;

    double last_coupon_time = d_coupons.back().first;
    bool is_last_coupon_at_expiry 
        = (std::fabs(last_coupon_time-d_expiry) < 1e-6);

    int number_of_coupons_before_expiry = get_number_of_coupons();
    if (is_last_coupon_at_expiry) {number_of_coupons_before_expiry--;} 

    for (int i=0; i<number_of_coupons_before_expiry; i++)
    {
        os << std::setw(spacing) << std::setprecision(num_digits) << std::fixed
           << d_coupons[i].first
           << std::setw(spacing) << std::setprecision(num_digits) << std::fixed
           << d_coupons[i].second
           << std::endl;
    }

    double last_payment = d_par_value;
    if (is_last_coupon_at_expiry) { last_payment += d_coupons.back().second; }

    os << std::setw(spacing) << std::setprecision(num_digits) << std::fixed
       << d_expiry 
       << std::setw(spacing) << std::setprecision(num_digits) << std::fixed
       << last_payment
       << std::endl;
}

double Bond::evaluate_from_zero_rate
    (const Time_dependent_parameters & zero_rate) const
{
    double bond_price = 0.0;
    int number_of_coupons = get_number_of_coupons();

    // coupons
    for (int i=0; i<number_of_coupons; i++)
    {
        double discount_factor 
            = std::exp( -zero_rate.time_product_diff(0, d_coupons[i].first) );
        bond_price += discount_factor * d_coupons[i].second;
    }

    // par value
    double discount_factor
        = std::exp( -zero_rate.time_product_diff(0, d_expiry) );
    bond_price += discount_factor * d_par_value;

    return bond_price;
}

double Bond::evaluate_from_instantaneous_rate
    (const Time_dependent_parameters & instantaneous_rate) const
{
    double bond_price = 0.0;
    int number_of_coupons = get_number_of_coupons();

    // coupons
    for (int i=0; i<number_of_coupons; i++)
    {
        double discount_factor 
            = std::exp( -instantaneous_rate.integral(0, d_coupons[i].first) );
        bond_price += discount_factor * d_coupons[i].second;
    }
    
    // par value
    double discount_factor
        = std::exp( -instantaneous_rate.integral(0, d_expiry) );
    bond_price += discount_factor * d_par_value;

    return bond_price;
}

double Bond::evaluate_from_yield( const Constant_parameters & yield,
                                  double * duration, 
                                  double * convexity ) const
{
    double bond_price = 0.0;
    (*duration) = 0.0;
    (*convexity) = 0.0;
    int number_of_coupons = get_number_of_coupons();

    // coupons
    for (int i=0; i<number_of_coupons; i++)
    {
        double discount_factor 
            = std::exp( -yield.integral(0, d_coupons[i].first) );
        bond_price += discount_factor * d_coupons[i].second;
        (*duration) += discount_factor * d_coupons[i].second * 
                       d_coupons[i].first;
        (*convexity) += discount_factor * d_coupons[i].second * 
                        d_coupons[i].first * d_coupons[i].first;
    }
    
    // par value
    double discount_factor
        = std::exp( -yield.integral(0, d_expiry) );
    bond_price += discount_factor * d_par_value;
    (*duration) += discount_factor * d_par_value * d_expiry;
    (*convexity) += discount_factor * d_par_value * d_expiry * d_expiry;

    (*duration) /= bond_price;
    (*convexity) /= bond_price;

    return bond_price;
}

#include <european_option.h>
#include <european_payoff.h>
#include <cumulative_normal.h>
#include <cmath>
#include <iostream>

double European_call::get_strike() const
{
    double spot = 0;
    double spot_incremental = 10;
    while ( d_vanilla_option.get_payoff(spot) <=0 )
    {
        spot += spot_incremental; 
    }

    return (spot-d_vanilla_option.get_payoff(spot));
}

double European_put::get_strike() const
{
    return d_vanilla_option.get_payoff(0.0);
}

double European_call::Black_Scholes_evaluate( double spot,
                                              const Parameters & vol,
                                              const Parameters & div,
                                              const Parameters & rate ) const
{
    double strike = get_strike();
    double expiry = d_vanilla_option.get_expiry();
    double vol_square = vol.square_integral(0, expiry);
    double vol_expiry_root = std::sqrt(vol_square);

    double d1 = (std::log(spot/strike) 
              + ( rate.integral(0, expiry)
                 -div.integral(0, expiry)
                 +0.5*vol_square )) / vol_expiry_root;
    double d2 = d1 - vol_expiry_root;

    return spot*std::exp(-div.integral(0, expiry)) * 
           Cumulative_normal::approximate_evaluate(d1) 
         - strike*std::exp(-rate.integral(0, expiry)) * 
           Cumulative_normal::approximate_evaluate(d2);
}

double European_put::Black_Scholes_evaluate( double spot,
                                             const Parameters & vol,
                                             const Parameters & div,
                                             const Parameters & rate ) const
{
    double strike = d_vanilla_option.get_payoff(0.0);
    double expiry = d_vanilla_option.get_expiry();
    double vol_square = vol.square_integral(0, expiry);
    double vol_expiry_root = std::sqrt(vol_square);

    double d1 = (std::log(spot/strike) 
              + ( rate.integral(0, expiry)
                 -div.integral(0, expiry)
                 +0.5*vol_square )) / vol_expiry_root;
    double d2 = d1 - vol_expiry_root;

    return strike*std::exp(-rate.integral(0, expiry)) *
           Cumulative_normal::approximate_evaluate(-d2)
         - spot*std::exp(-div.integral(0, expiry)) *
           Cumulative_normal::approximate_evaluate(-d1);
}

double European_call::Monte_Carlo_evaluate( double spot,
                                            const Parameters & vol,
                                            const Parameters & div,
                                            const Parameters & rate,
                                            long number_of_paths ) const
{
    return d_vanilla_option.Monte_Carlo_evaluate
        (spot, vol, div, rate, number_of_paths);
}

double European_put::Monte_Carlo_evaluate( double spot,
                                           const Parameters & vol,
                                           const Parameters & div,
                                           const Parameters & rate,
                                           long number_of_paths ) const
{
    return d_vanilla_option.Monte_Carlo_evaluate
        (spot, vol, div, rate, number_of_paths);
}


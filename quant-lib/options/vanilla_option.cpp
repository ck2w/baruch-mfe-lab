#include <vanilla_option.h>

double Vanilla_option::get_payoff(double spot) const
{
    return d_payoff(spot);
}

double Vanilla_option::get_expiry() const
{
    return d_expiry;
}
        
double Vanilla_option::Monte_Carlo_evaluate( double spot,
                                             const Parameters & vol,
                                             const Parameters & div,
                                             const Parameters & rate,
                                             long number_of_paths ) const
{
    return 0;
}


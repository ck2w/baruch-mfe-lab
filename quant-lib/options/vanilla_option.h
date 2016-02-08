#ifndef VANILLA_OPTION_H
#define VANILLA_OPTION_H
#include <option.h>
#include <payoff.h>
#include <parameters.h>

class Vanilla_option : public Option
{
    public:

        Vanilla_option() {}

        Vanilla_option( const Payoff_bridge & payoff, double expiry )
                      : d_payoff(payoff), d_expiry(expiry)
        {}

        double Monte_Carlo_evaluate( double spot,
                                     const Parameters & vol,
                                     const Parameters & div,
                                     const Parameters & rate,
                                     long number_of_paths ) const;

        double get_payoff(double spot) const;
        double get_expiry() const;

    private:

        Payoff_bridge d_payoff;
        double d_expiry;
};

#endif /* VANILLA_OPTION_H */

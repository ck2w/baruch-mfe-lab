#ifndef EUROPEAN_OPTION_H
#define EUROPEAN_OPTION_H
#include <option.h>
#include <vanilla_option.h>
#include <parameters.h>
#include <european_payoff.h>

class European_call : public Option
{
    public:

        European_call() {}

        European_call( double strike, double expiry )
                     : d_vanilla_option(European_call_payoff(strike), expiry)
        {}

        double get_strike() const;
        
        double Monte_Carlo_evaluate( double spot,
                                     const Parameters & vol,
                                     const Parameters & div,
                                     const Parameters & rate,
                                     long number_of_paths ) const;

        double Black_Scholes_evaluate( double spot,
                                       const Parameters & vol,
                                       const Parameters & div,
                                       const Parameters & rate ) const;
    
    private:

        Vanilla_option d_vanilla_option;
};

class European_put : public Option 
{
    public:

        European_put() {}

        European_put( double strike, double expiry )
                    : d_vanilla_option(European_put_payoff(strike), expiry)
        {}
        
        double get_strike() const;
        
        double Monte_Carlo_evaluate( double spot,
                                     const Parameters & vol,
                                     const Parameters & div,
                                     const Parameters & rate,
                                     long number_of_paths ) const;

        double Black_Scholes_evaluate( double spot,
                                       const Parameters & vol,
                                       const Parameters & div,
                                       const Parameters & rate ) const;
    
    private:

        Vanilla_option d_vanilla_option;
};

#endif /* EUROPEAN_OPTION_H */

#ifndef MC_H
#define MC_H 
#include <payoff.h>
#include <option_value.h>
#include <tuple>
#include <vector>

/******************************/
/*   Monte Carlo Simulation   */
/******************************/
class MonteCarlo 
{
    public:

        MonteCarlo( const Payoff & payoff,
                    double expiry,
                    double spot,
                    double rate,
                    double div,
                    double vol,
                    double (*ran)(int*, int*) );
        
        MonteCarlo( const Payoff & payoff,
                    double expiry,
                    double spot1, double spot2,
                    double rate,
                    double div1, double div2,
                    double vol1, double vol2,
                    double rho,
                    double (*ran)(int*, int*) );

        OptionValue BlackScholesValue() const;
        OptionValue evaluate(int N, 
                             bool controlVariate=false, 
                             bool useAntithetic=false,
                             bool momentMatching=false);
        OptionValue evaluateControlVariate(int N);
        OptionValue evaluateUseAntithetic(int N);
        OptionValue evaluateMomentMatching(int N);
        OptionValue evaluateControlVariateMomentMatching(int N);

        OptionValue evaluateBasket(int N);
        OptionValue evaluatePathDependentBasket(int N, int M);

    protected:
        
        int d_numOfUnderlyings;
        const Payoff* d_payoff;
        double d_expiry;
        double d_spot;
        double d_rate;
        double d_div;
        double d_vol;
        // if there are two underlyings
        double d_spot2;
        double d_div2;
        double d_vol2;
        double d_rho;
        // random number generator
        double (*d_ran)(int*, int*);
         // applicable when Black-Scholes
        double d_strike;
        OptionValue d_BlackScholes;

        double GetMeanEstimate(int N);
        double GetBetaEstimate(int N, bool momentMatching);
};

#endif /* MC_H */

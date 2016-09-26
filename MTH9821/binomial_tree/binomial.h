#ifndef BINOMIAL_TREE_H
#define BINOMIAL_TREE_H 
#include <payoff.h>
#include <option_value.h>
#include <tuple>
#include <vector>

/****************************/
/*   CRR Parameterization   */
/****************************/
class BinomialTree
{
    public:

        BinomialTree( const Payoff & payoff,
                      double expiry,
                      double spot,
                      double rate,
                      double div,
                      double vol );

        void print_config(int N);
        OptionValue BlackScholesValue() const;
        OptionValue evaluate(int N, 
                             bool isAmerican=false, 
                             bool useBS=false,
                             bool varReduction=false);
        OptionValue evaluateBBS(int N, bool isAmerican=false);
        // Var Reduction applicable to American options only
        OptionValue evaluateVarReduction(int N, bool useBS=false);
        OptionValue evaluateVarReductionBBS(int N);

    private:
        
        const Payoff* d_payoff;
        double d_expiry;
        double d_spot;
        double d_rate;
        double d_div;
        double d_vol;
         // applicable when Black-Scholes
        double d_strike;
        OptionValue d_BlackScholes;
};

#endif /* BINOMIAL_TREE_H */

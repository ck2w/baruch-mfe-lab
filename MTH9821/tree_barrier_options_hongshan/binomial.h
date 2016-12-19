#ifndef BINOMIAL_TREE_H
#define BINOMIAL_TREE_H 
#include "payoff.h"
#include "option_value.h"
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


		//divdend proportional to stock price, cannnot proceed when steps get large, memory inefficient
		OptionValue evaluate_prop_div(int N, double div, double freq, bool isAmerican = false);
		//modified version with less memory allocate
		OptionValue evaluate_prop_div2(int N, double div, double freq, bool isAmerican = false);
		
		//fix dividend
		OptionValue evaluate_fix_div(int N, double div, double freq,bool isAmerican = false);

		double spot_calculator_fix_div(double effective_spot,int k, int i, int N, double u, double div, double freq);

		//specific for problem(iii) of pricing discrete div paying option with binomial tree in hw10
		OptionValue evaluate_fix_prop_div_mix(int N, bool isAmerican = false);
		double spot_calculator_fix_prop_div_mix(double effective_spot, int k, int i, int N, double u);

		//right now down and out only:
		void barrierOption_exact(double strike,
						   double expiry,
						   double spot,
						   double rate,
						   double div,
						   double vol,
						   double barrier);
		//binomial tree for down and out call:
		OptionValue barrierOption_binomial_tree(int N, double barrier);
		

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

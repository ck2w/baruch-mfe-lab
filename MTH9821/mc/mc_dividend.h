#ifndef MC_DIVIDEND_H
#define MC_DIVIDEND_H 
#include <mc.h>
#include <payoff.h>
#include <option_value.h>
#include <dividend.h>
#include <tuple>
#include <vector>
#include <utility>


/******************************/
/*   Monte Carlo Simulation   */
/*   with Discrete Dividends  */
/******************************/
class MonteCarloDiscreteDividends : public MonteCarlo
{
    public:

        MonteCarloDiscreteDividends
            ( const Payoff & payoff,
              double expiry,
              double spot,
              double rate,
              double div,
              const std::vector<DiscreteDividend> & discreteDividends,
              double vol,
              double (*ran)(int*, int*) );

        OptionValue evaluate(int N, bool controlVariate=false);

    private:

        std::vector<DiscreteDividend> d_discreteDividends;
        std::pair<OptionValue,OptionValue> GetMeanEstimates(int N);
        std::pair<double,double> GetBetaEstimates(int N);
};

#endif /* MC_DIVIDEND_H */

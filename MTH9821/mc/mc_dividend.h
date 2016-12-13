#ifndef MC_DIVIDEND_H
#define MC_DIVIDEND_H 
#include <mc.h>
#include <payoff.h>
#include <option_value.h>
#include <dividend.h>
#include <tuple>
#include <vector>


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

        OptionValue evaluate(int N,
                             bool controlVariate=false, 
                             bool useAntithetic=false,
                             bool momentMatching=false);

    private:

        std::vector<DiscreteDividend> d_discreteDividends;

        //double d_t1; //time
        //double d_q1; //dividend
        //bool   d_p1; //proportional
};

#endif /* MC_DIVIDEND_H */

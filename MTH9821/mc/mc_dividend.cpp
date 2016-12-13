#include <mc_dividend.h>
#include <dividend.h>
#include <ran.h>
#include <payoff.h>
#include <option_value.h>
#include <tuple>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>

MonteCarloDiscreteDividends::MonteCarloDiscreteDividends
                   ( const Payoff & payoff, 
                     double expiry,
                     double spot,
                     double rate,
                     double div,
                     const std::vector<DiscreteDividend> & discreteDividends,
                     double vol,
                     double (*ran)(int*, int*) )
                   : MonteCarlo( payoff,
                                 expiry,
                                 spot,
                                 rate,
                                 div,
                                 vol,
                                 ran ),
                     d_discreteDividends(discreteDividends)
{
}

OptionValue MonteCarloDiscreteDividends::evaluate(int N,
                                                  bool controlVariate, 
                                                  bool useAntithetic,
                                                  bool momentMatching)
{
    int idum=1;
    int count=0;
    int sampleCount=0;
    double V = 0.0;
    double Var = 0.0;
    double Delta = 0.0;
    
    for (int k=0; k<N; k++) {
        double t=0;
        double s=d_spot;
        double sDivFactor=1;
        int numDividends = d_discreteDividends.size();

        // all discrete dividends
        for (int i=0; i<numDividends; i++) {
            DiscreteDividend dDiv = d_discreteDividends[i];
            double tDiv = dDiv.t;
            double vDiv = dDiv.v;
            double isFixed = dDiv.isFixed;

            double z = (*d_ran)(&idum, &count);
            double trend = (d_rate-d_div-0.5*d_vol*d_vol)*(tDiv-t);
            double random = d_vol*std::sqrt(tDiv-t)*z;
            s *= std::exp(trend+random);
            sDivFactor *= std::exp(trend+random);

            if (isFixed) { 
                s -= vDiv; 
            }
            else { 
                s *= (1-vDiv); 
                sDivFactor *= (1-vDiv);
            }

            t = tDiv;
        }

        // final payoff
        double z = (*d_ran)(&idum, &count);
        double trend = (d_rate-d_div-0.5*d_vol*d_vol)*(d_expiry-t);
        double random = d_vol*std::sqrt(d_expiry-t)*z;
        s *= std::exp(trend+random);
            
        double v = std::exp(-d_rate*d_expiry)*(*d_payoff)(s);
        V += v;
        Var += v*v;

        // Delta calculation
        double delta = 0.0;
        if ( (*d_payoff)(0) > 0 && d_strike > s ) {
            if ( d_strike > s ) {
                delta = -std::exp(-d_rate*d_expiry)*sDivFactor;
            }
        }
        else {
            if ( d_strike < s ) {
                delta = std::exp(-d_rate*d_expiry)*sDivFactor;
            }
        }

        Delta += delta;

        sampleCount++;
    }
    
    V /= sampleCount;
    Var /= (sampleCount-1);
    Var -= N*V*V/(sampleCount-1);
    Delta /= sampleCount;

    OptionValue optionValue;
    optionValue.price = V;
    optionValue.priceStd = std::sqrt(Var);
    optionValue.priceErr = optionValue.priceStd/std::sqrt(sampleCount);
    optionValue.delta = Delta;
    return optionValue;
}

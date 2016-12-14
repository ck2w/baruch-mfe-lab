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
#include <utility>

std::pair<OptionValue,OptionValue> 
    MonteCarloDiscreteDividends::GetMeanEstimates(int N)
{
    double V=0.0;
    double VNoDiv=0.0;
    double Delta=0.0;
    double DeltaNoDiv=0.0;

    int idum=1;
    int count=0;
    int sampleCount=0;
    for (int k=0; k<N; k++) {
        double t=0;
        double s=d_spot;
        double sNoDiv=d_spot;
        double divFactor=1;
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
            sNoDiv *= std::exp(trend+random);

            if (isFixed) { 
                s -= vDiv; 
            }
            else { 
                s *= (1-vDiv); 
            }

            t = tDiv;
        }

        // final payoff
        double z = (*d_ran)(&idum, &count);
        double trend = (d_rate-d_div-0.5*d_vol*d_vol)*(d_expiry-t);
        double random = d_vol*std::sqrt(d_expiry-t)*z;
        s *= std::exp(trend+random);
        sNoDiv *= std::exp(trend+random);
           
        // payoff with dividend 
        double v = std::exp(-d_rate*d_expiry)*(*d_payoff)(s);
        // payoff without dividend
        double vNoDiv = std::exp(-d_rate*d_expiry)*(*d_payoff)(sNoDiv);
        
        V += v;
        VNoDiv += vNoDiv;
        
        // Delta calculation
        double delta = 0.0;
        double deltaNoDiv = 0.0;
        if ( (*d_payoff)(0) > 0 ) {
            // delta with dividend
            if ( d_strike > s ) {
                delta = -std::exp(-d_rate*d_expiry)*sNoDiv*divFactor/d_spot;
            }
            // delta without dividend
            if ( d_strike > sNoDiv) {
                deltaNoDiv = -std::exp(-d_rate*d_expiry)*sNoDiv/d_spot;
            }
        }
        else {
            // delta with dividend
            if ( d_strike < s ) {
                delta = std::exp(-d_rate*d_expiry)*sNoDiv*divFactor/d_spot;
            }
            // delta without dividend
            if ( d_strike < sNoDiv) {
                deltaNoDiv = std::exp(-d_rate*d_expiry)*sNoDiv/d_spot;
            }
        }

        Delta += delta;
        DeltaNoDiv += deltaNoDiv;

        sampleCount++;
    }

    V /= sampleCount;
    VNoDiv /= sampleCount;
    Delta /= sampleCount;
    DeltaNoDiv /= sampleCount;

    OptionValue optionWithDividend;
    OptionValue optionNoDividend;
    optionWithDividend.price = V;
    optionWithDividend.delta = Delta;
    optionNoDividend.price = VNoDiv;
    optionNoDividend.delta = DeltaNoDiv;

    return std::make_pair(optionWithDividend, optionNoDividend);
}

std::pair<double,double> MonteCarloDiscreteDividends::GetBetaEstimates(int N)
{
    // get the sample mean for moment matching
    std::pair<OptionValue,OptionValue> mean = GetMeanEstimates(N);
    double vMean = mean.first.price;
    double deltaMean = mean.first.delta;
    double vNoDivMean = mean.second.price;
    double deltaNoDivMean = mean.second.delta;

    double vCovariance=0.0;
    double vVariance=0.0;
    double deltaCovariance=0.0;
    double deltaVariance=0.0;

    int idum=1;
    int count=0;
    for (int k=0; k<N; k++) {
        double t=0;
        double s=d_spot;
        double sNoDiv=d_spot;
        double divFactor=1;
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
            sNoDiv *= std::exp(trend+random);

            if (isFixed) { 
                s -= vDiv; 
            }
            else { 
                s *= (1-vDiv); 
            }

            t = tDiv;
        }

        // final payoff
        double z = (*d_ran)(&idum, &count);
        double trend = (d_rate-d_div-0.5*d_vol*d_vol)*(d_expiry-t);
        double random = d_vol*std::sqrt(d_expiry-t)*z;
        s *= std::exp(trend+random);
        sNoDiv *= std::exp(trend+random);
           
        // payoff with dividend 
        double v = std::exp(-d_rate*d_expiry)*(*d_payoff)(s);
        // payoff without dividend
        double vNoDiv = std::exp(-d_rate*d_expiry)*(*d_payoff)(sNoDiv);
        
        // Delta calculation
        double delta = 0.0;
        double deltaNoDiv = 0.0;
        if ( (*d_payoff)(0) > 0 ) {
            // delta with dividend
            if ( d_strike > s ) {
                delta = -std::exp(-d_rate*d_expiry)*sNoDiv*divFactor/d_spot;
            }
            // delta without dividend
            if ( d_strike > sNoDiv) {
                deltaNoDiv = -std::exp(-d_rate*d_expiry)*sNoDiv/d_spot;
            }
        }
        else {
            // delta with dividend
            if ( d_strike < s ) {
                delta = std::exp(-d_rate*d_expiry)*sNoDiv*divFactor/d_spot;
            }
            // delta without dividend
            if ( d_strike < sNoDiv) {
                deltaNoDiv = std::exp(-d_rate*d_expiry)*sNoDiv/d_spot;
            }
        }

        vCovariance += (v-vMean)*(vNoDiv-vNoDivMean);
        vVariance += (vNoDiv-vNoDivMean)*(vNoDiv-vNoDivMean);

        deltaCovariance 
            += (delta-deltaMean)*(deltaNoDiv-deltaNoDivMean);
        deltaVariance 
            += (deltaNoDiv-deltaNoDivMean)*(deltaNoDiv-deltaNoDivMean);
    }

    double betaPrice = vCovariance/vVariance;
    double betaDelta = deltaCovariance/deltaVariance;

    return std::make_pair(betaPrice, betaDelta);
}

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
    // Calculate the Black-Scholes value
    // in case of a vanilla call or put.
    d_strike = payoff.strike();
    std::tuple<OptionValue,OptionValue> res
        = BlackScholes(d_expiry,d_strike,d_spot,d_rate,d_div,d_vol);
    if ((*d_payoff)(0) > 0) { d_BlackScholes = std::get<1>(res);}
    else { d_BlackScholes = std::get<0>(res);}
}

OptionValue MonteCarloDiscreteDividends::evaluate(int N, bool controlVariate)
{
    double betaPrice=0.0;
    double betaDelta=0.0;
    if (controlVariate) { 
        std::pair<double,double> beta = GetBetaEstimates(N);
        betaPrice = beta.first;
        betaDelta = beta.second;
    }

    int idum=1;
    int count=0;
    int sampleCount=0;
    double V = 0.0;
    double Var = 0.0;
    double Delta = 0.0;
    
    for (int k=0; k<N; k++) {
        double t=0;
        double s=d_spot;
        double sNoDiv=d_spot;
        double divFactor=1;
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
            sNoDiv *= std::exp(trend+random);

            if (isFixed) { 
                s -= vDiv; 
            }
            else { 
                s *= (1-vDiv); 
                divFactor *= (1-vDiv);
            }

            t = tDiv;
        }

        // final payoff
        double z = (*d_ran)(&idum, &count);
        double trend = (d_rate-d_div-0.5*d_vol*d_vol)*(d_expiry-t);
        double random = d_vol*std::sqrt(d_expiry-t)*z;
        s *= std::exp(trend+random);
        sNoDiv *= std::exp(trend+random);
        double v = std::exp(-d_rate*d_expiry)*(*d_payoff)(s);

        if (controlVariate) { 
            // Use the same option on a non-dividend 
            // paying asset as a control variate.
            double vNoDiv = std::exp(-d_rate*d_expiry)*(*d_payoff)(sNoDiv);
            double vBS = d_BlackScholes.price;
            v -= betaPrice*(vNoDiv-vBS);
        }
        
        V += v;
        Var += v*v;

        // Delta calculation
        double delta = 0.0;
        if ( (*d_payoff)(0) > 0 ) {
            if ( d_strike > s ) {
                delta = -std::exp(-d_rate*d_expiry)*sNoDiv*divFactor/d_spot;
            }
        }
        else {
            if ( d_strike < s ) {
                delta = std::exp(-d_rate*d_expiry)*sNoDiv*divFactor/d_spot;
            }
        }

        if (controlVariate) {
            double deltaNoDiv = 0.0;
            double deltaBS = d_BlackScholes.delta;
            if ( (*d_payoff)(0) > 0 ) {
                if ( d_strike > sNoDiv ) {
                    deltaNoDiv = -std::exp(-d_rate*d_expiry)*sNoDiv/d_spot;
                }
            }
            else {
                if ( d_strike < sNoDiv ) {
                    deltaNoDiv = std::exp(-d_rate*d_expiry)*sNoDiv/d_spot;
                }
            }

            delta -= betaDelta*(deltaNoDiv-deltaBS);
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

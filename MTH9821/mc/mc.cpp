#include <mc.h>
#include <ran.h>
#include <payoff.h>
#include <option_value.h>
#include <black_scholes.h>
#include <tuple>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>

double MonteCarlo::GetMeanEstimate(int N)
{
    double S=0.0;
    int idum=1;
    for (int k=0; k<N; k++) {
        double z = (*d_ran)(&idum);
        double trend = (d_rate-d_div-0.5*d_vol*d_vol)*d_expiry;
        double random = d_vol*std::sqrt(d_expiry)*z;
        double s = d_spot*std::exp(trend+random);
        S += s;
    }

    return S/N;
}

double MonteCarlo::GetBetaEstimate(int N, bool momentMatching)
{
    // get the sample mean for moment matching
    double mean = GetMeanEstimate(N);

    double V=0.0;
    double SV=0.0;
    double S2=0.0;
    int idum=1;
    for (int k=0; k<N; k++) {
        double z = (*d_ran)(&idum);
        double trend = (d_rate-d_div-0.5*d_vol*d_vol)*d_expiry;
        double random = d_vol*std::sqrt(d_expiry)*z;

        double s = d_spot*std::exp(trend+random);
        if ( momentMatching ) {
            s *= d_spot*std::exp(d_rate*d_expiry)/mean;
        }

        double v = std::exp(-d_rate*d_expiry)*(*d_payoff)(s);
        SV += s*v;
        S2 += s*s;
        V += v;
    }

    double sExpectation=mean;
    if ( momentMatching ) {
        sExpectation = d_spot*std::exp(d_rate*d_expiry);
    }
    double numerator = SV - V*sExpectation;
    double denominator = S2 - N*sExpectation*sExpectation;
    return numerator/denominator;
}

MonteCarlo::MonteCarlo( const Payoff & payoff, 
                        double expiry,
                        double spot,
                        double rate,
                        double div,
                        double vol,
                        double (*ran)(int*) )
                      : d_payoff(&payoff),
                        d_expiry(expiry),
                        d_spot(spot),
                        d_rate(rate),
                        d_div(div),
                        d_vol(vol),
                        d_spot2(0),
                        d_div2(0),
                        d_vol2(0),
                        d_rho(0),
                        d_ran(ran)
{
    d_numOfUnderlyings = 1;
    // Calculate the Black-Scholes value
    // in case of a vanilla call or put.
    d_strike = payoff.strike();
    std::tuple<OptionValue,OptionValue> res
        = BlackScholes(d_expiry,d_strike,d_spot,d_rate,d_div,d_vol);
    if ((*d_payoff)(0) > 0) { d_BlackScholes = std::get<1>(res);}
    else { d_BlackScholes = std::get<0>(res);}
}
        
MonteCarlo::MonteCarlo( const Payoff & payoff,
                        double expiry,
                        double spot1, double spot2,
                        double rate,
                        double div1, double div2,
                        double vol1, double vol2,
                        double rho,
                        double (*ran)(int*) )
                      : d_payoff(&payoff),
                        d_expiry(expiry),
                        d_spot(spot1),
                        d_rate(rate),
                        d_div(div1),
                        d_vol(vol1),
                        d_spot2(spot2),
                        d_div2(div2),
                        d_vol2(vol2),
                        d_rho(rho),
                        d_ran(ran)
{
    d_numOfUnderlyings = 2;
}

OptionValue MonteCarlo::BlackScholesValue() const
{
    return d_BlackScholes; 
}

OptionValue MonteCarlo::evaluate(int N, 
                                 bool controlVariate, 
                                 bool useAntithetic,
                                 bool momentMatching)
{
    double b=0.0; // control variate parameter
    if ( controlVariate ) {
        b = GetBetaEstimate(N, momentMatching);
    }

    double mean=0.0; // 1st moment matching
    if ( momentMatching ) {
        mean = GetMeanEstimate(N);
    }

    int idum=1;
    double V = 0.0;
    double Var = 0.0;
    for (int k=0; k<N; k++) {
        double z = (*d_ran)(&idum);
        double trend = (d_rate-d_div-0.5*d_vol*d_vol)*d_expiry;
        double random = d_vol*std::sqrt(d_expiry)*z;

        double s = d_spot*std::exp(trend+random);
        if ( momentMatching ) {
            s *= d_spot*std::exp(d_rate*d_expiry)/mean;
        }
        
        double v = std::exp(-d_rate*d_expiry)*(*d_payoff)(s);
        
        if ( useAntithetic ) {
            double s2 = s*std::exp(-2*random);
            double v2 = std::exp(-d_rate*d_expiry)*(*d_payoff)(s2);
            s = (s+s2)/2.0;
            v = (v+v2)/2.0;
        }

        if ( controlVariate ) {
            v -= b*(s-std::exp(d_rate*d_expiry)*d_spot);
        }

        V += v;
        Var += v*v;
    }

    V /= N;
    Var /= (N-1);
    Var -= N*V*V/(N-1);

    OptionValue optionValue;
    optionValue.price = V;
    optionValue.priceStd = std::sqrt(Var);
    optionValue.priceErr = optionValue.priceStd/std::sqrt(N);
    return optionValue;
}

OptionValue MonteCarlo::evaluateControlVariate(int N)
{
    return evaluate(N, true, false, false);
}
        
OptionValue MonteCarlo::evaluateMomentMatching(int N)
{
    return evaluate(N, false, false, true);
}
        
OptionValue MonteCarlo::evaluateControlVariateMomentMatching(int N)
{
    return evaluate(N, true, false, true); 
}
        
OptionValue MonteCarlo::evaluateUseAntithetic(int N)
{
    return evaluate(N, false, true, false);
}

OptionValue MonteCarlo::evaluateBasket(int N)
{
    assert(d_numOfUnderlyings == 2);
    double rho2 = std::sqrt(1-d_rho*d_rho);
    
    int idum=1;
    double V = 0.0;
    double Var = 0.0;
    for (int k=0; k<N; k++) {
        double z1 = (*d_ran)(&idum);
        double z2 = d_rho*z1 + rho2*(*d_ran)(&idum);
        double trend1 = (d_rate-d_div-0.5*d_vol*d_vol)*d_expiry;
        double trend2 = (d_rate-d_div2-0.5*d_vol2*d_vol2)*d_expiry;
        double random1 = d_vol*std::sqrt(d_expiry)*z1;
        double random2 = d_vol2*std::sqrt(d_expiry)*z2;

        double s1 = d_spot*std::exp(trend1+random1);
        double s2 = d_spot2*std::exp(trend2+random2);
        double v = std::exp(-d_rate*d_expiry)*(*d_payoff)(s1,s2);
        
        V += v;
        Var += v*v;
    }

    V /= N;
    Var /= (N-1);
    Var -= N*V*V/(N-1);

    OptionValue optionValue;
    optionValue.price = V;
    optionValue.priceStd = std::sqrt(Var);
    optionValue.priceErr = optionValue.priceStd/std::sqrt(N);
    return optionValue;
}

OptionValue MonteCarlo::evaluatePathDependentBasket(int N, int M)
{
    assert(d_numOfUnderlyings == 2);
    double rho2 = std::sqrt(1-d_rho*d_rho);
    double dt = d_expiry/M;
    
    int idum=1;
    double V = 0.0;
    double Var = 0.0;
    for (int k=0; k<N; k++) {
        double s1 = d_spot;
        double s2 = d_spot2;
        double minS1S2 = s1+s2;
        for (int j=0; j<M; j++) {
            double z1 = (*d_ran)(&idum);
            double z2 = d_rho*z1 + rho2*(*d_ran)(&idum);
            double trend1 = (d_rate-d_div-0.5*d_vol*d_vol)*dt;
            double trend2 = (d_rate-d_div2-0.5*d_vol2*d_vol2)*dt;
            double random1 = d_vol*std::sqrt(dt)*z1;
            double random2 = d_vol2*std::sqrt(dt)*z2;
            s1 *= std::exp(trend1+random1);
            s2 *= std::exp(trend2+random2);
            if ( minS1S2 < s1+s2 ) {
                minS1S2 = s1+s2;
            }
        }

        double v = std::exp(-d_rate*d_expiry)*(*d_payoff)(minS1S2);
        
        V += v;
        Var += v*v;
    }

    V /= N;
    Var /= (N-1);
    Var -= N*V*V/(N-1);

    OptionValue optionValue;
    optionValue.price = V;
    optionValue.priceStd = std::sqrt(Var);
    optionValue.priceErr = optionValue.priceStd/std::sqrt(N);
    return optionValue;
}


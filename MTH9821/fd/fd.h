#ifndef FD_H
#define FD_H 
#include <option_value.h>
#include <evaluator.h>
#include <heat_pde.h>
#include <cmath>
#include <tuple>
#include <vector>
#include <iostream>

class VanillaPutTerminalCondition : public Evaluator
{
    public:

        VanillaPutTerminalCondition() {}
        VanillaPutTerminalCondition( double expiry,
                                     double spot,
                                     double strike,
                                     double rate,
                                     double div,
                                     double vol )
        {
            d_strike = strike;
            double c = (rate-div)/(vol*vol);
            d_a = c - 0.5;
            d_b = (c+0.5)*(c+0.5) + 2*div/(vol*vol);
        }

        virtual double operator()(double s) const
        {
            return d_strike*std::exp(d_a*s)*std::max(1-std::exp(s),0.0);
        }
        
    private:

        double d_strike;
        double d_a;
        double d_b;
};

class VanillaPutLeftBoundaryCondition : public Evaluator
{
    public:

        VanillaPutLeftBoundaryCondition() {}
        VanillaPutLeftBoundaryCondition( double expiry,
                                         double spot,
                                         double strike,
                                         double rate,
                                         double div,
                                         double vol )
                                       : d_expiry(expiry),
                                         d_spot(spot),
                                         d_strike(strike),
                                         d_rate(rate),
                                         d_div(div),
                                         d_vol(vol)
        {
            d_strike = strike;
            double c = (rate-div)/(vol*vol);
            d_a = c - 0.5;
            d_b = (c+0.5)*(c+0.5) + 2*div/(vol*vol);
            double xmid = std::log(spot/strike) + (rate-div-0.5*vol*vol)*expiry;
            d_xl = xmid - 3*vol*std::sqrt(expiry);
        }

        virtual double operator()(double s) const
        {
            double c1 = std::exp(-2*d_rate*s/(d_vol*d_vol));
            double c2 = std::exp(d_xl - 2*d_div*s/(d_vol*d_vol));
            return d_strike*std::exp(d_a*d_xl+d_b*s)*(c1-c2); 
        }
        
    private:

        double d_expiry;
        double d_spot;
        double d_strike;
        double d_rate;
        double d_div;
        double d_vol;
        double d_a;
        double d_b;
        double d_xl;
};

class VanillaPutRightBoundaryCondition : public Evaluator
{
    public:

        VanillaPutRightBoundaryCondition() {}
        VanillaPutRightBoundaryCondition( double expiry,
                                          double spot,
                                          double strike,
                                          double rate,
                                          double div,
                                          double vol )
        {}

        virtual double operator()(double s) const
        {
            return 0; 
        }
        
};

class AmericanPutLeftBoundaryCondition : public Evaluator
{
    public:

        AmericanPutLeftBoundaryCondition() {}
        AmericanPutLeftBoundaryCondition( double expiry,
                                          double spot,
                                          double strike,
                                          double rate,
                                          double div,
                                          double vol )
                                        : d_expiry(expiry),
                                          d_spot(spot),
                                          d_strike(strike),
                                          d_rate(rate),
                                          d_div(div),
                                          d_vol(vol)
        {
            d_strike = strike;
            double c = (rate-div)/(vol*vol);
            d_a = c - 0.5;
            d_b = (c+0.5)*(c+0.5) + 2*div/(vol*vol);
            double xmid = std::log(spot/strike) + (rate-div-0.5*vol*vol)*expiry;
            d_xl = xmid - 3*vol*std::sqrt(expiry);
            d_c = 1-std::exp(d_xl);
        }

        virtual double operator()(double s) const
        {
            return d_strike*d_c*std::exp(d_a*d_xl+d_b*s);
        }
        
    private:

        double d_expiry;
        double d_spot;
        double d_strike;
        double d_rate;
        double d_div;
        double d_vol;
        double d_a;
        double d_b;
        double d_c;
        double d_xl;
};

class AmericanPutEarlyExercisePremium : public Evaluator
{
    public:

        AmericanPutEarlyExercisePremium() {}
        AmericanPutEarlyExercisePremium( double expiry,
                                         double spot,
                                         double strike,
                                         double rate,
                                         double div,
                                         double vol )
                                       : d_expiry(expiry),
                                         d_spot(spot),
                                         d_strike(strike),
                                         d_rate(rate),
                                         d_div(div),
                                         d_vol(vol)
        {
            d_strike = strike;
            double c = (rate-div)/(vol*vol);
            d_a = c - 0.5;
            d_b = (c+0.5)*(c+0.5) + 2*div/(vol*vol);
        }

        virtual double operator()(double x, double t) const
        {
            return d_strike*std::exp(d_a*x+d_b*t)*std::max(1-std::exp(x), 0.0);
        }
        
    private:

        double d_expiry;
        double d_spot;
        double d_strike;
        double d_rate;
        double d_div;
        double d_vol;
        double d_a;
        double d_b;
};

enum FiniteDifferenceMethod 
{
    EulerForward,
    EulerBackwardByLU,
    EulerBackwardBySOR,
    CrankNicolsonByLU,
    CrankNicolsonBySOR,
    AmericanEulerForward,
    AmericanEulerBackwardByLU,
    AmericanEulerBackwardBySOR,
    AmericanCrankNicolsonByLU,
    AmericanCrankNicolsonBySOR,
};

class FiniteDifference 
{
    public:

        FiniteDifference( double expiry,
                          double spot,
                          double strike,
                          double rate,
                          double div,
                          double vol,
                          const Evaluator & terminalCondition,
                          const Evaluator & leftBoundaryCondition,
                          const Evaluator & rightBoundaryCondition,
                          const Evaluator & earlyExercisePremium=Evaluator() );
        
        OptionValue BlackScholesValue() const;
        void evaluate( int M, double alphaTemp, 
                       FiniteDifferenceMethod fdm=EulerForward,
                       double omega=1, double vExact=0 ); 

    private:

        double d_expiry;
        double d_spot;
        double d_strike;
        double d_rate;
        double d_div;
        double d_vol;

        const Evaluator * d_f;
        const Evaluator * d_gl;
        const Evaluator * d_gr;
        const Evaluator * d_prem;

        double d_xl;
        double d_xr;
        double d_tf;

        HeatPDE d_h;

        // applicable when Black-Scholes
        bool d_isPut;
        OptionValue d_BlackScholes;
};

#endif /* FD_H */

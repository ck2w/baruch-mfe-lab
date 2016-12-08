#ifndef FD_H
#define FD_H 
#include <option_value.h>
#include <evaluator.h>
#include <heat_pde.h>
#include <cmath>
#include <tuple>
#include <vector>
#include <iostream>

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
                          Evaluator & terminalCondition,
                          Evaluator & leftBoundaryCondition,
                          Evaluator & rightBoundaryCondition,
                          const Evaluator & earlyExercisePremium=Evaluator() );

        void setDefaultDomain();
        void setExpandedDomain(int M, double alpha, double T, double qDiv=0); 
        void setDomain(double xl, double xr, double tf, double ti);

        void discretizeDomainByTimeStepsAndAlphaTemp(int M, double alphaTemp);
        void discretizeDomainByIntervalsAndAlphaTemp(int N, double alphaTemp);
        void discretizeDomainByTimeStepsAndAlphaFixed(int M, double alpha);
        void discretizeDomain(int M, int N);
        void overrideTerminalCondition(const std::vector<double> & u);

        double getTi() const
        {
            return d_ti;
        }

        double getTf() const 
        {
            return d_tf;
        }

        double getXl() const
        {
            return d_xl;
        }

        double getXr() const
        {
            return d_xr;
        }

        int getTimeSteps() const
        {
            return d_M;
        }

        int getIntervals() const
        {
            return d_N;
        }

        double getAlpha() const
        {
            return d_alpha;
        }
        
        OptionValue BlackScholesValue() const;
       
        std::vector<double> evaluate( FiniteDifferenceMethod fdm=EulerForward, 
                                      double omega=1, 
                                      bool varReduction=false, 
                                      double vExact=0 ); 

    private:

        double d_expiry;
        double d_spot;
        double d_strike;
        double d_rate;
        double d_div;
        double d_vol;

        Evaluator * d_f;
        Evaluator * d_gl;
        Evaluator * d_gr;
        const Evaluator * d_prem;
        
        double d_xl;
        double d_xr;
        double d_tf;
        double d_ti;

        int d_M;
        int d_N;
        double d_alpha;

        HeatPDE d_h;

        bool d_terminalConditionOverriden;
        std::vector<double> d_u0;

        // applicable when Black-Scholes
        bool d_isPut;
        OptionValue d_BlackScholes;

        OptionValue getOptionValue( const std::vector<double> & u, 
                                    const std::vector<double> & uOld,
                                    double dt );
};

#endif /* FD_H */

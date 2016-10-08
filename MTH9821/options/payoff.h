#ifndef PAYOFF_H
#define PAYOFF_H 
#include <cmath>
#include <algorithm>
#include <black_scholes.h>

class Payoff
{
    public:

        Payoff() {}
        virtual double operator()(double s) const=0;
        virtual double strike() const=0;
};

class CallPayoff : public Payoff
{
    public:

        CallPayoff(double strike): d_strike(strike) {}

        virtual double operator()(double s) const
        {
            return std::max(s-d_strike, 0.0);
        }

        virtual double strike() const
        {
            return d_strike;
        }

    private:

        double d_strike;
};

class PutPayoff : public Payoff
{
    public:

        PutPayoff(double strike): d_strike(strike) {}

        virtual double operator()(double s) const
        {
            return std::max(d_strike-s, 0.0);
        }
        
        virtual double strike() const
        {
            return d_strike;
        }

    private:

        double d_strike;
};

#endif /* PAYOFF_H */

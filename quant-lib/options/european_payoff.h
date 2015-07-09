#ifndef EUROPEAN_PAYOFF_H
#define EUROPEAN_PAYOFF_H 
#include <payoff.h>

class European_call_payoff : public Payoff
{
    public:

        European_call_payoff(double strike): d_strike(strike) {}

        virtual double operator()(double spot) const;
        virtual Payoff* clone() const;
        virtual ~European_call_payoff() {}

    private:

        double d_strike;
};

class European_put_payoff : public Payoff
{
    public:

        European_put_payoff(double strike): d_strike(strike) {}
        
        virtual double operator()(double spot) const;
        virtual Payoff* clone() const;
        virtual ~European_put_payoff() {}

    private:

        double d_strike;
}; 

#endif /* EUROPEAN_PAYOFF_H */

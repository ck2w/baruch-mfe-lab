#ifndef PAYOFF_BRIDGE_H
#define PAYOFF_BRIDGE_H

class Payoff
{
    public:
        
        Payoff() {}

        virtual double operator()(double spot) const=0;
        virtual Payoff* clone() const=0;
        virtual ~Payoff() {}
};

class Payoff_bridge
{
    public:

        Payoff_bridge() {}
        Payoff_bridge( const Payoff_bridge & payoff_bridge );
        Payoff_bridge( const Payoff & payoff );

        Payoff_bridge& operator=( const Payoff_bridge & payoff_bridge );

        double operator()(double spot) const;
        ~Payoff_bridge();

    private:

        Payoff* d_payoff_ptr;
};

#endif /* PAYOFF_BRIDGE_H */

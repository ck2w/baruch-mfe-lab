#include <payoff.h>

Payoff_bridge::Payoff_bridge(const Payoff_bridge & payoff_bridge)
{
    d_payoff_ptr = payoff_bridge.d_payoff_ptr->clone();
}

Payoff_bridge::Payoff_bridge(const Payoff & payoff)
{
    d_payoff_ptr = payoff.clone();
}

Payoff_bridge& Payoff_bridge::operator=(const Payoff_bridge & payoff_bridge)
{
    if (this != &payoff_bridge)
    {
        delete d_payoff_ptr;
        d_payoff_ptr = payoff_bridge.d_payoff_ptr->clone();
    }

    return *this;
}

double Payoff_bridge::operator()(double spot) const
{
    return d_payoff_ptr->operator()(spot);
}

Payoff_bridge::~Payoff_bridge()
{
    delete d_payoff_ptr;
}


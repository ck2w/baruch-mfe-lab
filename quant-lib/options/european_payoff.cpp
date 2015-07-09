#include <european_payoff.h>
#include <cmath>
#include <algorithm>

double European_call_payoff::operator()(double spot) const
{
    return std::max(spot-d_strike, 0.0);
}

// virtual copy constructor
Payoff* European_call_payoff::clone() const
{
    return new European_call_payoff(*this);
}

double European_put_payoff::operator()(double spot) const
{
    return std::max(d_strike-spot, 0.0);
}

// virual copy constructor
Payoff* European_put_payoff::clone() const
{
    return new European_put_payoff(*this);
}


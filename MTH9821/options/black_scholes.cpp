#include <black_scholes.h>
#include <option_value.h>
#include <cmath>

std::tuple<OptionValue,OptionValue> BlackScholes(double expiry,
                                                 double strike,
                                                 double spot,
                                                 double rate,
                                                 double div,
                                                 double vol,
                                                 double qDiv)
{
    const double PI = 3.141592653589793238463;

    double vT1 = vol*vol*expiry; 
    double vT2 = std::sqrt(vT1);
    double d1 = (std::log(spot*(1-qDiv)/strike)+(rate-div)*expiry+0.5*vT1)/vT2;
    double d2 = d1-vT2;
    double Nd1 = 0.5*std::erfc(-d1/std::sqrt(2));
    double Nd2 = 0.5*std::erfc(-d2/std::sqrt(2));
    double SPvf = spot*std::exp(-div*expiry)*(1-qDiv);
    double KPvf = strike*std::exp(-rate*expiry);

    double call = SPvf*Nd1-KPvf*Nd2; 
    double put = KPvf*(1-Nd2)-SPvf*(1-Nd1);
    double callDelta = std::exp(-div*expiry)*Nd1;
    double putDelta = std::exp(-div*expiry)*(Nd1-1);
    double gamma 
        = std::exp(-div*expiry-0.5*d1*d1)/(vol*spot*std::sqrt(2*PI*expiry));
    double callTheta = div*spot*std::exp(-div*expiry)*Nd1
                     - rate*strike*std::exp(-rate*expiry)*Nd2
                     - 0.5*gamma*vol*vol*spot*spot;
    double putTheta = div*spot*std::exp(-div*expiry)*(Nd1-1)
                    - rate*strike*std::exp(-rate*expiry)*(Nd2-1)
                    - 0.5*gamma*vol*vol*spot*spot;

    OptionValue callOptionValue;
    callOptionValue.price = call;
    callOptionValue.delta = callDelta;
    callOptionValue.gamma = gamma;
    callOptionValue.theta = callTheta;

    OptionValue putOptionValue;
    putOptionValue.price = put;
    putOptionValue.delta = putDelta;
    putOptionValue.gamma = gamma;
    putOptionValue.theta = putTheta;

    return std::make_tuple(callOptionValue, putOptionValue);
}


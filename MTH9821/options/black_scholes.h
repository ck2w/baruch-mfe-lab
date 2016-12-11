#ifndef BLACK_SCHOLES_H
#define BLACK_SCHOLES_H 
#include <tuple>
#include <option_value.h>

// qDiv is a proportion of the stock price paid out as dividend
std::tuple<OptionValue, OptionValue> BlackScholes(double expiry,
                                                  double strike,
                                                  double spot,
                                                  double rate,
                                                  double div,
                                                  double vol,
                                                  double qDiv=0.0);

#endif /* BLACK_SCHOLES_H */

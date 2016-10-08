#ifndef BLACK_SCHOLES_H
#define BLACK_SCHOLES_H 
#include <tuple>
#include <option_value.h>

std::tuple<OptionValue, OptionValue> BlackScholes(double expiry,
                                                  double strike,
                                                  double spot,
                                                  double rate,
                                                  double div,
                                                  double vol);

#endif /* BLACK_SCHOLES_H */

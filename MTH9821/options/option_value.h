#ifndef OPTION_H
#define OPTION_H

/****************************/
/*  Option Valuation Result */
/****************************/
struct OptionValue 
{
    double price;
    double delta;
    double gamma;
    double theta;

    // used by Monte Carlo methods only
    double priceStd;
    double priceErr;
};

#endif /* OPTION_H */

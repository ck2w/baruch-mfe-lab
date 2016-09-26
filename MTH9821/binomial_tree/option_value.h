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
};

#endif /* OPTION_H */

#ifndef DIVIDEND_H
#define DIVIDEND_H 

/****************************/
/*  A discrete dividend     */
/****************************/
struct DiscreteDividend
{
    double t; // time of dividend
    double v; // fixed value or proportionality
    bool   isFixed; // fixed value or proportional
};

#endif /* DIVIDEND_H */

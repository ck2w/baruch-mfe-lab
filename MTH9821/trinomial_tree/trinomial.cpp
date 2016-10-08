#include <payoff.h>
#include <option_value.h>
#include <trinomial.h>
#include <black_scholes.h>
#include <tuple>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>

TrinomialTree::TrinomialTree( const Payoff & payoff, 
                              double expiry,
                              double spot,
                              double rate,
                              double div,
                              double vol )
                            : d_payoff(&payoff),
                              d_expiry(expiry),
                              d_spot(spot),
                              d_rate(rate),
                              d_div(div),
                              d_vol(vol)
{
    // Calculate the Black-Scholes value
    // in case of a vanilla call or put.
    d_strike = payoff.strike();
    std::tuple<OptionValue,OptionValue> res
        = BlackScholes(d_expiry,d_strike,d_spot,d_rate,d_div,d_vol);
    if ((*d_payoff)(0) > 0) { d_BlackScholes = std::get<1>(res);}
    else { d_BlackScholes = std::get<0>(res);}
}

void TrinomialTree::print_config(int N)
{
    double dt = d_expiry/N;
    double u = std::exp(d_vol*std::sqrt(3*dt));
    double d = 1/u;
    double pc = (d_rate-d_div-0.5*d_vol*d_vol)*std::sqrt(dt/(12*d_vol*d_vol));
    double pu = 1.0/6.0 + pc;
    double pd = 1.0/6.0 - pc;
    double pm = 2.0/3.0;
    std::cout << "dt = " << dt 
              << ", u = " << u
              << ", d = " << d
              << ", pu = " << pu
              << ", pd = " << pd
              << ", pm = " << pm
              << std::endl;
}

OptionValue TrinomialTree::BlackScholesValue() const
{
    return d_BlackScholes; 
}

OptionValue TrinomialTree::evaluate(int N, 
                                    bool isAmerican, 
                                    bool useBS, 
                                    bool varReduction)
{
    double dt = d_expiry/N;
    double u = std::exp(d_vol*std::sqrt(3*dt));
    double d = 1/u;
    double pc = (d_rate-d_div-0.5*d_vol*d_vol)*std::sqrt(dt/(12*d_vol*d_vol));
    double pu = 1.0/6.0 + pc;
    double pd = 1.0/6.0 - pc;
    double pm = 2.0/3.0;

    // If the payoff is vanilla, check call or put.
    // This is applicable only if usBS is true.
    bool isCall = true;
    if ( (*d_payoff)(0) > 0 ) { isCall = false; }

    /**********************************************
     * The storage is implemented as an array of  *
     * size M=(N+1)^2 where N is the number of    *
     * time steps. The i'th node at the k'th step *
     * is indexed at k^2 + i, i<=k.               *
     *                                            *
     *                        16                  *
     *                    09  17                  *
     *                 4  10  18                  *
     *              1  5  11  19                  *
     *           0  2  6  12  20                  *
     *              3  7  13  21                  *
     *                 8  14  22                  *
     *                    15  23                  *
     *                        24                  *
     *                                            *
     **********************************************/
    double M = (N+1)*(N+1);
    std::vector<double> tree(M,0);
    std::vector<double> eTree;
    // If use variance reduction, 
    // run a European binomial tree in parallel.
    if (varReduction) { eTree.resize(M,0); }

    // Terminal condition at time T, step N.
    int nn=N*N; //index base for step N
    int n2=2*N; //number of nodes for step N (minus 1)
    for (int i=0; i<=n2; i++) {
        double s = std::pow(u,N-i)*d_spot;
        tree[nn+i] = (*d_payoff)(s);
        if (varReduction) { eTree[nn+i] = (*d_payoff)(s); }
    }

    // Back-traverse the tree
    for (int k=(N-1); k>=0; k--) {
        int nn1=k*k;         //index base for step k
        int nn2=(k+1)*(k+1); //index base for step k+1
        int n2=2*k;          //number of nodes for step k (minus 1)
        for (int i=0; i<=n2; i++) {
            if ( useBS && k==(N-1) ) {
                // Binomial Black-Scholes
                double s = std::pow(u,k-i)*d_spot;
                std::tuple<OptionValue,OptionValue> res
                    = BlackScholes(dt,d_strike,s,d_rate,d_div,d_vol);
                if (isCall) { tree[nn1+i] = std::get<0>(res).price; }
                else { tree[nn1+i] = std::get<1>(res).price; }
                
                if (varReduction) { eTree[nn1+i] = tree[nn1+i]; }
            }
            else {
                double disc = std::exp(-d_rate*dt);
                tree[nn1+i] 
                    = disc*(pu*tree[nn2+i]+
                            pm*tree[nn2+i+1]+
                            pd*tree[nn2+i+2]);
                if (varReduction) { 
                    eTree[nn1+i] = disc*(pu*eTree[nn2+i]+
                                         pm*eTree[nn2+i+1]+
                                         pd*eTree[nn2+i+2]);
                }
            }

            if (isAmerican) {
                double s = std::pow(u,k-i)*d_spot;
                double payoff = (*d_payoff)(s);
                if (payoff > tree[nn1+i]) {
                    tree[nn1+i] = payoff;
                }
            }
        }
    }

    // record greeks
    double v00 = tree[0];
    double v10 = tree[1];
    double v11 = tree[2];
    double v12 = tree[3];
    double v20 = tree[4];
    double v21 = tree[5];
    double v22 = tree[6];
    double v23 = tree[7];
    double v24 = tree[8];
    double s00 = d_spot;
    double s10 = d_spot*u;
    double s11 = d_spot;
    double s12 = d_spot*d;
    double s20 = s10*u;
    double s21 = s10;
    double s22 = s10*d;
    double s23 = s12;
    double s24 = s12*d;

    // Control variate technique
    if (varReduction) {
        std::tuple<OptionValue,OptionValue> res00
            = BlackScholes(d_expiry,d_strike,s00,d_rate,d_div,d_vol);
        std::tuple<OptionValue,OptionValue> res10
            = BlackScholes(d_expiry-dt,d_strike,s10,d_rate,d_div,d_vol);
        std::tuple<OptionValue,OptionValue> res11
            = BlackScholes(d_expiry-dt,d_strike,s11,d_rate,d_div,d_vol);
        std::tuple<OptionValue,OptionValue> res12
            = BlackScholes(d_expiry-dt,d_strike,s12,d_rate,d_div,d_vol);
        std::tuple<OptionValue,OptionValue> res20
            = BlackScholes(d_expiry-2*dt,d_strike,s20,d_rate,d_div,d_vol);
        std::tuple<OptionValue,OptionValue> res21
            = BlackScholes(d_expiry-2*dt,d_strike,s21,d_rate,d_div,d_vol);
        std::tuple<OptionValue,OptionValue> res22
            = BlackScholes(d_expiry-2*dt,d_strike,s22,d_rate,d_div,d_vol);
        std::tuple<OptionValue,OptionValue> res23
            = BlackScholes(d_expiry-2*dt,d_strike,s23,d_rate,d_div,d_vol);
        std::tuple<OptionValue,OptionValue> res24
            = BlackScholes(d_expiry-2*dt,d_strike,s24,d_rate,d_div,d_vol);

        double vBS00, vBS10, vBS11, vBS12, vBS20, vBS21, vBS22, vBS23, vBS24;
        if (isCall) { 
            vBS00 = std::get<0>(res00).price; 
            vBS10 = std::get<0>(res10).price; 
            vBS11 = std::get<0>(res11).price; 
            vBS12 = std::get<0>(res12).price; 
            vBS20 = std::get<0>(res20).price; 
            vBS21 = std::get<0>(res21).price; 
            vBS22 = std::get<0>(res22).price; 
            vBS23 = std::get<0>(res23).price; 
            vBS24 = std::get<0>(res24).price; 
        }
        else { 
            vBS00 = std::get<1>(res00).price; 
            vBS10 = std::get<1>(res10).price; 
            vBS11 = std::get<1>(res11).price; 
            vBS12 = std::get<1>(res12).price; 
            vBS20 = std::get<1>(res20).price; 
            vBS21 = std::get<1>(res21).price; 
            vBS22 = std::get<1>(res22).price; 
            vBS23 = std::get<1>(res23).price; 
            vBS24 = std::get<1>(res24).price; 
        }

        double ev00 = eTree[0];
        double ev10 = eTree[1];
        double ev11 = eTree[2];
        double ev12 = eTree[3];
        double ev20 = eTree[4];
        double ev21 = eTree[5];
        double ev22 = eTree[6];
        double ev23 = eTree[7];
        double ev24 = eTree[8];

        v00 += (vBS00-ev00);
        v10 += (vBS10-ev10);
        v11 += (vBS11-ev11);
        v12 += (vBS12-ev12);
        v20 += (vBS20-ev20);
        v21 += (vBS21-ev21);
        v22 += (vBS22-ev22);
        v23 += (vBS23-ev23);
        v24 += (vBS24-ev24);
    }

    OptionValue optionValue;
    optionValue.price = v00;
    optionValue.delta = (v10-v12)/(s10-s12);
    optionValue.gamma = ((v20-v22)/(s20-s22)-(v22-v24)/(s22-s24))/(s10-s12);
    optionValue.theta = (v11-v00)/dt;
    return optionValue;
}

OptionValue TrinomialTree::evaluateTBS(int N, bool isAmerican)
{
    return evaluate(N, isAmerican, true);
}
        
OptionValue TrinomialTree::evaluateVarReduction(int N, bool useBS)
{
    return evaluate(N, true, useBS, true);
}

OptionValue TrinomialTree::evaluateVarReductionTBS(int N)
{
    return evaluateVarReduction(N, true);
}



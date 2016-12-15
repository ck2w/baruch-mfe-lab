#include <payoff.h>
#include <option_value.h>
#include <binomial.h>
#include <black_scholes.h>
#include <tuple>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>

BinomialTree::BinomialTree( const Payoff & payoff, 
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

void BinomialTree::print_config(int N)
{
    double dt = d_expiry/N;
    double u = std::exp(d_vol*std::sqrt(dt));
    double d = 1/u;
    double p = (std::exp((d_rate-d_div)*dt)-d)/(u-d);
    std::cout << "dt = " << dt 
              << ", u = " << u
              << ", d = " << d
              << ", p = " << p
              << std::endl;
}

OptionValue BinomialTree::BlackScholesValue() const
{
    return d_BlackScholes; 
}

OptionValue BinomialTree::evaluate(int N, 
                                   bool isAmerican, 
                                   bool useBS, 
                                   bool varReduction)
{
    double dt = d_expiry/N;
    double u = std::exp(d_vol*std::sqrt(dt));
    double d = 1/u;
    double p = (std::exp((d_rate-d_div)*dt)-d)/(u-d);

    // If the payoff is vanilla, check call or put.
    // This is applicable only if usBS is true.
    bool isCall = true;
    if ( (*d_payoff)(0) > 0 ) { isCall = false; }

    /*********************************************
     * The storage is implemented as an array of *
     * size M=(N+1)(N+2)/2 where N is the number *
     * of time steps. The i'th node at the k'th  *
     * step is indexed at k(k+1)/2 + i, i<=k.    *
     *                                           *
     *                   10                      *
     *                 6                         *
     *               3   11                      *
     *             1   7                         *
     *           0   4   12                      *
     *             2   8                         *
     *               5   13                      *
     *                 9                         *
     *                   14                      *
     *                                           *
     *********************************************/
    double M = (N+1)*(N+2)/2;
    std::vector<double> tree(M,0);
    std::vector<double> eTree;
    // If use variance reduction, 
    // run a European binomial tree in parallel.
    if (varReduction) { eTree.resize(M,0); }

    // Terminal condition at time T, step N.
    int n=N*(N+1)/2;
    for (int i=0; i<=N; i++) {
        double s = std::pow(u,N-2*i)*d_spot;
        tree[n+i] = (*d_payoff)(s);
        if (varReduction) { eTree[n+i] = (*d_payoff)(s); }
    }

    // Back-traverse the tree
    for (int k=(N-1); k>=0; k--) {
        int n1=k*(k+1)/2;
        int n2=(k+1)*(k+2)/2;
        for (int i=0; i<=k; i++) {
            if ( useBS && k==(N-1) ) {
                // Binomial Black-Scholes
                double s = std::pow(u,k-2*i)*d_spot;
                std::tuple<OptionValue,OptionValue> res
                    = BlackScholes(dt,d_strike,s,d_rate,d_div,d_vol);
                if (isCall) { tree[n1+i] = std::get<0>(res).price; }
                else { tree[n1+i] = std::get<1>(res).price; }
                
                if (varReduction) { eTree[n1+i] = tree[n1+i]; }
            }
            else {
                double disc = std::exp(-d_rate*dt);
                tree[n1+i] = disc*(p*tree[n2+i]+(1-p)*tree[n2+i+1]);
                if (varReduction) { 
                    eTree[n1+i] = disc*(p*eTree[n2+i]+(1-p)*eTree[n2+i+1]);
                }
            }

            if (isAmerican) {
                double s = std::pow(u,k-2*i)*d_spot;
                double payoff = (*d_payoff)(s);
                if (payoff > tree[n1+i]) {
                    tree[n1+i] = payoff;
                }
            }
        }
    }
    
    // record greeks
    double v00 = tree[0];
    double v10 = tree[1];
    double v11 = tree[2];
    double v20 = tree[3];
    double v21 = tree[4];
    double v22 = tree[5];
    double s00 = d_spot;
    double s10 = d_spot*u;
    double s11 = d_spot*d;
    double s20 = s10*u;
    double s21 = s10*d;
    double s22 = s11*d;

    // Control variate technique
    OptionValue varCorrection; 
    if (varReduction) {
        std::tuple<OptionValue,OptionValue> res00
            = BlackScholes(d_expiry,d_strike,s00,d_rate,d_div,d_vol);

        double vBS;
        double deltaBS;
        double gammaBS;
        double thetaBS;
        if (isCall) { 
            vBS = std::get<0>(res00).price; 
            deltaBS = std::get<0>(res00).delta;
            gammaBS = std::get<0>(res00).gamma;
            thetaBS = std::get<0>(res00).theta;
        }
        else { 
            vBS = std::get<1>(res00).price; 
            deltaBS = std::get<1>(res00).delta;
            gammaBS = std::get<1>(res00).gamma;
            thetaBS = std::get<1>(res00).theta;
        }

        double ev00 = eTree[0];
        double ev10 = eTree[1];
        double ev11 = eTree[2];
        double ev20 = eTree[3];
        double ev21 = eTree[4];
        double ev22 = eTree[5];
        
        double eDelta = (ev10 - ev11)/(s10-s11);
        double eGamma = 2*((ev20-ev21)/(s20-s21)-(ev21-ev22)/(s21-s22))/(s20-s22);
        double eTheta = (ev21-ev00)/(2*dt);

        varCorrection.price = vBS-ev00;
        varCorrection.delta = deltaBS-eDelta;
        varCorrection.gamma = gammaBS-eGamma;
        varCorrection.theta = thetaBS-eTheta;
    }
   
    if (isAmerican && varReduction) {
        std::cout << std::endl
                  << "======== First three time steps in European tree ========" 
                  << std::endl;
        std::cout << "First time step:" 
                  << std::fixed 
                  << std::setprecision(9) 
                  << eTree[0]
                  << std::endl;
        std::cout << "Second time step:" 
                  << std::fixed 
                  << std::setprecision(9) 
                  << eTree[1] << "," << eTree[2]
                  << std::endl;
        std::cout << "Third time step:" 
                  << std::fixed 
                  << std::setprecision(9) 
                  << eTree[3] << "," << eTree[4] << "," << eTree[5]
                  << std::endl;
        std::cout << "---------------------------------------------------------" 
                  << std::endl;
        
    }

    // print out option values on some nodes
    std::cout << std::endl
              << "======== First three time steps in Amercian tree ========" 
              << std::endl;
    std::cout << "First time step:" 
              << std::fixed 
              << std::setprecision(9) 
              << v00 
              << std::endl;
    std::cout << "Second time step:" 
              << std::fixed 
              << std::setprecision(9) 
              << v10 << "," << v11
              << std::endl;
    std::cout << "Third time step:" 
              << std::fixed 
              << std::setprecision(9) 
              << v20 << "," << v21 << "," << v22
              << std::endl;
    std::cout << "---------------------------------------------------------" 
              << std::endl;

    OptionValue optionValue;
    optionValue.price = v00;
    optionValue.delta = (v10-v11)/(s10-s11);
    optionValue.gamma = 2*((v20-v21)/(s20-s21)-(v21-v22)/(s21-s22))/(s20-s22);
    optionValue.theta = (v21-v00)/(2*dt);

    if (isAmerican && varReduction) {
        optionValue.price += varCorrection.price;
        optionValue.delta += varCorrection.delta;
        optionValue.gamma += varCorrection.gamma;
        optionValue.theta += varCorrection.theta;
    }

    return optionValue;
}

OptionValue BinomialTree::evaluateBBS(int N, bool isAmerican)
{
    return evaluate(N, isAmerican, true);
}
        
OptionValue BinomialTree::evaluateVarReduction(int N, bool useBS)
{
    return evaluate(N, true, useBS, true);
}

OptionValue BinomialTree::evaluateVarReductionBBS(int N)
{
    return evaluateVarReduction(N, true);
}



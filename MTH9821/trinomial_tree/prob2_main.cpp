#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <payoff.h>
#include <option_value.h>
#include <trinomial.h>

/**********************************
 *  Run Trinomial American option  *
 **********************************/
int main(int argc, char* argv[])
{
    // parameters
    double T = 1;
    double K = 40;
    double S = 41;
    double r = 0.03;
    double q = 0.01;
    double vol = 0.3;
    bool isAmerican = true;
    int N=8;

    // print precision
    int p=9;

    PutPayoff vanillaPut(K);
    TrinomialTree b1(vanillaPut, T, S, r, q, vol);
    double bsPrice = 3.970043924;
    double bsDelta = -0.388620465;
    double bsGamma = 0.032101174;
    double bsTheta = -1.990533559;

    std::cout << "Exact solution:,"
              << std::fixed
              << std::setprecision(p)
              << "," << bsPrice
              << "," << bsDelta
              << "," << bsGamma
              << "," << bsTheta
              << std::endl;

    int n=10;
    for (int k=0; k<N; k++) {
        std::cout << n << std::fixed << std::setprecision(p);
        OptionValue v = b1.evaluate(n,isAmerican);
        std::cout << "," << v.price
                  << "," << std::fabs(v.price-bsPrice)
                  << "," << n*std::fabs(v.price-bsPrice)
                  << "," << n*n*std::fabs(v.price-bsPrice)
                  << "," << v.delta << "," << std::fabs(v.delta-bsDelta)
                  << "," << v.gamma << "," << std::fabs(v.gamma-bsGamma)
                  << "," << v.theta << "," << std::fabs(v.theta-bsTheta); 

        // Trinomial Black-Scholes (TBS)
        OptionValue vTBS = b1.evaluateTBS(n,isAmerican);
        std::cout << "," << vTBS.price
                  << "," << std::fabs(vTBS.price-bsPrice)
                  << "," << n*std::fabs(vTBS.price-bsPrice)
                  << "," << n*n*std::fabs(vTBS.price-bsPrice)
                  << "," << vTBS.delta << "," << std::fabs(vTBS.delta-bsDelta)
                  << "," << vTBS.gamma << "," << std::fabs(vTBS.gamma-bsGamma)
                  << "," << vTBS.theta << "," << std::fabs(vTBS.theta-bsTheta); 

        // Trinomial Black-Scholes 
        // with Richardson extrapolation (TBSR)
        if (n/2>0) {
            OptionValue vBbsOld = b1.evaluateTBS(n/2,isAmerican);
            double vTBSR = 2*vTBS.price-vBbsOld.price;
            double tbsrDelta = 2*vTBS.delta-vBbsOld.delta;
            double tbsrGamma = 2*vTBS.gamma-vBbsOld.gamma;
            double tbsrTheta = 2*vTBS.theta-vBbsOld.theta;
            std::cout << "," << vTBSR
                  << "," << std::fabs(vTBSR-bsPrice)
                  << "," << n*std::fabs(vTBSR-bsPrice)
                  << "," << n*n*std::fabs(vTBSR-bsPrice)
                  << "," << tbsrDelta << "," << std::fabs(tbsrDelta-bsDelta)
                  << "," << tbsrGamma << "," << std::fabs(tbsrGamma-bsGamma)
                  << "," << tbsrTheta << "," << std::fabs(tbsrTheta-bsTheta); 
        }

        std::cout << std::endl;
        
        n *= 2; // double the tree size
    }

    return 0;
}


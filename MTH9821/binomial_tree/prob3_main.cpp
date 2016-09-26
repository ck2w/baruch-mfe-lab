#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <payoff.h>
#include <option_value.h>
#include <binomial.h>

/**********************************
 *  Run binomial American option  *
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
    BinomialTree b1(vanillaPut, T, S, r, q, vol);
    double bsPrice = 3.969988406;
    double bsDelta = -0.388621062;
    double bsGamma = 0.032101612;
    double bsTheta = -1.990567871;

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

        // Average of odd/even
        OptionValue vNext = b1.evaluate(n+1,isAmerican);
        double vAverage = (v.price+vNext.price)/2;
        double aveDelta = (v.delta+vNext.delta)/2;
        double aveGamma = (v.gamma+vNext.gamma)/2;
        double aveTheta = (v.theta+vNext.theta)/2;
        std::cout << "," << vAverage
                  << "," << std::fabs(vAverage-bsPrice)
                  << "," << n*std::fabs(vAverage-bsPrice)
                  << "," << n*n*std::fabs(vAverage-bsPrice)
                  << "," << aveDelta << "," << std::fabs(aveDelta-bsDelta)
                  << "," << aveGamma << "," << std::fabs(aveGamma-bsGamma)
                  << "," << aveTheta << "," << std::fabs(aveTheta-bsTheta); 

        // Binomial Black-Scholes (BBS)
        OptionValue vBBS = b1.evaluateBBS(n,isAmerican);
        std::cout << "," << vBBS.price
                  << "," << std::fabs(vBBS.price-bsPrice)
                  << "," << n*std::fabs(vBBS.price-bsPrice)
                  << "," << n*n*std::fabs(vBBS.price-bsPrice)
                  << "," << vBBS.delta << "," << std::fabs(vBBS.delta-bsDelta)
                  << "," << vBBS.gamma << "," << std::fabs(vBBS.gamma-bsGamma)
                  << "," << vBBS.theta << "," << std::fabs(vBBS.theta-bsTheta); 

        // Binomial Black-Scholes 
        // with Richardson extrapolation (BBSR)
        if (n/2>0) {
            OptionValue vBbsOld = b1.evaluateBBS(n/2,isAmerican);
            double vBBSR = 2*vBBS.price-vBbsOld.price;
            double bbsrDelta = 2*vBBS.delta-vBbsOld.delta;
            double bbsrGamma = 2*vBBS.gamma-vBbsOld.gamma;
            double bbsrTheta = 2*vBBS.theta-vBbsOld.theta;
            std::cout << "," << vBBSR
                  << "," << std::fabs(vBBSR-bsPrice)
                  << "," << n*std::fabs(vBBSR-bsPrice)
                  << "," << n*n*std::fabs(vBBSR-bsPrice)
                  << "," << bbsrDelta << "," << std::fabs(bbsrDelta-bsDelta)
                  << "," << bbsrGamma << "," << std::fabs(bbsrGamma-bsGamma)
                  << "," << bbsrTheta << "," << std::fabs(bbsrTheta-bsTheta); 
        }

        std::cout << std::endl;
        
        n *= 2; // double the tree size
    }

    return 0;
}


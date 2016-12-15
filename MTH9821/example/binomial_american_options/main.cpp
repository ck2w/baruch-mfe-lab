#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <payoff.h>
#include <option_value.h>
#include <binomial.h>

/**********************************
 *  On binomial American options  *
 **********************************/
int main(int argc, char* argv[])
{
    // parameters
    double T = 7.0/12.0;
    double K = 50;
    double S = 47;
    double r = 0.03;
    double q = 0.01;
    double vol = 0.25;
    PutPayoff vanillaPut(K);
    // maximum number of cases
    int N=10;
    // number steps in the first case
    int n=200;
    // do variance reduction?
    bool varReduction = true;
    // stop until convergence?
    bool stopUntilConverge = true;
    double tol=1e-4;
    // do average odd/even?
    bool doAverageOddEven=false;
    // do binomial Black-Scholes?
    bool doBBS=false;
    // do binomial Black-Scholes with Richardson?
    bool doBBSR=false;
    // have known exact solution?
    bool knowExact = false;
    double exactPrice = 0; 
    double exactDelta = 0;
    double exactGamma = 0; 
    double exactTheta = 0; 
    // print precision
    int p=9;

    BinomialTree b1(vanillaPut, T, S, r, q, vol);

    if (knowExact) {
        std::cout << "Exact solution:,"
                  << std::fixed
                  << std::setprecision(p)
                  << "," << exactPrice
                  << "," << exactDelta
                  << "," << exactGamma
                  << "," << exactTheta
                  << std::endl;
    }

    double vPrevious=0.0;
    for (int k=0; k<N; k++) {
        OptionValue v = b1.evaluate(n,true,false,varReduction);
        std::cout << n << std::fixed << std::setprecision(p);
        std::cout << "," << v.price;
        if (knowExact) {
            std::cout << "," << std::fabs(v.price-exactPrice)
                      << "," << n*std::fabs(v.price-exactPrice)
                      << "," << n*n*std::fabs(v.price-exactPrice);
        }
        std::cout << "," << v.delta;
        if (knowExact) { std::cout << "," << std::fabs(v.delta-exactDelta); } 
        std::cout << "," << v.gamma;
        if (knowExact) { std::cout << "," << std::fabs(v.gamma-exactGamma); }
        std::cout << "," << v.theta;
        if (knowExact) { std::cout << "," << std::fabs(v.theta-exactTheta); }

        // Average of odd/even
        if (doAverageOddEven) {
            OptionValue vNext = b1.evaluate(n+1,true,false,varReduction);
            double vAverage = (v.price+vNext.price)/2;
            double aveDelta = (v.delta+vNext.delta)/2;
            double aveGamma = (v.gamma+vNext.gamma)/2;
            double aveTheta = (v.theta+vNext.theta)/2;
            std::cout << "," << vAverage;
            if (knowExact) {
                std::cout << "," << std::fabs(vAverage-exactPrice)
                          << "," << n*std::fabs(vAverage-exactPrice)
                          << "," << n*n*std::fabs(vAverage-exactPrice);
            }
            std::cout << "," << aveDelta;
            if (knowExact) { std::cout << "," << std::fabs(aveDelta-exactDelta); } 
            std::cout << "," << aveGamma;
            if (knowExact) { std::cout << "," << std::fabs(aveGamma-exactGamma); }
            std::cout << "," << aveTheta;
            if (knowExact) { std::cout << "," << std::fabs(aveTheta-exactTheta); }
        }

        // Binomial Black-Scholes (BBS)
        if (doBBS) {
            OptionValue vBBS = b1.evaluate(n,true,true,varReduction);
            std::cout << "," << vBBS.price;
            if (knowExact) {
                std::cout << "," << std::fabs(vBBS.price-exactPrice)
                          << "," << n*std::fabs(vBBS.price-exactPrice)
                          << "," << n*n*std::fabs(vBBS.price-exactPrice);
            }
            std::cout << "," << vBBS.delta;
            if (knowExact) { std::cout << "," << std::fabs(vBBS.delta-exactDelta); } 
            std::cout << "," << vBBS.gamma;
            if (knowExact) { std::cout << "," << std::fabs(vBBS.gamma-exactGamma); }
            std::cout << "," << vBBS.theta;
            if (knowExact) { std::cout << "," << std::fabs(vBBS.theta-exactTheta); }

            // Binomial Black-Scholes 
            // with Richardson extrapolation (BBSR)
            if (doBBSR) {
                if (n/2>0) {
                    OptionValue vBbsOld = b1.evaluate(n/2,true,true,varReduction);
                    double vBBSR = 2*vBBS.price-vBbsOld.price;
                    double bbsrDelta = 2*vBBS.delta-vBbsOld.delta;
                    double bbsrGamma = 2*vBBS.gamma-vBbsOld.gamma;
                    double bbsrTheta = 2*vBBS.theta-vBbsOld.theta;
                    std::cout << "," << vBBSR;
                    if (knowExact) {
                        std::cout << "," << std::fabs(vBBSR-exactPrice)
                                  << "," << n*std::fabs(vBBSR-exactPrice)
                                  << "," << n*n*std::fabs(vBBSR-exactPrice);
                    }
                    std::cout << "," << bbsrDelta;
                    if (knowExact) { std::cout << "," << std::fabs(bbsrDelta-exactDelta); } 
                    std::cout << "," << bbsrGamma;
                    if (knowExact) { std::cout << "," << std::fabs(bbsrGamma-exactGamma); }
                    std::cout << "," << bbsrTheta;
                    if (knowExact) { std::cout << "," << std::fabs(bbsrTheta-exactTheta); }
                }
            }
        }

        std::cout << std::endl;

        if ( stopUntilConverge && std::fabs(v.price-vPrevious)<tol ) {
            std::cout << "|vCurent-vPrevious| = " 
                      << std::fabs(v.price-vPrevious)
                      << std::endl;
            break;
        }

        vPrevious = v.price;
        
        n *= 2; // double the tree size
    }

    return 0;
}


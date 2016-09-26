#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <payoff.h>
#include <binomial.h>

int main(int argc, char* argv[])
{
    // parameters
    double T = 1;
    double K = 40;
    double S = 41;
    double r = 0.03;
    double q = 0.01;
    double vol = 0.3;
    bool isAmerican = false;
    int N=100;

    // print precision
    int p=9;

    PutPayoff vanillaPut(K);
    BinomialTree b1(vanillaPut, T, S, r, q, vol);
    double vBS = b1.BlackScholesValue().price;

    // output column header
    std::cout << "N,Binomial,ABT,BBS,BBSR,Shanks";
    if (!isAmerican) {
        std::cout << ",Black-Scholes";
        std::cout << ",Binomial Error,ABT Error";
        std::cout << ",BBS Error,BBSR Error,Shanks Error";
    }
    std::cout << std::endl;

    double vOld = vanillaPut(K);
    double vOld2 = vOld;
    for (int n=1; n<=N; n++) {
        double v = b1.evaluate(n,isAmerican).price;
        // Average of odd/even
        double vAverage = (v+vOld)/2;
        // Binomial Black-Scholes (BBS)
        double vBBS = b1.evaluateBBS(n,isAmerican).price;

        // Record convergence
        std::cout << n 
                  << std::fixed
                  << std::setprecision(p)
                  << "," << v
                  << "," << vAverage
                  << "," << vBBS;

        // Binomial Black-Scholes 
        // with Richardson extrapolation (BBSR)
        double vBBSR = 0;
        if (n/2>0) {
            double vBbsOld = b1.evaluateBBS(n/2,isAmerican).price;
            vBBSR = 2*vBBS-vBbsOld;
            std::cout << "," << vBBSR;
        }

        // Shanks transformation
        double vShanks = 0;
        if (n>2) {
            vShanks = (v*vOld2-vOld*vOld)/(v-2*vOld+vOld2);
            std::cout << "," << vShanks;
        }

        // If European option, compare with Black-Scholes
        if (!isAmerican) {
            std::cout << "," << vBS
                      << "," << std::fabs(v-vBS)
                      << "," << std::fabs(vAverage-vBS)
                      << "," << std::fabs(vBBS-vBS);
            if (n/2>0) {std::cout << "," << std::fabs(vBBSR-vBS);}
            if (n>2) {std::cout << "," << std::fabs(vShanks-vBS);}
        }
        std::cout << std::endl;

        // Update old values
        vOld2 = vOld;
        vOld = v;
    }

    return 0;
}


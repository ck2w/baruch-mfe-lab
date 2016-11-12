#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <payoff.h>
#include <ran.h>
#include <mc.h>

int main(int argc, char* argv[])
{
    // parameters
    double T = 0.5;
    double K = 55;
    double S = 50;
    double r = 0.04;
    double q = 0.0;
    double vol = 0.3;
    int n=10;

    // print precision
    int p=9;

    PutPayoff vanillaPut(K);
    MonteCarlo mc(vanillaPut, T, S, r, q, vol, &bmnormal);
    double vBS = mc.BlackScholesValue().price;

    std::cout << "N,V,|V-vBS|,Std,Err,vBS" << std::endl;

    int nk=1;
    for (int k=0; k<n; k++) {
        int N = 10000*nk;
        OptionValue optionValue = mc.evaluateControlVariate(N);
        double v = optionValue.price;
        double std = optionValue.priceStd;
        double err = optionValue.priceErr;
        std::cout << N 
                  << std::fixed 
                  << std::setprecision(p)
                  << "," << v
                  << "," << std::fabs(v-vBS)
                  << "," << std
                  << "," << err
                  << "," << vBS
                  << std::endl;
        
        nk *= 2;
    }

    return 0;
}


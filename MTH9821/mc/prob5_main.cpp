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
    double K = 50;
    double S1 = 25;
    double S2 = 30;
    double r = 0.05;
    double q = 0.0;
    double vol1 = 0.3;
    double vol2 = 0.2;
    double rho = 0.35;
    int n=9;

    // print precision
    int p=9;

    CallPayoff basketCall(K);
    MonteCarlo mc(basketCall, T, S1, S2, r, q, q, vol1, vol2, rho, &bmnormal);

    std::cout << "N,V,Std,Err" << std::endl;

    int nk=1;
    for (int k=0; k<n; k++) {
        int N = 10000*nk;
        OptionValue optionValue = mc.evaluateBasket(N);
        double v = optionValue.price;
        double std = optionValue.priceStd;
        double err = optionValue.priceErr;
        std::cout << N 
                  << std::fixed 
                  << std::setprecision(p)
                  << "," << v
                  << "," << std
                  << "," << err
                  << std::endl;
        
        nk *= 2;
    }

    return 0;
}


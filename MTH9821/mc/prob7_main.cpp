#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <payoff.h>
#include <ran.h>
#include <mc_dividend.h>
#include <dividend.h>

int main(int argc, char* argv[])
{
    // parameters
    double T = 7.0/12.0;
    double K = 55.55;
    double S = 50;
    double r = 0.02;
    double q = 0.0;
    double vol = 0.3;
    int n=9;
    bool controlVariate = true;

    std::vector<DiscreteDividend> div;
    div.push_back({2.0/12.0, 0.50, true});
    div.push_back({4.0/12.0, 0.01, false});
    div.push_back({6.0/12.0, 0.75, true});

    // print precision
    int p=9;

    PutPayoff vanillaPut(K);
    MonteCarloDiscreteDividends mc(vanillaPut, T, S, r, q, div, vol, &bmnormal);

    int nk=1;
    for (int k=0; k<n; k++) {
        int N = 10000*nk;
        OptionValue optionValue = mc.evaluate(N, controlVariate);
        double v = optionValue.price;
        double std = optionValue.priceStd;
        double err = optionValue.priceErr;
        double delta = optionValue.delta;
        std::cout << N 
                  << std::fixed 
                  << std::setprecision(p)
                  << "," << v
                  << "," << delta
                  << "," << std
                  << "," << err
                  << std::endl;
        
        nk *= 2;
    }

    return 0;
}


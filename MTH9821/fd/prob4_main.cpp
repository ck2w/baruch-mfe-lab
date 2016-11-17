#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <payoff.h>
#include <fd.h>

int main(int argc, char* argv[])
{
    // parameters
    double T = 0.75;
    double S = 41;
    double K = 40;
    double r = 0.04;
    double q = 0.02;
    double vol = 0.35;

    int M = 4;
    double alphaTemp = 0.45;
    double omega = 1.2;

    FiniteDifferenceMethod fdm = CrankNicolsonBySOR;

    PutPayoff vanillaPut(K);
    FiniteDifference fd(vanillaPut, T, S, r, q, vol);
    M=4;   fd.evaluate(M, alphaTemp, fdm, omega);
    M=16;  fd.evaluate(M, alphaTemp, fdm, omega);
    M=64;  fd.evaluate(M, alphaTemp, fdm, omega);
    M=256; fd.evaluate(M, alphaTemp, fdm, omega);

    return 0;
}


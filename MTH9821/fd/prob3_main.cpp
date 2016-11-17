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

    PutPayoff vanillaPut(K);
    FiniteDifference fd(vanillaPut, T, S, r, q, vol);
    std::cout << "Euler forward" << std::endl;
    fd.evaluate(M, alphaTemp, EulerForward);
    std::cout << "Euler backward By LU " << std::endl;
    fd.evaluate(M, alphaTemp, EulerBackwardByLU);
    std::cout << "Euler backward By SOR " << std::endl;
    fd.evaluate(M, alphaTemp, EulerBackwardBySOR, omega);
    std::cout << "Crank Nicolson By LU " << std::endl;
    fd.evaluate(M, alphaTemp, CrankNicolsonByLU);
    std::cout << "Crank Nicolson By SOR " << std::endl;
    fd.evaluate(M, alphaTemp, CrankNicolsonBySOR, omega);

    return 0;
}


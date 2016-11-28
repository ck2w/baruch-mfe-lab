#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
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

    // Exact solution of American put
    double vExact = 4.083817051176386;

    int M = 4;
    double alphaTemp = 0.45;
    double omega = 1.2;

    //FiniteDifferenceMethod fdm = AmericanEulerForward; 
    //FiniteDifferenceMethod fdm = AmericanEulerBackwardByLU; 
    //FiniteDifferenceMethod fdm = AmericanEulerBackwardBySOR; 
    //FiniteDifferenceMethod fdm = AmericanCrankNicolsonByLU; 
    FiniteDifferenceMethod fdm = AmericanCrankNicolsonBySOR; 

    VanillaPutTerminalCondition f(T, S, K, r, q, vol);
    AmericanPutLeftBoundaryCondition gl(T, S, K, r, q, vol);
    VanillaPutRightBoundaryCondition gr(T, S, K, r, q, vol);
    AmericanPutEarlyExercisePremium prem(T, S, K, r, q, vol);

    FiniteDifference fd(T, S, K, r, q, vol, f, gl, gr, prem);
    M=4;   fd.evaluate(M, alphaTemp, fdm, omega, vExact);
    M=16;  fd.evaluate(M, alphaTemp, fdm, omega, vExact);
    M=64;  fd.evaluate(M, alphaTemp, fdm, omega, vExact);
    M=256; fd.evaluate(M, alphaTemp, fdm, omega, vExact);
    //M=512; fd.evaluate(M, alphaTemp, fdm, omega, vExact);

    return 0;
}


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

    int M = 4;
    double alphaTemp = 0.45;
    double omega = 1.2;

    FiniteDifferenceMethod fdm = CrankNicolsonBySOR;

    VanillaPutTerminalCondition f(T, S, K, r, q, vol);
    VanillaPutLeftBoundaryCondition gl(T, S, K, r, q, vol);
    VanillaPutRightBoundaryCondition gr(T, S, K, r, q, vol);

    FiniteDifference fd(T, S, K, r, q, vol, f, gl, gr);
    M=4;   fd.evaluate1(M, alphaTemp, fdm, omega);
    M=16;  fd.evaluate1(M, alphaTemp, fdm, omega);
    M=64;  fd.evaluate1(M, alphaTemp, fdm, omega);
    M=256; fd.evaluate1(M, alphaTemp, fdm, omega);

    return 0;
}


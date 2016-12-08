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

    VanillaPutTerminalCondition f(T, S, K, r, q, vol);
    VanillaPutLeftBoundaryCondition gl(T, S, K, r, q, vol);
    VanillaPutRightBoundaryCondition gr(T, S, K, r, q, vol);

    FiniteDifference fd(T, S, K, r, q, vol, f, gl, gr);
    std::cout << "Euler forward" << std::endl;
    fd.evaluate1(M, alphaTemp, EulerForward);
    std::cout << "Euler backward By LU " << std::endl;
    fd.evaluate1(M, alphaTemp, EulerBackwardByLU);
    std::cout << "Euler backward By SOR " << std::endl;
    fd.evaluate1(M, alphaTemp, EulerBackwardBySOR, omega);
    std::cout << "Crank Nicolson By LU " << std::endl;
    fd.evaluate1(M, alphaTemp, CrankNicolsonByLU);
    std::cout << "Crank Nicolson By SOR " << std::endl;
    fd.evaluate1(M, alphaTemp, CrankNicolsonBySOR, omega);

    return 0;
}


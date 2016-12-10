#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <fd.h>
#include <boundary_conditions.h>

int main(int argc, char* argv[])
{
    // parameters
    double T = 0.5;
    double S = 42;
    double K = 40;
    double B = 35;
    double r = 0.05;
    double q = 0.03;
    double vol = 0.3;

    int M = 16;
    double alphaTemp = 0.4;
    double omega = 1.2;

    FiniteDifferenceMethod fdm = EulerBackwardByLU; 

    VanillaCallTerminalCondition f(T, S, K, r, q, vol);
    VanillaCallLeftBoundaryCondition gl(T, S, K, r, q, vol);
    VanillaCallRightBoundaryCondition gr(T, S, K, r, q, vol);

    FiniteDifference fd(T, S, K, r, q, vol, f, gl, gr);
    fd.setToBarrierOption(B);

    // Computational domain
    double xCompute = std::log(S/K);
    double xLeft = std::log(B/K);
    double tFinal = 0.5*T*vol*vol;
    double dt = tFinal/M;
    double dxTemp = std::sqrt(dt/alphaTemp);
    int Nleft = static_cast<int>(std::floor((xCompute-xLeft)/dxTemp)); 
    double dx = (xCompute-xLeft)/Nleft;
    double alpha = dt/(dx*dx);
    double xRight = xCompute + (r-q-0.5*vol*vol)*T+3*vol*std::sqrt(T);
    int Nright = static_cast<int>(std::ceil((xRight-xCompute)/dx));
    xRight = xCompute + Nright*dx;
    int N = Nleft+Nright;

    std::cout << "alpha = " << alpha << std::endl;
    std::cout << "xl = " << xLeft << ", xr = " << xRight << std::endl;
    std::cout << "Nl = " << Nleft << ", Nr = " << Nright << std::endl;
    std::cout << "M = " << M << ", N = " << N << std::endl;
    std::cout << "dx = " << dx << ", dt = " << dt << std::endl;

    fd.setDomain(xLeft, xRight, tFinal, 0);
    fd.discretizeDomain(M,N);

    fd.evaluate(fdm, omega);

    return 0;
}


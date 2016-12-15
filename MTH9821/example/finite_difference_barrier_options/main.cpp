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
    double T = 5.0/12.0;
    double S = 52;
    double K = 55;
    double r = 0.02;
    double q = 0.005;
    double vol = 0.25;
    
    bool doubleBarrier=true;
    double B1 = 45;
    double B2 = 62;

    int M = 4;
    double alphaTemp = 0.5;
    double omega = 1.2;

    FiniteDifferenceMethod fdm = EulerForward; 

    VanillaPutTerminalCondition f(T, S, K, r, q, vol);
    VanillaPutLeftBoundaryCondition gl(T, S, K, r, q, vol);
    VanillaPutRightBoundaryCondition gr(T, S, K, r, q, vol);

    // If it's a single barrier option, do the following:
    //FiniteDifference fd(T, S, K, r, q, vol, f, gl, gr);

    // If it's double barrier option, do the following:
    // simple replacing the right boundary by zero, too.
    FiniteDifference fd(T, S, K, r, q, vol, f, gr, gr);

    // setToBarrierOption if it's a single barrier
    // setToDoubleBarrierOption does not do anything
    // because there is no exact solution for double barrier.
    fd.setToDoubleBarrierOption(B1,B2);

    // Computational domain
    double xCompute = std::log(S/K);

    // Terminal time
    double tFinal = 0.5*T*vol*vol;
    double dt = tFinal/M;
    double dxTemp = std::sqrt(dt/alphaTemp);

    // Left boundary
    // Single barrier
    /*
    double xLeft = std::log(B1/K);
    int Nleft = static_cast<int>(std::floor((xCompute-xLeft)/dxTemp)); 
    double dx = (xCompute-xLeft)/Nleft;
    double alpha = dt/(dx*dx);
    double xRight = xCompute + (r-q-0.5*vol*vol)*T+3*vol*std::sqrt(T);
    int Nright = static_cast<int>(std::ceil((xRight-xCompute)/dx));
    xRight = xCompute + Nright*dx;
    int N = Nleft+Nright;
    */

    // Double barrier
    double xLeft = std::log(B1/K);
    double xRight = std::log(B2/K);
    int N = static_cast<int>(std::floor((xRight-xLeft)/dxTemp));
    double dx = (xRight-xLeft)/N;
    double alpha = dt/(dx*dx);


    std::cout << std::fixed << std::setprecision(9) 
              << "x Compute = " << xCompute << std::endl;
    std::cout << std::fixed << std::setprecision(9) 
              << "alpha = " << alpha << std::endl;
    std::cout << std::fixed << std::setprecision(9) 
              << "xl = " << xLeft << ", xr = " << xRight << std::endl;
    std::cout << std::fixed << std::setprecision(9) 
              << "t final = " << tFinal << std::endl;
    //std::cout << "Nl = " << Nleft << ", Nr = " << Nright << std::endl;
    std::cout << "M = " << M << ", N = " << N << std::endl;
    std::cout << std::fixed << std::setprecision(9) 
              << "dx = " << dx << ", dt = " << dt << std::endl;

    fd.setDomain(xLeft, xRight, tFinal, 0);
    fd.discretizeDomain(M,N);

    fd.evaluate(fdm, omega);

    return 0;
}


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <fd.h>

double f(double x) 
{
    return std::exp(x);
}

double gl(double t)
{
    return std::exp(t-2);
}

double gr(double t)
{
    return std::exp(t+2);
}

double uExact(double x)
{
    return std::exp(x+1);
}

int main(int argc, char* argv[])
{
    double xl = -2;
    double xr = 2;
    double tf = 1;

    HeatPDE h(xl,xr,tf,&f,&gl,&gr);

    int M,N,dM,dN;

    // Case I
    M=8;
    N=8;
    dM=1;
    dN=1;
    std::vector<double> u_8_4(N+1,0);
    std::cout << " --- Forward Euler --- " << std::endl;
    h.fdSolveForwardEuler(M, N, &u_8_4, dM, dN);
    std::cout << " --- Backward Euler --- " << std::endl;
    h.fdSolveBackwardEulerByLU(M, N, &u_8_4, dM, dN);
    std::cout << " --- Backward Euler by SOR --- " << std::endl;
    h.fdSolveBackwardEulerBySOR(M, N, 1.2, &u_8_4, dM, dN);
    std::cout << " --- Crank-Nicolson --- " << std::endl;
    h.fdSolveCrankNicolsonByLU(M, N, &u_8_4, dM, dN);
    std::cout << " --- Crank-Nicolson by SOR --- " << std::endl;
    h.fdSolveCrankNicolsonBySOR(M, N, 1.2, &u_8_4, dM, dN);
  
    std::cout << " --- exact --- " << std::endl; 
    double dx = (xr-xl)/N;
    for (int n=0; n<N+1; n++) {
        double x = xl + n*dx;
        if ( n%dN == 0 ) {
            std::cout << "," << uExact(x);
        }
    }
    std::cout << std::endl;

    return 0;
}


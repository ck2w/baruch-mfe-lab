#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <utility>
#include <heat_pde.h>

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

void computeApproximationError( const std::vector<double> & u, 
                                double (*uExact)(double),
                                double xl, double xr )
{
    int N = u.size()-1;
    double dx = (xr-xl)/N;

    double maxPointwiseError = 0.0;
    double rmsError = 0.0;

    for (int i=1; i<N; i++) {

        double x = xl + i*dx;
        double uExactValue = (*uExact)(x);

        double diff = std::fabs(u[i]-uExactValue);
        if ( diff > maxPointwiseError ) { maxPointwiseError = diff; }
        rmsError += (diff*diff)/(uExactValue*uExactValue);
    }

    rmsError /= double(N+1);
    rmsError = std::sqrt(rmsError);

    std::cout << std::setprecision(9) 
              << maxPointwiseError << ","
              << rmsError
              << std::endl;
}

int main(int argc, char* argv[])
{
    double xl = -2;
    double xr = 2;
    double tf = 1;

    HeatPDE h(xl,xr,tf,&f,&gl,&gr);

    int M,N;

    // Available Methods:
    // - fdSolve-ForwardEuler
    // - fdSolve-BackwardEulerByLU
    // - fdSolve-BackwardEulerBySOR
    // - fdSolve-CrankNicolsonByLU
    // - fdSolve-CrankNicolsonBySOR
    
    std::vector<double> u;
    std::cout << " --- CrankNicolsonBySOR --- " << std::endl;
    
    M=8; N=4; u.resize(N+1,0);
    h.fdSolveCrankNicolsonBySOR(M, N, 1.2, &u); 
    computeApproximationError(u, &uExact, xl, xr);
    
    M=32; N=8; u.resize(N+1,0); 
    h.fdSolveCrankNicolsonBySOR(M, N, 1.2, &u);
    computeApproximationError(u, &uExact, xl, xr);
    
    M=128; N=16; u.resize(N+1,0);
    h.fdSolveCrankNicolsonBySOR(M, N, 1.2, &u);
    computeApproximationError(u, &uExact, xl, xr);
    
    M=512; N=32; u.resize(N+1,0);
    h.fdSolveCrankNicolsonBySOR(M, N, 1.2, &u);
    computeApproximationError(u, &uExact, xl, xr);
    
    M=8; N=8; u.resize(N+1,0);
    h.fdSolveCrankNicolsonBySOR(M, N, 1.2, &u);
    computeApproximationError(u, &uExact, xl, xr);
    
    M=32; N=16; u.resize(N+1,0);
    h.fdSolveCrankNicolsonBySOR(M, N, 1.2, &u);
    computeApproximationError(u, &uExact, xl, xr);
    
    M=128; N=32; u.resize(N+1,0);
    h.fdSolveCrankNicolsonBySOR(M, N, 1.2, &u);
    computeApproximationError(u, &uExact, xl, xr);
    
    M=512; N=64; u.resize(N+1,0);
    h.fdSolveCrankNicolsonBySOR(M, N, 1.2, &u);
    computeApproximationError(u, &uExact, xl, xr);
    
    M=8; N=16; u.resize(N+1,0);
    h.fdSolveCrankNicolsonBySOR(M, N, 1.2, &u);
    computeApproximationError(u, &uExact, xl, xr);
    
    M=32; N=32; u.resize(N+1,0);
    h.fdSolveCrankNicolsonBySOR(M, N, 1.2, &u);
    computeApproximationError(u, &uExact, xl, xr);
    
    M=128; N=64; u.resize(N+1,0);
    h.fdSolveCrankNicolsonBySOR(M, N, 1.2, &u);
    computeApproximationError(u, &uExact, xl, xr);
    
    M=512; N=128; u.resize(N+1,0);
    h.fdSolveCrankNicolsonBySOR(M, N, 1.2, &u);
    computeApproximationError(u, &uExact, xl, xr);

    return 0;
}


#include <fd.h>
#include <evaluator.h>
#include <payoff.h>
#include <option_value.h>
#include <black_scholes.h>
#include <heat_pde.h>
#include <tuple>
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>

FiniteDifference::FiniteDifference( const Payoff & payoff, 
                                    double expiry,
                                    double spot,
                                    double rate,
                                    double div,
                                    double vol )
                                  : d_payoff(&payoff),
                                    d_expiry(expiry),
                                    d_spot(spot),
                                    d_rate(rate),
                                    d_div(div),
                                    d_vol(vol)
{
    // Calculate the Black-Scholes value
    // in case of a vanilla call or put.
    d_strike = payoff.strike();
    std::tuple<OptionValue,OptionValue> res
        = BlackScholes(d_expiry,d_strike,d_spot,d_rate,d_div,d_vol);
    if ((*d_payoff)(0) > 0) { d_BlackScholes = std::get<1>(res);}
    else { d_BlackScholes = std::get<0>(res);}

    // Set up computational domain 
    d_tf = 0.5*expiry*vol*vol;
    double xmid = std::log(spot/d_strike) + (rate-div-0.5*vol*vol)*expiry;
    d_xl = xmid - 3*vol*std::sqrt(expiry);
    d_xr = xmid + 3*vol*std::sqrt(expiry);

    // Set up boudary condition
    d_f = VanillaPutTerminalCondition(expiry, spot, d_strike, rate, div, vol);
    d_gl = VanillaPutLeftBoundaryCondition
               (expiry, spot, d_strike, rate, div, vol);
    d_gr = VanillaPutRightBoundaryCondition
               (expiry, spot, d_strike, rate, div, vol);

    // Set up PDE
    d_h = HeatPDE(d_xl, d_xr, d_tf, d_f, d_gl, d_gr);
}

OptionValue FiniteDifference::BlackScholesValue() const
{
    return d_BlackScholes;
}

void FiniteDifference::evaluate( int M, double alphaTemp,
                                 FiniteDifferenceMethod fdm,
                                 double omega )
{
    // Discretize the computational domain
    double dt = d_tf/M;
    double dxTemp = std::sqrt(dt/alphaTemp);
    int N = static_cast<int>(std::floor((d_xr - d_xl)/dxTemp));
    double dx = (d_xr-d_xl)/N;

    // Discretize and solve the PDE
    std::vector<double> u(N+1,0);
    int dM=0;
    int dN=0;
    switch (fdm)
    {
        case EulerForward: 
            d_h.fdSolveForwardEuler(M, N, &u, dM, dN);
            break;
        case EulerBackwardByLU:
            d_h.fdSolveBackwardEulerByLU(M, N, &u, dM, dN);
            break;
        case EulerBackwardBySOR:
            d_h.fdSolveBackwardEulerBySOR(M, N, omega, &u, dM, dN);
            break;
        case CrankNicolsonByLU:
            d_h.fdSolveCrankNicolsonByLU(M, N, &u, dM, dN);
            break;
        case CrankNicolsonBySOR:
            d_h.fdSolveCrankNicolsonBySOR(M, N, omega, &u, dM, dN);
            break;
    }

    // Use the numerical solution to find option values and greeks
    double c = (d_rate-d_div)/(d_vol*d_vol);
    double a = c-0.5;
    double b = (c+0.5)*(c+0.5) + 2*d_div/(d_vol*d_vol);
    
    double xCompute = std::log(d_spot/d_strike);
    int spotIndex = (xCompute-d_xl)/dx;
    double xi = d_xl + spotIndex*dx;
    double xj = d_xl + (spotIndex+1)*dx;
    double si = d_strike*std::exp(xi);
    double sj = d_strike*std::exp(xj);
    double ui = u[spotIndex];
    double uj = u[spotIndex+1];
    double vi = std::exp(-a*xi-b*d_tf)*ui;
    double vj = std::exp(-a*xj-b*d_tf)*uj;
    //std::cout << "xi = " << xi << ", xj = " << xj << std::endl; 
    //std::cout << "si = " << si << ", sj = " << sj << std::endl; 
    //std::cout << "ui = " << ui << ", uj = " << uj << std::endl; 
    //std::cout << "vi = " << vi << ", vj = " << vj << std::endl; 

    // Approximation I: interpolation of v
    double vApproximate1 = ((sj-d_spot)*vi + (d_spot-si)*vj)/(sj-si);

    // Approximation II: interpolation of u
    double uApproximate = ((xj-xCompute)*ui + (xCompute-xi)*uj)/(xj-xi);
    double vApproximate2 = std::exp(-a*xCompute-b*d_tf)*uApproximate;

    // Point-wise Error
    double bsPrice = d_BlackScholes.price;
    //double bsDelta = d_BlackScholes.delta;
    //double bsGamma = d_BlackScholes.gamma;
    //double bsTheta = d_BlackScholes.theta;
    double error1 = std::fabs(vApproximate1-bsPrice);
    double error2 = std::fabs(vApproximate2-bsPrice);

    // Root-Mean-Squared (RMS) Error
    double error3 = 0;
    int count = 0;
    for (int n=0; n<N+1; n++) {
        double x = d_xl + n*dx;
        double s = d_strike*std::exp(x);
        // Black-Scholes value
        OptionValue vBlackScholes;
        std::tuple<OptionValue,OptionValue> res
            = BlackScholes(d_expiry,d_strike,s,d_rate,d_div,d_vol);
        if ((*d_payoff)(0) > 0) { vBlackScholes = std::get<1>(res);}
        else { vBlackScholes = std::get<0>(res);}

        double bsPrice = vBlackScholes.price;
        if (bsPrice > 0.00001*s) {
            double vApproximate = std::exp(-a*x-b*d_tf)*u[n];
            double diff = (vApproximate - bsPrice)/bsPrice;
            double diff2 = diff*diff;
            error3 += diff2;
            count++;
        }
    }

    error3 /= count;
    error3 = std::sqrt(error3);

    // Compute Greeks I: Delta and Gamma
    double xii = d_xl + (spotIndex-1)*dx;
    double xjj = d_xl + (spotIndex+2)*dx;
    double sii = d_strike*std::exp(xii);
    double sjj = d_strike*std::exp(xjj);
    double uii = u[spotIndex-1];
    double ujj = u[spotIndex+2];
    double vii = std::exp(-a*xii-b*d_tf)*uii;
    double vjj = std::exp(-a*xjj-b*d_tf)*ujj;

    double delta = (vj-vi)/(sj-si);
    double deltaLeft = (vi-vii)/(si-sii);
    double deltaRight = (vjj-vj)/(sjj-sj);
    double gamma = 2*(deltaRight-deltaLeft)/(sjj+sj-si-sii);

    // Compute Greeks II: Theta
    std::vector<double> uOld = d_h.getPrevSolution();

    double dT = 2*dt/(d_vol*d_vol);
    double uOldi = uOld[spotIndex];
    double uOldj = uOld[spotIndex+1];
    double vOldi = std::exp(-a*xi-b*(d_tf-dt))*uOldi;
    double vOldj = std::exp(-a*xj-b*(d_tf-dt))*uOldj;

    double vOldApproximate1 = ((sj-d_spot)*vOldi + (d_spot-si)*vOldj)/(sj-si);
    // minus sign: backward B-S equation converted to forward heat equation
    double theta = -(vApproximate1-vOldApproximate1)/dT;

    std::cout << std::fixed
              << std::setprecision(9)
              << error1 << ",,"
              << error2 << ",,"
              << error3 << ",,"
              << delta << "," //<< bsDelta << ","
              << gamma << "," //<< bsGamma << ","
              << theta << "," //<< bsTheta
              << std::endl;
    
} 

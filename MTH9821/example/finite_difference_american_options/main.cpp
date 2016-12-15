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
    double S = 27;
    double K = 30;
    double r = 0.03;
    double q = 0.01;
    double vol = 0.35;

    // Exact solution of American put
    double vExact = 4.083817051176386;

    double alphaTemp = 0.5;
    double omega = 1.2;
    bool useVarReduction = true;

    FiniteDifferenceMethod fdm = AmericanEulerForward; 
    //FiniteDifferenceMethod fdm = AmericanEulerBackwardByLU; 
    //FiniteDifferenceMethod fdm = AmericanEulerBackwardBySOR; 
    //FiniteDifferenceMethod fdm = AmericanCrankNicolsonByLU; 
    
    // Important!!! Check iterative_solve/linear_iterate.cpp
    // for stopping conditions.
    //FiniteDifferenceMethod fdm = AmericanCrankNicolsonBySOR; 

    AmericanPutTerminalCondition f(T, S, K, r, q, vol);
    AmericanPutLeftBoundaryCondition gl(T, S, K, r, q, vol);
    AmericanPutRightBoundaryCondition gr(T, S, K, r, q, vol);
    AmericanPutEarlyExercisePremium prem(T, S, K, r, q, vol);

    FiniteDifference fd(T, S, K, r, q, vol, f, gl, gr, prem);
    fd.setDefaultDomain();

    int M=4;
    int n=4; // number of cases
    for (int i=0; i<n; i++) {
        fd.discretizeDomainByTimeStepsAndAlphaTemp(M, alphaTemp);
        
        //std::cout << "M =, " << fd.getTimeSteps() << std::endl;
        //std::cout << "N =, " << fd.getIntervals() << std::endl;
        //std::cout << "Alpha =, " << fd.getAlpha() << std::endl;
        //std::cout << "xl =, " << std::fixed << std::setprecision(9) << fd.getXl() << std::endl;
        //std::cout << "xr =, " << std::fixed << std::setprecision(9) << fd.getXr() << std::endl;
        //std::cout << "tf =, " << std::fixed << std::setprecision(9) << fd.getTf() << std::endl;
        //std::cout << "ti =, " << std::fixed << std::setprecision(9) << fd.getTi() << std::endl;
        //std::cout << "dt =, " << std::fixed << std::setprecision(9) 
        //          << (fd.getTf()-fd.getTi())/fd.getTimeSteps() << std::endl;
        //std::cout << "dx =, " << std::fixed << std::setprecision(9) 
        //          << (fd.getXr()-fd.getXl())/fd.getIntervals() << std::endl;

        std::cout << "Parameters:"
                  << std::fixed
                  << std::setprecision(9)
                  << "," << fd.getTimeSteps()
                  << "," << fd.getAlpha()
                  << "," << fd.getIntervals()
                  << "," << fd.getXl()
                  << "," << fd.getXr()
                  << "," << std::log(S/K)
                  << "," << fd.getTf()
                  << "," << (fd.getTf()-fd.getTi())/fd.getTimeSteps()
                  << "," << (fd.getXr()-fd.getXl())/fd.getIntervals()
                  << std::endl;

        fd.evaluate(fdm, omega, useVarReduction, vExact);
        M *= 4;
    }

    return 0;
}


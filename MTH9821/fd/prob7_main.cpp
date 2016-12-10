#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <fd.h>
#include <dividend.h>
#include <boundary_conditions.h>

int main(int argc, char* argv[])
{
    // parameters
    double T = 1.0;
    double S = 52;
    double K = 50;
    double r = 0.03;
    double q = 0.0;
    double vol = 0.3;
    DiscreteDividend div = {5.0/12.0, 0.02, false};

    double tDiv = div.t;
    double vDiv = div.v;

    int M = 4;
    double alpha = 0.4;
    double omega = 1.2;

    FiniteDifferenceMethod fdm = EulerForward; 

    VanillaCallTerminalCondition f(T, S, K, r, q, vol);
    VanillaCallLeftBoundaryCondition gl(T, S, K, r, q, vol);
    VanillaCallRightBoundaryCondition gr(T, S, K, r, q, vol);

    FiniteDifference fd1(T-tDiv, S, K, r, q, vol, f, gl, gr);
    fd1.setExpandedDomain(M, alpha, T, vDiv);
    fd1.discretizeDomainByTimeStepsAndAlphaFixed(M, alpha);

    //std::cout << "M1 =, " << fd1.getTimeSteps() << std::endl;
    //std::cout << "N1 =, " << fd1.getIntervals() << std::endl;
    //std::cout << "xl1 =, " << std::fixed << std::setprecision(9) << fd1.getXl() << std::endl;
    //std::cout << "xr1 =, " << std::fixed << std::setprecision(9) << fd1.getXr() << std::endl;
    //std::cout << "tf1 =, " << std::fixed << std::setprecision(9) << fd1.getTf() << std::endl;
    //std::cout << "ti1 =, " << std::fixed << std::setprecision(9) << fd1.getTi() << std::endl;
    //std::cout << "dt1 =, " << std::fixed << std::setprecision(9) 
    //          << (fd1.getTf()-fd1.getTi())/fd1.getTimeSteps() << std::endl;
    //std::cout << "dx1 =, " << std::fixed << std::setprecision(9) 
    //          << (fd1.getXr()-fd1.getXl())/fd1.getIntervals() << std::endl;

    std::vector<double> u1 = fd1.evaluate(fdm, omega);

    //for (int i=0; i<u1.size(); i++) { std::cout << u1[i] << ", "; }
    //std::cout << std::endl;

    double xl = fd1.getXl() - std::log(1-div.v);
    double xr = fd1.getXr() - std::log(1-div.v);
    int N1 = fd1.getIntervals();
   
    double tau_f = 0.5*T*vol*vol;
    double tau_i = 0.5*(T-tDiv)*vol*vol;
    FiniteDifference fd2(tDiv, S, K, r, q, vol, f, gl, gr);
    fd2.setDomain(xl, xr, tau_f, tau_i);
    fd2.discretizeDomainByIntervalsAndAlphaTemp(N1, alpha);

    //std::cout << "M2 =, " << fd2.getTimeSteps() << std::endl;
    //std::cout << "Alpha2 =, " << fd2.getAlpha() << std::endl;
    //std::cout << "N2 =, " << fd1.getIntervals() << std::endl;
    //std::cout << "xl2 =, " << std::fixed << std::setprecision(9) << fd2.getXl() << std::endl;
    //std::cout << "xr2 =, " << std::fixed << std::setprecision(9) << fd2.getXr() << std::endl;
    //std::cout << "tf2 =, " << std::fixed << std::setprecision(9) << fd2.getTf() << std::endl;
    //std::cout << "ti2 =, " << std::fixed << std::setprecision(9) << fd2.getTi() << std::endl;
    //std::cout << "dt2 =, " << std::fixed << std::setprecision(9) 
    //          << (fd2.getTf()-fd2.getTi())/fd2.getTimeSteps() << std::endl;
    //std::cout << "dx2 =, " << std::fixed << std::setprecision(9) 
    //          << (fd2.getXr()-fd2.getXl())/fd2.getIntervals() << std::endl;

    //std::cout << "tauDiv =, " << tau_i << std::endl;

    //std::cout << std::fixed << std::setprecision(9)
    //          << fd1.getTimeSteps() << ","
    //          << fd2.getTimeSteps() << ","
    //          << fd2.getAlpha() << ","
    //          << fd1.getIntervals() << ","
    //          << fd1.getXl() << ","
    //          << fd1.getXr() << ","
    //          << fd2.getXl() << ","
    //          << fd2.getXr() << ","
    //          << tau_i << ","
    //          << (fd1.getTf()-fd1.getTi())/fd1.getTimeSteps() << ","
    //          << (fd2.getTf()-fd2.getTi())/fd2.getTimeSteps() << ","
    //          << (fd1.getXr()-fd1.getXl())/fd1.getIntervals()
    //          << std::endl;

    fd2.overrideTerminalCondition(u1);
    fd2.evaluate(fdm, omega);

    return 0;
}


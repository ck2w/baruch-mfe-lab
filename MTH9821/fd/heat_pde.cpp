#include <heat_pde.h>
#include <updater.h>
#include <evaluator.h>
#include <cassert>
#include <iostream>
#include <iomanip>

HeatPDE::HeatPDE( double xl,
                  double xr,
                  double tf,
                  double ti,
                  double (*f)(double),
                  double (*gl)(double),
                  double (*gr)(double) )
                : d_xl(xl), d_xr(xr), d_tf(tf), d_ti(ti)
{
    d_f = new FunctionEvaluator(f);
    d_gl = new FunctionEvaluator(gl);
    d_gr = new FunctionEvaluator(gr);
    d_functionalTerminal = true;
}

HeatPDE::HeatPDE( double xl,
                  double xr,
                  double tf,
                  double ti,
                  const std::vector<double> & u0,
                  const Evaluator & gl,
                  const Evaluator & gr,
                  const Evaluator & prem )
                : d_xl(xl), d_xr(xr), d_tf(tf), d_ti(ti),
                  d_gl(&gl), d_gr(&gr), d_prem(&prem)
{
    d_functionalTerminal = false;
    d_u0 = u0;
}

HeatPDE::HeatPDE( double xl,
                  double xr,
                  double tf,
                  double ti,
                  const Evaluator & f,
                  const Evaluator & gl,
                  const Evaluator & gr,
                  const Evaluator & prem )
                : d_xl(xl), d_xr(xr), d_tf(tf), d_ti(ti),
                  d_f(&f), d_gl(&gl), d_gr(&gr), d_prem(&prem)
{
    d_functionalTerminal = true;
}

void HeatPDE::print(double t, const std::vector<double> & u, int step)
{
    int uLen = u.size();
    std::cout << std::fixed << std::setprecision(9) << t;
    for (int n=0; n<uLen; n++) {
        if ( n%step == 0 ) {
            std::cout << "," << u[n];
        }
    }
    std::cout << std::endl;
}

std::vector<double> HeatPDE::earlyExercisePremiumAtGivemTime(int N, double t,
                                                             bool xReversed)
{
    std::vector<double> earlyExercisePremium(N+1);
    double dx = (d_xr-d_xl)/N;
    for (int n=0; n<N+1; n++) {
        double x = d_xl + n*dx;
        if (xReversed) { x = d_xr - n*dx; }
        earlyExercisePremium[n] = (*d_prem)(x,t);
    }

    return earlyExercisePremium;
}

int HeatPDE::earlyExerciseBoundaryAtGivenTime(const std::vector<double> & u,
                                              const std::vector<double> & p,
                                              bool xReversed)
{
    int uLen = u.size();
    if (!xReversed) {
        for (int n=0; n<uLen-1; n++) {
            if (u[n] == p[n] && u[n+1] > p[n+1]) {
                return n;
            }
        }
    
        return uLen;
    }
    else {
        for (int n=uLen; n>0;) {
            n--;
            if (u[n] == p[n] && u[n-1] > p[n-1]) {
                return n;
            }
        }

        return 0;
    }
}
        
void HeatPDE::fdSolve( int M, int N, std::vector<double>* u, Updater* up, 
                       bool isAmerican, bool xReversed, int dM, int dN )
{
    assert(u->size() == (unsigned int)(N+1));
    d_earlyExerciseBoundary.resize(0);
    bool doPrint = false;
    if ( dM>0 && dN>0 ) { doPrint = true; }

    // grid parameters 
    double dt = (d_tf-d_ti)/M;
    double dx = (d_xr-d_xl)/N;
    double c = dt/(dx*dx);

    up->config(c,N);

    std::vector<double> uNew(N+1,0);
    // initial condition
    for (int n=0; n<N+1; n++) {
        double x = d_xl + n*dx;
        if (xReversed) { 
            x = d_xr - n*dx; 
            std::reverse(std::begin(*u), std::end(*u));
        }

        // initialize if functional terminal condition is given
        if (d_functionalTerminal) { 
            (*u)[n] = (*d_f)(x); 
        }

        // print x-coordinates
        if (doPrint && n%dN == 0) { std::cout << "," << x; }
    }

    //std::cout << "==============" << std::endl;
    //std::cout << "u initial: " << std::endl;
    //for (int i=0; i<(*u).size(); i++) { std::cout << std::setprecision(3) << (*u)[i] << ", "; }
    //std::cout << std::endl;
    //std::cout << "--------------" << std::endl;

    d_earlyExerciseBoundary.push_back(0);

    if (doPrint) { 
        std::cout << std::endl;
        print(0, *u, dN); 
    }

    for (int m=0; m<M; m++) {
        // boundary conditions at the next time step
        double t = d_ti + (m+1)*dt;
        uNew[0] = (*d_gl)(t);
        uNew[N] = (*d_gr)(t);
        // Some methods like Backward Euler by LU.
        // works in reversed x-coordinates
        if (xReversed) {
            uNew[N] = (*d_gl)(t);
            uNew[0] = (*d_gr)(t);
        }

        // early exercise boundary at the current time step
        std::vector<double> prem;
        if (isAmerican) {
            prem = earlyExercisePremiumAtGivemTime(N, t, xReversed);
        }

        // time evolution
        up->update((*u), &uNew, prem);

        // store the solution at T-dt
        if ( m == M-1 ) { d_uOld = (*u); }
        // write back to output
        for (int n=0; n<=N; n++) { (*u)[n] = uNew[n]; }

        if (isAmerican) {
            int nb = earlyExerciseBoundaryAtGivenTime((*u), prem, xReversed);
            if (!xReversed) {
                d_earlyExerciseBoundary.push_back(d_xl+nb*dx);
            }
            else {
                d_earlyExerciseBoundary.push_back(d_xr-nb*dx);
            }
        }
    
        if ( doPrint && m%dM == 0 ) { print(t, *u, dN); }
    }

    if (xReversed) {
        std::reverse(std::begin(*u), std::end(*u));
    }
    
    //std::cout << "u final: " << std::endl;
    //for (int i=0; i<(*u).size(); i++) { std::cout << std::setprecision(3) << (*u)[i] << ", "; }
    //std::cout << std::endl;
    //std::cout << "=============" << std::endl;

}
        
void HeatPDE::fdSolveForwardEuler(int M, int N, 
                                  std::vector<double>* u,
                                  int dM, int dN)
{
    ForwardEulerUpdater up;
    fdSolve(M,N,u,&up,false,false,dM,dN);
}

void HeatPDE::fdSolveBackwardEulerByLU(int M, int N, 
                                       std::vector<double>* u,
                                       int dM, int dN)
{
    BackwardEulerUpdater up;
    fdSolve(M,N,u,&up,false,false,dM,dN);
}

void HeatPDE::fdSolveBackwardEulerBySOR(int M, int N, double w, 
                                        std::vector<double>* u,
                                        int dM, int dN)
{
    BackwardEulerSorUpdater up(w);
    fdSolve(M,N,u,&up,false,false,dM,dN);
}

void HeatPDE::fdSolveCrankNicolsonByLU(int M, int N, 
                                       std::vector<double>* u,
                                       int dM, int dN)
{
    CrankNicolsonUpdater up;
    fdSolve(M,N,u,&up,false,false,dM,dN);
}

void HeatPDE::fdSolveCrankNicolsonBySOR(int M, int N, double w, 
                                        std::vector<double>* u,
                                        int dM, int dN)
{
    CrankNicolsonSorUpdater up(w);
    fdSolve(M,N,u,&up,false,false,dM,dN);
}

void HeatPDE::fdSolveAmericanForwardEuler(int M, int N, 
                                          std::vector<double>* u,
                                          int dM, int dN)
{
    ForwardEulerUpdater up;
    fdSolve(M,N,u,&up,true,false,dM,dN);
}

void HeatPDE::fdSolveAmericanBackwardEulerByLU(int M, int N, 
                                               std::vector<double>* u,
                                               int dM, int dN)
{
    BackwardEulerUpdater up;
    fdSolve(M,N,u,&up,true,true,dM,dN);
}

void HeatPDE::fdSolveAmericanBackwardEulerBySOR(int M, int N, double w,
                                                std::vector<double>* u,
                                                int dM, int dN)
{
    BackwardEulerSorUpdater up(w);
    fdSolve(M,N,u,&up,true,false,dM,dN);
}

void HeatPDE::fdSolveAmericanCrankNicolsonByLU(int M, int N,
                                               std::vector<double>* u,
                                               int dM, int dN)
{
    CrankNicolsonUpdater up;
    fdSolve(M,N,u,&up,true,true,dM,dN);
}

void HeatPDE::fdSolveAmericanCrankNicolsonBySOR(int M, int N, double w, 
                                                std::vector<double>* u,
                                                int dM, int dN)
{
    CrankNicolsonSorUpdater up(w);
    fdSolve(M,N,u,&up,true,false,dM,dN);
}

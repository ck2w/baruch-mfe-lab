#include <fd.h>
#include <updater.h>
#include <cassert>
#include <iostream>
#include <iomanip>

void HeatPDE::print(double t, const std::vector<double> & u, int step)
{
    int uLen = u.size();
    std::cout << std::setprecision(9) << t;
    for (int n=0; n<uLen; n++) {
        if ( n%step == 0 ) {
            std::cout << "," << u[n];
        }
    }
    std::cout << std::endl;
}
        
void HeatPDE::fdSolve( int M, int N, std::vector<double>* u, Updater* up,
                       int dM, int dN )
{
    assert(u->size() == (unsigned int)(N+1));
    bool doPrint = false;
    if ( dM>0 && dN>0 ) { doPrint = true; }

    // grid parameters 
    double dt = d_tf/M;
    double dx = (d_xr-d_xl)/N;
    double c = dt/(dx*dx);

    up->config(c,N);

    std::vector<double> uNew(N+1,0);
    // initial condition
    for (int n=0; n<N+1; n++) {
        double x = d_xl + n*dx;
        (*u)[n] = d_f(x);
    }

    if (doPrint) { print(0, *u, dN); }

    for (int m=0; m<M; m++) {
        // boundary conditions
        double t = (m+1)*dt;
        uNew[0] = d_gl(t);
        uNew[N] = d_gr(t);

        // time evolution
        up->update((*u), &uNew);

        // write back
        for (int n=0; n<=N; n++) { (*u)[n] = uNew[n]; }
    
        if ( doPrint && m%dM == 0 ) { print(t, *u, dN); }
    }

}
        
void HeatPDE::fdSolveForwardEuler(int M, int N, 
                                  std::vector<double>* u,
                                  int dM, int dN)
{
    ForwardEulerUpdater up;
    fdSolve(M,N,u,&up,dM,dN);
}

void HeatPDE::fdSolveBackwardEulerByLU(int M, int N, 
                                       std::vector<double>* u,
                                       int dM, int dN)
{
    BackwardEulerUpdater up;
    fdSolve(M,N,u,&up,dM,dN);
}

void HeatPDE::fdSolveBackwardEulerBySOR(int M, int N, double w, 
                                        std::vector<double>* u,
                                        int dM, int dN)
{
    BackwardEulerSorUpdater up(w);
    fdSolve(M,N,u,&up,dM,dN);
}

void HeatPDE::fdSolveCrankNicolsonByLU(int M, int N, 
                                       std::vector<double>* u,
                                       int dM, int dN)
{
    CrankNicolsonUpdater up;
    fdSolve(M,N,u,&up,dM,dN);
}

void HeatPDE::fdSolveCrankNicolsonBySOR(int M, int N, double w, 
                                        std::vector<double>* u,
                                        int dM, int dN)
{
    CrankNicolsonSorUpdater up(w);
    fdSolve(M,N,u,&up,dM,dN);
}


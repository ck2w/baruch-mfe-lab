#ifndef HEAT_PDE_H
#define HEAT_PDE_H 
#include <vector>
#include <Eigen/Dense>
#include <tuple>
#include <updater.h>
#include <evaluator.h>

/********************************************
 *  Finite Different Method: Heat Equation  *
 *                                          *
 *  Input:                                  *
 *  - xl: x-left                            *
 *  - xr: x-right                           *
 *  - tf: t-final                           *
 *  which define the computatioal domain    *
 *  - f(x)                                  *
 *  - gl(t): g-left                         *
 *  - gr(t): g-right                        *
 *  which define the boundary conditions    *
 *  - M: # of intervals on t-axis           *
 *  - N: # of intervals on x-axis           *
 *                                          *
 ********************************************/
class HeatPDE
{
    public:

        HeatPDE( double xl,
                 double xr,
                 double tf,
                 double (*f)(double),
                 double (*gl)(double),
                 double (*gr)(double) );
        
        HeatPDE( double xl,
                 double xr,
                 double tf,
                 const Evaluator & f,
                 const Evaluator & gl,
                 const Evaluator & gf );

        void fdSolveForwardEuler(int M, int N, 
                                 std::vector<double>* u,
                                 int dM=0, int dN=0);
        void fdSolveBackwardEulerByLU(int M, int N, 
                                      std::vector<double>* u,
                                      int dM=0, int dN=0);
        void fdSolveBackwardEulerBySOR(int M, int N, double w, 
                                       std::vector<double>* u,
                                       int dM=0, int dN=0);
        void fdSolveCrankNicolsonByLU(int M, int N, 
                                      std::vector<double>* u,
                                      int dM=0, int dN=0);
        void fdSolveCrankNicolsonBySOR(int M, int N, double w, 
                                       std::vector<double>* u,
                                       int dM=0, int dN=0);

    private:

        double d_xl;
        double d_xr;
        double d_tf;

        const Evaluator * d_f;
        const Evaluator * d_gl;
        const Evaluator * d_gr;

        void print(double t, const std::vector<double> & u, int step);

        void fdSolve( int M, int N, std::vector<double>* u, Updater* up,
                      int dM=0, int dN=0 );
};


#endif /* HEAT_PDE_H */

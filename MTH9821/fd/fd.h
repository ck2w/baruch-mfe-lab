#ifndef FD_H
#define FD_H 
#include <vector>
#include <Eigen/Dense>
#include <tuple>
#include <updater.h>

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
                 double (*gr)(double) )
               : d_xl(xl), d_xr(xr), d_tf(tf),
                 d_f(f), d_gl(gl), d_gr(gr)
        {}

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

        double (*d_f)(double);
        double (*d_gl)(double);
        double (*d_gr)(double);

        void print(double t, const std::vector<double> & u, int step);

        void fdSolve( int M, int N, std::vector<double>* u, Updater* up,
                      int dM=0, int dN=0 );
};


#endif /* FD_H */

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

        HeatPDE() {}

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
                 const Evaluator & gf,
                 const Evaluator & prem=Evaluator() );

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
        void fdSolveAmericanForwardEuler(int M, int N, 
                                         std::vector<double>* u,
                                         int dM=0, int dN=0);
        void fdSolveAmericanBackwardEulerByLU(int M, int N, 
                                              std::vector<double>* u,
                                              int dM=0, int dN=0);
        void fdSolveAmericanBackwardEulerBySOR(int M, int N, double w,
                                               std::vector<double>* u,
                                               int dM=0, int dN=0);
        void fdSolveAmericanCrankNicolsonByLU(int M, int N,
                                              std::vector<double>* u,
                                              int dM=0, int dN=0);
        void fdSolveAmericanCrankNicolsonBySOR(int M, int N, double w, 
                                               std::vector<double>* u,
                                               int dM=0, int dN=0);

        std::vector<double> getPrevSolution(bool xReversed) const
        {
            std::vector<double> uOld = d_uOld;
            if (xReversed) { std::reverse(std::begin(uOld), std::end(uOld)); }
            return uOld;
        }

        const std::vector<double> & getEarlyExerciseBoundary() const
        {
            return d_earlyExerciseBoundary;
        }

    private:

        double d_xl;
        double d_xr;
        double d_tf;

        const Evaluator * d_f;
        const Evaluator * d_gl;
        const Evaluator * d_gr;
        const Evaluator * d_prem;

        std::vector<double> d_uOld;
        std::vector<double> d_earlyExerciseBoundary;

        void print(double t, const std::vector<double> & u, int step);

        void fdSolve( int M, int N, std::vector<double>* u, Updater* up,
                      bool isAmerican=false, bool xReversed=false, 
                      int dM=0, int dN=0 );

        std::vector<double> earlyExercisePremiumAtGivemTime(int N, double t,
                                                            bool xReversed);

        int earlyExerciseBoundaryAtGivenTime(const std::vector<double> & u,
                                             const std::vector<double> & p,
                                             bool xReversed);
};


#endif /* HEAT_PDE_H */

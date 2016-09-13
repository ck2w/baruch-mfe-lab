#include <laplace2d.h>
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <tuple>
#include <iostream>

std::tuple<int, double, double> laplace2d(const Eigen::ArrayXXd & a,
                                          const Eigen::ArrayXXd & b,
                                          const Eigen::ArrayXXd & c,
                                          const Eigen::ArrayXXd & d,
                                          const Eigen::ArrayXXd & e,
                                          const Eigen::ArrayXXd & f,
                                          Eigen::ArrayXXd * u,
                                          double omega,
                                          double tol,
                                          bool chebyshev)
{
    assert(omega<2);  // required by convergence
    assert(omega>=1); // enforce over-relaxation
    // estimated jacobi spectral radius
    double rjac = std::sqrt(1.0-(2.0/omega-1.0)*(2.0/omega-1.0));

    const int MAXITS=1000000;

    double anorm=0.0;
    double anormf=0.0;
    double resid=0.0;

    int imax=a.rows();
    int jmax=a.cols();
    assert(imax == b.rows());
    assert(imax == c.rows());
    assert(imax == d.rows());
    assert(imax == e.rows());
    assert(imax == f.rows());
    assert(imax == u->rows());
    assert(jmax == b.cols());
    assert(jmax == c.cols());
    assert(jmax == d.cols());
    assert(jmax == e.cols());
    assert(jmax == f.cols());
    assert(jmax == u->cols());

    // Compute the initial norm of residual and terminate 
    // iteration when norm has been reduced by a tolerance
    // factor. Assume the initial values of u are all zero.
    for (int i=1; i<imax-1; i++) {
        for (int j=1; j<jmax-1; j++) {
            anormf += f(i,j)*f(i,j);
        }
    }
    anormf = std::sqrt(anormf);
    if (anormf < 1e-16) {anormf = 1.0;}

    for (int n=0; n<MAXITS; n++) {
        anorm = 0.0;
        int isw=1;
        //red-black checkerboard ordering
        for (int ipass=0; ipass<2; ipass++) {
            int jsw=isw;
            for (int i=1; i<imax-1; i++) {
                for (int j=jsw; j<jmax-1; j+=2) {
                    resid = a(i,j)*(*u)(i+1,j)
                          + b(i,j)*(*u)(i-1,j)
                          + c(i,j)*(*u)(i,j+1)
                          + d(i,j)*(*u)(i,j-1)
                          + e(i,j)*(*u)(i,j)
                          - f(i,j);
                    anorm += resid*resid;
                    (*u)(i,j) -= omega*resid/e(i,j);
                }
                jsw=3-jsw; //switch between 1 and 2
            }
            isw=3-isw; //switch between 1 and 2

            // chebyshev acceleration
            if (chebyshev) {
                if (n==0 && ipass==0) {
                    omega = 1.0/(1.0-0.5*rjac*rjac);
                }
                else {
                    omega = 1.0/(1.0-0.25*rjac*rjac*omega);
                }
            }
        }

        anorm = std::sqrt(anorm);
        if (anorm < tol*anormf) {
            return std::make_tuple(n,anorm, anorm/anormf);
        }
    }

    bool MAXITS_EXCEEDED = false;
    assert(!MAXITS_EXCEEDED);

    return std::make_tuple(MAXITS,anormf,1);
}



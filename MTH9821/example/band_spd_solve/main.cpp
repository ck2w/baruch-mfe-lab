#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <Eigen/Dense>
#include <linear_solve.h>

int main(int argc, char* argv[])
{
    double w=1.0;
    int n=1024;

    int N=n*n;
    double dx = w/(n+1);
    double dx2 = dx*dx;

    Eigen::ArrayXXd mat = Eigen::ArrayXXd::Zero(N,n+1);
    //diagonal
    for (int i=0; i<N; i++) {mat(i,0) = 4;}
    //upper diagonal
    for (int i=0; i<N-1; i++) 
    {
        if ((i+1)%n != 0) {
            mat(i,1) = -1;
        }
    }

    for (int i=0; i<N-n; i++) {
        mat(i,n) = -1;
    }

    Eigen::MatrixXd vec = Eigen::VectorXd::Zero(N);
    for (int i=0; i<n; i++) {
        double x = (i+1)*dx;
        for (int j=0; j<n; j++) {
            int index = i*n+j;
            double y = (j+1)*dx;
            // source
            double f = (x*x+y*y-2)*std::sin(x)*std::sin(y)
                       -2*x*std::cos(x)*std::sin(y)
                       -2*y*std::sin(x)*std::cos(y);
            vec(index) = f*dx2;
            
            // boundary
            if (j == 0) {
            }
            if (j == n-1) {
                double ub = 0.5*(x*x+1)*std::sin(x)*std::sin(w);
                vec(index) += ub;
            }
            if (i == 0) {
            }
            if (i == n-1) {
                double ub = 0.5*(y*y+1)*std::sin(y)*std::sin(w);
                vec(index) += ub;
            }
        }
    }

    Eigen::VectorXd x = banded_spd_solve(mat, n, vec);
    
    /*************** exact solution *******************/
    Eigen::VectorXd u_exact = Eigen::VectorXd::Zero(N);
    for (int i=0; i<n; i++) {
        double x = (i+1)*dx;
        for (int j=0; j<n; j++) {
            int index = i*n+j;
            double y = (j+1)*dx;
            u_exact(index) = 0.5*(x*x+y*y)*std::sin(x)*std::sin(y);
        }
    }
    double error=0.0;
    for (int i=0; i<N; i++) {
        double diff = std::fabs(x(i)-u_exact(i));
        if (diff>error) {
            error = diff;
        }
    }
    std::cout << "Approximation error:  " << error << std::endl;
    return 0;
}


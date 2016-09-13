#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <cmath>
#include <tuple>
#include <load.h>
#include <laplace2d.h>
#include <Eigen/Dense>

int main(int argc, char* argv[])
{
    /******************* grid **********************/
    double wx = 1.0;      // x-width
    double wy = 1.0;      // y-width 
    int nx = 64;           // # of interior x-points
    int ny = 64;           // # of interior y-points
    /**************** relaxation *******************/
    double omega = 1.0;   // relaxation parameter
    double tol   = 1e-6;  // tolerance factor
    bool cheby   = false; // chebyshev acceleration
    /***********************************************/
    
    double dx  = wx/(nx+1); // x-pixel size
    double dy  = wy/(ny+1); // y-pixel size
    double dx2 = dx*dx;     // x-pixel size square
    double dy2 = dy*dy;     // y-pixel size square
    double r2  = dx2/dy2;   // aspect ratio square
    double eij = 2*(1+r2);  // diagonal term

    int Nx=nx+2;            // # of x-points
    int Ny=ny+2;            // # of y-points
    Eigen::ArrayXXd a = Eigen::ArrayXXd::Zero(Nx,Ny);
    Eigen::ArrayXXd b = Eigen::ArrayXXd::Zero(Nx,Ny);
    Eigen::ArrayXXd c = Eigen::ArrayXXd::Zero(Nx,Ny);
    Eigen::ArrayXXd d = Eigen::ArrayXXd::Zero(Nx,Ny);
    Eigen::ArrayXXd e = Eigen::ArrayXXd::Zero(Nx,Ny);
    Eigen::ArrayXXd f = Eigen::ArrayXXd::Zero(Nx,Ny);
    Eigen::ArrayXXd u = Eigen::ArrayXXd::Zero(Nx,Ny);

    for (int i=0; i<Nx; i++) {
        for (int j=0; j<Ny; j++) {
            // laplacian
            a(i,j) = -1;
            b(i,j) = -1;
            c(i,j) = -r2;
            d(i,j) = -r2;
            e(i,j) = eij;
        }
    }
    
    /****************** source ************************/
    for (int i=0; i<Nx; i++) {
        double x=i*dx;
        for (int j=0; j<Ny; j++) {
            double y=j*dy;
            f(i,j) = (x*x+y*y-2)*std::sin(x)*std::sin(y)
                     -2*x*std::cos(x)*std::sin(y)
                     -2*y*std::sin(x)*std::cos(y);
            f(i,j) *= dx2;
        }
    }

    /****************** boundary **********************/
    for (int i=0; i<Nx; i++) {
        double x=i*dx;
        u(i,0)    = 0;
        u(i,Ny-1) = 0.5*(x*x+1)*std::sin(x)*std::sin(wx); 
    }
    
    for (int j=0; j<Ny; j++) {
        double y=j*dy;
        u(0,j)    = 0;
        u(Ny-1,j) = 0.5*(y*y+1)*std::sin(y)*std::sin(wy); 
    }
    /**************************************************/

    std::tuple<int, double, double> res 
        = laplace2d(a,b,c,d,e,f,&u,omega,tol,cheby);

    std::cout << "X-dimension:          " << wx << std::endl;
    std::cout << "Y-dimension:          " << wy << std::endl;
    std::cout << "X-grid size:          " << Nx << std::endl;
    std::cout << "Y-grid size:          " << Ny << std::endl;
    std::cout << "x-pixel:              " << dx << std::endl;
    std::cout << "y-pixel:              " << dy << std::endl;
    std::cout << "Relaxation parameter: " << omega << std::endl;
    std::cout << "Tolerance factor:     " << tol << std::endl;
    std::cout << "Chebyshev accelerate: " << (cheby ? "YES":"NO") << std::endl;
    std::cout << "Number of iterations: " << std::get<0>(res) << std::endl;
    std::cout << "Absolute residual:    " << std::get<1>(res) << std::endl;
    std::cout << "Relative residual:    " << std::get<2>(res) << std::endl;
    //std::cout << u << std::endl;
    
    /*************** exact solution *******************/
    Eigen::ArrayXXd u_exact = Eigen::ArrayXXd::Zero(Nx,Ny);
    for (int i=0; i<Nx; i++) {
        double x=i*dx;
        for (int j=0; j<Ny; j++) {
            double y=j*dy;
            u_exact(i,j) = 0.5*(x*x+y*y)*std::sin(x)*std::sin(y);
        }
    }
    //std::cout << "------ exact solution ------" << std::endl;
    //std::cout << u_exact << std::endl;
    /**************************************************/
    double error=0.0;
    for (int i=1; i<Nx-1; i++) {
        for (int j=1; j<Ny-1; j++) {
            double diff = u(i,j)-u_exact(i,j);
            error += diff*diff*dx*dy;
        }
    }
    error = std::sqrt(error);
    std::cout << "Approximation error:  " << error << std::endl;
   
    return 0;
}


#include <updater.h>
#include <lu.h>
#include <triangular_solve.h>
#include <sor.h>
#include <vector>
#include <tuple>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

void ForwardEulerUpdater::config( double c, int N )
{
    Updater::setAlpha(c);
} 

void ForwardEulerUpdater::update( const std::vector<double> & uOld,
                                  std::vector<double>* uNew,
                                  const std::vector<double> & prem ) const
{
    bool isAmerican = (prem.size()>0);
    int N = uOld.size()-1;
    double c = getAlpha();
    for (int n=1; n<N; n++) {
        (*uNew)[n] = c*uOld[n+1]+(1-2*c)*uOld[n]+c*uOld[n-1];
        // overriding with early exercise premium
        if (isAmerican && prem[n]>(*uNew)[n]) { (*uNew)[n] = prem[n]; }
    }
}

void BackwardEulerUpdater::config( double c, int N )
{
    Updater::setAlpha(c);
    // Construct A-matrix LU-decomposition
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N-1,N-1); 
    for (int i=0; i<N-1; i++) {
        A(i,i) = 1+2*c;
        if (i>=1) A(i,i-1) = -c;
        if (i+1<=N-2) A(i,i+1) = -c;
    }
    
    std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> res = lu_no_pivoting(A);
    d_L = std::get<0>(res);
    d_U = std::get<1>(res);
} 

void BackwardEulerUpdater::update( const std::vector<double> & uOld,
                                   std::vector<double>* uNew,
                                   const std::vector<double> & prem ) const
{
    int N = uOld.size()-1;
    double c = getAlpha();

    // Construct b-vector
    Eigen::VectorXd b = Eigen::VectorXd::Zero(N-1);
    for (int i=0; i<N-1; i++) { b(i) = uOld[i+1]; }

    b(0) += c*(*uNew)[0];
    b(N-2) += c*(*uNew)[N];
    
    // Construct early exercise premium vector
    Eigen::VectorXd p = Eigen::VectorXd::Zero(0);
    if ( prem.size()>0 ) {
        p = Eigen::VectorXd::Zero(N-1);
        for (int i=0; i<N-1; i++) { p(i) = prem[i+1]; }
    }
   
    // update 
    Eigen::VectorXd y = forward_subst(d_L, b);
    Eigen::VectorXd x = backward_subst(d_U, y, p);
    
    for (int i=0; i<N-1; i++) { (*uNew)[i+1] = x(i); }
}

void BackwardEulerSorUpdater::config( double c, int N )
{
    Updater::setAlpha(c);
    // Construct A-matrix LU-decomposition
    d_A = Eigen::MatrixXd::Zero(N-1,N-1); 
    for (int i=0; i<N-1; i++) {
        d_A(i,i) = 1+2*c;
        if (i>=1) d_A(i,i-1) = -c;
        if (i+1<=N-2) d_A(i,i+1) = -c;
    }
} 

void BackwardEulerSorUpdater::update( const std::vector<double> & uOld,
                                      std::vector<double>* uNew,
                                      const std::vector<double> & prem ) const
{
    int N = uOld.size()-1;
    double c = getAlpha();

    // Construct b-vector
    Eigen::VectorXd b = Eigen::VectorXd::Zero(N-1);
    for (int i=0; i<N-1; i++) { b(i) = uOld[i+1]; }
    b(0) += c*(*uNew)[0];
    b(N-2) += c*(*uNew)[N];

    // Use the state vector at the previous time step as initial estimate
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(N-1);
    for (int i=0; i<N-1; i++) { x0(i) = uOld[i+1]; }
    
    // Construct early exercise premium vector
    Eigen::VectorXd p = Eigen::VectorXd::Zero(0);
    if ( prem.size()>0 ) {
        p = Eigen::VectorXd::Zero(N-1);
        for (int i=0; i<N-1; i++) { 
            p(i) = prem[i+1]; 
            // If early exercise premium is given, the initial 
            // estimate x0 is set to early exercise premium.
            // This is a reasonable choice because:
            // (1) At the 1st time step, i.e., m=0,
            //     early_exercise_premium = K*exp(ax)*max(1-exp(x),0),
            //     which is exactly the same as the terminal condition;
            // (2) At the successive time steps, i.e., m>0,
            //     early_exercise_premium serves as a good approximation
            //     of the value of state vectors. 
            x0(i) = prem[i+1];
        }
    }

    // update 
    std::tuple<Eigen::VectorXd, int> res = sor(d_omega, d_A, b, d_tol, x0, p);
    Eigen::VectorXd x = std::get<0>(res);
    
    for (int i=0; i<N-1; i++) { (*uNew)[i+1] = x(i); }
}

void CrankNicolsonUpdater::config( double c, int N )
{
    Updater::setAlpha(c);
    // Construct A-matrix LU-decomposition
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N-1,N-1); 
    for (int i=0; i<N-1; i++) {
        A(i,i) = 1+c;
        if (i>=1) A(i,i-1) = -0.5*c;
        if (i+1<=N-2) A(i,i+1) = -0.5*c;
    }
    
    std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> res = lu_no_pivoting(A);
    d_L = std::get<0>(res);
    d_U = std::get<1>(res);
    
    // Construct B-matrix
    d_B = Eigen::MatrixXd::Zero(N-1,N-1); 
    for (int i=0; i<N-1; i++) {
        d_B(i,i) = 1-c;
        if (i>=1) d_B(i,i-1) = 0.5*c;
        if (i+1<=N-2) d_B(i,i+1) = 0.5*c;
    }
} 

void CrankNicolsonUpdater::update( const std::vector<double> & uOld,
                                   std::vector<double>* uNew,
                                   const std::vector<double> & prem ) const
{
    int N = uOld.size()-1;
    double c = getAlpha();

    // uOld restricted in the interior region
    Eigen::VectorXd u = Eigen::VectorXd::Zero(N-1);
    for (int i=0; i<N-1; i++) { u(i) = uOld[i+1]; }

    // Construct b-vector
    Eigen::VectorXd b = d_B*u;

    b(0) += 0.5*c*((*uNew)[0]+uOld[0]);
    b(N-2) += 0.5*c*((*uNew)[N]+uOld[N]);
    
    // Construct early exercise premium vector
    Eigen::VectorXd p = Eigen::VectorXd::Zero(0);
    if ( prem.size()>0 ) {
        p = Eigen::VectorXd::Zero(N-1);
        for (int i=0; i<N-1; i++) { p(i) = prem[i+1]; }
    }
    
    // update 
    Eigen::VectorXd y = forward_subst(d_L, b);
    Eigen::VectorXd x = backward_subst(d_U, y, p);
    
    for (int i=0; i<N-1; i++) { (*uNew)[i+1] = x(i); }
}

void CrankNicolsonSorUpdater::config( double c, int N )
{
    Updater::setAlpha(c);
    // Construct A-matrix LU-decomposition
    d_A = Eigen::MatrixXd::Zero(N-1,N-1); 
    for (int i=0; i<N-1; i++) {
        d_A(i,i) = 1+c;
        if (i>=1) d_A(i,i-1) = -0.5*c;
        if (i+1<=N-2) d_A(i,i+1) = -0.5*c;
    }
    
    // Construct B-matrix
    d_B = Eigen::MatrixXd::Zero(N-1,N-1); 
    for (int i=0; i<N-1; i++) {
        d_B(i,i) = 1-c;
        if (i>=1) d_B(i,i-1) = 0.5*c;
        if (i+1<=N-2) d_B(i,i+1) = 0.5*c;
    }
} 

void CrankNicolsonSorUpdater::update( const std::vector<double> & uOld,
                                      std::vector<double>* uNew,
                                      const std::vector<double> & prem ) const
{
    int N = uOld.size()-1;
    double c = getAlpha();

    // uOld restricted in the interior region
    Eigen::VectorXd u = Eigen::VectorXd::Zero(N-1);
    for (int i=0; i<N-1; i++) { u(i) = uOld[i+1]; }

    // Construct b-vector
    Eigen::VectorXd b = d_B*u;

    b(0) += 0.5*c*((*uNew)[0]+uOld[0]);
    b(N-2) += 0.5*c*((*uNew)[N]+uOld[N]);

    // Use the state vector at the previous time step as initial estimate
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(N-1);
    for (int i=0; i<N-1; i++) { x0(i) = uOld[i+1]; }
    
    // Construct early exercise premium vector
    Eigen::VectorXd p = Eigen::VectorXd::Zero(0);
    if ( prem.size()>0 ) {
        p = Eigen::VectorXd::Zero(N-1);
        for (int i=0; i<N-1; i++) { 
            p(i) = prem[i+1]; 
            // If early exercise premium is given, the initial 
            // estimate x0 is set to early exercise premium.
            // This is a reasonable choice because:
            // (1) At the 1st time step, i.e., m=0,
            //     early_exercise_premium = K*exp(ax)*max(1-exp(x),0),
            //     which is exactly the same as the terminal condition;
            // (2) At the successive time steps, i.e., m>0,
            //     early_exercise_premium serves as a good approximation
            //     of the value of state vectors. 
            x0(i) = prem[i+1];
        }
    }
    
    // update 
    std::tuple<Eigen::VectorXd, int> res = sor(d_omega, d_A, b, d_tol, x0, p);
    Eigen::VectorXd x = std::get<0>(res);
    
    for (int i=0; i<N-1; i++) { (*uNew)[i+1] = x(i); }
}

#include <ran.h>
#include <norminv.h>
#include <iostream>
#include <cmath>
#include <random>

double std_normal(int* idum)
{
    static std::default_random_engine generator(*idum);
    std::normal_distribution<double> distribution(0.0,1.0);
    return distribution(generator);
}

double ran(int* idum)
{
    int a = 39373;
    int k = 2147483647;

    int q = int(k/a);
    int r = k%a;
    int t = int((*idum)/q);
    (*idum) = a*((*idum)-q*t)-t*r;
    if ((*idum)<0) {(*idum) += k;}
    return double((*idum))/k;

    /*
    int a = 39373;
    int k = 2147483647;
    double ki = 4.656612875245797e-10;
    // Schrange's algorithm
    // k = a*q + r
    int q = 54542;
    int r = 1481;

    int j = (*idum)/q;
    (*idum) = a*((*idum)-j*q)-r*j;
    if ((*idum)<0) { (*idum) += k; }

    return (*idum)*ki;
    */
}

double normal(int* idum)
{
    double u = ran(idum);
    return norminv(u);
}

double arnormal(int* idum)
{
    double u1=0;
    double u2=0;
    double u3=0;
    double x=0;

    do {
        u1 = ran(idum);
        u2 = ran(idum);
        u3 = ran(idum);
        x = -std::log(u1);
    } while ( u2 > std::exp(-0.5*(x-1)*(x-1)) );

    if ( u3 <= 0.5 ) {
        x = -x;
    }

    return x;
}

double bmnormal(int* idum)
{
    static bool z2Available = false;
    static double z1 = 0;
    static double z2 = 0;

    if (z2Available) {
        z2Available = false;
        return z2;
    }
    else {
        double u1=0;
        double u2=0;
        double x=0;

        do {
            u1 = ran(idum);
            u2 = ran(idum);
            u1 = 2*u1-1;
            u2 = 2*u2-1;
            x = u1*u1+u2*u2;
        } while ( x>1 );

        double y = std::sqrt(-2*std::log(x)/x);
        z1 = u1*y;
        z2 = u2*y;

        z2Available = true;

        return z1;
    }
}


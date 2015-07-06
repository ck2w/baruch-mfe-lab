#include <cumulative_normal.h>
#include <cmath>

double Cumulative_normal::approximate_evaluate(double x)
{
    double z = std::fabs(x);
    double y = 1.0/(1.0+a0*z);
    double y2 = y*y;
    double y3 = y2*y;
    double y4 = y2*y2;
    double y5 = y3*y2;
    double v = 1.0
        -std::exp(-x*x/2)*(a1*y+a2*y2+a3*y3+a4*y4+a5*y5)/std::sqrt(2*M_PI);

    if (x>0)
    {
        return v;
    }
    else
    {
        return (1-v);
    }
}

double Cumulative_normal::exact_evaluate(double x)
{
    return 0.5 * erfc(-x/std::sqrt(2));
}

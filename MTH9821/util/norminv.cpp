#include <norminv.h>
#include <limits>
#include <cmath>

double norminv(double u)
{
    if (u<=0 || u>=1) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double a0 = 2.50662823884;
    double a1 = -18.61500062529;
    double a2 = 41.39119773534;
    double a3 = -25.44106049637;
    double b0 = -8.47351093090;
    double b1 = 23.08336743743;
    double b2 = -21.06224101826;
    double b3 = 3.13082909833;
    double c0 = 0.3374754822726147;
    double c1 = 0.9761690190917186;
    double c2 = 0.1607979714918209;
    double c3 = 0.0276438810333863;
    double c4 = 0.0038405729373609;
    double c5 = 0.0003951896511919;
    double c6 = 0.0000321767881768;
    double c7 = 0.0000002888167364;
    double c8 = 0.0000003960315187;

    double x = 0;
    double y = u-0.5;
    if ( std::fabs(y) < 0.42 ) {
        double r = y*y;
        x = y*(((a3*r+a2)*r+a1)*r+a0)/((((b3*r+b2)*r+b1)*r+b0)*r+1);
    }
    else {
        double r = u;
        if (y>0) {
            r = 1-u;
        }
        r = std::log(-std::log(r));
        x = c0+r*(c1+r*(c2+r*(c3+r*(c4+r*(c5+r*(c6+r*(c7+r*c8)))))));
        if (y>0) {
            x = -x;
        }
    }

    return x;
}


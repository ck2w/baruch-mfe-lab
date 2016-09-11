#include <band_mult.h>
#include <Eigen/Dense>
#include <cassert>
#include <algorithm>

Eigen::VectorXd band_mult(const Eigen::ArrayXXd & a, int m1, int m2,
                          const Eigen::VectorXd & x)
{
    int width=m1+m2+1;
    int n = a.rows();
    Eigen::VectorXd y = Eigen::VectorXd::Zero(n);

    for (int i=0; i<n; i++) {
        int k=i-m1;
        int tmploop = std::min(width,n-k);
        for (int j=std::max(0,-k); j<tmploop; j++) {
            y(i) += a(i,j)*x(j+k);
        }
    }

    return y;
}

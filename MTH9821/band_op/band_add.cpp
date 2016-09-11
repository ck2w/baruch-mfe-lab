#include <band_add.h>
#include <Eigen/Dense>
#include <cassert>
#include <algorithm>

Eigen::ArrayXXd band_add(const Eigen::ArrayXXd & a, int m1, int m2,
                         const Eigen::ArrayXXd & b, int n1, int n2)
{
    int n=a.rows();
    assert(n == b.rows());

    int l1 = std::max(m1,n1);
    int l2 = std::max(m2,n2);
    int l=l1+l2+1;

    Eigen::ArrayXXd c = Eigen::ArrayXXd::Zero(n,l);

    c.col(l1) = a.col(m1) + b.col(n1);

    for (int i=1; i<=l1; i++) {
        if (i <= m1) {
            c.col(l1-i) += a.col(m1-i);
        }
        if (i <= n1) {
            c.col(l1-i) += b.col(n1-i);
        }
    }
    
    for (int i=1; i<=l2; i++) {
        if (i <= m2) {
            c.col(l1+i) += a.col(m1+i);
        }
        if (i <= n2) {
            c.col(l1+i) += b.col(n1+i);
        }
    }

    return c;
}

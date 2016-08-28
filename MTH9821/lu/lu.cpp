#include <lu.h>
#include <Eigen/Dense>
#include <cassert>
#include <tuple>

static void lu_helper(int k, int n,
                      Eigen::MatrixXd * A,
                      Eigen::MatrixXd * L,
                      Eigen::MatrixXd * U)
{
    for (int i=k; i<n; i++) {
        (*U)(k,i) = (*A)(k,i);
        (*L)(i,k) = (*A)(i,k)/(*U)(k,k);
    }

    for (int i=k+1; i<n; i++) {
        for (int j=k+1; j<n; j++) {
            (*A)(i,j) -= (*L)(i,k)*(*U)(k,j);
        }
    }
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> 
    lu_no_pivoting(const Eigen::MatrixXd & A)
{
    Eigen::MatrixXd Acopy = A;
    int n = Acopy.rows();
    assert(n == Acopy.cols());

    Eigen::MatrixXd L(n,n);
    Eigen::MatrixXd U(n,n);

    L.triangularView<Eigen::StrictlyUpper>().setZero();
    U.triangularView<Eigen::StrictlyLower>().setZero();

    for (int k=0; k<n-1; k++) {
        lu_helper(k, n, &Acopy, &L, &U);
    }

    L(n-1,n-1) = 1;
    U(n-1,n-1) = Acopy(n-1,n-1);

    return std::make_tuple(L,U);
}

std::tuple<Eigen::VectorXi, Eigen::MatrixXd, Eigen::MatrixXd> 
    lu_row_pivoting(const Eigen::MatrixXd & A)
{
    Eigen::MatrixXd Acopy = A;
    int n = Acopy.rows();
    assert(n == Acopy.cols());

    Eigen::VectorXi p = Eigen::VectorXi::LinSpaced(n,1,n);
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n,n);;
    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(n,n);;

    for (int k=0; k<n-1; k++) {
        int maxRow, maxCol;
        Acopy.block(k,k,n-k,1).array().abs().maxCoeff(&maxRow, &maxCol);
        Acopy.row(k).swap(Acopy.row(maxRow+k));
        p.row(k).swap(p.row(maxRow+k));
        L.row(k).swap(L.row(maxRow+k));
        lu_helper(k, n, &Acopy, &L, &U);
    }
    
    L(n-1,n-1) = 1;
    U(n-1,n-1) = Acopy(n-1,n-1);

    return std::make_tuple(p,L,U);
}


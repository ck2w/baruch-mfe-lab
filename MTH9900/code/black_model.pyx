import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt, log, exp, fabs, erfc, erf, fmax, fmin

np.import_array()

DEF SQRT2PI = 2.5066282746310002
DEF ISQRT2PI = 0.3989422804014327

cdef double norm_cdf(double x):
    return 0.5 * erfc(-0.70710678118654746 * x)

#cdef double norm_pdf(double x):
#    return ISQRT2PI*exp(-0.5 * x**2)

#cdef double lognormal_pdf(double x, double mu, double sigma):
#    return exp(-(log(x)-mu)**2 / (2 * sigma**2)) / (x * sigma * SQRT2PI)

# forward is assumed to be 1, opttype = 1 for call, -1 for put, 0 for timevalue
cdef double black_impv_cython( double K, double T, double V, int opttype,
                               double tol, int maxiter):

    if K <= 0 or T <= 0:
        return float('NaN')

    V -= fmax(opttype * (1-K), 0)

    if V < 0 or V >= fmin(K,1):
        return float('NaN')

    if V == 0:
        return 0

    cdef int j = 1
    cdef double x0, x1, d1, p, sqrtT
    sqrtT = sqrt(T)
    tol = tol*sqrtT
    p = log(K)

    if K >= 1:
        x0 = sqrt(2*p)
        x1 = x0 - (0.5 - K*norm_cdf(-x0) - V) * SQRT2PI
        while fabs(x0-x1) > tol and j < maxiter:
            x0 = x1
            d1 = -p/x1 + 0.5*x1
            x1 = x1 - (norm_cdf(d1) - K*norm_cdf(d1-x1) - V) * SQRT2PI * exp(0.5 * d1**2)
            j += 1
        return x1/sqrtT
    else:
        x0 = sqrt(-2*p)
        x1 = x0 - (0.5*K - norm_cdf(-x0) - V)*SQRT2PI/K
        while fabs(x0-x1) > tol and j < maxiter:
            x0 = x1
            d1 = -p/x1 + 0.5*x1
            x1 = x1 - (K*norm_cdf(x1-d1) - norm_cdf(-d1) - V) * SQRT2PI * exp(0.5 * d1**2)
            j += 1
        return x1/sqrtT


@cython.boundscheck(False)
@cython.wraparound(False)
def black_impv(K, T, F, V, opttype=1, double tol=1e-6, int maxiter=500):

    F = np.asfarray(F)
    K = np.asfarray(K)/F
    T = np.asfarray(T)
    V = np.asfarray(V)/F
    opttype = np.asarray(opttype, dtype=int)

    out = np.empty(np.broadcast(K,T,V,opttype).shape, dtype=float)

    cdef np.broadcast it = np.broadcast(K,T,V,opttype,out)

    while np.PyArray_MultiIter_NOTDONE(it):
        (<double*>np.PyArray_MultiIter_DATA(it,4))[0] = black_impv_cython(
            (<double*>np.PyArray_MultiIter_DATA(it,0))[0],
            (<double*>np.PyArray_MultiIter_DATA(it,1))[0],
            (<double*>np.PyArray_MultiIter_DATA(it,2))[0],
            (<int*>np.PyArray_MultiIter_DATA(it,3))[0],
            tol, maxiter)
        np.PyArray_MultiIter_NEXT(it)

    return out+0



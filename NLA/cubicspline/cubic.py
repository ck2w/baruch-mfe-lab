#!/usr/bin/anaconda3/bin/python

import math
import numpy as np

# input
#x = np.array([0.0, 1.0/3.0, 5.0/6.0, 4.0/3.0, 11.0/6.0])
#y = np.array([0.0421751, 0.0421751, 0.0200164, 0.0223065, 0.0239350])

#x = np.array([0.0, 1.0/6.0, 5.0/12.0, 11.0/12.0, 5.0/4.0])
#y = np.array([0.0075, 0.012012, 0.015651, 0.019815, 0.018206])

x = np.array([2,5,10])
Y = np.array([[4.69,4.57,4.63],
              [4.81,4.69,4.73],
              [4.81,4.70,4.74],
              [4.79,4.77,4.81],
              [4.79,4.77,4.80],
              [4.83,4.73,4.79],
              [4.81,4.72,4.76],
              [4.81,4.74,4.77],
              [4.83,4.77,4.80],
              [4.81,4.75,4.77],
              [4.82,4.77,4.80],
              [4.82,4.76,4.80],
              [4.80,4.75,4.78],
              [4.78,4.72,4.73],
              [4.79,4.71,4.73]])

for i in range(15):
    y = Y[i]

    n = len(x)

    z = np.zeros(n)

    for i in range(1,n-1):
        z[i] = 6.0 * ( (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]) )

    z = z[1:n-1]

    M = np.array([[0.0]*n for i in range(n)])

    for i in range(1,n-1):
        M[i][i] = 2.0 * (x[i+1]-x[i-1])

    for i in range(1,n-2):
        M[i][i+1] = x[i+1]-x[i]

    for i in range(2,n-1):
        M[i][i-1] = x[i]-x[i-1]

    M = np.array([M[i][1:n-1] for i in range(1,n-1)])

    w = np.linalg.solve(M, z)
    w = np.concatenate(([0], w, [0]))

    c = np.zeros(n)
    d = np.zeros(n)
    for i in range(1,n):
        c[i] = 0.5*(w[i-1]*x[i]-w[i]*x[i-1])/(x[i]-x[i-1])
        d[i] = (w[i]-w[i-1])/(6.0*(x[i]-x[i-1]))

    q = np.zeros(n)
    r = np.zeros(n)
    for i in range(1,n):
        q[i-1] = y[i-1] - c[i]*x[i-1]*x[i-1] - d[i]*x[i-1]*x[i-1]*x[i-1]
        r[i] = y[i] - c[i]*x[i]*x[i] - d[i]*x[i]*x[i]*x[i]

    a = np.zeros(n)
    b = np.zeros(n)
    for i in range(1,n):
        a[i] = (q[i-1]*x[i]-r[i]*x[i-1])/(x[i]-x[i-1])
        b[i] = (r[i]-q[i-1])/(x[i]-x[i-1])

    print(w)
    print(a)
    print(b)
    print(c)
    print(d)
    t=3
    print(a[1]+b[1]*t+c[1]*t*t+d[1]*t*t*t)


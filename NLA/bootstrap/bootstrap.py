#!/usr/bin/anaconda3/bin/python

import math
import numpy as np

# input
#x = np.array([0.0, 1.0/3.0, 5.0/6.0, 4.0/3.0, 11.0/6.0])
#y = np.array([0.0421751, 0.0421751, 0.0200164, 0.0223065, 0.0239350])

x = np.array([0.0, 1.0/12.0, 4.0/12.0, 10.0/12.0, 15.0/12.0, 21.0/12.0])
y = np.array([0.007, 0.0, 0,0, 0.0, 0.0, 0.0])
disc = np.array([None, 0.9980, 0.9935, 0.9820, 0.9775, 0.9620])
size = len(disc)
for i in range(size):
    if disc[i] is not None:
        y[i] = -math.log(disc[i])/x[i]

print("Interest Rates:")
for i in range(size):
    print("{0:.10f}".format(y[i]))

Y = [y]
numOfCurves = len(Y)
for i in range(numOfCurves):
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
    print("Tridiagonal system to solve:")
    print(M)
    print("More precision form:")
    for i in range(4):
        for j in range(4):
            print("{0:.10f}".format(M[i][j]), end=",")
        print("")



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

    print("Cubic spline coefficients:")
    print(w)
    print(a)
    print(b)
    print(c)
    print(d)

    print("Pricing a bond:")
    face=100.0
    value = 0.0
    coupon = face*(0.035/2.0)
    t = [1.0/12.0, 7.0/12.0, 13.0/12.0, 19.0/12.0]
    indices = [1,3,4,5]
    numOfPayments = len(indices)
    for i in range(numOfPayments):
        time = t[i]
        index = indices[i]
        rate = a[index] + b[index]*time + c[index]*time*time + d[index]*time*time*time
        print(time, rate)
        disc = math.exp(-rate*time)
        value += coupon*disc
    value += face*disc
    print(value)






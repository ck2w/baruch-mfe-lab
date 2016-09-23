#!/usr/bin/anaconda3/bin/python

import numpy as np
from scipy.optimize import fsolve
from scipy.stats import norm
import math

from io import StringIO

a = np.genfromtxt('data.csv',delimiter=',')

ns = len(a)
A = np.zeros((ns,2))

k_column=0
c_column=1
p_column=2
c_minus_p_column=3
A[:,0] = np.ones(ns)
A[:,1] = -a[:,0]
b = a[:, c_minus_p_column]

#m, c = np.linalg.lstsq(A, y)[0]
#print(m, c)

print("=== A ===")
for i in range(14):
    for j in range(2):
        print(A[i][j], end=",")
    print("")

print("=== b ===")
for i in range(14):
    print(b[i])

ATA = A.T.dot(A)
ATb = A.T.dot(b)

#print(ATA)
#print(ATb)

x = np.linalg.solve(ATA,ATb)

PVF = x[0]
disc = x[1]

print("=== PVF ===")
print(PVF)
print("=== disc ===")
print(disc)

T = 201/252
K = 0 
C = 0
P = 0

def call(x):
    d1 = math.log(PVF/(K*disc))/(x*math.sqrt(T))+0.5*x*math.sqrt(T)
    d2 = math.log(PVF/(K*disc))/(x*math.sqrt(T))-0.5*x*math.sqrt(T)
    cum1 = norm.cdf(d1)
    cum2 = norm.cdf(d2) 
    return PVF*cum1-K*disc*cum2-C

def put(x):
    d1 = math.log(PVF/(K*disc))/(x*math.sqrt(T))+0.5*x*math.sqrt(T)
    d2 = math.log(PVF/(K*disc))/(x*math.sqrt(T))-0.5*x*math.sqrt(T)
    cum1 = norm.cdf(-d1)
    cum2 = norm.cdf(-d2) 
    return -PVF*cum1+K*disc*cum2-P

karray = a[:,k_column]
carray = a[:,c_column]
parray = a[:,p_column]
for i in range(ns):
    K = karray[i]
    C = carray[i]
    P = parray[i]
    ivc = fsolve(call, 0.1)
    ivp = fsolve(put, 0.1)
    print(ivc, ivp)



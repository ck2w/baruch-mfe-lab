#!/usr/bin/anaconda3/bin/python

import numpy as np
from scipy.optimize import fsolve
from scipy.stats import norm
import math

from io import StringIO

a = np.genfromtxt('data.csv',delimiter=',')

ns = len(a)
A = np.zeros((ns,2))

A[:,0] = np.ones(ns)
A[:,1] = -a[:,0]
b = a[:, 7]

#m, c = np.linalg.lstsq(A, y)[0]
#print(m, c)

ATA = A.T.dot(A)
ATb = A.T.dot(b)

#print(ATA)
#print(ATb)

x = np.linalg.solve(ATA,ATb)

PVF = x[0]
disc = x[1]

print(PVF, disc)

T = 166/252
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

karray = a[:,0]
carray = a[:,3]
parray = a[:,6]
for i in range(ns):
    K = karray[i]
    C = carray[i]
    P = parray[i]
    ivc = fsolve(call, 0.1)
    ivp = fsolve(put, 0.1)
    print(ivc, ivp)



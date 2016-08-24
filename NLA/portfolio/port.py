#!/usr/bin/anaconda3/bin/python

import math
import numpy as np

mu = np.array([0.05, 0.09, 0.1])
sigma = np.array([[ 0.15*0.15,       -0.15*0.20*0.25, 0.15*0.25*0.50],
                  [-0.15*0.20*0.25,        0.20*0.20, 0.20*0.25*0.25],
                  [ 0.15*0.25*0.50,   0.20*0.25*0.25,      0.25*0.25]])
rf = 0.02
muP = 0.08
sigmaP = 0.27

# Tangency portfolio
ones = np.ones(len(mu))
nu = mu - rf
x = np.linalg.solve(sigma, nu)
den = np.dot(ones, x)
wT = x / den

print("Asset allocation of tangency portfolio:")
print(wT)

muT = rf + np.dot(nu, wT)
print("Expected return of tangency portfolio:", muT)

sigmaT = math.sqrt(np.dot(wT, np.dot(sigma, wT)))
print("Standard deviation of return of tangency portfolio:", sigmaT)

sharpeT = (muT-rf)/sigmaT
print("Sharpe ratio of tangency portfolio:", sharpeT)

# Min-variance portfolio
wMinCash = 1.0 - (muP - rf)/np.dot(nu, wT)
wMin = (1.0 - wMinCash)*wT

print("Asset allocation of min variance portfolio:")
print(wMin)
print("Cash allocation of min variance portfolio:", wMinCash)

sigmaMin = math.sqrt(np.dot(wMin, np.dot(sigma, wMin)))
print("Standard deviation of return of min variance portfolio:", sigmaMin)

sharpeMin = (muP-rf)/sigmaMin
print("Sharpe ratio of min variance portfolio:", sharpeMin)

# Max-return portfolio
sign = np.dot(ones, x)
wMaxCash = 0
if sign > 0:
    wMaxCash = 1.0 - sigmaP/sigmaT
else:
    wMaxCash = 1.0 + sigmaP/sigmaT

wMax = (1.0 - wMaxCash)*wT
muMax = rf + np.dot(nu, wMax)

print("Asset allocation of max return portfolio:")
print(wMax)
print("Cash allocation of max return portfolio:", wMaxCash)

print("Expected return of max return portfolio:", muMax)

sharpeMax = (muMax-rf)/sigmaP
print("Sharpe ratio of max return portfolio:", sharpeMax)

# Min-variance portfolio fully invested in assets
y = np.linalg.solve(sigma, ones)
denF = np.dot(ones, y)
wF = y/denF

print("Asset allocation of min variance portfolio fully invested in assets:")
print(wF)
muF = rf + np.dot(nu, wF)
print("Expected return of min variance portfolio fully invested in assets:", muF)
sigmaF = 1.0 / math.sqrt(denF)
print("Standard deviation of min variance portfolio fully invested in assets:", sigmaF)
sharpeF = (muF-rf)/sigmaF
print("Sharpe ratio of min variance portfolio fully invested in assets:", sharpeF)

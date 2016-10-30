#!/home/quan/anaconda3/bin/python

import math
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

f, axarr = plt.subplots(1)

# input parameters
R = 0.4 # recovery rate
r = 0.015 # riskless interest rate
timeStepsPerYear = 12
rollDatesPerYear = 4

# constants
bps = 0.0001

# data
name = 'GE'
CdsMaturities = np.array([1,2,3,4,5,7,10])
CdsParSpreads = np.array([19.35,25.45,31.85,38.0,47.9,71.0,93.6])*bps

# derived paramters
rollDatesInterval = int(timeStepsPerYear/rollDatesPerYear)
delta = 1.0/rollDatesPerYear # interval between roll dates
dt = 1.0/timeStepsPerYear # time step width for Stieltjes integral
c = (1-R)/delta
numOfMaturities = len(CdsMaturities)
maxMaturity = CdsMaturities[numOfMaturities-1]
maxNumOfTimeSteps = maxMaturity*timeStepsPerYear

def CalcParSpread(x, minUnknownIndex, maxUnknownIndex, knownIntensities):
    numerator = 0.0
    denominator = 0.0
    for i in range(maxUnknownIndex):
        l = 0.0
        if i >= minUnknownIndex: l = x
        else: l = knownIntensities[i]

        P = math.exp(-r*i*dt) + math.exp(-r*(i+1)*dt)
        S = math.exp(-l*i*dt) - math.exp(-l*(i+1)*dt)
        numerator += P*S
        
    for i in range(0,maxUnknownIndex,rollDatesInterval):
        l = 0.0
        if i >= minUnknownIndex: l = x
        else: l = knownIntensities[i]

        P = math.exp(-r*(i+rollDatesInterval)*dt)
        S = math.exp(-l*i*dt) + math.exp(-l*(i+rollDatesInterval)*dt)
        denominator += P*S
        
    return c*numerator/denominator

def CalcForwardParSpread(startMonth, endMonth, knownIntensities):
    numerator = 0.0
    denominator = 0.0
    for i in range(startMonth, endMonth):
        l = knownIntensities[i]
        P = math.exp(-r*i*dt) + math.exp(-r*(i+1)*dt)
        S = math.exp(-l*i*dt) - math.exp(-l*(i+1)*dt)
        numerator += P*S
        
    for i in range(startMonth, endMonth, rollDatesInterval):
        l = knownIntensities[i]
        P = math.exp(-r*(i+rollDatesInterval)*dt)
        S = math.exp(-l*i*dt) + math.exp(-l*(i+rollDatesInterval)*dt)
        denominator += P*S
        
    return c*numerator/denominator


def func(x, *data):
    return CalcParSpread(x, data[0], data[1], data[2])-data[3]

def CalcIntensity(parSpreads, knownIntensities):
    minUnknownIndex = None
    maxUnknownIndex = None
    currentParSpread = None
    for i in range(numOfMaturities):
        currentParSpread = parSpreads[i]
        minUnknownIndex = 0
        if i>0: minUnknownIndex = CdsMaturities[i-1]*timeStepsPerYear
        maxUnknownIndex = CdsMaturities[i]*timeStepsPerYear

        lambda0 = currentParSpread/(1-R)
        data = (minUnknownIndex, maxUnknownIndex, knownIntensities, currentParSpread)
        y = fsolve(func, lambda0, args=data)
        for j in range(minUnknownIndex,maxUnknownIndex):
            knownIntensities[j] = y[0]

    survivalRates = []
    for i in range(maxNumOfTimeSteps):
        survivalRates.append(math.exp(-(i+1)*dt*knownIntensities[i]))
  
    return survivalRates

timeAxis = []
for i in range(maxNumOfTimeSteps):
    timeAxis.append((i+1)*dt)

CdsParSpreads = np.array([19.35,25.45,31.85,38.0,47.9,71.0,93.6])*bps
knownIntensities = [None]*maxNumOfTimeSteps
survivalRates = CalcIntensity(CdsParSpreads, knownIntensities)
axarr.plot(timeAxis, survivalRates, label="GE")
print(CalcForwardParSpread(24,60,knownIntensities))

CdsParSpreads = np.array([38.8,51.5,61.45,72.5,89.6,111.05,131.7])*bps
knownIntensities = [None]*maxNumOfTimeSteps
survivalRates = CalcIntensity(CdsParSpreads, knownIntensities)
axarr.plot(timeAxis, survivalRates, label="JPM")
print(CalcForwardParSpread(24,60,knownIntensities))

CdsParSpreads = np.array([152.5,200.75,235.75,241.3,264.2,269.8,272.05])*bps
knownIntensities = [None]*maxNumOfTimeSteps
survivalRates = CalcIntensity(CdsParSpreads, knownIntensities)
axarr.plot(timeAxis, survivalRates, label="Axis")
print(CalcForwardParSpread(24,60,knownIntensities))

CdsParSpreads = np.array([393.35,463.2,563.25,696.4,727.2,749.4,735.5])*bps
knownIntensities = [None]*maxNumOfTimeSteps
survivalRates = CalcIntensity(CdsParSpreads, knownIntensities)
axarr.plot(timeAxis, survivalRates, label="MBIA")
print(CalcForwardParSpread(24,60,knownIntensities))

plt.legend(bbox_to_anchor=(0.05,0.1,0.2,0.2), loc=2,
        ncol=1, mode="expand", borderaxespad=0.)
axarr.set_xlabel('Time (years)', fontsize=20)
axarr.set_ylabel('Survival Probabiliy', fontsize=20)
plt.tight_layout()
plt.show()

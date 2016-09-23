#!/usr/bin/anaconda3/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt

f, axarr = plt.subplots(2)

print("-------------------------")
print("Simulation of Cox Process")
print("-------------------------")

lambda0 = 0.03 # initial intensity
theta = 0.005  # intensity drift
sigma = 0.008  # intensity diffusion
t = 10.0       # time horizon
ns = 1000000   # number of samples

# time points to record counts
timer = np.array([1,2,3,4,5,6,7,8,9,10])

dt = 0.01      # time step width
n = int(t/dt)  # time step number
dt_root = math.sqrt(dt)
theta_dt = theta*dt

indexer = (timer/dt).astype(int);
m = len(timer)
stats = np.zeros(m)

for k in range(ns):

    if k>1 and k%1000 == 0:
        print(k, "iterations...")
        for i in range(m):
            ti = timer[i]
            s = math.exp(-lambda0*ti-0.5*theta*ti*ti+sigma*sigma*ti*ti*ti/6.0)
            print(timer[i],",",s,",",stats[i],",",k-stats[i],",",1-stats[i]/k)

    timesteps = np.zeros(n+1)
    timesteps[0] = 0.0
    lambdas = np.zeros(n+1)
    lambdas[0] = lambda0
    counts = np.zeros(n+1)
    counts[0] = 0

    for i in range(n):
        # march forward in time
        timesteps[i+1] = timesteps[i] + dt
        # generate intensity process 
        normal = np.random.normal(0, dt_root)
        lambdas[i+1] = lambdas[i] + theta_dt + sigma*normal
        # generate counting process
        lambda_dt = lambdas[i+1]*dt
        uniform = np.random.random()
        if uniform <= lambda_dt:
            counts[i+1] = counts[i]+1
        else:
            counts[i+1] = counts[i]

    for i in range(m):
        if counts[indexer[i]] > 0:
            stats[i] += 1


#axarr[0].plot(timesteps, lambdas)
#axarr[1].plot(timesteps, counts)
#plt.tight_layout()
#plt.show()

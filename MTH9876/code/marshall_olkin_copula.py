#!/home/quan/anaconda3/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt

lambda1=1.7
lambda2=0.7
lambda12=0.3

ns=2000

T1  = np.random.exponential(1.0/lambda1,  ns)
T2  = np.random.exponential(1.0/lambda2,  ns)
T12 = np.random.exponential(1.0/lambda12, ns)

tau1 = np.minimum(T1, T12)
tau2 = np.minimum(T2, T12)

u1 = np.exp(-(lambda1+lambda12)*tau1)
u2 = np.exp(-(lambda2+lambda12)*tau2)

plt.figure(figsize=(8,8))
plt.plot(u1,u2,'o', markersize=3)
plt.xlabel("u1")
plt.ylabel("u2")
plt.show()

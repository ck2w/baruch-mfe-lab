#!/usr/bin/anaconda3/bin/python

import math
import numpy as np
import matplotlib.pyplot as plt

f, axarr = plt.subplots(1)

F=100.0
T=10.0
sigma=0.3
r=0.001
n=1000

T_sqrt = math.sqrt(T)
TwoPi_sqrt = math.sqrt(2*np.pi)

def func(x, V0):
    d = (math.log(V0/(F-x)) + (r-0.5*sigma*sigma)*T)/(sigma*T_sqrt)
    return math.exp(-0.5*d*d)/(sigma*TwoPi_sqrt*T_sqrt*(F-x))
    
vfunc = np.vectorize(func)
x = np.linspace(0, F, n) 
x = x[:n-1]

v0 = [95,99,120,140]
for V0 in v0:
    y = vfunc(x, V0)
    label_str = '$V_0 = ' + str(V0) + '$'
    axarr.plot(x, y, label=label_str)

plt.legend(bbox_to_anchor=(0.05,0.75,0.3,0.2), loc=2,
        ncol=1, mode="expand", borderaxespad=0.)
axarr.set_title('Loss Given Default ($\sigma=0.3, r=0.1\%$)', fontsize=20)
axarr.set_xlabel('x', fontsize=20)
axarr.set_ylabel('LGD', fontsize=20)
plt.tight_layout()
plt.show()

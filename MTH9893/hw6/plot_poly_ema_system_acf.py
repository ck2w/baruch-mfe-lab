#!/usr/bin/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt

from poly_ema import poly_ema
from compute_system_acf import compute_system_acf

Nwindow = 1000
Neff = 50

# Problem 2 (i)

f, axarr = plt.subplots(2, 1)

h1 = poly_ema(1, Neff, Nwindow)
h2 = poly_ema(2, Neff, Nwindow)
h3 = poly_ema(3, Neff, Nwindow)
h4 = poly_ema(4, Neff, Nwindow)
h5 = poly_ema(5, Neff, Nwindow)
    
L = np.arange(-Nwindow,Nwindow+1)

kh1 = compute_system_acf(h1)
kh2 = compute_system_acf(h2)
kh3 = compute_system_acf(h3)
kh4 = compute_system_acf(h4)
kh5 = compute_system_acf(h5)

khn1 = kh1 / kh1[Nwindow]
khn2 = kh2 / kh2[Nwindow]
khn3 = kh3 / kh3[Nwindow]
khn4 = kh4 / kh4[Nwindow]
khn5 = kh5 / kh5[Nwindow]

axarr[0].set_title("System Auto-correlation function - m = 1,2,3,4,5")
axarr[0].plot(L, kh1, label='m=1')
axarr[0].plot(L, kh2, label='m=2')
axarr[0].plot(L, kh3, label='m=3')
axarr[0].plot(L, kh4, label='m=4')
axarr[0].plot(L, kh5, label='m=5')

constants = np.ones(2*Nwindow+1)
constants = 0.05*constants

axarr[1].set_title("Normalized System Auto-correlation function - m = 1,2,3,4,5")
axarr[1].plot(L, khn1, label='m=1')
axarr[1].plot(L, khn2, label='m=2')
axarr[1].plot(L, khn3, label='m=3')
axarr[1].plot(L, khn4, label='m=4')
axarr[1].plot(L, khn5, label='m=5')
axarr[1].plot(L, constants, 'b--' )

for i in np.arange(Nwindow, 2*Nwindow+1):
    if khn1[i] < 0.05: 
        print("m = 1: l* = ", str(i-Nwindow))
        break

for i in np.arange(Nwindow, 2*Nwindow+1):
    if khn2[i] < 0.05: 
        print("m = 2: l* = ", str(i-Nwindow))
        break

for i in np.arange(Nwindow, 2*Nwindow+1):
    if khn3[i] < 0.05: 
        print("m = 3: l* = ", str(i-Nwindow))
        break

for i in np.arange(Nwindow, 2*Nwindow+1):
    if khn4[i] < 0.05: 
        print("m = 4: l* = ", str(i-Nwindow))
        break

for i in np.arange(Nwindow, 2*Nwindow+1):
    if khn5[i] < 0.05: 
        print("m = 5: l* = ", str(i-Nwindow))
        break

plt.legend(bbox_to_anchor=(0.6, 0.5, 0.2, 0.3), loc=1,
           ncol=1, mode="expand", borderaxespad=0.)

plt.tight_layout()
plt.show()

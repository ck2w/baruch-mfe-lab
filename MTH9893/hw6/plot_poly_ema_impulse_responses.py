#!/usr/bin/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt

from poly_ema import poly_ema
from compute_M0_M1_RHW import compute_M0_M1_RHW

from apply_poly_ema_m1_filter import apply_poly_ema_m1_filter
from apply_poly_ema_m2_filter import apply_poly_ema_m2_filter
from apply_poly_ema_m3_filter import apply_poly_ema_m3_filter
from apply_poly_ema_m4_filter import apply_poly_ema_m4_filter
from apply_poly_ema_m5_filter import apply_poly_ema_m5_filter

Nwindow = 1000
Neff = 50

# Problem 2 (e)

f, axarr = plt.subplots(3, 2)

h1 = poly_ema(1, Neff, Nwindow)
h2 = poly_ema(2, Neff, Nwindow)
h3 = poly_ema(3, Neff, Nwindow)
h4 = poly_ema(4, Neff, Nwindow)
h5 = poly_ema(5, Neff, Nwindow)

param_h1 = compute_M0_M1_RHW(h1)
param_h2 = compute_M0_M1_RHW(h2)
param_h3 = compute_M0_M1_RHW(h3)
param_h4 = compute_M0_M1_RHW(h4)
param_h5 = compute_M0_M1_RHW(h5)

axarr[0,0].set_title('m=1,2,3,4,5')
axarr[0,0].plot(h1)
axarr[0,0].plot(h2)
axarr[0,0].plot(h3)
axarr[0,0].plot(h3)
axarr[0,0].plot(h4)

axarr[0,1].set_title('m=1. open circle - fde')
axarr[0,1].plot(h1, label='m=1')
axarr[1,0].set_title('m=2; open circle - fde')
axarr[1,0].plot(h2, label='m=2')
axarr[1,1].set_title('m=3; open circle - fde')
axarr[1,1].plot(h3, label='m=3')
axarr[2,0].set_title('m=4; open circle - fde')
axarr[2,0].plot(h4, label='m=4')
axarr[2,1].set_title('m=5; open circle - fde')
axarr[2,1].plot(h5, label='m=5')

# Problem 2 (f)

print("M0 --- M1 --- RHW")
print(round(param_h1[0],6), round(param_h1[1],6), round(param_h1[2],6))
print(round(param_h2[0],6), round(param_h2[1],6), round(param_h2[2],6))
print(round(param_h3[0],6), round(param_h3[1],6), round(param_h3[2],6))
print(round(param_h4[0],6), round(param_h4[1],6), round(param_h4[2],6))
print(round(param_h5[0],6), round(param_h5[1],6), round(param_h5[2],6))

# Problem 2 (h)

lag = 1
impulse = np.zeros(Nwindow)
impulse[lag] = 1
candidate = apply_poly_ema_m1_filter(impulse, Neff)
impulse_response_fde = candidate[lag:]
axarr[0,1].plot(impulse_response_fde, 'o' ,markerfacecolor='none')

lag = 2
impulse = np.zeros(Nwindow)
impulse[lag] = 1
candidate = apply_poly_ema_m2_filter(impulse, Neff)
impulse_response_fde = candidate[lag:]
axarr[1,0].plot(impulse_response_fde, 'o' ,markerfacecolor='none')

lag = 3
impulse = np.zeros(Nwindow)
impulse[lag] = 1
candidate = apply_poly_ema_m3_filter(impulse, Neff)
impulse_response_fde = candidate[lag:]
axarr[1,1].plot(impulse_response_fde, 'o' ,markerfacecolor='none')

lag = 4
impulse = np.zeros(Nwindow)
impulse[lag] = 1
candidate = apply_poly_ema_m4_filter(impulse, Neff)
impulse_response_fde = candidate[lag:]
axarr[2,0].plot(impulse_response_fde, 'o' ,markerfacecolor='none')

lag = 5
impulse = np.zeros(Nwindow)
impulse[lag] = 1
candidate = apply_poly_ema_m5_filter(impulse, Neff)
impulse_response_fde = candidate[lag:]
axarr[2,1].plot(impulse_response_fde, 'o' ,markerfacecolor='none')

plt.tight_layout()
plt.show()



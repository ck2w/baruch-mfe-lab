#!/usr/bin/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage

from ema import ema
from delta import delta

def circular_convolve(series, impulse):
    "circular convolution"
    return np.fft.ifft(np.fft.fft(impulse)*np.fft.fft(series)).real

f, axarr = plt.subplots(3, 1)

Nwindow = 256
Neff = 32
Ndelta = 192

print("Problem 11 (a): Ema and delta impulse responses")

h_ema = ema(Neff, Nwindow)
h_delta = delta(Ndelta, Nwindow)

axarr[0].plot(h_ema)
axarr[0].set_title('ema with Neff = 32')
axarr[1].plot(h_delta)
axarr[1].set_title('delta at k = 192')

print("Problem 11 (b): Causal convolution")
candidate1 = np.convolve(h_delta, h_ema)
truncated1 = candidate1[:Nwindow]

candidate2 = circular_convolve(h_delta, h_ema)

axarr[2].plot(truncated1)
axarr[2].plot(candidate2)
axarr[2].set_title('Causal(blue) and Non-causal(red) convolution')

plt.tight_layout()
plt.show()




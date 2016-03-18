#!/usr/bin/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt

from ema import ema
from box import box
from comb import comb

Nwindow = 1024
Nperiod = 4

# (a) Ema impulse response
Neff = 128
h_ema = ema(Neff, Nwindow)

# (b) Box impulse response
Nbox = Neff / (1-np.exp(-1))
h_box = box(Nbox, Nwindow)

# (c) Comb impulse response
h_comb = comb(Nperiod, Nwindow)

# (d) Sampled ema
h_ema_sampled = np.multiply(h_comb, h_ema)

# (e) Sampled box
h_box_sampled = np.multiply(h_comb, h_box)

# (f-h) Amplitude spectra 
H_ema = np.absolute(np.fft.fft(h_ema))
H_box = np.absolute(np.fft.fft(h_box))
H_ema_sampled = Nperiod * np.absolute(np.fft.fft(h_ema_sampled))
H_box_sampled = Nperiod * np.absolute(np.fft.fft(h_box_sampled))

f, axarr = plt.subplots(2, 2)

axarr[0, 0].plot(h_ema)
axarr[0, 0].plot(h_ema_sampled)
axarr[0, 0].set_title('Ema and sampled ema')

axarr[0, 1].plot(h_box)
axarr[0, 1].plot(h_box_sampled)
axarr[0, 1].set_title('Box and sampled box')

axarr[1, 0].plot(H_ema)
axarr[1, 0].plot(H_ema_sampled)
axarr[1, 0].set_title('Ema and sampled ema spectra')

axarr[1, 1].plot(H_box)
axarr[1, 1].plot(H_box_sampled)
axarr[1, 1].set_title('Box and sampled box spectra')

plt.tight_layout()
plt.show()

# (i) Remarks on the sample period and the period of spectral replications
# The sample period is 4.
# The period of spectral replication is 1024/4 = 256.

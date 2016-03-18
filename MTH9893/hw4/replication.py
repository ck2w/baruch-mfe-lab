#!/usr/bin/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt

from ema import ema
from box import box
from comb import comb

Nwindow = 1024

# (a) Ema and box
Neff = 32
h_ema = ema(Neff, Nwindow)

Nbox = Neff / (1-np.exp(-1))
h_box = box(Nbox, Nwindow)

# (b) Comb with period = 256
Nperiod = 256
h_comb = comb(Nperiod, Nwindow)

# (c-f) Convolution
candidate1 = np.convolve(h_comb, h_ema)
candidate2 = np.convolve(h_comb, h_box)

h_ema_replicated = candidate1[:Nwindow]
h_box_replicated = candidate2[:Nwindow]

H_ema = np.absolute(np.fft.fft(h_ema))
H_box = np.absolute(np.fft.fft(h_box))
H_ema_replicated = np.absolute(np.fft.fft(h_ema_replicated)) / Nperiod
H_box_replicated = np.absolute(np.fft.fft(h_box_replicated)) / Nperiod

f, axarr = plt.subplots(2, 2)

axarr[0, 0].plot(h_ema)
axarr[0, 0].plot(h_ema_replicated)
axarr[0, 0].set_title('Ema and replicated ema')

axarr[0, 1].plot(h_box)
axarr[0, 1].plot(h_box_replicated)
axarr[0, 1].set_title('Box and replicated box')

axarr[1, 0].plot(H_ema)
axarr[1, 0].plot(H_ema_replicated)
axarr[1, 0].set_title('Ema and replicated ema spectra')

axarr[1, 1].plot(H_box)
axarr[1, 1].plot(H_box_replicated)
axarr[1, 1].set_title('Box and replicated box spectra')

plt.tight_layout()
plt.show()

# (g) Comment on the duality between sampling and replication
# Sampling in time domain corresponds to replication in frequency domain.
# Replication in time domain corresponds to sampling in frequency domain.

#!/usr/bin/anaconda3/bin/python

import csv
import numpy as np
import matplotlib.pyplot as plt

from scipy import optimize

from box import box
from ema import ema
from calc_box_ema_spectra_mse import calc_box_ema_spectra_mse

Nwindow = 1024
Nbox = 16

f, axarr = plt.subplots(2, 2)

# a) Indicative responses

Neff = Nbox / (1-np.exp(-1))

h_box = box(Nbox, Nwindow)
h_ema = ema(Neff, Nwindow)

axarr[0, 0].plot(h_box[:150])
axarr[0, 0].plot(h_ema[:150])
axarr[0, 0].set_title('Box and ema impulse response')

# b) Spectra

h_box_fft = np.fft.fft(h_box)
h_ema_fft = np.fft.fft(h_ema)

h_box_fft_abs = np.absolute(h_box_fft)
h_ema_fft_abs = np.absolute(h_ema_fft)

axarr[0, 1].plot(h_box_fft_abs)
axarr[0, 1].plot(h_ema_fft_abs, 'r--')

# c) Cumulative amplitude spectra

h_box_fft_abs_sum = np.cumsum(h_box_fft_abs[:Nwindow/2])
h_ema_fft_abs_sum = np.cumsum(h_ema_fft_abs[:Nwindow/2])

axarr[1, 0].plot(h_box_fft_abs_sum)
axarr[1, 0].plot(h_ema_fft_abs_sum, 'r--')

# d) See file calc_box_ema_spectra_mse.py

# e) Minimize the MSE f) Plot optimized cumulative spectrum

mse_func = lambda x: calc_box_ema_spectra_mse(x, Nbox, Nwindow)

res = optimize.minimize(mse_func, Neff)
Neff_optimal = res['x'][0]

h_ema_optimal = ema(Neff_optimal, Nwindow)
h_ema_optimal_fft = np.fft.fft(h_ema_optimal)
h_ema_optimal_fft_abs = np.absolute(h_ema_optimal_fft)
h_ema_optimal_fft_abs_sum = np.cumsum(h_ema_optimal_fft_abs[:Nwindow/2])

axarr[0, 1].plot(h_ema_optimal_fft_abs)
axarr[0, 1].set_title('Box and ema spectra, Neff*/Nbox = ' + str(round(Neff_optimal/Nbox, 2))) 

axarr[1, 0].plot(h_ema_optimal_fft_abs_sum)
axarr[1, 0].set_title('Cumulative spetra, Neff*/Nbox = ' + str(round(Neff_optimal/Nbox, 2))) 

# g) Convolve the price series

with open('jpm_trades.csv', 'r') as trades_csv:
    trades = np.array(list(csv.reader(trades_csv))[1:]).astype('double')
    prices = trades[:,1]

    candidate1 = np.convolve(h_box, prices-prices[0]) + prices[0]
    box_prices = candidate1[:len(prices)]

    candidate2 = np.convolve(h_ema_optimal, prices-prices[0]) + prices[0]
    ema_optimal_prices = candidate2[:len(prices)]

    axarr[1, 1].plot(box_prices)
    axarr[1, 1].plot(ema_optimal_prices)
    axarr[1, 1].set_title('Box and ema smoothed prices.') 

plt.tight_layout()
plt.show()



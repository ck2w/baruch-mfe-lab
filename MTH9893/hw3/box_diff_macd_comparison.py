#!/usr/bin/anaconda3/bin/python

import csv
import numpy as np
import matplotlib.pyplot as plt

from scipy import optimize

from box_diff import box_diff
from macd import macd
from switch_to_composite_unit_gauge import switch_to_composite_unit_gauge
from calc_box_diff_macd_spectra_mse import calc_box_diff_macd_spectra_mse

Nwindow = 1024
Nbox = 16
Neff_neg = 16
Neff_pos = Neff_neg/3.0


f, axarr = plt.subplots(2, 2)

# a) Indicative responses

h_box_diff = box_diff(Nbox, Nwindow)
h_macd = switch_to_composite_unit_gauge(macd(Neff_pos, Neff_neg, Nwindow))

axarr[0, 0].plot(h_box_diff[:150])
axarr[0, 0].plot(h_macd[:150])
axarr[0, 0].set_title('Box and macd differencer impulse response')

# b) Spectra

h_box_diff_fft = np.fft.fft(h_box_diff)
h_macd_fft = np.fft.fft(h_macd)

h_box_diff_fft_abs = np.absolute(h_box_diff_fft)
h_macd_fft_abs = np.absolute(h_macd_fft)

axarr[0, 1].plot(h_box_diff_fft_abs)
axarr[0, 1].plot(h_macd_fft_abs, 'r--')

# c) Cumulative amplitude spectra

h_box_diff_fft_abs_sum = np.cumsum(h_box_diff_fft_abs[:Nwindow/2])
h_macd_fft_abs_sum = np.cumsum(h_macd_fft_abs[:Nwindow/2])

axarr[1, 0].plot(h_box_diff_fft_abs_sum)
axarr[1, 0].plot(h_macd_fft_abs_sum, 'r--')

# d) See file calc_box_diff_macd_spectra_mse.py

# e) Minimize the MSE f) Plot optimized cumulative spectrum

mse_func = lambda x: calc_box_diff_macd_spectra_mse(x, Nbox, Nwindow)

res = optimize.minimize(mse_func, Neff_neg)
Neff_optimal = res['x'][0]

h_macd_optimal = switch_to_composite_unit_gauge(macd(Neff_optimal/3.0, Neff_optimal, Nwindow))
h_macd_optimal_fft = np.fft.fft(h_macd_optimal)
h_macd_optimal_fft_abs = np.absolute(h_macd_optimal_fft)
h_macd_optimal_fft_abs_sum = np.cumsum(h_macd_optimal_fft_abs[:Nwindow/2])

axarr[0, 1].plot(h_macd_optimal_fft_abs)
axarr[0, 1].set_title('Box diff and macd spectra, Nbox/Neff* = ' + str(round(Nbox/Neff_optimal, 3))) 

axarr[1, 0].plot(h_macd_optimal_fft_abs_sum)
axarr[1, 0].set_title('Cumulative spetra, Nbox/Neff* = ' + str(round(Nbox/Neff_optimal, 3))) 

# g) Convolve the price series

with open('jpm_trades.csv', 'r') as trades_csv:
    trades = np.array(list(csv.reader(trades_csv))[1:]).astype('double')
    prices = trades[:,1]

    candidate1 = np.convolve(h_box_diff, prices-prices[0])
    box_diff_prices = candidate1[:len(prices)]

    candidate2 = np.convolve(h_macd_optimal, prices-prices[0])
    macd_optimal_prices = candidate2[:len(prices)]

    axarr[1, 1].plot(box_diff_prices)
    axarr[1, 1].plot(macd_optimal_prices)
    axarr[1, 1].set_title('Box and macd differenced prices.') 


plt.tight_layout()
plt.show()



#!/usr/bin/anaconda3/bin/python

import csv
import numpy as np
import matplotlib.pyplot as plt

from scipy import optimize

from ema import ema
from int_macd_poly import int_macd_poly
from calc_ema_imacd_spectra_mse import calc_ema_imacd_spectra_mse

Nwindow = 1024
Neff_ema = 16
Neff_imacd = 16

f, axarr = plt.subplots(2, 2)

# a) Indicative responses

h_ema = ema(Neff_ema, Nwindow)
h_imacd = int_macd_poly(Neff_imacd, Nwindow)

axarr[0, 0].plot(h_ema[:150])
axarr[0, 0].plot(h_imacd[:150])
axarr[0, 0].set_title('Ema and integrated macd poly.')

# b) Spectra

h_ema_fft = np.fft.fft(h_ema)
h_imacd_fft = np.fft.fft(h_imacd)

h_ema_fft_abs = np.absolute(h_ema_fft)
h_imacd_fft_abs = np.absolute(h_imacd_fft)

axarr[0, 1].plot(h_ema_fft_abs)
axarr[0, 1].plot(h_imacd_fft_abs, 'r--')

# c) Cumulative amplitude spectra

h_ema_fft_abs_sum = np.cumsum(h_ema_fft_abs[:Nwindow/2])
h_imacd_fft_abs_sum = np.cumsum(h_imacd_fft_abs[:Nwindow/2])

axarr[1, 0].plot(h_ema_fft_abs_sum)
axarr[1, 0].plot(h_imacd_fft_abs_sum, 'r--')

# d) See file calc_ema_imacd_spectra_mse.py

# e) Minimize the MSE f) Plot optimized cumulative spectrum

print(calc_ema_imacd_spectra_mse(Neff_imacd, Neff_ema, Nwindow))

mse_func = lambda x: calc_ema_imacd_spectra_mse(x, Neff_ema, Nwindow)

res = optimize.minimize(mse_func, Neff_imacd)
Neff_imacd_optimal = res['x'][0]

h_imacd_optimal = int_macd_poly(Neff_imacd_optimal, Nwindow)
h_imacd_optimal_fft = np.fft.fft(h_imacd_optimal)
h_imacd_optimal_fft_abs = np.absolute(h_imacd_optimal_fft)
h_imacd_optimal_fft_abs_sum = np.cumsum(h_imacd_optimal_fft_abs[:Nwindow/2])

axarr[0, 1].plot(h_imacd_optimal_fft_abs)
axarr[0, 1].set_title('Ema and int-macd spectra, Neff_imacd/Neff_ema = ' + str(round(Neff_imacd_optimal/Neff_ema, 2))) 

axarr[1, 0].plot(h_imacd_optimal_fft_abs_sum)
axarr[1, 0].set_title('Cumulative spetra, Neff_imacd/Neff_ema = ' + str(round(Neff_imacd_optimal/Neff_ema, 2))) 

# g) Convolve the price series

with open('jpm_trades.csv', 'r') as trades_csv:
    trades = np.array(list(csv.reader(trades_csv))[1:]).astype('double')
    prices = trades[:,1]

    candidate1 = np.convolve(h_ema, prices-prices[0]) + prices[0]
    ema_prices = candidate1[:len(prices)]

    candidate2 = np.convolve(h_imacd_optimal, prices-prices[0]) + prices[0]
    imacd_optimal_prices = candidate2[:len(prices)]

    axarr[1, 1].plot(ema_prices)
    axarr[1, 1].plot(imacd_optimal_prices)
    axarr[1, 1].set_title('Ema and int-macd smoothed prices.') 

plt.tight_layout()
plt.show()



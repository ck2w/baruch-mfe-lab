#!/usr/bin/anaconda3/bin/python

import csv
import numpy as np
import matplotlib.pyplot as plt

from ema import ema
from delta import delta

f, axarr = plt.subplots(3, 2)

print("Problem 7: Ema and delta delay on price series")
    
Nwindow = 500

h_ema_1 = ema(2, Nwindow)
h_ema_2 = ema(5, Nwindow)
h_ema_3 = ema(10, Nwindow)
h_ema_4 = ema(20, Nwindow)
h_ema_5 = ema(50, Nwindow)
h_ema_6 = ema(100, Nwindow)

h_delta_1 = delta(2, Nwindow)
h_delta_2 = delta(5, Nwindow)
h_delta_3 = delta(10, Nwindow)
h_delta_4 = delta(20, Nwindow)
h_delta_5 = delta(50, Nwindow)
h_delta_6 = delta(100, Nwindow)

with open('jpm_trades.csv', 'r') as trades_csv:
    trades = np.array(list(csv.reader(trades_csv))[1:]).astype('double')
    prices = trades[:,1]

    candidate1 = np.convolve(prices-prices[0], h_ema_1) + prices[0]
    y1 = candidate1[:len(prices)]
    candidate11 = np.convolve(prices-prices[0], h_delta_1) + prices[0]
    y11 = candidate11[:len(prices)]
    axarr[0, 0].plot(y1)
    axarr[0, 0].plot(y11)
    axarr[0, 0].set_title('Ema smoothed and delta delay with Neff = 2')
    
    candidate2 = np.convolve(prices-prices[0], h_ema_2) + prices[0]
    y2 = candidate2[:len(prices)]
    candidate22 = np.convolve(prices-prices[0], h_delta_2) + prices[0]
    y22 = candidate22[:len(prices)]
    axarr[0, 1].plot(y2)
    axarr[0, 1].plot(y22)
    axarr[0, 1].set_title('Ema smoothed and delta delay with Neff = 5')

    candidate3 = np.convolve(prices-prices[0], h_ema_3) + prices[0]
    y3 = candidate3[:len(prices)]
    candidate33 = np.convolve(prices-prices[0], h_delta_3) + prices[0]
    y33 = candidate33[:len(prices)]
    axarr[1, 0].plot(y3)
    axarr[1, 0].plot(y33)
    axarr[1, 0].set_title('Ema smoothed and delta delay with Neff = 10')

    candidate4 = np.convolve(prices-prices[0], h_ema_4) + prices[0]
    y4 = candidate4[:len(prices)]
    candidate44 = np.convolve(prices-prices[0], h_delta_4) + prices[0]
    y44 = candidate44[:len(prices)]
    axarr[1, 1].plot(y4)
    axarr[1, 1].plot(y44)
    axarr[1, 1].set_title('Ema smoothed and delta delay with Neff = 20')

    candidate5 = np.convolve(prices-prices[0], h_ema_5) + prices[0]
    y5 = candidate5[:len(prices)]
    candidate55 = np.convolve(prices-prices[0], h_delta_5) + prices[0]
    y55 = candidate55[:len(prices)]
    axarr[2, 0].plot(y5)
    axarr[2, 0].plot(y55)
    axarr[2, 0].set_title('Ema smoothed and delta delay with Neff = 50')

    candidate6 = np.convolve(prices-prices[0], h_ema_6) + prices[0]
    y6 = candidate6[:len(prices)]
    candidate66 = np.convolve(prices-prices[0], h_delta_6) + prices[0]
    y66 = candidate66[:len(prices)]
    axarr[2, 1].plot(y6)
    axarr[2, 1].plot(y66)
    axarr[2, 1].set_title('Ema smoothed and delta delay with Neff = 100')

plt.tight_layout()
plt.show()

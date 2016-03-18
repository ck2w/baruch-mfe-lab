#!/usr/bin/anaconda3/bin/python

import csv
import numpy as np
import matplotlib.pyplot as plt

from ema import ema

f, axarr = plt.subplots(4, 2)

print("Problem 6: Ema on price series")
    
print("(a) Ema impulse response.")

Nwindow = 500

h_ema_1 = ema(2, Nwindow)
h_ema_2 = ema(5, Nwindow)
h_ema_3 = ema(10, Nwindow)
h_ema_4 = ema(20, Nwindow)
h_ema_5 = ema(50, Nwindow)
h_ema_6 = ema(100, Nwindow)

axarr[0, 0].plot(h_ema_1)
axarr[0, 0].plot(h_ema_2)
axarr[0, 0].plot(h_ema_3)
axarr[0, 0].plot(h_ema_4)
axarr[0, 0].plot(h_ema_5)
axarr[0, 0].plot(h_ema_6)
axarr[0, 0].set_title('Ema impulse responses with Neff = 2,5,10,20,50,100')

with open('jpm_trades.csv', 'r') as trades_csv:
    trades = np.array(list(csv.reader(trades_csv))[1:]).astype('double')
    prices = trades[:,1]

    candidate1 = np.convolve(prices-prices[0], h_ema_1) + prices[0]
    y1 = candidate1[:len(prices)]
    axarr[1, 0].plot(y1)
    axarr[1, 0].plot(prices)
    axarr[1, 0].set_title('Ema smoothed with Neff = 2')
    
    candidate2 = np.convolve(prices-prices[0], h_ema_2) + prices[0]
    y2 = candidate2[:len(prices)]
    axarr[1, 1].plot(y2)
    axarr[1, 1].plot(prices)
    axarr[1, 1].set_title('Ema smoothed with Neff = 5')

    candidate3 = np.convolve(prices-prices[0], h_ema_3) + prices[0]
    y3 = candidate3[:len(prices)]
    axarr[2, 0].plot(y3)
    axarr[2, 0].plot(prices)
    axarr[2, 0].set_title('Ema smoothed with Neff = 10')

    candidate4 = np.convolve(prices-prices[0], h_ema_4) + prices[0]
    y4 = candidate4[:len(prices)]
    axarr[2, 1].plot(y4)
    axarr[2, 1].plot(prices)
    axarr[2, 1].set_title('Ema smoothed with Neff = 20')

    candidate5 = np.convolve(prices-prices[0], h_ema_5) + prices[0]
    y5 = candidate5[:len(prices)]
    axarr[3, 0].plot(y5)
    axarr[3, 0].plot(prices)
    axarr[3, 0].set_title('Ema smoothed with Neff = 50')

    candidate6 = np.convolve(prices-prices[0], h_ema_6) + prices[0]
    y6 = candidate6[:len(prices)]
    axarr[3, 1].plot(y6)
    axarr[3, 1].plot(prices)
    axarr[3, 1].set_title('Ema smoothed with Neff = 100')

plt.tight_layout()
plt.show()

#!/usr/bin/anaconda3/bin/python

import csv
import numpy as np
import matplotlib.pyplot as plt

from diff import diff

f, axarr = plt.subplots(5, 1)

Nwindow=500

h_diff_1 = diff(2, 4, Nwindow)
h_diff_2 = diff(5, 10, Nwindow)
h_diff_3 = diff(10, 20, Nwindow)
h_diff_4 = diff(20, 40, Nwindow)
h_diff_5 = diff(50, 100, Nwindow)

with open('jpm_trades.csv', 'r') as trades_csv:
    trades = np.array(list(csv.reader(trades_csv))[1:]).astype('double')
    prices = trades[:,1]
    log_prices = np.log(prices)

    candidate1 = np.convolve(log_prices-log_prices[0], h_diff_1)
    y1 = candidate1[:len(prices)]
    axarr[0].plot(y1)
    axarr[0].set_title('Macd-M1 of log px series: Neff+ = 2')
    
    candidate2 = np.convolve(log_prices-log_prices[0], h_diff_2)
    y2 = candidate2[:len(prices)]
    axarr[1].plot(y2)
    axarr[1].set_title('Macd-M1 of log px series: Neff+ = 5')

    candidate3 = np.convolve(log_prices-log_prices[0], h_diff_3)
    y3 = candidate3[:len(prices)]
    axarr[2].plot(y3)
    axarr[2].set_title('Macd-M1 of log px series: Neff+ = 10')
    
    candidate4 = np.convolve(log_prices-log_prices[0], h_diff_4)
    y4 = candidate4[:len(prices)]
    axarr[3].plot(y4)
    axarr[3].set_title('Macd-M1 of log px series: Neff+ = 20')
    
    candidate5 = np.convolve(log_prices-log_prices[0], h_diff_5)
    y5 = candidate5[:len(prices)]
    axarr[4].plot(y5)
    axarr[4].set_title('Macd-M1 of log px series: Neff+ = 50')

plt.tight_layout()
plt.show()

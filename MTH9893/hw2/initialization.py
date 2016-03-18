#!/usr/bin/anaconda3/bin/python

import csv
import numpy as np
import matplotlib.pyplot as plt

from ema import ema
from macd import macd

f, axarr = plt.subplots(5, 1)

print("Problem 4: Initialization")

Neff = 20
Nwindow = 100
h_ema = ema(Neff, Nwindow)
    
Neff_pos = 20
Neff_neg = 40
h_macd = macd(Neff_pos, Neff_neg, Nwindow)

with open('jpm_trades.csv', 'r') as trades_csv:
    trades = np.array(list(csv.reader(trades_csv))[1:]).astype('double')
    prices = trades[:,1]
    candidate = np.convolve(prices, h_ema)
    truncated = candidate[:len(prices)]

    print("(a) Incorrect Initialization")
    axarr[0].plot(prices)
    axarr[0].plot(truncated)
    axarr[0].set_title('length of price series:' + str(len(prices)) + \
                       '; length of truncated series:' + str(len(truncated)))
    
    print("(b) Correct Initialization")
    candidate = np.convolve(prices-prices[0], h_ema)
    truncated = candidate[:len(prices)] 
    conditioned = truncated + prices[0]
    axarr[1].plot(prices)
    axarr[1].plot(conditioned)
    axarr[1].set_title('length of price series:' + str(len(prices)) + \
                       '; length of conditioned series:' + str(len(conditioned)))
    
    print("(c) Initialization for differencer")
    candidate = np.convolve(prices, h_macd)
    y1 = candidate[:len(prices)]

    candidate = np.convolve(prices-prices[0], h_macd)
    y3 = candidate[:len(prices)]
    y2 = y3 + prices[0]

    axarr[2].plot(y1)
    axarr[2].set_title('Macd, no initial conditioning')    
    
    axarr[3].plot(y2)
    axarr[3].set_title('Macd, initial conditioning with adding initial value back')    
    
    axarr[4].plot(y3)
    axarr[4].set_title('Macd, initial conditioning without adding initial value back')    

plt.tight_layout()
plt.show()

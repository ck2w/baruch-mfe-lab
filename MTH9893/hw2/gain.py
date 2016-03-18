#!/usr/bin/anaconda3/bin/python

import csv
import numpy as np
import matplotlib.pyplot as plt

from ema import ema
from macd import macd

f, axarr = plt.subplots(2, 3)

print("Problem 5: Gain")

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

    print("(a) Ema gain")
    # The gain = 1 is the correct choice.
    # The factor g must be multiplied with the offset prices[0]
    # because 
    # conv(g h_ema, px[n]-px[0]) = conv(g h_ema, px[n]) - g px[0]
    # produces a term g px[0] because the gain of the ema impulse
    # response itself is one.

    g = 0.2
    candidate1 = np.convolve(prices-prices[0], g*h_ema) + g*prices[0]
    y1 = candidate1[:len(prices)] 
    axarr[0, 0].plot(prices)
    axarr[0, 0].plot(y1)
    axarr[0, 0].set_title('Ema gain: g=0.2')
    
    g = 1.0
    candidate2 = np.convolve(prices-prices[0], g*h_ema) + g*prices[0]
    y2 = candidate2[:len(prices)] 
    axarr[0, 1].plot(prices)
    axarr[0, 1].plot(y2)
    axarr[0, 1].set_title('Ema gain: g=1.0')
    
    g = 5.0
    candidate3 = np.convolve(prices-prices[0], g*h_ema) + g*prices[0]
    y3 = candidate3[:len(prices)] 
    axarr[0, 2].plot(prices)
    axarr[0, 2].plot(y3)
    axarr[0, 2].set_title('Ema gain: g=5')
    
    print("(b) Macd gain")
    # The gain of the macd impulse response itself is zero.
    # There is no reason to add the term g px[0] back.

    g = 0.2
    candidate4 = np.convolve(prices-prices[0], g*h_macd) + g*prices[0]
    y4 = candidate4[:len(prices)] 
    axarr[1, 0].plot(y4)
    axarr[1, 0].set_title('Macd gain: g=0.2')
    
    g = 1.0
    candidate5 = np.convolve(prices-prices[0], g*h_macd) + g*prices[0]
    y5 = candidate5[:len(prices)] 
    axarr[1, 1].plot(y5)
    axarr[1, 1].set_title('Macd gain: g=1.0')
    
    g = 5.0
    candidate6 = np.convolve(prices-prices[0], g*h_macd) + g*prices[0]
    y6 = candidate6[:len(prices)] 
    axarr[1, 2].plot(y6)
    axarr[1, 2].set_title('Macd gain: g=5')

plt.tight_layout()
plt.show()

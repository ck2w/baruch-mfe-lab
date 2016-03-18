#!/usr/bin/anaconda3/bin/python

import csv
import numpy as np
import matplotlib.pyplot as plt

from ema import ema

f, axarr = plt.subplots(2, 1)

print("Problem 3: Trunction")

Neff = 20
Nwindow = 100
h_ema = ema(Neff, Nwindow)

with open('jpm_trades.csv', 'r') as trades_csv:
    trades = np.array(list(csv.reader(trades_csv))[1:]).astype('double')
    prices = trades[:,1]
    candidate = np.convolve(prices, h_ema)

    print("(a) Incorrect Truncation")
    axarr[0].plot(prices)
    axarr[0].plot(candidate)
    axarr[0].set_title('length of price series:' + str(len(prices)) + \
                       '; length of convolved series:' + str(len(candidate)))

    print("(b) Correct Truncation")
    # only strictly valid region is used
    truncated = candidate[:len(prices)]
    axarr[1].plot(prices)
    axarr[1].plot(truncated)
    axarr[1].set_title('length of price series:' + str(len(prices)) + \
                       '; length of truncated series:' + str(len(truncated)))

plt.tight_layout()
plt.show()

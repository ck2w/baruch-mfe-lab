#!/usr/bin/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
import csv

from compute_ab_ema import compute_ab_ema

f, axarr = plt.subplots(3, 2)

with open('jpm_quotes.csv', 'r') as quotes_csv:
    quotes = np.array(list(csv.reader(quotes_csv))[1:]).astype('double')
    prices = quotes[:,1]
    counts = quotes[:,2]

    Neff = 32
    results = compute_ab_ema(Neff, prices, counts)
    
    ema_prices = results[0]
    vwap_prices = results[2]

    ema_counts = results[1]

    Neff_star = results[3]
    Neff_star_sorted = np.sort(Neff_star)
    prob_axis = np.arange(len(Neff_star)) / (len(Neff_star) - 1)

    axarr[0, 0].set_title('Prices, Ema(dashed) and VWAP')
    axarr[0, 0].plot(prices)
    axarr[0, 0].plot(ema_prices, '--')
    axarr[0, 0].plot(vwap_prices)

    axarr[0, 1].set_title('Intensity series (quote counts)')
    axarr[0, 1].plot(counts)

    axarr[1, 0].set_title('Neff*')
    axarr[1, 0].plot(Neff_star)
    
    axarr[1, 1].set_title('CDF of Neff*')
    axarr[1, 1].plot(Neff_star_sorted, prob_axis)

    axarr[2, 0].set_title('Neff* .vs. Quote Intensity')
    axarr[2, 0].scatter(counts[1:], Neff_star[1:])
    
    axarr[2, 1].set_title('Neff* .vs. Prices')
    axarr[2, 1].scatter(prices[1:], Neff_star[1:])

    

plt.tight_layout()
plt.show()

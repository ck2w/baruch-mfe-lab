#!/usr/bin/anaconda3/bin/python

import csv
import numpy as np
import matplotlib.pyplot as plt

f, axarr = plt.subplots(2, 2)

with open('jpm_quotes.csv', 'r') as quotes_csv:
    quotes = np.array(list(csv.reader(quotes_csv))[1:]).astype('double')
    
    axarr[0, 0].plot(quotes[:,1])
    axarr[0, 0].set_title('quote mid')
    
    axarr[0, 1].plot(quotes[:,2])
    axarr[0, 1].set_title('quote cnt')

with open('jpm_trades.csv', 'r') as trades_csv:
    trades = np.array(list(csv.reader(trades_csv))[1:]).astype('double')
    
    axarr[1, 0].plot(trades[:,1])
    axarr[1, 0].set_title('trade price')
    
    axarr[1, 1].plot(trades[:,4])
    axarr[1, 1].set_title('trade volume')
    
plt.tight_layout()
plt.show()

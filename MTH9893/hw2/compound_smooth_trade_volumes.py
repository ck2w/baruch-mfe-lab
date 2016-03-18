#!/usr/bin/anaconda3/bin/python

import csv
import numpy as np
import matplotlib.pyplot as plt

from ema import ema
from step import step
from macd import macd
from compound_smooth_impulse_response import compound_smooth_impulse_response

f, axarr = plt.subplots(3, 1)

Nwindow = 500
Neff = 21
Neff_pos = 7
Neff_neg = 14
h_ema = ema(Neff, Nwindow)
h_step = step(Nwindow)
h_macd = macd(Neff_pos, Neff_neg, Nwindow)

with open('jpm_trades.csv', 'r') as trades_csv:
    trades = np.array(list(csv.reader(trades_csv))[1:]).astype('double')
    volumes = trades[:,4]

    #ema   
    print("Problem 10 (a): Ema on trade volume series")
    candidate1 = np.convolve(volumes-volumes[0], h_ema) + volumes[0]
    truncated1 = candidate1[:len(volumes)]
    axarr[0].plot(truncated1)
    axarr[0].plot(volumes)
    axarr[0].set_title('Ema applied to volume series with Neff = 21')
    
    print("Problem 10 (b): Compound smoothing on volume series with Neff+ = 7, Neff- = 14")
    # compound step 1
    candidate2 = np.convolve(volumes-volumes[0], h_step) + Nwindow * volumes[0]
    truncated2 = candidate2[:len(volumes)]

    # compound step 2
    gauge = 1.0 / (Neff_neg - Neff_pos)
    candidate3 = np.convolve(truncated2-truncated2[0], gauge * h_macd) + volumes[0] 
    truncated3 = candidate3[:len(volumes)]
    axarr[1].plot(truncated1)
    axarr[1].plot(truncated3)
    axarr[1].set_title('Ema and compound smoothing: Neff+ = 7, Neff- = 14.')
    
    print("Problem 10 (c): New compound smoothing impulse response on volume series with Neff+ = 7, Neff- = 14")
    # the result is exactly the same as the compound smoothing obtained in part (b)
    h_cmp = compound_smooth_impulse_response(Neff_pos, Neff_neg, Nwindow)
    candidate4 = np.convolve(volumes-volumes[0], h_cmp) + volumes[0]
    truncated4 = candidate4[:len(volumes)]
    axarr[2].plot(truncated3)
    axarr[2].plot(truncated4)
    axarr[2].set_title('New compound IR smoothing (conv of step and macd) and previous.')

plt.tight_layout()
plt.show()

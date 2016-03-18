#!/usr/bin/anaconda3/bin/python

import numpy as np

from ema import ema

def macd(Neff_pos, Neff_neg, Nwindow):
    "moving average convergence/divergence"

    h_pos = ema(Neff_pos, Nwindow)
    h_neg = ema(Neff_neg, Nwindow)
    h = h_pos - h_neg

    return h

if __name__ == "__main__":
    print(macd(4, 2, 10))

    

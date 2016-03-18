#!/usr/bin/anaconda3/bin/python

import numpy as np

def ema(Neff, Nwindow):
    "an exponential moving average function"

    p = Neff/ (Neff + 1.0)
    n = np.arange(Nwindow)
    h = (1-p) * p**n

    return h

if __name__ == "__main__":
    print(ema(4, 10))
    print(ema(2, 10))

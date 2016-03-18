#!/usr/bin/anaconda3/bin/python

import numpy as np

from ema import ema
from ema_poly1 import ema_poly1

def macd_poly(Neff, Nwindow):
    "Macd-Poly differencer"

    if Neff > 0:
        h_pos = ema(Neff/3.0, Nwindow)
        h_neg = ema_poly1(Neff, Nwindow)
        h = h_pos - h_neg
    else:
        h = np.zeros(Nwindow)
        h[0]=1
        h[1]=-1

    return h
    
if __name__ == "__main__":
    print(macd_poly(4, 10))


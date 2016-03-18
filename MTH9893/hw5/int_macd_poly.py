#!/usr/bin/anaconda3/bin/python

import numpy as np

from macd_poly import macd_poly

def int_macd_poly(Neff, Nwindow):
    "Integrated Macd-Poly function"

    h_macd_poly = macd_poly(Neff, Nwindow)

    # integrate
    g = 3.0 / (2*Neff)
    h = g * np.cumsum(h_macd_poly) 

    return h
    
if __name__ == "__main__":
    print(int_macd_poly(4, 10))


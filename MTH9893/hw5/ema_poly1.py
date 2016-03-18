#!/usr/bin/anaconda3/bin/python

import numpy as np

def ema_poly1(Neff, Nwindow):
    "Ema-Poly1 function"

    if Neff > 0:
        p = ( Neff / 2.0 ) / ( (Neff / 2.0) + 1)
        n = np.arange(Nwindow)
        h = (1-p)*(1-p)*(n+1)*np.exp(n*np.log(p))
    else:
        h = np.zeros(Nwindow)
        h[0] = 1
    
    return h
    
if __name__ == "__main__":
    print(ema_poly1(4, 10))


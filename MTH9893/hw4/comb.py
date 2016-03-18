#!/usr/bin/anaconda3/bin/python

import numpy as np

def comb(Nperiod, Nwindow):
    "a finite comb with period Nperiod"

    h = np.zeros(Nwindow)
    if Nperiod > 0:
        h[0:Nwindow:Nperiod] = 1
    else:
        h[0] = 1

    return h

if __name__ == "__main__":
    print(comb(4, 20))
        

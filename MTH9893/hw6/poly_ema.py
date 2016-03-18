#!/usr/bin/anaconda3/bin/python

import numpy as np

def poly_ema(m, Neff, Nwindow):
    "polynomial ema impulse response"

    p = Neff / (Neff/m + 1)
    pm = p/m
    n = np.arange(Nwindow)

    if m == 1:
        h = (1-pm) * pm**n
    elif m == 2:
        h = (1-pm) * (1-pm) * (n+1) * pm**n
    elif m == 3:
        h = (1-pm) * (1-pm) * (1-pm) * (n+1) * (n+2) * pm**n / 2
    elif m == 4:
        h = (1-pm) * (1-pm) * (1-pm) * (1-pm) * (n+1) * (n+2) * (n+3) * pm**n / 6
    elif m == 5:
        h = (1-pm) * (1-pm) * (1-pm) * (1-pm) * (1-pm) * (n+1) * (n+2) * (n+3) * (n+4) * pm**n / 24
    else:
        h = np.zeros(Nwindow)

    return h
    

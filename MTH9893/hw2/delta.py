#!/usr/bin/anaconda3/bin/python

import numpy as np

def delta(k, Nwindow):
    "a kronecker delta function"

    h = np.zeros(Nwindow)

    if k >=0 and k<Nwindow:
        h[k] = 1

    return h

if __name__ == "__main__":
    print(delta(3, 10))

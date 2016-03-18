#!/usr/bin/anaconda3/bin/python

import numpy as np

def step(Nwindow):
    "a heaviside step funcion"

    h = np.ones(Nwindow)
    return h

if __name__ == "__main__":
    print(step(10))

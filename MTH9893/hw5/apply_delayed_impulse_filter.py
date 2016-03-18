#!/usr/bin/anaconda3/bin/python

import numpy as np

def apply_delayed_impulse_filter(x, N):
    "Delayed impulse: y[n] = x[n-N]"

    NT = len(x)
    y = np.zeros(NT)
   
    for i in range(NT):
        if i >= N:
            y[i] = x[i-N]

    return y

if __name__ == "__main__":
    x = np.array([1,2,3,4,5,6,7,8,9,10])
    print(apply_delayed_impulse_filter(x, 6))

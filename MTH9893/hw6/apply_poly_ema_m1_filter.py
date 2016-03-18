#!/usr/bin/anaconda3/bin/python

import numpy as np

def apply_poly_ema_m1_filter(x, N):
    "Ema filter: y[n] = p y[n-1] + (1-p) x[n]"

    p = N/(N+1)

    lag = 1
    
    x_expanded = np.concatenate((np.zeros(lag), x-x[0]))

    expanded_length = len(x_expanded)

    y_expanded = np.zeros(expanded_length)

    for i in range(lag, expanded_length):
        y_expanded[i] = p*y_expanded[i-1] + (1-p)*x_expanded[i]

    y = y_expanded[lag:] + x[0]

    return y

if __name__ == "__main__":
    x = np.array([11,12,13,14,15,16,17,18,19,20])
    print(apply_ema_filter(x, 6))

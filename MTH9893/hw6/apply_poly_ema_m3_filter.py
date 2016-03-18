#!/usr/bin/anaconda3/bin/python

import numpy as np

def apply_poly_ema_m3_filter(x, N):
    "y[n] = p y[n-1] - 3 * (p/3)^2 y[n-2] + (p/3)^3 y[n-3] + (1-p/3)^3 x[n] "

    p = N / (1+N/3) 

    lag = 3

    x_expanded = np.concatenate((np.zeros(lag), x-x[0]))

    expanded_length = len(x_expanded)

    y_expanded = np.zeros(expanded_length)

    for i in range(lag, expanded_length):
        y_expanded[i] = p*y_expanded[i-1] - 3*(p/3)**2 * y_expanded[i-2] + (p/3)**3 * y_expanded[i-3] + (1-p/3)**3 * x_expanded[i]

    y = y_expanded[lag:] + x[0]

    return y

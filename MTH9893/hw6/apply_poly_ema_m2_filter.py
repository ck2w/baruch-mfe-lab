#!/usr/bin/anaconda3/bin/python

import numpy as np

def apply_poly_ema_m2_filter(x, N):
    "y[n] = p y[n-1] -(p/2)^2 y[n-2] + (1-p/2)^2 x[n]"

    p = N / (1+N/2)

    lag = 2

    x_expanded = np.concatenate((np.zeros(lag), x-x[0]))

    expanded_length = len(x_expanded)

    y_expanded = np.zeros(expanded_length)
    
    for i in range(lag, expanded_length):
        y_expanded[i] = p*y_expanded[i-1] - (p/2)*(p/2)*y_expanded[i-2] + (1-p/2)*(1-p/2)*x_expanded[i]

    y = y_expanded[lag:] + x[0]

    return y

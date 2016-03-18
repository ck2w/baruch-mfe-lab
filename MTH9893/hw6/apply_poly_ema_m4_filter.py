#!/usr/bin/anaconda3/bin/python

import numpy as np

def apply_poly_ema_m4_filter(x, N):
    "y[n] = p y[n-1] - 6 * (p/4)^2 y[n-2] + 4 * (p/4)^3 y[n-3] - (p/4)^4 y[n-4] + (1-p/4)^4 x[n] "

    p = N / (1+N/4) 

    lag = 4 

    x_expanded = np.concatenate((np.zeros(lag), x-x[0]))

    expanded_length = len(x_expanded)

    y_expanded = np.zeros(expanded_length)
    
    for i in range(lag, expanded_length):
        y_expanded[i] = p*y_expanded[i-1] - 6*(p/4)**2 * y_expanded[i-2] + 4*(p/4)**3 * y_expanded[i-3] - (p/4)**4 * y_expanded[i-4] + (1-p/4)**4 * x_expanded[i]

    y = y_expanded[lag:] + x[0]

    return y

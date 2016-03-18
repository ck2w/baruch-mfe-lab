#!/usr/bin/anaconda3/bin/python

import numpy as np

def apply_poly_ema_m5_filter(x, N):
    "y[n] = p y[n-1] - 10 * (p/5)^2 y[n-2] + 10 * (p/5)^3 y[n-3] - 5 * (p/5)^4 y[n-4] + (p/5)^5 y[n-5] + (1-p/4)^4 x[n] "

    p = N / (1+N/5) 

    lag = 5

    x_expanded = np.concatenate((np.zeros(lag), x-x[0]))

    expanded_length = len(x_expanded)

    y_expanded = np.zeros(expanded_length)
    
    for i in range(lag, expanded_length):
        y_expanded[i] = p*y_expanded[i-1] - 10*(p/5)**2 * y_expanded[i-2] + 10*(p/5)**3 * y_expanded[i-3] - 5*(p/5)**4 * y_expanded[i-4] + (p/5)**5 * y_expanded[i-5] + (1-p/5)**5 * x_expanded[i]

    y = y_expanded[lag:] + x[0]

    return y

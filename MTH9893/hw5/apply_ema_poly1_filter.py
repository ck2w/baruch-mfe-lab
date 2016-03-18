#!/usr/bin/anaconda3/bin/python

import numpy as np

def apply_ema_poly1_filter(x, N):
    "Ema poly1 filter: y[n] = 2p y[n-1] -p^2 y[n-2] + (1-p)^2 x[n]"

    p = (N/2)/((N/2)+1)

    lag = 2

    x_expanded = np.concatenate((np.zeros(lag), x-x[0]))

    expanded_length = len(x_expanded)

    y_expanded = np.zeros(expanded_length)
    
    for i in range(lag, expanded_length):
        y_expanded[i] = 2*p*y_expanded[i-1] - p*p*y_expanded[i-2] + (1-p)*(1-p)*x_expanded[i]

    y = y_expanded[lag:] + x[0]

    return y

if __name__ == "__main__":
    x = np.array([11,12,13,14,15,16,17,18,19,20])
    print(apply_ema_poly1_filter(x, 6))

    


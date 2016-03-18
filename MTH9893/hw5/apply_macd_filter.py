#!/usr/bin/anaconda3/bin/python

import numpy as np

def apply_macd_filter(x, Neff_pos, Neff_neg):
    "Macd filter: y[n] = (p+ + p-) y[n-1] - p+p- y[n-2] - (p+ - p-) (x[n] - x[n-1])"

    p_pos = Neff_pos / (Neff_pos + 1)
    p_neg = Neff_neg / (Neff_neg + 1)

    lag = 2

    x_expanded = np.concatenate((np.zeros(lag), x-x[0]))

    expanded_length = len(x_expanded)
    
    y_expanded = np.zeros(expanded_length)
    
    a = p_pos + p_neg
    b = p_pos - p_neg
    c = p_pos * p_neg
    for i in range(lag, expanded_length):
        y_expanded[i] = a * y_expanded[i-1] - c * y_expanded[i-2] - b * (x_expanded[i] - x_expanded[i-1])

    y = y_expanded[lag:]

    return y

if __name__ == "__main__":
    x = np.array([11,12,13,14,15,16,17,18,19,20])
    print(apply_macd_filter(x, 2, 4))

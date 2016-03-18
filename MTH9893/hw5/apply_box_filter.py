#!/usr/bin/anaconda3/bin/python

import numpy as np

def apply_box_filter(x, N):
    "Box filter: y[n] = y[n-1] + (x[n] - x[n-N]) / N"

    lag = max(1,N)

    x_expanded = np.concatenate((np.zeros(lag), x - x[0]))

    expanded_length = len(x_expanded)

    y_expanded = np.zeros(expanded_length)

    for i in range(lag, expanded_length):
        y_expanded[i] = y_expanded[i-1] + (x_expanded[i] - x_expanded[i-N])/N

    y = y_expanded[lag:] + x[0]

    return y


if __name__ == "__main__":
    x = np.array([11,12,13,14,15,16,17,18,19,20])
    print(apply_box_filter(x, 6))

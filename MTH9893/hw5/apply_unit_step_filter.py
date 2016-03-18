#!/usr/bin/anaconda3/bin/python

import numpy as np

def apply_unit_step_filter(x):
    "Unit step: y[n] = y[n-1] + x[n]"

    lag = 1

    x_expanded = np.concatenate((np.zeros(lag), x))

    expanded_length = len(x_expanded)

    y_expanded = np.zeros(expanded_length)

    for i in range(lag, expanded_length):
        y_expanded[i] = y_expanded[i-1] + x_expanded[i]
   
    y = y_expanded[lag:] 
    
    return y
    

if __name__ == "__main__":
    x = np.array([1,2,3,4,5,6,7,8,9,10])
    print(apply_unit_step_filter(x))

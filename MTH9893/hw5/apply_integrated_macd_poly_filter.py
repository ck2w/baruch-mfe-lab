#!/usr/bin/anaconda3/bin/python

import numpy as np

def apply_integrated_macd_poly_filter(x, N):
    "Integrated macd poly filter"

    p1 = (N/3) / ((N/3)+1)
    p2 = (N/2) / ((N/2)+1)

    g = 2*N/3
    
    x0 = ( 2*p2 - p2*p2 - p1 ) / g
    x1 = -p2*p2*(1-p1) / g

    y1 = 2*p2+p1
    y2 = -p2*(p2+2*p1)
    y3 = p2*p2*p1 

    lag = 3

    x_expanded = np.concatenate((np.zeros(lag), x-x[0]))

    expanded_length = len(x_expanded)

    y_expanded = np.zeros(expanded_length)

    for i in range(lag, expanded_length):
        y_expanded[i] = y1*y_expanded[i-1] +y2*y_expanded[i-2] + y3*y_expanded[i-3] + x0*x_expanded[i] + x1*x_expanded[i-1]  

    y = y_expanded[lag:] + x[0]
    
    return y

if __name__ == "__main__":
    x = np.array([11,12,13,14,15,16,17,18,19,20])
    print(apply_integrated_macd_poly_filter(x, 6))

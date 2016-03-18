#!/usr/bin/anaconda3/bin/python

import numpy as np

def compute_M0_M1_RHW(h):
    "Compute M0, M1, and RHW"
    
    n = np.arange(len(h))

    M0 = np.sum(h)
    M1 = np.sum(n*h)
    M2 = np.sum(n*n*h)
    RHW = np.sqrt(M2-M1*M1)/M1

    return (M0, M1, RHW)

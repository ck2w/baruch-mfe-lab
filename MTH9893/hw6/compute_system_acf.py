#!/usr/bin/anaconda3/bin/python

import numpy as np

def compute_system_acf(h):
    "Compute system auto-correlation function"

    N = len(h)
    L = np.arange(-N,N+1)

    kh = np.zeros(len(L))
    for l in L:
        for i in np.arange(N):
            if i-l >= 0 and i-l < N:
                kh[l+N] += h[i]*h[i-l]

    return kh

if __name__ == "__main__":
    print(compute_system_acf([3,4,5]))

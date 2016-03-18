#!/usr/bin/anaconda3/bin/python

import numpy as np

from scipy import signal

def partial_fraction_expansion(b, a):
    "Return the pole, residue, and polynomial of a partial fraction expansion"

    return signal.residue(b, a)

if __name__ == "__main__":
    b = [1,0,0,0]
    a = [1,-0.35,-0.18,0.14]
    results = partial_fraction_expansion(b, a)
    print(results[0])
    print(results[1])
    print(results[2])

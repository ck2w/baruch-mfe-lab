#!/usr/bin/anaconda3/bin/python

import numpy as np

def box(Nbox, Nwindow):
    "a box function"

    h = np.zeros(Nwindow)
    h[:Nbox] = 1.0 / Nbox

    return h

if __name__ == "__main__":
    print(box(4, 10))

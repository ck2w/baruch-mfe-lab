#!/usr/bin/anaconda3/bin/python

import numpy as np

def box_diff(Nbox, Nwindow):
    "a box differencer"

    h = np.zeros(Nwindow)
    h[:Nbox] = 1.0 / Nbox;
    h[Nbox:2*Nbox] = -1.0 / Nbox;
    
    return h
    
if __name__ == "__main__":
    print(box_diff(4, 10))

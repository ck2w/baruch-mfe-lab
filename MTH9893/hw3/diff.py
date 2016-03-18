#!/usr/bin/anaconda3/bin/python

import numpy as np
from delta import delta

def diff(Neff_pos, Neff_neg, Nwindow):
    "a delta differencer with two arms - called Macd M1 in homework"

    h_pos = delta(Neff_pos, Nwindow)
    h_neg = delta(Neff_neg, Nwindow)
    h = h_pos - h_neg

    return h

if __name__ == "__main__":
    print(diff(3,7,10))



#!/usr/bin/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt

from step import step
from macd import macd

def compound_smooth_impulse_response(Neff_pos, Neff_neg, Nwindow):
    "a compound smoothing impulse response"

    h_step = step(Nwindow)
    h_macd = macd(Neff_pos, Neff_neg, Nwindow)

    gauge = 1.0 / (Neff_neg - Neff_pos)

    candidate = np.convolve(gauge*h_macd, h_step)
    return candidate[:Nwindow]

if __name__ == "__main__":
    
    f, axarr = plt.subplots(1, 1)
    Nwindow = 500
    Neff_pos = 7
    Neff_neg = 14
    h_step = step(Nwindow)
    h_macd = macd(Neff_pos, Neff_neg, Nwindow)
    h_csir = compound_smooth_impulse_response(Neff_pos, Neff_neg, Nwindow)
    axarr.plot(h_step)
    axarr.plot(h_macd)
    axarr.plot(h_csir)
    
    axarr.set_title('Unit step, macd, and compound impulse responses')

    plt.tight_layout()
    plt.show()

     

    

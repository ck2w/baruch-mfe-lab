#!/usr/bin/anaconda3/bin/python

import numpy as np

from ema import ema
from box import box

def calc_box_ema_spectra_mse(Neff, Nbox, Nwindow):
    "Mean-Squared Error\
     1) generates h_ema and h_box given the input parameters;\
     2) takes the amplitude of the FFT for each response\
     3) computes the cumulative sum of each gain spectrum\
     4) computes and retuns the MSE of the two cumulative gain spectra"

    h_box = box(Nbox, Nwindow)
    h_ema = ema(Neff, Nwindow)

    h_box_fft = np.fft.fft(h_box)
    h_ema_fft = np.fft.fft(h_ema)

    h_box_fft_abs = np.absolute(h_box_fft) 
    h_ema_fft_abs = np.absolute(h_ema_fft) 

    h_box_fft_abs_sum = np.cumsum(h_box_fft_abs)
    h_ema_fft_abs_sum = np.cumsum(h_ema_fft_abs)

    mse = np.sum(np.square(h_box_fft_abs_sum - h_ema_fft_abs_sum)) / Nwindow
   
    return mse 

    

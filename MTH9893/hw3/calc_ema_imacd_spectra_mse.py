#!/usr/bin/anaconda3/bin/python

import numpy as np

from ema import ema
from int_macd_poly import int_macd_poly

def calc_ema_imacd_spectra_mse(Neff_imacd, Neff_ema, Nwindow):
    "Mean-Squared Error\
     1) generates h_ema and h_imacd given the input parameters;\
     2) takes the amplitude of the FFT for each response\
     3) computes the cumulative sum of each gain spectrum\
     4) computes and retuns the MSE of the two cumulative gain spectra"

    h_ema = ema(Neff_ema, Nwindow)
    h_imacd = int_macd_poly(Neff_imacd, Nwindow)

    h_ema_fft = np.fft.fft(h_ema)
    h_imacd_fft = np.fft.fft(h_imacd)

    h_ema_fft_abs = np.absolute(h_ema_fft) 
    h_imacd_fft_abs = np.absolute(h_imacd_fft) 

    h_ema_fft_abs_sum = np.cumsum(h_ema_fft_abs)
    h_imacd_fft_abs_sum = np.cumsum(h_imacd_fft_abs)

    mse = np.sum(np.square(h_ema_fft_abs_sum - h_imacd_fft_abs_sum)) / Nwindow
   
    return mse 

    

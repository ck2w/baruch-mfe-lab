#!/usr/bin/anaconda3/bin/python

import numpy as np

from box_diff import box_diff
from macd import macd
from switch_to_composite_unit_gauge import switch_to_composite_unit_gauge

def calc_box_diff_macd_spectra_mse(Neff, Nbox, Nwindow):
    "Mean-Squared Error\
     1) generates h_box_diff and h_macd given the input parameters;\
     2) takes the amplitude of the FFT for each response\
     3) computes the cumulative sum of each gain spectrum\
     4) computes and retuns the MSE of the two cumulative gain spectra"

    h_box_diff = box_diff(Nbox, Nwindow)
    h_macd = switch_to_composite_unit_gauge(macd(Neff/3.0, Neff, Nwindow))
    #h_macd = macd(Neff, Nwindow)

    h_box_diff_fft = np.fft.fft(h_box_diff)
    h_macd_fft = np.fft.fft(h_macd)

    h_box_diff_fft_abs = np.absolute(h_box_diff_fft) 
    h_macd_fft_abs = np.absolute(h_macd_fft) 

    h_box_diff_fft_abs_sum = np.cumsum(h_box_diff_fft_abs)
    h_macd_fft_abs_sum = np.cumsum(h_macd_fft_abs)

    mse = np.sum(np.square(h_box_diff_fft_abs_sum - h_macd_fft_abs_sum)) / Nwindow
   
    return mse 

    

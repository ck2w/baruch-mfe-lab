#!/usr/bin/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
import csv

from delta import delta
from step import step
from box import box
from ema import ema
from ema_poly1 import ema_poly1
from macd import macd
from int_macd_poly import int_macd_poly

from apply_macd_filter import apply_macd_filter
from apply_ema_poly1_filter import apply_ema_poly1_filter
from apply_ema_filter import apply_ema_filter
from apply_box_filter import apply_box_filter
from apply_unit_step_filter import apply_unit_step_filter
from apply_delayed_impulse_filter import apply_delayed_impulse_filter 
from apply_integrated_macd_poly_filter import apply_integrated_macd_poly_filter


f, axarr = plt.subplots(3, 2)

with open('jpm_trades.csv', 'r') as trades_csv:
    trades = np.array(list(csv.reader(trades_csv))[1:]).astype('double')
    prices = trades[:,1]

    Nwindow = len(prices)
    
    # a) Delayed impulse

    Ndelay = 32
    lag = 0

    y_fde = apply_delayed_impulse_filter(prices, Ndelay)
    h_delta = delta(Ndelay, Nwindow)
    candidate = np.convolve(prices, h_delta)
    y_convolution = candidate[:len(prices)]

    x_axis = np.arange(Ndelay, len(prices))

    axarr[0, 0].set_title("Ideal delay")
    axarr[0, 0].plot(x_axis, y_fde[Ndelay:])
    axarr[0, 0].plot(x_axis, y_convolution[Ndelay:], 'o', markerfacecolor='none')
    
    # b) Box

    Nbox = 32
    lag = 1

    y_fde = apply_box_filter(prices, Nbox)
    h_box = box(Nbox, Nwindow)
    candidate = np.convolve(prices - prices[0], h_box) + prices[0]
    y_convolution = candidate[:len(prices)]

    axarr[0, 1].set_title("Box")
    axarr[0, 1].plot(y_fde)
    axarr[0, 1].plot(y_convolution, 'o', markerfacecolor='none')
    
    # c) Ema

    Neff = 32
    lag = 1

    y_fde = apply_ema_filter(prices, Neff)
    h_ema = ema(Neff, Nwindow)
    candidate = np.convolve(prices - prices[0], h_ema) + prices[0]
    y_convolution = candidate[:len(prices)]

    axarr[1, 0].set_title("Ema")
    axarr[1, 0].plot(y_fde)
    axarr[1, 0].plot(y_convolution, 'o', markerfacecolor='none')
    
    # d) Ema poly1

    Neff = 32
    lag = 2

    y_fde = apply_ema_poly1_filter(prices, Neff)
    h_ema_poly1 = ema_poly1(Neff, Nwindow)
    candidate = np.convolve(prices - prices[0], h_ema_poly1) + prices[0]
    y_convolution = candidate[:len(prices)]

    axarr[1, 1].set_title("Ema poly1")
    axarr[1, 1].plot(y_fde)
    axarr[1, 1].plot(y_convolution, 'o', markerfacecolor='none')

    # e) Integrated Macd-Poly    
    
    Neff = 32
    lag = 3

    y_fde = apply_integrated_macd_poly_filter(prices, Neff)
    h_int_macd_poly = int_macd_poly(Neff, Nwindow)
    candidate = np.convolve(prices - prices[0], h_int_macd_poly) + prices[0]
    y_convolution = candidate[:len(prices)]

    axarr[2, 0].set_title("Integrated macd poly")
    axarr[2, 0].plot(y_fde)
    axarr[2, 0].plot(y_convolution, 'o', markerfacecolor='none')

    # f) Macd

    Neff_pos = 16
    Neff_neg = 32
    lag = 2

    y_fde = apply_macd_filter(prices, Neff_pos, Neff_neg)
    h_macd = macd(Neff_pos, Neff_neg, Nwindow)
    candidate = np.convolve(prices - prices[0], h_macd)
    y_convolution = candidate[:len(prices)]

    axarr[2, 1].set_title("Macd")
    axarr[2, 1].plot(y_fde)
    axarr[2, 1].plot(y_convolution, 'o', markerfacecolor='none')
    
    
    
plt.tight_layout()
plt.show()
    

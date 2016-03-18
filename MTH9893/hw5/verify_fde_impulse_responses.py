#!/usr/bin/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt

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

Nwindow = 256

# a) Delayed impulse

Ndelay = 32
lag = 0

impulse = np.zeros(Nwindow)
impulse[lag] = 1

candidate = apply_delayed_impulse_filter(impulse, Ndelay)
impulse_response_fde = candidate[lag:]
impulse_response_direct = delta(Ndelay, Nwindow)

axarr[0, 0].set_ylim(0, 1.1)
axarr[0, 0].set_title("Ideal delay")
axarr[0, 0].plot(impulse_response_fde)
axarr[0, 0].plot(impulse_response_direct, 'o', markerfacecolor='none')

# b) Box

Nbox = 32
lag = 1

impulse = np.zeros(Nwindow)
impulse[lag] = 1

candidate = apply_box_filter(impulse, Nbox)
impulse_response_fde = candidate[lag:]
impulse_response_direct = box(Nbox, Nwindow)

axarr[0, 1].set_title("Box")
axarr[0, 1].plot(impulse_response_fde)
axarr[0, 1].plot(impulse_response_direct, 'o',markerfacecolor='none')

# c) Ema

Neff = 32
lag = 1

impulse = np.zeros(Nwindow)
impulse[lag] = 1

candidate = apply_ema_filter(impulse, Neff)
impulse_response_fde = candidate[lag:]
impulse_response_direct = ema(Neff, Nwindow)

axarr[1, 0].set_title("Ema")
axarr[1, 0].plot(impulse_response_fde)
axarr[1, 0].plot(impulse_response_direct, 'o', markerfacecolor='none')

# d) Ema poly1

Neff = 32
lag = 2 

impulse = np.zeros(Nwindow)
impulse[lag] = 1

candidate = apply_ema_poly1_filter(impulse, Neff)
impulse_response_fde = candidate[lag:]
impulse_response_direct = ema_poly1(Neff, Nwindow)

axarr[1, 1].set_title("Ema poly1")
axarr[1, 1].plot(impulse_response_fde)
axarr[1, 1].plot(impulse_response_direct, 'o', markerfacecolor='none')

# e) Integrated Macd-Poly 

Neff = 32
Neff_modified = 32
lag = 3 

impulse = np.zeros(Nwindow)
impulse[lag] = 1

candidate = apply_integrated_macd_poly_filter(impulse, Neff)
impulse_response_fde = candidate[lag:]
impulse_response_direct = int_macd_poly(Neff_modified, Nwindow)

axarr[2, 0].set_title("Integrated Macd-Poly")
axarr[2, 0].plot(impulse_response_fde)
axarr[2, 0].plot(impulse_response_direct, 'o', markerfacecolor='none')

# f) Macd

Neff_pos = 16
Neff_neg = 32
lag = 2 

impulse = np.zeros(Nwindow)
impulse[lag] = 1

candidate = apply_macd_filter(impulse, Neff_pos, Neff_neg)
impulse_response_fde = candidate[lag:]
impulse_response_direct = macd(Neff_pos, Neff_neg, Nwindow)

axarr[2, 1].set_title("Macd")
axarr[2, 1].plot(impulse_response_fde)
axarr[2, 1].plot(impulse_response_direct, 'o', markerfacecolor='none')


plt.tight_layout()
plt.show()

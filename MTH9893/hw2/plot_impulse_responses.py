#!/usr/bin/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt

from step import step
from box import box
from delta import delta
from ema import ema
from macd import macd
from diff import diff

Nwindow = 400
Neff_pos = 20
Neff_neg = 40
Neff = 20
Nbox = 20

h_step = step(Nwindow);
h_box = box(Nbox, Nwindow);
h_delta = delta(Neff_pos, Nwindow);
h_ema = ema(Neff_pos, Nwindow);
h_macd = macd(Neff_pos, Neff_neg, Nwindow);
h_diff = diff(Neff_pos, Neff_neg, Nwindow);

f, axarr = plt.subplots(3, 2)

axarr[0, 0].plot(h_step)
axarr[0, 0].set_title('step')

axarr[0, 1].plot(h_box)
axarr[0, 1].set_title('box')

axarr[1, 0].plot(h_delta)
axarr[1, 0].set_title('delta')

axarr[1, 1].plot(h_diff)
axarr[1, 1].set_title('diff - macd m1')

axarr[2, 0].plot(h_ema)
axarr[2, 0].set_title('ema')

axarr[2, 1].plot(h_macd)
axarr[2, 1].set_title('macd')

plt.tight_layout()
plt.show()



#!/usr/bin/anaconda3/bin/python

import numpy as np

def switch_to_composite_unit_gauge(h_in):
    "Convert input array to unit gauge"

    s = h_in[h_in>0].sum()
    return h_in/s

#!/usr/bin/anaconda3/bin/python

import numpy as np

from apply_ema_filter import apply_ema_filter

def compute_ab_ema(Neff, a, b):
    "Compute vwap with a and b series"

    ab = a * b
    
    a_fde = apply_ema_filter(a, Neff)
    b_fde = apply_ema_filter(b, Neff)
    ab_fde = apply_ema_filter(ab, Neff)

    a_star = ab_fde / b_fde

    b_fde_lag = np.concatenate(( np.zeros(1), b_fde ))[:len(b_fde)]
    
    p = Neff/(Neff+1)
    Neff_star = (p * b_fde_lag) / (b_fde - p * b_fde_lag)

    return (a_fde, b_fde, a_star, Neff_star)

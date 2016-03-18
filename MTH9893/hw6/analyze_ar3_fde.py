#!/usr/bin/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
import cmath

from partial_fraction_expansion import partial_fraction_expansion

# This problem is about analyzing an AR(3) equation
#
# y[n] = 0.0047 + 0.35 y[n-1] + 0.18 y[n-2] -0.14 y[n-3] + e[n]
#
# in the framework of finite difference equation
#
# The transformed equation is
#
# Q(z) Y(z) = c U(z) + E(z)
#
# where U(z) = 1 / (1-z^-1).
#

# (d) find poles

c = 0.0047
b = [1,0,0,0]
a = [1,-0.35,-0.18,0.14]
results = partial_fraction_expansion(b, a)
residues = results[0]
poles = results[1]

print("Poles:", poles)
print("Residues:", residues)

# (f-g) impulse reponses

f, axarr = plt.subplots(3, 1)

n = np.arange(-3, 20)

h_u = np.zeros(len(n))

for i in range(3):
    h_u = h_u + residues[i] * poles[i]**n

h_u = np.real(h_u)
axarr[0].set_title("h[n]")
axarr[0].plot(h_u)

u_convolved = np.convolve(h_u, np.ones(len(n)))[:len(h_u)]
axarr[1].set_title("h[n]*u[n]")
axarr[1].plot(u_convolved)

h_ar3 = h_u + c * u_convolved
axarr[2].plot(h_ar3)

# (h) solve finite difference equation
lag = 3
x_expanded = np.concatenate((np.zeros(lag), np.zeros(len(n))))
x_expanded[lag] = 1

expanded_length = len(x_expanded)

y_expanded = np.zeros(expanded_length)

y_expanded[0] = 0
y_expanded[1] = 0
y_expanded[2] = 0

for i in range(lag, expanded_length):
    y_expanded[i] = 0.0047 + 0.35 * y_expanded[i-1] + 0.18 * y_expanded[i-2] - 0.14 * y_expanded[i-3] + x_expanded[i]

y = y_expanded[1:len(h_u)]

axarr[2].plot(y, 'o', markerfacecolor='none')

# (i) compute gain: analytic and numerical

gain_numeric = np.sum(h_u)
gain_analytic = 1 / np.sum(a)

axarr[2].set_title("Impulse response - line, compared with FDE - open circle,"
                    "\n Analytic gain: " + str(round(gain_analytic, 6)) + ", numerical gain: " + str(round(gain_numeric, 6)))

plt.tight_layout()
plt.show()

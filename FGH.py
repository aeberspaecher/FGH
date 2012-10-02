#!/usr/bin/env python
#-*- coding:utf-8 -*-

# Use a minimal FGH implementation to compute eigenfumctions for bound 1d Hamiltonians

import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import h as H, hcont as HCont  # TODO: fix lowercase name
from scipy.linalg import eigh
from scipy.integrate import simps

L = 15  # length of x-interval to use
N = 2**10  # number of samples
NplotMin = 0  # quantum number of first eigenfunction to plot
Nplot = 7  # number of eigenfunctions to plot

# define a potential:
pot = lambda x: +0.5 * 3.0**2 * (x-L/2.)**2  # harmonic oscillator potential, centered
#pot = lambda x: -0.5 * 10.0**2 * (x-L/2.)**2 + (x-L/2.)**4  # double well potential

# define a grid:
x = np.linspace(0, L, N)

# sample potential:
Vsampled = pot(x)
Hsampled = H(L, Vsampled)

# diagonalize the Hamiltonian matrix:
E, psi = eigh(Hsampled)

WFscaleFactor = (np.max(Vsampled) - np.min(Vsampled))/Nplot
plt.plot(x, Vsampled, ls="-", c="k", lw=2)
for i in range(NplotMin, NplotMin+Nplot):
    # physically normalize WF
    WFnorm = np.sqrt(simps(np.abs(psi[:,i])**2, x=x))
    psi[:,i] /= WFnorm
    plt.plot(x, WFscaleFactor*np.abs(psi[:,i])**2 + E[i], ls="-", lw=1)
    print("E[%s] = %s"%(i, E[i]))

plt.show()

# TODO: plot potential *and* WF
# TODO: perform physical normalization of the WFs

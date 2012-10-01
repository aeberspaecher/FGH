#!/usr/bin/env python
#-*- coding:utf-8 -*-

# Use a minimal FGH implementation to compute eigenfumctions for bound 1d Hamiltonians

import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian import h as H  # TODO: fix lowercase name
from scipy.linalg import eigh

L = 10  # length of x-interval to use
N = 2**11  # number of samples
NplotMin = 0  # quantum number of first eigenfunction to plot
Nplot = 3  # number of eigenfunctions to plot

# define a potential:
pot = lambda x: 0.5 * 3.0**2 * (x-L/2.)**2 # harmonic oscillator potential, centered

# define a grid:
x = np.linspace(0, L, N)

# sample potential:
Vsampled = pot(x)
Hsampled = H(L, Vsampled)

# diagonalize the Hamiltonian
E, psi = eigh(Hsampled)

#plt.plot(x, Vsampled, c="k", ls="-", lw=2)
for i in range(NplotMin, NplotMin+Nplot):
    plt.plot(x, np.abs(psi[:,i])**2, ls="--", lw=2)
    print("E[%s] = %s"%(i, E[i]))

#plt.plot(x, np.abs((3/np.pi)**0.25 * np.exp(-3*((x-L/2.)**2)/2))**2)  # groundstate WF for omega=3
plt.show()

# TODO: plot potential and WF
# TODO: perform physical normalization of the WFs

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 11:46:39 2017

@author: alexandra
"""

import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import interp2d
import sys


plt.rcParams['axes.labelsize'] = 9
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.serif'] = ['Arial']

WIDTH = 510  # the number latex spits out
FACTOR = 0.6  # the fraction of the width you'd like the figure to occupy
fig_width_pt  = WIDTH * FACTOR
inches_per_pt = 1.0 / 72.27
golden_ratio  = (np.sqrt(5) - 1.0) / 2  # because it looks good
fig_width_in  = fig_width_pt * inches_per_pt  # figure width in inches
fig_height_in = fig_width_in * golden_ratio   # figure height in inches
fig_dims      = [fig_width_in, fig_height_in] # fig dims as a list


def bm_velocity(P, K0, K1, dx):
    U = np.copy(P)
    n = P.shape[0]
    for i in range(n):
        dp = np.gradient(P[i,:], dx)
        k = np.copy(dp)
        k[dp>=0] = K1
        k[dp<0] = K0
        U[i,:] = -k*dp
    return U


# Units
cm = 1e-2
um = 1e-4 * cm
dyn = 1
pa = 10 * dyn/cm**2
s = 1

# Scaled variables
r0 = 20*um
R = r0/r0
Le = 10*R
We = 0.2*R
Disp = 4*um/r0
w_ast = 10*um/r0
ast0 = Le/2-w_ast/2
gap = 1*um/r0
tf = 50*s
Y = 1.0e6 * pa
nu = 0.49
lam = Y*nu/((1+nu)*(1-2*nu))
mu = Y/(2*(1+nu))
Mu = mu/mu
Lam = lam/mu

sim = "oneast"

disp = np.loadtxt("./data/disp.csv", delimiter=',')

# import data
Nx = 50
Ny = 500
Nt = 200
ipad = np.zeros((Nt, Ny))
for i in range(0, Nt): 
    srr = np.loadtxt("./data/%s/stress.%d.csv" % (sim, i), delimiter=',',
                          skiprows=1)[:,0].reshape((Nx, Ny)).T

    # IPAD
    Pbm = -srr
    K0 = 1.0 * (1e-11 * cm**2)/(1.5e-3 * pa*s)
    K1 = 1.0 * (1e-11 * cm**2)/(1.5e-3 * pa*s)
    K0 = K0 / r0**2 * mu
    K1 = K1 / r0**2 * mu
    
    dx = Le/Pbm.shape[0]
    Ubm = bm_velocity(Pbm, K0, K1, dx) 
    
    ipad[i,:] = np.mean(Ubm, axis=1)

np.savetxt("./data/%s/ipad.csv" % (sim), ipad[:,30:-30]*r0, delimiter=",")

print(np.mean(ipad[30:-30])*r0/um*s)
    
# Mean vel over time
times = np.linspace(0, tf, 11)
x = np.linspace(0, Le*r0/um, ipad.shape[1])
f, axarr = plt.subplots(1, 1)
f.set_size_inches(w=fig_dims[0], h=fig_dims[1])
for t in range(0, len(times)):
    axarr.plot(x[30:-30], ipad[t,30:-30]*r0/um*s, label='%.1f s' % (times[t]))
#    axarr.plot(x, ipad[t,:]*r0/um**s, label='%.1f s' % (times[t]))
axarr.set_xlabel("z (μm)")
axarr.set_ylabel("IPAD (μm/s)")
axarr.legend(bbox_to_anchor=(1.0, 1.03))
#plt.savefig("./figures/ipad_oneast.png", dpi=600, bbox_inches='tight')
#plt.show()

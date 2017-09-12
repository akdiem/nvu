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


def initial_conditions(Ri):
    U = Ri * (0.1*um)
    return U


def F(U, ubm):
    return U*ubm
    

def S(U):
    return 0
    
    
def lax_wendroff(U_prev, F_prev, S_prev, u_bm, dt, dx):
    # u_prev = [U[m-1], U[m], U[m+1]], a_prev, p_prev analogously
    U_np_mp = (U_prev[2]+U_prev[1])/2 + dt/2 * (-(F_prev[2]-F_prev[1])/dx +\
                (S_prev[2]+S_prev[1])/2)
    U_np_mm = (U_prev[1]+U_prev[0])/2 + dt/2 * (-(F_prev[1]-F_prev[0])/dx +\
                (S_prev[1]+S_prev[0])/2)
                
    F_np_mp = F(U_np_mp, u_bm)
    F_np_mm = F(U_np_mm, u_bm)
    S_np_mp = S(U_np_mp)
    S_np_mm = S(U_np_mm)
    
    U_np = U_prev[1] - dt/dx * (F_np_mp-F_np_mm) + dt/2 * (S_np_mp+S_np_mm)
    return U_np


def numerical(U, ubm, time, dt, dx, x, L):
    for i in range(1,len(time)):
        # test cfl condition
        v = (max(U[i-1,:])) 

        # inlet boundary condition
        U[i,0] = U[0,0]
                
        for j in range(1,len(x)-1):
            u_prev = U[i-1,j-1:j+2]
            f_prev = u_prev * ubm[i-1,j-1:j+2]
            s_prev = np.array([0,0,0])
            if len(u_prev) == 2: # at the end of the array
                u_prev = U[i-1,j-1:]
                f_prev = u_prev * ubm[i-1,j-1:]
                s_prev = np.array([0,0,0])
            U[i,j] = lax_wendroff(u_prev, f_prev, s_prev, ubm[i,j], dt, dx)                                                                    
        # outlet boundary condition
        U[i,-1] = U[0,-1]
    
    return U


# Units
cm = 1e-2
um = 1e-4 * cm
dyn = 1
pa = 10 * dyn/cm**2
s = 1
mi = 60*s
ul = 1e9 * um**3

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

ipad = np.loadtxt("./data/%s/ipad.csv" % (sim), delimiter=',') / r0

Nt, Nx = ipad.shape
Ri = np.zeros((Nt, Nx))
for i in range(0, Nt):
    Ri[i,:] = np.loadtxt("./data/%s/disp_line.%d.csv" % (sim, i),
      delimiter=',', skiprows=1)[30:-31,0]
Ri = Ri/r0 + 0.5*We # location plus displacement
U = initial_conditions(Ri)
u0 = U[0,:]

x = np.linspace(0, Le, 500)[30:-30]
y = np.linspace(0, tf, Nt)
dt = x[1] - x[0]
dx = y[1] - y[0]
U = numerical(U, ipad, y, dt, dx, x, Le)

q = np.zeros((Nt, Nx))
for i in range(U.shape[0]):
    q[i,:] = 2 * np.pi * ipad[i,:] * U[i,:]
    
print(np.mean(q * r0**3/ul*mi))

# Mean vel over time
times = np.linspace(0, tf, 11)
x = np.linspace(0, Le*r0/um, ipad.shape[1])
f, axarr = plt.subplots(1, 1)
f.set_size_inches(w=fig_dims[0], h=fig_dims[1])
for t in range(0, len(times)):
    axarr.plot(x[30:-30], q[t,30:-30]*r0**3/ul*mi, label='%.1f s' % (times[t]))
axarr.set_xlabel("z (μm3)")
axarr.set_ylabel("IPAD (μl/min)")
#plt.savefig("./figures/flux_oneast.png", dpi=600, bbox_inches='tight')

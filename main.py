# -*- coding: utf-8 -*-

from nvu import nvu, utils

import numpy as np
from sys import argv


def init(r0):
    potassium_s = 2.92655044308714e-8
    ip3 = 5.37611796987610e-10
    calcium_a = 1.47220569018281e-07
    h = 0.404507631346124
    ss = 0.0161921297424289
    eet = 4.78801348065449e-07
    nbk = 6.24930194169376e-5
    Vk = -0.0814061063457068
    potassium_p = 0.00353809145071707
    calcium_p = 4.60269585230500e-06
    k = 8.01125818473096e-09
    Vm = 8.33004194103223e-05
    n = 0.283859572906570
    x = 2*np.pi*r0
    calcium_smc = 3.41385670857693e-07
    omega = 0.536911672725179
    yy = 0.000115089683436595
    amyloid = 3.75e-8
    return [potassium_s, ip3, calcium_a, h, ss, eet, nbk, Vk, potassium_p,
            calcium_p, k, Vm, n, x, calcium_smc, omega, yy, amyloid]
    
    
def K_glut_release(t1, t2, uM=0, s=0, **kwargs):
    sizeJrho = 1600
    sec = sizeJrho/(t2-t1)
    Max_neural_Kplus = 0.55*uM/s
    Max_neural_glut = 0.5
    Jrho_IN = np.zeros((sizeJrho,3))
    Jrho_IN[:,0] = np.linspace(t1, t2, sizeJrho)
    it1 = int(2.5*sec)
    it2 = int(sec)
    it3 = int(15*sec)
    it4 = int(0.25*sec)
    pos = it1
    Jrho_IN[pos+1:pos+it2+1,1] = Max_neural_Kplus * np.linspace(0, 1, it2)
    Jrho_IN[pos+1:pos+it2+1,2] = Max_neural_glut * np.linspace(0, 1, it2)
    pos += it2
    Jrho_IN[pos+1:pos+it3+1,1] = Max_neural_Kplus * np.ones(it3)
    Jrho_IN[pos+1:pos+it3+1,2] = Max_neural_glut * np.ones(it3)
    pos += it3
    Jrho_IN[pos+1:pos+it4+1,1] = Max_neural_Kplus * np.linspace(1, 0, it4)
    Jrho_IN[pos+1:pos+it4+1,2] = Max_neural_glut * np.linspace(1, 0, it4)
    return Jrho_IN


def main(fparam, fig_dims):
    units, param = utils.read_config(fparam)
    
    r0 = 20 * units['um']
    y0 = init(r0)
    x_rel = y0[13]
    sol = np.zeros(len(y0))

    # Equilibration
    t1 = -20
    t2 = 0
    nt = 100
    Jrho_IN = np.zeros((nt,3))
    Jrho_IN[:,0] = np.linspace(t1, t2, nt)
    t = np.linspace(t1, t2, nt)
    sol = nvu.run_simulation(t, y0, Jrho_IN, x_rel, units, param)
    y0 = sol[-1,:]
    
    # Plot solution
    nvu.plot_solution(t, sol, fig_dims, **units)
    
    # Simulation
    t1 = 0
    t2 = 50 
    nt = 200
    Jrho_IN = K_glut_release(t1, t2, **units)
    t = np.linspace(t1, t2, nt)    
    sol = nvu.run_simulation(t, y0, Jrho_IN, x_rel, units, param)
    
#    plt.figure(figsize=fig_dims)
#    plt.plot(t, sol[:,14]/uM, label="", lw=2)
#    plt.ylabel("Ca2+ smc (uM)")
#    plt.show()
    
    # Plot solution
    nvu.plot_solution(t, sol, fig_dims, **units)
    
    # Export radius data
#    r = sol[:,13]/(2*np.pi)
#    r_diff = (param['Sx']/param['Lx'])/2
#    Ra = r - r_diff
#    Rb = r + r_diff
#    np.savetxt('data/Ra.csv', Ra/units['um'], delimiter=',')
#    np.savetxt('data/Rb.csv', Rb/units['um'], delimiter=',')


if __name__ == "__main__":
    script, param = argv
    
    WIDTH = 510  # the number latex snp.pits out
    FACTOR = 1.0  # the fraction of the width you'd like the figure to occupy
    fig_width_pt  = WIDTH * FACTOR
    inches_per_pt = 1.0 / 72.27
    golden_ratio  = (np.sqrt(5) - 1.0) / 3  # because it looks good
    fig_width_in  = fig_width_pt * inches_per_pt  # figure width in inches
    fig_height_in = fig_width_in * golden_ratio   # figure height in inches
    fig_dims    = [fig_width_in, fig_height_in] # fig dims as a list
    
    main(param, fig_dims)
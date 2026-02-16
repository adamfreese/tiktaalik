"""
benchmarks.py

part of the tiktaalik package for GPD evolution
by Adam Freese

Methods to benchmark the accuracy of the finite element method.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from . import continuum, matrices, model, pars

mpl.rc('font',size=30,family='cmr10',weight='normal')
mpl.rc('text',usetex=True)
mpl.rc('text.latex', preamble=r"\usepackage{bm,amsmath,amssymb,amsfonts,mathrsfs}")
plt.rcParams["axes.formatter.use_mathtext"] = True

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Evolution shift benchmark

def shift_benchmark(key='NS', xi=0.1, nx=81, nlo=False, ns_type=1, grid_type=2):
    ''' Benchmark of the shift induced by a kernel on the right-hand side of the
    GPD evolution equation.
    ------
    Input:
        - key (string)
            'NS', 'qq', 'qg', 'gq', or 'gg'
            What kind of kernel is being benchmarked
        - xi (float)
            skewness value
            For xi<=0.1, grid_type=2 is recommended
            For xi>=0.1, grid_type=1 is recommended
        - nx (int)
            number of x points
            if grid_type=2, then using one more than a multiple of four will
            give x=xi and x=-xi on the grid
        - nlo (bool)
            True or NLO evolution, False for LO evolution
        - ns_type (int)
            +1 or -1
            relevant only for NLO evolution; whether it's a q+qbar or
            q-qbar non-singlet mixture that's being evolved.
        - grid_type (int)
            1 : linear x spacing
            2 : log-linear hybrid spacing
    '''
    # Set the grids
    matrices.set_x_xi_grids(nx, xi, grid_type)
    matrices.set_Q2_grid(np.array([4,5]))
    # Retrieve x grid
    x_pixel = matrices.get_x_grid()
    # Retrieve kernel matrix (interpixel method)
    K  = _get_kernel(key=key, ns_type=ns_type, nlo=nlo)
    # Interpixel shift
    H0 = _get_gpd(x_pixel, xi=xi, key=key)
    dH_pixel = np.einsum('ij,j->i', K[:,:,0], H0)
    # Continuum shift ("ground truth")
    x_truth = _make_continuum_x(xi)
    dH_truth = _get_continuum_shift(x_truth, key=key, xi=xi, nlo=nlo, ns_type=ns_type)
    # Error
    dH_truth_2 = _get_continuum_shift(x_pixel, key=key, xi=xi, nlo=nlo, ns_type=ns_type)
    error = 100*abs(dH_pixel - dH_truth_2) / abs(dH_truth_2)
    # Set up the plot
    nrows, ncols = 2, 1
    fig, (ax1, ax2) = plt.subplots(
            nrows, ncols,
            gridspec_kw={'height_ratios': [3,1]},
            figsize=(8,8),
            layout = 'constrained'
            )
    ax1.plot(x_truth, dH_truth, '-', label=r'Truth',     color='xkcd:lavender')
    ax1.plot(x_pixel, dH_pixel, '+', label=r'tiktaalik', color='xkcd:forest green')
    # Error
    ax2.plot(x_pixel, error, '+', label=r'tiktaalik', color='xkcd:forest green')
    # Plot labels etc
    ax2.set_ylim((1e-6,1e3))
    ax2.plot(x_truth, x_truth*0+1, linewidth=1, color='tab:gray')
    for ax in [ax1, ax2]:
        _plot_xi_lines(ax, xi)
        ax.set_xlim((-1,1))
        if(grid_type==2):
            ax.set_xscale('symlog', linthresh=xi)
    ax2.set_yscale('log')
    ax2.set_xlabel(r'$x$')
    ax1.get_xaxis().set_visible(False)
    ax1.set_ylabel(r'$\int \mathrm{d} y \, K(x,y,\xi) H(y)$')
    ax2.set_ylabel(r'percent error')
    # Finish
    fig.show()
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wilson coefficient benchmark

def wilson_benchmark(
        key='q',
        nx = 101,
        t = 0,
        grid_type = 2,
        nlo = False
        ):
    # Set the grids
    xi = np.geomspace(1e-4, 1, 100)
    matrices.set_Q2_grid(np.array([4,5]))
    matrices.set_x_xi_grids(nx, xi, grid_type)
    # Retrieve x grid
    x = matrices.get_x_grid()
    # Retrieve coefficient matrix (interpixel method)
    C = _get_dvcs_coefficient(key=key, nlo=nlo)
    # Interpiixel CFF
    Hq = _dvcs_quark_combo(x, xi, t)
    cff_pixel = np.einsum('ij,ji...->i...', C, Hq)[:,0]
    # Continuum CFF ("ground truth")
    cff_truth = _get_continuum_cff(xi, key=key, Q2=4, nlo=nlo)
    # Set up plot
    nrows, ncols = 2, 1
    fig, (ax1, ax2) = plt.subplots(
            nrows, ncols,
            gridspec_kw={'height_ratios': [3,1]},
            figsize=(8,8),
            layout = 'constrained'
            )
    ax1.plot(xi, xi*np.real(cff_truth), '-',  label=r'Truth     (real)', color='xkcd:electric blue')
    ax1.plot(xi, xi*np.real(cff_pixel), '+',  label=r'tiktaalik (real)', color='xkcd:forest green')
    ax1.plot(xi, xi*np.imag(cff_truth), '--', label=r'Truth     (imag)', color='xkcd:ochre')
    ax1.plot(xi, xi*np.imag(cff_pixel), 'x',  label=r'tiktaalik (imag)', color='xkcd:rich purple')
    # Error
    ImErr = 100*abs(np.imag(cff_truth-cff_pixel) / np.imag(cff_truth))
    ReErr = 100*abs(np.real(cff_truth-cff_pixel) / np.real(cff_truth))
    ax2.plot(xi, ReErr, '+', label=r'Error (real)', color='xkcd:forest green')
    ax2.plot(xi, ImErr, 'x', label=r'Error (imag)', color='xkcd:rich purple')
    # Finish
    for ax in [ax1, ax2]:
        ax.set_xscale('log')
    ax2.set_xlabel(r'$\xi$')
    ax1.set_ylabel(r'$\mathcal{H}(\xi)$')
    ax2.set_ylabel(r'percent error')
    _ = ax1.legend(prop = { 'size' : 17 }, loc=1)
    fig.show()
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Auxiliary functions

def _get_kernel(key='NS', nfl=4, ns_type=1, nlo=False):
    if(key=='NS'):
        K = matrices.kernel_VQQ(nfl=nfl, nlo=nlo, ns_type=ns_type)
    elif(key=='qq'):
        K = matrices.kernel_VQQ(nfl=nfl, nlo=nlo, ns_type=0)
    elif(key=='qg'):
        K = matrices.kernel_VQG(nfl=nfl, nlo=nlo)
    elif(key=='gq'):
        K = matrices.kernel_VGQ(nfl=nfl, nlo=nlo)
    elif(key=='gg'):
        K = matrices.kernel_VGG(nfl=nfl, nlo=nlo)
    else:
        raise ValueError("Key "+key+" unrecognized.")
    return K

def _get_dvcs_coefficient(key='q', nlo=False):
    if(key=='q'):
        C = matrices.dvcs_Cq(nlo=nlo)[:,:,0]
    elif(key=='g'):
        C = matrices.dvcs_Cg(nlo=nlo)[:,:,0]
    else:
        raise ValueError("Key "+key+" unrecognized.")
    return C

def _get_gpd(x, xi=0.1, key='NS'):
    H = np.zeros(x.shape)
    if(key=='NS'):
        H = _ns_gpd(x, xi)
    elif(key[1]=='g'):
        H = _gluon_gpd(x, xi)
    elif(key[1]=='q'):
        H = _singlet_gpd(x, xi)
    H[np.isnan(H)] = 0
    H[np.isinf(H)] = 0
    return H

def _get_continuum_shift(x, key='NS', xi=0.1, Q2=pars.mc2, nlo=False, ns_type=1):
    if(key=='NS'):
        continuum_shift = continuum.shift_cNS(x, xi, Q2, nlo, ns_type)[:,0]
    elif(key=='qq'):
        continuum_shift = continuum.shift_cQQ(x, xi, Q2, nlo)[:,0]
    elif(key=='qg'):
        continuum_shift = continuum.shift_cQG(x, xi, Q2, nlo)[:,0]
    elif(key=='gq'):
        continuum_shift = continuum.shift_cGQ(x, xi, Q2, nlo)[:,0]
    elif(key=='gg'):
        continuum_shift = continuum.shift_cGG(x, xi, Q2, nlo)[:,0]
    else:
        raise ValueError("Key "+key+" unrecognized.")
    return continuum_shift

def _get_continuum_cff(xi, key='q', Q2=4, nlo=False):
    if(key=='q'):
        continuum_cff = continuum.cff_q(xi, Q2, nlo=nlo)
    elif(key=='g'):
        continuum_cff = continuum.cff_g(xi, Q2, nlo=nlo)
    else:
        raise ValueError("Key "+key+" unrecognized.")
    return continuum_cff

def _ns_gpd(x, xi):
    T3 = model.Hu(x,xi,0) - model.Hd(x,xi,0)
    T3 = T3 - np.flip(T3, axis=0)
    return T3[:,0,0]

def _singlet_gpd(x, xi):
    return model.H_singlet(x,xi,0)[:,0,0]

def _gluon_gpd(x, xi):
    return model.Hg(x,xi,0)[:,0,0]

def _dvcs_quark_combo(x, xi, t):
    H = (
            4/9*(model.Hu(x,xi,t) - model.Hu(-x,xi,t))
            +
            1/9*(model.Hd(x,xi,t) - model.Hd(-x,xi,t))
            +
            1/9*(model.Hs(x,xi,t) - model.Hs(-x,xi,t))
            )
    return H

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot utilities

def _make_continuum_x(xi):
    N = 200
    x1 = np.geomspace( -1, -xi, N)
    x2 = np.linspace( -xi,  xi, N)
    x3 = np.geomspace( xi,   1, N)
    x = np.concatenate((x1[:N-1], x2[0:N], x3[1:]))
    return x

def _plot_xi_lines(ax, xi):
    ymin, ymax = ax.get_ylim()
    ax.vlines( xi, ymin, ymax, color='tab:gray', linewidth=1)
    ax.vlines(-xi, ymin, ymax, color='tab:gray', linewidth=1)
    ax.set_ylim((ymin,ymax))
    return

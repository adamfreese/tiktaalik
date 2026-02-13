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
    # Define a ground truth function
    gt_fun = _get_gt_fun(key=key, xi=xi, nlo=nlo, ns_type=ns_type)
    # Make the benchmark plot object
    symlog = False
    if(grid_type==2):
        symlog = True
    plot = shiftplot(gt_fun, xi=xi, symlog = symlog)
    # Get the shifts from tiktaalik
    matrices.set_x_xi_grids(nx, xi, grid_type)
    K  = _get_kernel(key=key, ns_type=ns_type, nlo=nlo)
    x = matrices.get_x_grid()[:,0]
    H0 = _gpd_key(x, xi=xi, key=key)
    dH = np.einsum('ij,j->i', K[:,:,0], H0)
    # Compare
    plot.plot_data(x, dH, label=r'tiktaalik')
    plot.show()
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Auxiliary functions

def _get_gt_fun(key='NS', xi=0.1, Q2=pars.mc2, nlo=False, ns_type=1):
    def gt_fun_NS(x):
        return continuum.shift_cNS(x, xi, Q2, nlo, ns_type)[:,0]
    def gt_fun_QQ(x):
        return continuum.shift_cQQ(x, xi, Q2, nlo)[:,0]
    def gt_fun_QG(x):
        return continuum.shift_cQG(x, xi, Q2, nlo)[:,0]
    def gt_fun_GQ(x):
        return continuum.shift_cGQ(x, xi, Q2, nlo)[:,0]
    def gt_fun_GG(x):
        return continuum.shift_cGG(x, xi, Q2, nlo)[:,0]
    if(key=='NS'):
        gt_fun = gt_fun_NS
    elif(key=='qq'):
        gt_fun = gt_fun_QQ
    elif(key=='qg'):
        gt_fun = gt_fun_QG
    elif(key=='gq'):
        gt_fun = gt_fun_GQ
    elif(key=='gg'):
        gt_fun = gt_fun_GG
    else:
        raise ValueError("Key "+key+" unrecognized.")
    return gt_fun

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

def _gpd_key(x, xi=0.1, key='NS'):
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

def _ns_gpd(x, xi):
    T3 = model.Hu(x,xi,0) - model.Hd(x,xi,0)
    T3 = T3 - np.flip(T3, axis=0)
    return T3[:,0,0]

def _singlet_gpd(x, xi):
    return model.H_singlet(x,xi,0)[:,0,0]

def _gluon_gpd(x, xi):
    return model.Hg(x,xi,0)[:,0,0]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A plot object for accuracy demos

class shiftplot:
    ''' Class for making benchmark plots of evolution shifts. '''

    def __init__(self, gt_fun, xi=0.1, symlog=False):
        nrows,ncols=2,1
        self.fig, (self.ax1, self.ax2) = plt.subplots(
                nrows, ncols,
                gridspec_kw={'height_ratios': [3,1]},
                figsize=(8,8),
                layout = 'constrained'
                )
        self.gt_fun = gt_fun
        #
        self.xi = xi
        self.symlog = symlog
        self.plot_baselines()
        self.plot_gt()
        return

    def plot_gt(self):
        N = 200
        x1 = np.geomspace(      -1, -self.xi, N)
        x2 = np.linspace( -self.xi,  self.xi, N)
        x3 = np.geomspace( self.xi,        1, N)
        x = np.concatenate((x1[:N-1], x2[0:N], x3[1:]))
        H = self.gt_fun(x)
        self.ax1.plot(x, H, '-', linewidth=2, label=r'Ground truth', color='xkcd:true green')
        return

    def plot_baselines(self):
        x = np.linspace(-1, 1, 100)
        zero = np.zeros(x.shape)
        one  = zero + 1
        self.ax1.plot(x, zero, '-', color='tab:gray', linewidth=1)
        self.ax2.plot(x, one,  '-', color='tab:gray', linewidth=1)
        return

    def plot_xi_lines(self):
        for ax in [self.ax1, self.ax2]:
            ymin, ymax = ax.get_ylim()
            ax.vlines( self.xi, ymin, ymax, color='tab:gray', linewidth=1)
            ax.vlines(-self.xi, ymin, ymax, color='tab:gray', linewidth=1)
            ax.set_ylim((ymin,ymax))
        return

    def plot_data(self, x, H, label=None):
        self.ax1.plot(x, H, '.', label='tiktaalik', color='xkcd:rich purple')
        truth = self.gt_fun(x)
        error = 100*abs(H - truth) / (abs(truth))
        self.ax2.plot(x, error, '.', color='xkcd:rich purple')
        return

    def show(self):
        # Log scale for error plot
        self.ax2.set_yscale('log')
        self.ax2.yaxis.get_major_locator().numticks = 3
        self.ax2.set_ylim((1e-6,1e3))
        ymin, ymax = self.ax2.get_ylim()
        if(ymax > 1e3):
            self.ax2.set_ylim((ymin,1e3))
        # xi lines
        self.plot_xi_lines()
        # x and y labels, and title
        self.ax1.get_xaxis().set_ticklabels([])
        self.ax1.get_xaxis().set_visible(False)
        self.ax2.set_xlabel(r'$x$', fontsize=30)
        shift_text = r'$\int \mathrm{d} y \, K(x,y,\xi) H(y)$'
        self.ax1.set_ylabel(shift_text, fontsize=30)
        self.ax2.set_ylabel(r'error (\%)', fontsize=34)
        # legend
        _ = self.ax1.legend(prop = { 'size' : 24 })
        # Do we scale it...?
        if(self.symlog):
            self.ax1.set_xscale('symlog', linthresh=self.xi)
            self.ax2.set_xscale('symlog', linthresh=self.xi)
        self.ax1.set_xlim((-1,1))
        self.ax2.set_xlim((-1,1))
        # the end
        self.fig.show()
        return

import numpy as np
import tiktaalik as tk

# Plotting things
import matplotlib as mpl
import matplotlib.pyplot as py
import matplotlib.ticker as mticker
mpl.rc('font',size=26,family='cmr10',weight='normal')
mpl.rc('text',usetex=True)
mpl.rc('text.latex', preamble=r"\usepackage{bm,amsmath,amssymb,amsfonts,mathrsfs}")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot demos

def plot_LO(nx=100, full_xi_list=False):
    if(full_xi_list):
        xi_array = np.linspace(0.5e-5, 5e-5, 10)
    else:
        xi_array = np.array([5e-6, 4e-5, 5e-5])
    #
    #
    nrows, ncols = 1, 1
    fig = py.figure(figsize=(ncols*8,nrows*6),layout='constrained')
    ax = py.subplot(nrows,ncols,1)
    #
    xmin = 1
    xmax = 1
    ymin = 1e-2
    ymax = ymin
    #ymin = 0
    for xi in xi_array:
        x = tk.matrices.pixelspace(nx, xi=xi, grid_type=2)
        dH = get_xi_shift(xi, nx=nx, nlo=False)
        ax.plot(x, dH, 'x', label=r'$\xi='+'{:.1f}'.format(xi*1e5)+r'\cdot10^{-5}$')
        # Update plot limits
        xmin = min(xmin, np.abs(x).min())
        #ymin = min(ymin, np.abs(dH).min())
        ymax = max(ymax, np.abs(dH).max())
    ymax *= 1.1
    xmin /= 1.1
    # Limit view to positive x, use log scale, to see what's going on better
    ax.set_xlim((xmin, xmax))
    ax.set_ylim((ymin, ymax))
    ax.set_xscale('log')
    ax.set_yscale('log')
    #
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'Single shift')
    _ = ax.legend(prop = { 'size' : 17 }, loc=3, ncol=2)
    fig.savefig('derp.pdf')
    fig.savefig('small_xi_LO_tk.pdf')
    fig.savefig('small_xi_LO_tk.png')
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Now add in continuum stuff
    x = np.geomspace(xmin, xmax/1.1, 20)
    for xi in xi_array:
        dH = get_continuum_shift(x, xi, nlo=False)
        ax.plot(x, dH, '--')
    fig.savefig('derp.pdf')
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Methods to generate the needed data

def get_xi_shift(xi, nx=100, nlo=False):
    # Fixed
    key = 'NS'
    ns_type = 1
    grid_type = 2
    # Initialize kernels, make shift
    tk.matrices.initialize_kernels(nx, xi, grid_type=grid_type)
    K  = get_kernel(key=key, ns_type=ns_type, nlo=nlo)
    x = tk.matrices.pixelspace(nx, xi=xi, grid_type=grid_type)
    H0 = gpd_key(x, xi=xi, key=key)
    dH = np.einsum('ij,j->i', K[:,:,0], H0)
    return dH

def get_continuum_shift(x, xi, key='NS', Q2=tk.pars.mc2, nlo=False, ns_type=1):
    if(key=='NS'):
        shift = tk.testing.test_shift_cNS(x, xi, Q2, nlo, ns_type)[:,0]
    elif(key=='qq'):
        shift = tk.testing.test_shift_cQQ(x, xi, Q2, nlo)[:,0]
    elif(key=='qg'):
        shift = tk.testing.test_shift_cQG(x, xi, Q2, nlo)[:,0]
    elif(key=='gq'):
        shift = tk.testing.test_shift_cGQ(x, xi, Q2, nlo)[:,0]
    elif(key=='gg'):
        shift = tk.testing.test_shift_cGG(x, xi, Q2, nlo)[:,0]
    else:
        raise ValueError("Key "+key+" unrecognized.")
    return shift

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Auxiliary functions

def get_kernel(key='NS', nfl=4, ns_type=1, nlo=False):
    if(key=='NS'):
        K = tk.matrices.kernel_VQQ(nfl=nfl, nlo=nlo, ns_type=ns_type)
    elif(key=='qq'):
        K = tk.matrices.kernel_VQQ(nfl=nfl, nlo=nlo, ns_type=0)
    elif(key=='qg'):
        K = tk.matrices.kernel_VQG(nfl=nfl, nlo=nlo)
    elif(key=='gq'):
        K = tk.matrices.kernel_VGQ(nfl=nfl, nlo=nlo)
    elif(key=='gg'):
        K = tk.matrices.kernel_VGG(nfl=nfl, nlo=nlo)
    else:
        raise ValueError("Key "+key+" unrecognized.")
    return K

def gpd_key(x, xi=0.5, key='NS'):
    H = np.zeros(x.shape)
    if(key=='NS'):
        H = ns_gpd(x, xi)
    elif(key[1]=='g'):
        H = gluon_gpd(x, xi)
    elif(key[1]=='q'):
        H = singlet_gpd(x, xi)
    H[np.isnan(H)] = 0
    H[np.isinf(H)] = 0
    # May need to interpolate to fix x=xi and x=-xi if using Marco's grid
    xidist = np.min(abs(x-xi))
    if(xidist < 1e-6):
        iA = (np.abs(x+xi)).argmin()
        iB = (np.abs(x-xi)).argmin()
        H[iA] = 0.5*(H[iA-1] + H[iA+1])
        H[iB] = 0.5*(H[iB-1] + H[iB+1])
    return H

def ns_gpd(x, xi):
    T3 = tk.model.Hu(x,xi,0) - tk.model.Hd(x,xi,0)
    T3 = T3 - np.flip(T3, axis=0)
    return T3[:,0,0]

def singlet_gpd(x, xi):
    return tk.model.H_singlet(x,xi,0)[:,0,0]

def gluon_gpd(x, xi):
    return tk.model.Hg(x,xi,0)[:,0,0]

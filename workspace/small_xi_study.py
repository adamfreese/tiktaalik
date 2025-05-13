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
# Global plot properties

lins = ['-','--', '-.', ':']
dots   = ['*', 'x', '+', '.']
cols = ['xkcd:forest green', 'xkcd:rich purple', 'xkcd:ochre', 'xkcd:electric blue']

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot demos

def test_gk_plot(
        xi_array=np.array([1e-1, 1e-2, 1e-3, 1e-4]),
        nx=100,
        key='NS'
        ):
    N = xi_array.shape[0]
    # Set up plot canvas
    nrows, ncols = 1, 1
    fig = py.figure(figsize=(ncols*8,nrows*6),layout='constrained')
    ax = py.subplot(nrows,ncols,1)
    for n in range(N):
        # Compute things to be plotted
        xi = xi_array[n]
        x = tk.matrices.pixelspace(nx, xi=xi, grid_type=2)
        H = gpd_key(x, xi=xi, key=key)
        ax.plot(x, H, dots[n],
                color=cols[n],
                label=r'$\xi='+'{:.1e}$'.format(xi)
                )
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'Single shift')
    _ = ax.legend(prop = { 'size' : 17 }, loc=1)
    fig.savefig('derp.pdf')
    return

def plot_multi_xi(
        key='NS',
        xi_array=np.array([1e-3, 1e-4, 1e-5, 3e-6]), # LO limit
        #xi_array=np.array([1e-3, 1e-4, 2e-5]), # NLO limit
        nx=100,
        nlo=False, continuum=False
                  ):
    N = xi_array.shape[0]
    # Set up plot canvas
    nrows, ncols = 1, 1
    fig = py.figure(figsize=(ncols*8,nrows*6),layout='constrained')
    ax = py.subplot(nrows,ncols,1)
    # Set up plot limits
    xmin = 1
    xmax = 1
    ymin = 1e-2
    ymax = ymin
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # The finite element stuff
    for n in range(N):
        # Compute things to be plotted
        xi = xi_array[n]
        x = tk.matrices.pixelspace(nx, xi=xi, grid_type=2)
        dH = get_xi_shift(xi, nx=nx, nlo=nlo, key=key)
        # Help with QQ...?
        if(key[0]=='q'):
            dH *= x
        ax.plot(x, dH, dots[n],
                color=cols[n],
                label=r'$\xi='+'{:.1e}$'.format(xi)
                )
        # Update plot limits
        xmin = min(xmin, np.abs(x).min())
        ymax = max(ymax, np.abs(dH).max())
    ymax *= 1.1
    xmin /= 1.1
    # Limit view to positive x, use log scale, to see what's going on better
    # (Only do this for NS I guess?)
    ax.set_xlim((xmin, xmax))
    ax.set_xscale('log')
    ##ax.set_ylim((ymin, ymax))
    ##ax.set_yscale('log')
    # 
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'Single shift')
    _ = ax.legend(prop = { 'size' : 17 }, loc=3)
    # Vertical line test
    #if(nlo):
    #    ymin, ymax = ax.get_ylim()
    #    ax.vlines(2e-3, ymin, ymax, color='tab:gray', linewidth=1)
    #    ax.set_ylim((ymin,ymax))
    # Save
    fig.savefig('derp.pdf')
    if(nlo):
        fig.savefig('small_xi_NLO_tk.pdf')
        fig.savefig('small_xi_NLO_tk.png')
    else:
        fig.savefig('small_xi_LO_tk.pdf')
        fig.savefig('small_xi_LO_tk.png')
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Now add in continuum stuff
    if(continuum):
        x = np.geomspace(xmin, xmax/1.1, 666)
        for n in range(N):
            xi = xi_array[n]
            dH = get_continuum_shift(x, xi, nlo=nlo)
            ax.plot(x, dH, lins[n],
                    color=cols[n],
                    )
        fig.savefig('derp.pdf')
        if(nlo):
            fig.savefig('small_xi_NLO_ct.pdf')
            fig.savefig('small_xi_NLO_ct.png')
        else:
            fig.savefig('small_xi_LO_ct.pdf')
            fig.savefig('small_xi_LO_ct.png')
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Methods to generate the needed data

def get_xi_shift(xi, nx=100, nlo=False, key='NS'):
    # Fixed
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
    # My GK model code becomes unstable for sufficiently small xi
    # Just replace with 1e-3 if it's smaller than this
    if(xi < 5e-3):
        xi = 5e-3
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

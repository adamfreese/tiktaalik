import numpy as np
import tiktaalik as tk

# Plotting things
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cmasher as cmr

mpl.rc('font',size=30,family='cmr10',weight='normal')
mpl.rc('text',usetex=True)
mpl.rc('text.latex', preamble=r"\usepackage{bm,amsmath,amssymb,amsfonts,mathrsfs}")
plt.rcParams["axes.formatter.use_mathtext"] = True

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Interpixel demo

def interpixel_demo_plots():
    n_pixels = 41
    i_highl = 10
    nx = 666
    xi = 0.3
    grid_type = 2
    # Set the grids used by tiktaalik
    tk.matrices.set_x_xi_grids(n_pixels, xi, grid_type=grid_type, lagrange_order=5)
    # Define a ground truth function using the GK model
    x_truth = np.linspace(-1, 1, nx)
    y_truth = tk.model.Hu(x_truth, xi, 0)[:,0,0]
    # Pixelation of the ground truth
    x_pixel = tk.matrices.get_x_grid()[:,0]
    y_pixel = tk.model.Hu(x_pixel, xi, 0)[:,0,0]
    # Build up interpixels
    y_inter = np.zeros((n_pixels, nx))
    for i_pixel in range(n_pixels):
        y_inter[i_pixel,:] = tk.matrices.interpixel(i_pixel, x_truth, xi)
        y_inter[i_pixel,:] *= y_pixel[i_pixel]
    y_recon = np.sum(y_inter, axis=0)
    # Multiple plot approach
    nrows,ncols=1,1
    figs = []
    axes = []
    for n in range(8):
        figs += [ plt.figure(figsize=(ncols*8,nrows*6), layout='constrained') ]
        axes += [ plt.subplot(nrows,ncols,1) ]
    # 1. Ground truth only
    axes[0].plot(x_truth, y_truth, '-', linewidth=1.9, color='xkcd:rich purple')
    # 2. Pixels on top (ground truth + pixels)
    axes[1].plot(x_truth, y_truth, '-', linewidth=1.9, color='xkcd:rich purple')
    axes[1].plot(x_pixel, y_pixel, 'o', color='xkcd:forest green')
    # 3. Remove ground truth (pixels only)
    axes[2].plot(x_pixel, y_pixel, 'o', color='xkcd:forest green')
    # 4. Interpixels on top (pixels + interpixels)
    for i_pixel in range(n_pixels):
        axes[3].plot(x_truth, y_inter[i_pixel,:], ':', linewidth=1.3, color='xkcd:light purple')
    axes[3].plot(x_truth, y_inter[i_highl,:], '-', linewidth=1.6, color='xkcd:blue')
    axes[3].plot(x_pixel, y_pixel, 'o', color='xkcd:forest green')
    # 5. Remove pixels (interpixels only)
    for i_pixel in range(n_pixels):
        axes[4].plot(x_truth, y_inter[i_pixel,:], ':', linewidth=1.3, color='xkcd:light purple')
    axes[4].plot(x_truth, y_inter[i_highl,:], '-', linewidth=1.6, color='xkcd:blue')
    # 6. Reconstruction (interpixels + reconstruction)
    for i_pixel in range(n_pixels):
        axes[5].plot(x_truth, y_inter[i_pixel,:], ':', linewidth=1.3, color='xkcd:light purple')
    axes[5].plot(x_truth, y_inter[i_highl,:], '-', linewidth=1.6, color='xkcd:blue')
    axes[5].plot(x_truth, y_recon, '--', linewidth=1.3, color='xkcd:lawn green')
    # 7. Remove interpixels (reconstruction only)
    axes[6].plot(x_truth, y_recon, '--', linewidth=1.3, color='xkcd:lawn green')
    # 8. Add ground truth (reconstruction + ground truth)
    axes[7].plot(x_truth, y_truth, '-', linewidth=1.9, color='xkcd:rich purple')
    axes[7].plot(x_truth, y_recon, '--', linewidth=1.3, color='xkcd:lawn green')
    ymax = 6.7
    # Post-processing for all plots
    for n in range(8):
        # Lines for zero and skewness
        axes[n].vlines( xi, -ymax, ymax, color='tab:gray', linewidth=1)
        axes[n].vlines(-xi, -ymax, ymax, color='tab:gray', linewidth=1)
        axes[n].plot(x_truth, x_truth*0, linewidth=1, color='tab:gray')
        # Labels
        axes[n].set_xlabel(r'$x$')
        # Consistent x and y ranges
        axes[n].set_xlim((-1,1))
        axes[n].set_ylim((-ymax,ymax))
        # Save figures
        figs[n].patch.set_alpha(0)
        figs[n].savefig('interpixeldemo{:d}.pdf'.format(n+1))
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Visualization of matrix evolution

def evomatrix_demo():
    n_pixels = 17
    xi = 0.3
    grid_type = 2
    # Set the grids used by tiktaalik, and evolution order
    tk.matrices.set_x_xi_grids(n_pixels, xi, grid_type=grid_type, lagrange_order=5)
    tk.matrices.set_Q2_grid(np.geomspace(4,17,11))
    tk.matrices.do_nlo_evolution()
    # Create initial and final GPD, along with matrix
    x = tk.matrices.get_x_grid()
    dx = x[1] - x[0]
    x_extra = np.append(x, x[-1] + dx)
    H0 = tk.model.Hu(x, xi, 0)[:,:,0] + tk.model.Hu(-x, xi, 0)[:,:,0]
    M = tk.matrices.matrix_VNS(ns_type=-1)
    HQ = np.einsum('xyzQ,yz->xzQ', M, H0)[:,:,-1]
    # Three figures
    fig0 = plt.figure(figsize=(1,n_pixels), layout='constrained')
    ax0 = plt.subplot(1,1,1)
    figQ = plt.figure(figsize=(1,n_pixels), layout='constrained')
    axQ = plt.subplot(1,1,1)
    figM =plt.figure(figsize=(n_pixels,n_pixels), layout='constrained')
    axM = plt.subplot(1,1,1)
    # Normalization for color maps
    Hmax = max(abs(H0).max(), abs(HQ).max())
    normM = mpl.colors.SymLogNorm(linthresh=1e-2, linscale=1, vmin=-1, vmax=1)
    # Plot initial and final GPDs as vertical pixel plots
    c0 = ax0.pcolormesh(np.array([0,1]), x_extra, H0, vmin=-Hmax, vmax=Hmax, cmap=cmr.iceburn, shading='flat')
    cQ = axQ.pcolormesh(np.array([0,1]), x_extra, HQ, vmin=-Hmax, vmax=Hmax, cmap=cmr.iceburn, shading='flat')
    cM = axM.pcolormesh(x_extra, x_extra, M[:,:,0,-1], norm=normM, cmap=cmr.iceburn, shading='flat')
    # Post-processing
    for ax in [ax0, axQ, axM]:
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
    for fig in [fig0, figQ, figM]:
        fig.patch.set_alpha(0)
    # Save plots
    fig0.savefig('evomatdemo0.pdf')
    figQ.savefig('evomatdemoQ.pdf')
    figM.savefig('evomatdemoM.pdf')
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Benchmark plots

def accuracy_benchmarks():
    xi_lin = 0.5
    xi_log = 1e-3
    nx = 81
    # NS
    fig_ns_lin = tk.benchmarks.shift_benchmark(
            key = 'NS',
            xi = xi_lin,
            nx = nx,
            nlo = True,
            ns_type = 1,
            grid_type = 1
            )
    fig_ns_lin.savefig('benchmark_ns_lin.pdf')
    fig_ns_log = tk.benchmarks.shift_benchmark(
            key = 'NS',
            xi = xi_log,
            nx = nx,
            nlo = True,
            ns_type = 1,
            grid_type = 2
            )
    fig_ns_log.savefig('benchmark_ns_log.pdf')
    # Quark from quark
    fig_qq_lin = tk.benchmarks.shift_benchmark(
            key = 'qq',
            xi = xi_lin,
            nx = nx,
            nlo = True,
            ns_type = 0,
            grid_type = 1
            )
    fig_qq_lin.savefig('benchmark_qq_lin.pdf')
    fig_qq_log = tk.benchmarks.shift_benchmark(
            key = 'qq',
            xi = xi_log,
            nx = nx,
            nlo = True,
            ns_type = 0,
            grid_type = 2
            )
    fig_qq_log.savefig('benchmark_qq_log.pdf')
    # Quark from gluon
    fig_qg_lin = tk.benchmarks.shift_benchmark(
            key = 'qg',
            xi = xi_lin,
            nx = nx,
            nlo = True,
            ns_type = 0,
            grid_type = 1
            )
    fig_qg_lin.savefig('benchmark_qg_lin.pdf')
    fig_qg_log = tk.benchmarks.shift_benchmark(
            key = 'qg',
            xi = xi_log,
            nx = nx,
            nlo = True,
            ns_type = 0,
            grid_type = 2
            )
    fig_qg_log.savefig('benchmark_qg_log.pdf')
    # Quark from gluon
    fig_gq_lin = tk.benchmarks.shift_benchmark(
            key = 'gq',
            xi = xi_lin,
            nx = nx,
            nlo = True,
            ns_type = 0,
            grid_type = 1
            )
    fig_gq_lin.savefig('benchmark_gq_lin.pdf')
    fig_gq_log = tk.benchmarks.shift_benchmark(
            key = 'gq',
            xi = xi_log,
            nx = nx,
            nlo = True,
            ns_type = 0,
            grid_type = 2
            )
    fig_gq_log.savefig('benchmark_gq_log.pdf')
    # Gluon from gluon
    fig_gg_lin = tk.benchmarks.shift_benchmark(
            key = 'gg',
            xi = xi_lin,
            nx = nx,
            nlo = True,
            ns_type = 0,
            grid_type = 1
            )
    fig_gg_lin.savefig('benchmark_gg_lin.pdf')
    fig_gg_log = tk.benchmarks.shift_benchmark(
            key = 'gg',
            xi = xi_log,
            nx = nx,
            nlo = True,
            ns_type = 0,
            grid_type = 2
            )
    fig_gg_log.savefig('benchmark_gg_log.pdf')
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Evolution comparisons

def make_nonlo_plots():
    ns_evolution_plots(xi=0.5, grid_type=1)
    ns_evolution_plots(xi=1e-3, grid_type=2)
    singlet_evolution_plots(xi=0.5, grid_type=1)
    singlet_evolution_plots(xi=1e-3, grid_type=2)
    return

def ns_evolution_plots(xi=1e-3, grid_type=2):
    # Create x grid
    nx = 81
    tk.matrices.set_x_xi_grids(nx, xi, grid_type)
    x = tk.matrices.get_x_grid()[:,0]
    # Create Q2 grid
    Q2 = np.geomspace(4, 17, 11)
    tk.matrices.set_Q2_grid(Q2)
    # Model scale GPD
    H0 = tk.model.Hu(x, xi, 0)[:,0,0] + tk.model.Hu(-x, xi, 0)[:,0,0]
    # Retrieve singlet evolution matrices
    tk.matrices.do_lo_evolution()
    M1 = tk.matrices.matrix_VNS(ns_type=-1)[:,:,0,:]
    tk.matrices.do_nlo_evolution()
    M2 = tk.matrices.matrix_VNS(ns_type=-1)[:,:,0,:]
    # Get evolved GPDs
    H1 = np.einsum('xyq,y->xq', M1, H0)[:,-1]
    H2 = np.einsum('xyq,y->xq', M2, H0)[:,-1]
    # Plots
    nrows, ncols = 2, 1
    figQ, (axQ1, axQ2) = plt.subplots(
            nrows, ncols,
            gridspec_kw={'height_ratios': [3,1]},
            figsize=(8,8),
            layout = 'constrained'
            )
    # Plot the GPD results
    axQ1.plot(x, H0, 'o', label=r'Initial', color='xkcd:forest green')
    axQ1.plot(x, H1, 'x', label=r'LO', color='xkcd:rich purple')
    axQ1.plot(x, H2, '+', label=r'NLO', color='xkcd:ochre')
    # Plot differences from initial
    axQ2.plot(x, H0-H0, 'o', color='xkcd:forest green')
    axQ2.plot(x, H1-H0, 'x', color='xkcd:rich purple')
    axQ2.plot(x, H2-H0, '+', color='xkcd:ochre')
    # Post-processing
    for ax in [axQ1, axQ2]:
        _plot_xi_lines(ax, xi)
        ax.set_xlim((-1,1))
        ax.plot(x, 0*x, color='tab:gray', linewidth=1)
        if(grid_type==2):
            ax.set_xscale('symlog', linthresh=xi)
    axQ1.set_ylabel(r'$H_u^-(x,\xi)$')
    for ax in [axQ2]:
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'Change')
    for ax in [axQ1]:
        ax.get_xaxis().set_visible(False)
        legend = ax.legend(prop = { 'size' : 26 })
        legend.get_frame().set_facecolor('#f8f8f8')
    for fig in [figQ]:
        fig.patch.set_alpha(0)
    # Save
    figQ.savefig('lonlo_ns_{:d}.pdf'.format(grid_type))
    return

def singlet_evolution_plots(xi=1e-3, grid_type=2):
    # Create x grid
    nx = 81
    tk.matrices.set_x_xi_grids(nx, xi, grid_type)
    x = tk.matrices.get_x_grid()[:,0]
    # Create Q2 grid
    Q2 = np.geomspace(4, 17, 11)
    tk.matrices.set_Q2_grid(Q2)
    # Model scale GPD
    Hq0 = tk.model.H_singlet(x, xi, 0)[:,0,0]
    Hg0 = tk.model.Hg(x, xi, 0)[:,0,0]
    HS0 = np.concat((Hq0, Hg0))
    # Retrieve singlet evolution matrices
    tk.matrices.do_lo_evolution()
    M1 = tk.matrices.matrix_VSG()[:,:,0,:]
    tk.matrices.do_nlo_evolution()
    M2 = tk.matrices.matrix_VSG()[:,:,0,:]
    # Get evolved GPDs
    HS1 = np.einsum('xyq,y->xq', M1, HS0)[:,-1]
    HS2 = np.einsum('xyq,y->xq', M2, HS0)[:,-1]
    Hq1 = HS1[0:nx]
    Hg1 = HS1[nx:2*nx]
    Hq2 = HS2[0:nx]
    Hg2 = HS2[nx:2*nx]
    # Quark and gluon plots
    nrows, ncols = 2, 1
    figQ, (axQ1, axQ2) = plt.subplots(
            nrows, ncols,
            gridspec_kw={'height_ratios': [3,1]},
            figsize=(8,8),
            layout = 'constrained'
            )
    figG, (axG1, axG2) = plt.subplots(
            nrows, ncols,
            gridspec_kw={'height_ratios': [3,1]},
            figsize=(8,8),
            layout = 'constrained'
            )
    # Plot the GPD results
    axQ1.plot(x, Hq0, 'o', label=r'Initial', color='xkcd:forest green')
    axQ1.plot(x, Hq1, 'x', label=r'LO', color='xkcd:rich purple')
    axQ1.plot(x, Hq2, '+', label=r'NLO', color='xkcd:ochre')
    axG1.plot(x, Hg0, 'o', label=r'Initial', color='xkcd:forest green')
    axG1.plot(x, Hg1, 'x', label=r'LO', color='xkcd:rich purple')
    axG1.plot(x, Hg2, '+', label=r'NLO', color='xkcd:ochre')
    # Plot differences from initial
    axQ2.plot(x, Hq0-Hq0, 'o', color='xkcd:forest green')
    axQ2.plot(x, Hq1-Hq0, 'x', color='xkcd:rich purple')
    axQ2.plot(x, Hq2-Hq0, '+', color='xkcd:ochre')
    axG2.plot(x, Hg0-Hg0, 'o', color='xkcd:forest green')
    axG2.plot(x, Hg1-Hg0, 'x', color='xkcd:rich purple')
    axG2.plot(x, Hg2-Hg0, '+', color='xkcd:ochre')
    # Post-processing
    for ax in [axQ1, axQ2, axG1, axG2]:
        _plot_xi_lines(ax, xi)
        ax.set_xlim((-1,1))
        ax.plot(x, 0*x, color='tab:gray', linewidth=1)
        if(grid_type==2):
            ax.set_xscale('symlog', linthresh=xi)
    axQ1.set_ylabel(r'$H_S(x,\xi)$')
    axG1.set_ylabel(r'$H_g(x,\xi)$')
    for ax in [axQ2, axG2]:
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'Change')
    for ax in [axQ1, axG1]:
        ax.get_xaxis().set_visible(False)
        legend = ax.legend(prop = { 'size' : 26 })
        legend.get_frame().set_facecolor('#f8f8f8')
    for fig in [figQ, figG]:
        fig.patch.set_alpha(0)
    # Save
    figQ.savefig('lonlo_q_{:d}.pdf'.format(grid_type))
    figG.savefig('lonlo_g_{:d}.pdf'.format(grid_type))
    return


def generic_NS_evolution(xi=0.1, nx=81, nlo=False, ns_type=1, grid_type=2):
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
    ax1.plot(x_truth, dH_truth, '-', label=r'Truth',     color='tab:orange')
    ax1.plot(x_pixel, dH_pixel, '+', label=r'tiktaalik', color='tab:blue')
    # Error
    ax2.plot(x_pixel, error, '+', label=r'tiktaalik', color='tab:blue')
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
    _ = ax1.legend(prop = { 'size' : 26 })
    fig.patch.set_alpha(0)
    return fig



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot utilities

def _plot_xi_lines(ax, xi):
    ymin, ymax = ax.get_ylim()
    ax.vlines( xi, ymin, ymax, color='tab:gray', linewidth=1)
    ax.vlines(-xi, ymin, ymax, color='tab:gray', linewidth=1)
    ax.set_ylim((ymin,ymax))
    return

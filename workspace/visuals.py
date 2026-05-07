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
    n_pixels = 17
    i_highl = 11
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

def _color(n):
    colors = [
            'tab:blue',
            'tab:orange',
            'tab:green',
            'tab:red',
            'tab:purple',
            'tab:brown',
            'tab:pink',
            'tab:olive',
            'tab:cyan'
            ]
    Nmax = len(colors) - 1
    if(n > Nmax):
        return _color(n-Nmax)
    return colors[n]

def interpixel_demo_alt():
    n_pixels = 13
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
    axes[0].plot(x_truth, y_truth, '-', linewidth=2, color='xkcd:rich purple')
    # 2. Pixels on top (ground truth + pixels)
    axes[1].plot(x_truth, y_truth, '-', linewidth=2, color='xkcd:rich purple')
    for i_pixel in range(n_pixels):
        axes[1].plot(x_pixel[i_pixel], y_pixel[i_pixel], 'o', color=_color(i_pixel), markeredgecolor='black')
    # 3. Remove ground truth (pixels only)
    for i_pixel in range(n_pixels):
        axes[2].plot(x_pixel[i_pixel], y_pixel[i_pixel], 'o', color=_color(i_pixel), markeredgecolor='black')
    # 4. Interpixels on top (pixels + interpixels)
    for i_pixel in range(n_pixels):
        axes[3].plot(x_truth, y_inter[i_pixel,:], '-', linewidth=1, color=_color(i_pixel))
    for i_pixel in range(n_pixels):
        axes[3].plot(x_pixel[i_pixel], y_pixel[i_pixel], 'o', color=_color(i_pixel), markeredgecolor='black')
    # 5. Remove pixels (interpixels only)
    for i_pixel in range(n_pixels):
        axes[4].plot(x_truth, y_inter[i_pixel,:], '-', linewidth=1, color=_color(i_pixel))
    # 6. Reconstruction (interpixels + reconstruction)
    for i_pixel in range(n_pixels):
        axes[5].plot(x_truth, y_inter[i_pixel,:], '-', linewidth=1, color=_color(i_pixel))
    axes[5].plot(x_truth, y_recon, '--', linewidth=2, color='black')
    # 7. Remove interpixels (reconstruction only)
    axes[6].plot(x_truth, y_recon, '--', linewidth=2, color='xkcd:lawn green')
    # 8. Add ground truth (reconstruction + ground truth)
    axes[7].plot(x_truth, y_truth, '-', linewidth=2, color='xkcd:rich purple')
    axes[7].plot(x_truth, y_recon, '--', linewidth=2, color='xkcd:lawn green')
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
    n_pixels = 33
    xi = 0.5
    grid_type = 1
    # Set the grids used by tiktaalik, and evolution order
    tk.matrices.set_x_xi_grids(n_pixels, xi, grid_type=grid_type, lagrange_order=5)
    tk.matrices.set_Q2_grid(np.geomspace(4,17,11))
    tk.matrices.do_nlo_evolution()
    # Create initial and final GPD, along with matrix
    x = tk.matrices.get_x_grid()
    dx = x[1] - x[0]
    x_extra = np.append(x, x[-1] + dx)
    H0 = tk.model.H3plus(x, xi, 0)[:,:,0]
    M = tk.matrices.matrix_VNS(ns_type=1)
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
    c0 = ax0.pcolormesh(np.array([0,1]), x_extra[::-1], H0, vmin=-Hmax, vmax=Hmax, cmap=cmr.fusion_r, shading='flat')
    cQ = axQ.pcolormesh(np.array([0,1]), x_extra[::-1], HQ, vmin=-Hmax, vmax=Hmax, cmap=cmr.fusion_r, shading='flat')
    cM = axM.pcolormesh(x_extra, x_extra[::-1], M[:,:,0,-1], norm=normM, cmap=cmr.fusion_r, shading='flat')
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
    ns_evolution_plots(xi=1e-4, grid_type=2)
    singlet_evolution_plots(xi=0.5, grid_type=1)
    singlet_evolution_plots(xi=1e-4, grid_type=2)
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
    H0 = tk.model.H3plus(x, xi, 0)[:,0,0]
    # Retrieve singlet evolution matrices
    tk.matrices.do_lo_evolution()
    M1 = tk.matrices.matrix_VNS(ns_type=1)[:,:,0,:]
    tk.matrices.do_nlo_evolution()
    M2 = tk.matrices.matrix_VNS(ns_type=1)[:,:,0,:]
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
    axQ1.plot(x, H0, 'o', label=r'Initial', color='tab:blue')
    axQ1.plot(x, H1, 'x', label=r'LO',      color='tab:orange')
    axQ1.plot(x, H2, '+', label=r'NLO',     color='tab:green')
    # Plot differences from initial
    axQ2.plot(x, H0-H0, 'o', color='tab:blue')
    axQ2.plot(x, H1-H0, 'x', color='tab:orange')
    axQ2.plot(x, H2-H0, '+', color='tab:green')
    # Post-processing
    for ax in [axQ1, axQ2]:
        _plot_xi_lines(ax, xi)
        ax.set_xlim((-1,1))
        ax.plot(x, 0*x, color='tab:gray', linewidth=1)
        if(grid_type==2):
            ax.set_xscale('symlog', linthresh=xi)
    axQ1.set_ylabel(r'$H_3^+(x,\xi)$')
    for ax in [axQ2]:
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'Change')
        if(grid_type==2):
            ax.set_xticks(ax.get_xticks()[::2])
    for ax in [axQ1]:
        ax.get_xaxis().set_visible(False)
        legend = ax.legend(prop = { 'size' : 26 })
        legend.get_frame().set_facecolor('#f8f8f8')
    for fig in [figQ]:
        fig.patch.set_alpha(0)
    bbox = dict(facecolor='#f8f8f8', alpha=0.76, edgecolor='gray', boxstyle='round,pad=0.2')
    axQ1.annotate(
            r'\textbf{NS}', xy=(0.03,0.05), xycoords='axes fraction',
            bbox=bbox
            )
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
    axQ1.plot(x, Hq0, 'o', label=r'Initial', color='tab:blue')
    axQ1.plot(x, Hq1, 'x', label=r'LO',      color='tab:orange')
    axQ1.plot(x, Hq2, '+', label=r'NLO',     color='tab:green')
    axG1.plot(x, Hg0, 'o', label=r'Initial', color='tab:blue')
    axG1.plot(x, Hg1, 'x', label=r'LO',      color='tab:orange')
    axG1.plot(x, Hg2, '+', label=r'NLO',     color='tab:green')
    # Plot differences from initial
    axQ2.plot(x, Hq0-Hq0, 'o', color='tab:blue')
    axQ2.plot(x, Hq1-Hq0, 'x', color='tab:orange')
    axQ2.plot(x, Hq2-Hq0, '+', color='tab:green')
    axG2.plot(x, Hg0-Hg0, 'o', color='tab:blue')
    axG2.plot(x, Hg1-Hg0, 'x', color='tab:orange')
    axG2.plot(x, Hg2-Hg0, '+', color='tab:green')
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
        if(grid_type==2):
            ax.set_xticks(ax.get_xticks()[::2])
    for ax in [axQ1, axG1]:
        ax.get_xaxis().set_visible(False)
        legend = ax.legend(prop = { 'size' : 26 })
        legend.get_frame().set_facecolor('#f8f8f8')
    for fig in [figQ, figG]:
        fig.patch.set_alpha(0)
    bbox = dict(facecolor='#f8f8f8', alpha=0.76, edgecolor='gray', boxstyle='round,pad=0.2')
    axQ1.annotate(
            r'\textbf{quark}', xy=(0.03,0.05), xycoords='axes fraction',
            bbox=bbox
            )
    axG1.annotate(
            r'\textbf{gluon}', xy=(0.03,0.05), xycoords='axes fraction',
            bbox=bbox
            )
    # Save
    figQ.savefig('lonlo_q_{:d}.pdf'.format(grid_type))
    figG.savefig('lonlo_g_{:d}.pdf'.format(grid_type))
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wilson demos

def wilson_demo():
    nx = 401
    fig_q_lo = tk.benchmarks.wilson_benchmark(
        key='q',
        nx = nx,
        t = 0,
        grid_type = 2,
        nlo = False
        )
    fig_q_lo.savefig('dvcs_q_lo.pdf')
    fig_q_nlo = tk.benchmarks.wilson_benchmark(
        key='q',
        nx = nx,
        t = 0,
        grid_type = 2,
        nlo = True
        )
    fig_q_nlo.savefig('dvcs_q_nlo.pdf')
    fig_g_nlo = tk.benchmarks.wilson_benchmark(
        key='g',
        nx = nx,
        t = 0,
        grid_type = 2,
        nlo = True
        )
    fig_g_nlo.savefig('dvcs_g_nlo.pdf')
    return



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot utilities

def _plot_xi_lines(ax, xi):
    ymin, ymax = ax.get_ylim()
    ax.vlines( xi, ymin, ymax, color='tab:gray', linewidth=1)
    ax.vlines(-xi, ymin, ymax, color='tab:gray', linewidth=1)
    ax.set_ylim((ymin,ymax))
    return

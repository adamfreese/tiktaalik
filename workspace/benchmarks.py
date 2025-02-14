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
# The tests

def shift_benchmark(key='NS', xi=0.5, nx=40, nlo=False, nstype=1):
    # Define a ground truth function
    gt_fun = get_gt_fun(key=key, xi=xi, nlo=nlo, nstype=nstype)
    # Make the benchmark plot object
    bm = bmplot(
            gt_fun,
            xi=xi,
            title = r'\textbf{'+key+r'}',
            ylabel = None,
            legend = True,
            filename = 'benchmark_' + key
            )
    tk.matrices.initialize_kernels(nx, xi)
    K  = get_kernel(key=key, nstype=nstype, nlo=nlo)
    x = tk.matrices.pixelspace(nx)
    H0 = gpd_key(x, xi=xi, key=key)
    dH = np.einsum('ij,j->i', K[:,:,0], H0)
    bm.plot_data(x, dH, label=r'$n_x={:d}$'.format(nx))
    # Finish
    bm.finish()
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Auxiliary functions

def get_gt_fun(key='NS', xi=0.5, Q2=tk.pars.mc2, nlo=False, nstype=1):
    def gt_fun_NS(x):
        return tk.testing.test_shift_cNS(x, xi, Q2, nlo, nstype)[:,0]
    def gt_fun_QQ(x):
        return tk.testing.test_shift_cQQ(x, xi, Q2, nlo)[:,0]
    def gt_fun_QG(x):
        return tk.testing.test_shift_cQG(x, xi, Q2, nlo)[:,0]
    def gt_fun_GQ(x):
        return tk.testing.test_shift_cGQ(x, xi, Q2, nlo)[:,0]
    def gt_fun_GG(x):
        return tk.testing.test_shift_cGG(x, xi, Q2, nlo)[:,0]
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

def get_kernel(key='NS', nfl=4, nstype=1, nlo=False):
    if(key=='NS'):
        K = tk.matrices.kernel_VQQ(nfl=nfl, nlo=nlo, nstype=nstype)
    elif(key=='qq'):
        K = tk.matrices.kernel_VQQ(nfl=nfl, nlo=nlo, nstype=0)
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A plot object for accuracy demos

class bmplot:

    def __init__(self,
            gt_fun, xi=0.5, title=None, ylabel=None, legend=False,
            axy=((0.04,0.04)),
            axy2=((0.04,0.92)),
            key = None, showkey = False,
            filename='derp'):
        nrows,ncols=2,1
        self.fig, (self.ax1, self.ax2) = py.subplots(
                nrows, ncols,
                gridspec_kw={'height_ratios': [3,1]},
                figsize=(8,8),
                layout = 'constrained'
                )
        self.gt_fun = gt_fun
        #
        self.n = 0
        self.xi = xi
        self.title = title
        self.ylabel = ylabel
        self.filename = filename
        self.legend = legend
        self.axy = axy
        self.axy2 = axy2
        self.key = key
        self.showkey = showkey
        #
        self.color = ('xkcd:rich purple', 'xkcd:forest green', 'mediumblue', 'tab:orange')
        self.fmt   = ('x', '+', '.', 's')
        self.size  = (6, 8, 8, 6)
        #
        self.plot_baselines()
        self.plot_gt()
        return

    def plot_gt(self):
        #x = np.linspace(-1, 1, 666)
        # TEMPORARY
        x1 = np.geomspace(      -1, -self.xi, 100)
        x2 = np.linspace( -self.xi,  self.xi, 100)
        x3 = np.geomspace( self.xi,        1, 100)
        x = np.concatenate((x1[:99], x2[1:99], x3[1:]))
        H = self.gt_fun(x)
        self.ax1.plot(x, H, '-', color='black', linewidth=2, label=r'Ground truth')
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
        self.ax1.plot(x, H, self.fmt[self.n], color=self.color[self.n], label=label, markersize=self.size[self.n])
        truth = self.gt_fun(x)
        error = 100*abs(H - truth) / (abs(truth) + 1e-15)
        self.ax2.plot(x, error, self.fmt[self.n], color=self.color[self.n], markersize=self.size[self.n])
        self.n += 1
        return

    def finish(self):
        # Log scale for error plot
        self.ax2.set_yscale('log')
        self.ax2.yaxis.get_major_locator().numticks = 3
        ##self.ax2.set_ylim((1e-6,1e3))
        ymin, ymax = self.ax2.get_ylim()
        if(ymax > 1e3):
            self.ax2.set_ylim((ymin,1e3))
        # xi lines
        self.plot_xi_lines()
        # x and y labels, and title
        ##self.ax1.set_xlabel(r'$x$', fontsize=27)
        #self.ax1.get_xaxis().set_visible(False)
        self.ax1.get_xaxis().set_ticklabels([])
        self.ax2.set_xlabel(r'$x$', fontsize=30)
        if(self.ylabel is not None):
            self.ax1.set_ylabel(self.ylabel, fontsize=30)
            self.ax2.set_ylabel(r'error (\%)', fontsize=34)
        if(self.axy is not None):
            self.ax1.annotate(self.title, xy=self.axy, xycoords="axes fraction",
                    bbox=dict(
                        facecolor='white', edgecolor='0.8',
                        boxstyle='round,pad=0.1',
                        alpha=0.8
                        ),
                    fontsize=30
                    )
        if(self.showkey):
            self.ax1.annotate(self.key, xy=self.axy2, xycoords="axes fraction",
                    bbox=dict(
                        facecolor='white', edgecolor='0.8',
                        boxstyle='round,pad=0.1',
                        alpha=0.8
                        ),
                    fontsize=30
                    )
        # legend
        if(self.legend):
            loc = 0
            if(self.xi<0.1):
                loc = 4
            _ = self.ax1.legend(prop = { 'size' : 24 }, loc=loc)
            markersizes = [1, 12, 14, 12, 14]
            for count, legend_handle in enumerate(self.ax1.get_legend().legend_handles):
                legend_handle.set(markersize = markersizes[count])
        # Temporary
        #self.ax1.set_xscale('symlog', linthresh=self.xi)
        #self.ax2.set_xscale('symlog', linthresh=self.xi)
        # the end
        #self.fig.patch.set_alpha(0)
        self.fig.savefig('derp.pdf')
        self.fig.savefig('derp.png')
        self.fig.savefig(self.filename+'.pdf')
        self.fig.savefig(self.filename+'.png')
        return

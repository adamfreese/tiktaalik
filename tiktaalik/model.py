"""
model.py

part of the tiktaalik package for GPD evolution
by Adam Freese

This implements a model due to Goloshkokov and Kroll.
The specific reference consulted is was:
  P. Kroll, H. Moutarde, F. Sabatie
  European Physical Journal C (2013) 73:2278
  arxiv:1210.6975
  Kroll:2012sm
see /f90src/model/gk.f90 for more details.
"""

import numpy as np
from .f90wrap.model import dummy as f90src

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model GPDs

def Hu(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    if(x.ndim==2):
        H = f90src.hu_wrap_2d(x,xi,t)
    else:
        H = f90src.hu_wrap(x,xi,t)
    return H

def Hd(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    if(x.ndim==2):
        H = f90src.hd_wrap_2d(x,xi,t)
    else:
        H = f90src.hd_wrap(x,xi,t)
    return H

def Hs(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    if(x.ndim==2):
        H = f90src.hs_wrap_2d(x,xi,t)
    else:
        H = f90src.hs_wrap(x,xi,t)
    return H

def Hg(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    if(x.ndim==2):
        H = f90src.hg_wrap_2d(x,xi,t)
    else:
        H = f90src.hg_wrap(x,xi,t)
    return H

def Hu_tilde(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    if(x.ndim==2):
        H = f90src.htu_wrap_2d(x,xi,t)
    else:
        H = f90src.htu_wrap(x,xi,t)
    return H

def Hd_tilde(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    if(x.ndim==2):
        H = f90src.htd_wrap_2d(x,xi,t)
    else:
        H = f90src.htd_wrap(x,xi,t)
    return H

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Helpful quark combinations

def H_singlet(x, xi, t):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    if(np.isscalar(t)):
        t = np.array([t])
    f = (
              Hu(x,xi,t) - Hu(-x,xi,t)
            + Hd(x,xi,t) - Hd(-x,xi,t)
            + Hs(x,xi,t) - Hs(-x,xi,t)
            )
    return f

def H3plus(x, xi, t):
    Hup = Hu(x, xi, t) - Hu(-x, xi, t)
    Hdp = Hd(x, xi, t) - Hd(-x, xi, t)
    return Hup - Hdp

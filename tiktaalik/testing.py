"""
testing.py

part of the tiktaalik package for GPD evolution
by Adam Freese

This file contains methods for assisting validation and making sure the code
works right.
"""

import numpy as np
import pandas
from .f90wrap.testing import dummy as f90src

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Continuum shifts

def test_shift_cNS(x, xi, Q2, nlo=False, nstype=1):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cns(x,xi,Q2,nlo,nstype)
    return V

def test_shift_cQQ(x, xi, Q2, nlo=False):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cqq(x,xi,Q2,nlo)
    return V

def test_shift_cQG(x, xi, Q2, nlo=False):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cqg(x,xi,Q2,nlo)
    return V

def test_shift_cGQ(x, xi, Q2, nlo=False):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cgq(x,xi,Q2,nlo)
    return V

def test_shift_cGG(x, xi, Q2, nlo=False):
    if(np.isscalar(x)):
        x = np.array([x])
    if(np.isscalar(xi)):
        xi = np.array([xi])
    V = f90src.test_shift_cgg(x,xi,Q2,nlo)
    return V

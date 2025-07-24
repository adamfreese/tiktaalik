"""
matrices.py

part of the tiktaalik package for GPD evolution
by Adam Freese

This file contains methods to construct evolution matrices.
"""

import numpy as np
from .f90wrap.matrices import dummy as f90src

from . import pars

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A global dictionary containing information about what's been initialized
# It's initialized with some sane default values

matrix_dict = {
        'nQ2' : 2,
        'Q2'  : np.array([pars.mc2, pars.mb2]),
        'nxi' : 5,
        'xi'  : np.linspace(0.1, 0.5, 5),
        'nx'  : 80,
        'grid_type' : 1,
        'lagrange_order' : 5,
        'kernel_init' : False,
        'evomat_init' : False,
        'wilson_init' : False,
        'nlo_evo' : False
        }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The autmatically-run initialization routine, called when tiktaalik starts

def first_initialization():
    # TODO: docstring
    f90src.initialize_x_xi_wrap(matrix_dict['nx'],
                                matrix_dict['xi'],
                                matrix_dict['grid_type'],
                                matrix_dict['lagrange_order']
                                )
    f90src.initialize_q2_wrap(matrix_dict['Q2'])
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Methods for the user to retrieve the x, xi and Q2 grids in use

def get_x_grid():
    ''' Returns the x grid in use.
    The returned grid is two-dimensional, with dimensions (nx,nxi).
    The reason for this is that, depending on the grid type,
    the x spacing may be xi-dependent.
    Thus, this is effectively an array of nxi x arrays.
    '''
    nx  = matrix_dict['nx']
    nxi = matrix_dict['nxi']
    x = f90src.get_x_wrap(nx, nxi)
    return x

def get_xi_array():
    nxi = matrix_dict['nxi']
    xi = f90src.get_xi_wrap(nxi)
    return xi

def get_Q2_array():
    nQ2 = matrix_dict['nQ2']
    Q2 = f90src.get_q2_wrap(nQ2)
    return Q2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Methods for the user to set x, xi and Q2 grids

def set_x_xi_grids(nx, xi, grid_type, lagrange_order=5):
    # TODO: docstring
    ''' TODO: dicstring '''
    # Assert requirements
    assert(nx==(nx//2)*2)
    assert(nx>=6)
    # Cache the data sent by the user
    matrix_dict['nx']  = nx
    if(np.isscalar(xi)):
        xi = np.array([xi])
    matrix_dict['nxi'] = xi.shape[0]
    matrix_dict['xi']  = xi
    matrix_dict['grid_type'] = grid_type
    # Set the x and xi grids in the Fortran code
    f90src.initialize_x_xi_wrap(nx, xi, grid_type, lagrange_order)
    # Mark all matrices as uninitialized, since their shape may now be invalid
    matrix_dict['kernel_init'] = False
    matrix_dict['evomat_init'] = False
    matrix_dict['wilson_init'] = False
    return

def set_Q2_grid(Q2):
    # TODO: docstring
    ''' TODO: dicstring '''
    # Assert requirements
    assert(Q2.shape[0] >= 2)
    # Cache the data sent by the user
    matrix_dict['nQ2'] = Q2.shape[0]
    matrix_dict['Q2']  = Q2
    # Set the Q2 grid in the Fortran code
    f90src.initialize_q2_wrap(Q2)
    # Mark matrices as uninitialized that depend on the Q2 shape/range
    matrix_dict['evomat_init'] = False
    matrix_dict['wilson_init'] = False
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Routines to obtain evolution kernel matrices

def kernel_VQQ(Q2=pars.mc2, nfl=4, nlo=False, ns_type=1):
    ''' Evolution kernel for Q->Q, helicity-independent (V-type).
    The ns_type parameter can be -1, 0, or 1. The +1 and -1 values are for
    plus-type non-singlet distributions such as (u+ubar)-(d+dbar),
    and minus-type non-singlet distributions such as u-ubar, respectively.
    The 0 value is for the singlet Q->Q kernel.
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO too.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    if(matrix_dict['kernel_init']==False):
        f90src.make_kernels_wrap()
        matrix_dict['kernel_init'] = True
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vqq_wrap(Q2, nx, nxi, nfl, nlo, ns_type)
    return K

def kernel_VQG(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for G->Q, helicity-independent (V-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    if(matrix_dict['kernel_init']==False):
        f90src.make_kernels_wrap()
        matrix_dict['kernel_init'] = True
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vqg_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_VGQ(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for Q->G, helicity-independent (V-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    if(matrix_dict['kernel_init']==False):
        f90src.make_kernels_wrap()
        matrix_dict['kernel_init'] = True
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vgq_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_VGG(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for G->G, helicity-independent (V-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    if(matrix_dict['kernel_init']==False):
        f90src.make_kernels_wrap()
        matrix_dict['kernel_init'] = True
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_vgg_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_AQQ(Q2=pars.mc2, nfl=4, nlo=False, ns_type=1):
    ''' Evolution kernel for Q->Q, helicity-dependent (A-type).
    The ns_type parameter can be -1, 0, or 1. The +1 and -1 values are for
    plus-type non-singlet distributions such as (u+ubar)-(d+dbar),
    and minus-type non-singlet distributions such as u-ubar, respectively.
    The 0 value is for the singlet Q->Q kernel.
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO too.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    if(matrix_dict['kernel_init']==False):
        f90src.make_kernels_wrap()
        matrix_dict['kernel_init'] = True
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_aqq_wrap(Q2, nx, nxi, nfl, nlo, ns_type)
    return K

def kernel_AQG(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for G->Q, helicity-dependent (A-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    if(matrix_dict['kernel_init']==False):
        f90src.make_kernels_wrap()
        matrix_dict['kernel_init'] = True
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_aqg_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_AGQ(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for Q->G, helicity-dependent (A-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    if(matrix_dict['kernel_init']==False):
        f90src.make_kernels_wrap()
        matrix_dict['kernel_init'] = True
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_agq_wrap(Q2, nx, nxi, nfl, nlo)
    return K

def kernel_AGG(Q2=pars.mc2, nfl=4, nlo=False):
    ''' Evolution kernel for G->G, helicity-dependent (A-type).
    The kernel depends on Q2 through alphaQCD, nfl, and whether it's NLO.
    '''
    assert(nfl==3 or nfl==4 or nfl==5)
    if(matrix_dict['kernel_init']==False):
        f90src.make_kernels_wrap()
        matrix_dict['kernel_init'] = True
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    K = f90src.evokernel_agg_wrap(Q2, nx, nxi, nfl, nlo)
    return K

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Routines to obtain evolution matrices
# TODO: initialize if uninitialized, or nlo truth value contradicts cached value

def matrix_VNS(ns_type=1):
    ''' Evolution matrix for *non-singlet* Q->Q, helicity-independent (V-type).
    The ns_type parameter must be 1 or -1, and distinguishes between
    plus-type non-singlet distributions such as (u+ubar)-(d+dbar),
    and minus-type non-singlet distributions such as u-ubar.
    The ns_type parameter only matters for NLO.
    '''
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    M = f90src.evomatrix_vns_wrap(nx, nxi, nQ2, ns_type)
    M[np.isnan(M)] = 0
    M[np.isinf(M)] = 0
    return M

def matrix_VSG():
    ''' Evolution matrix for singlet sector, helicity-independent (V-type).
    '''
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    M = f90src.evomatrix_vsg_wrap(nx, nxi, nQ2)
    M[np.isnan(M)] = 0
    M[np.isinf(M)] = 0
    return M

def matrix_ANS(ns_type=1):
    ''' Evolution matrix for *non-singlet* Q->Q, helicity-dependent (A-type).
    The ns_type parameter must be 1 or -1, and distinguishes between
    plus-type non-singlet distributions such as (u+ubar)-(d+dbar),
    and minus-type non-singlet distributions such as u-ubar.
    The ns_type parameter only matters for NLO.
    '''
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    M = f90src.evomatrix_ans_wrap(nx, nxi, nQ2, ns_type)
    M[np.isnan(M)] = 0
    M[np.isinf(M)] = 0
    return M

def matrix_ASG():
    ''' Evolution matrix for singlet sector, helicity-dependent (A-type).
    '''
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    M = f90src.evomatrix_asg_wrap(nx, nxi, nQ2)
    M[np.isnan(M)] = 0
    M[np.isinf(M)] = 0
    return M

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Routines to obtain Wilson coefficient matrices

def dvcs_Cq(nlo=False):
    # TODO: docstring
    # Get the matrix dimensions
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    # Make sure the Wilson coefficient matrices are initialized
    if(matrix_dict['wilson_init']==False):
        f90src.make_wilson_wrap(nQ2)
        matrix_dict['wilson_init'] = True
    # Get the coefficients
    C = f90src.dvcs_cq_wrap(nx, nxi, nQ2, nlo)
    return C

def dvcs_Cg(nlo=False):
    # TODO: docstring
    # Get the matrix dimensions
    nx  = f90src.get_nx_wrap()
    nxi = f90src.get_nxi_wrap()
    nQ2 = f90src.get_nq2_wrap()
    # Make sure the Wilson coefficient matrices are initialized
    if(matrix_dict['wilson_init']==False):
        f90src.make_wilson_wrap(nQ2)
        matrix_dict['wilson_init'] = True
    # Get the coefficients
    C = f90src.dvcs_cg_wrap(nx, nxi, nQ2, nlo)
    return C

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Q2 related methods

def Q2space(Q2i, Q2f, nQ2):
    ''' Creates something that's mostly a geomspace,
    but injects any quark mass thresholds if they're present.
    '''
    nThresh = int(Q2i < pars.mc2 and Q2f > pars.mc2) + int(Q2i < pars.mb2 and Q2f > pars.mb2)
    Q2 = np.geomspace(Q2i, Q2f, nQ2-nThresh)
    if(Q2i < pars.mc2 and Q2f > pars.mc2):
        Q2 = np.append(Q2, pars.mc2)
    if(Q2i < pars.mb2 and Q2f > pars.mb2):
        Q2 = np.append(Q2, pars.mb2)
    Q2 = np.sort(Q2)
    return Q2

def passes(Q2_array, m2):
    ''' Routine to check if m2 is located between the minimum and maximum
    values in Q2_array. Assumes Q2_array is monotonically ordered.
    Used to check whether a mass threshold is passed during evolution.
    '''
    if(Q2_array[0] < m2-pars.epsilon and Q2_array[-1] > m2+pars.epsilon):
        return True
    return False

def get_nfl(Q2):
    ''' Gets the effective number of flavors. Returns either 3, 4 or 5,
    under the assumption that Q2 below the strange threshold or above the
    truth threshold is never passed.
    This routine  is used by the Evolver class.
    '''
    if(Q2 < pars.mc2):
        return 3
    elif(Q2 < pars.mb2):
        return 4
    return 5

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Deprecated methods

def initialize_kernels(nx, xi, grid_type=1):
    ''' Deprecated routine. '''
    #This **MUST** be called before any other methods:
    #- Before methods to get kernels.
    #- Before methods to get evolution matrices.
    #- Before methods to initialize evolution matrices.
    #INPUT:
    #- nx .......... integer, number of x points.
    #                Must be even, and at least 6.
    #                Spacing depends on grid_type
    #- xi .......... numpy.array, with the xi points to be used; if a scalar
    #                is passed, it will be turned into a one-component array.
    #- grid_type ... integer, specifying the grid type.
    #                1 : nx linearly spaced points from -1+1/nx to 1-1/nx
    #                2 : nx/4 log-spaced points in each of the DGLAP regions,
    #                    and nx/2 linearly-spaced points in the ERBL region
    #OUTPUT:
    #- None
    #NOTES:
    #- If the user wishes to change either nx or xi, just call this method again.
    #  Calling this method deallocates any existing evolution matrices.
    print("The routine initialize_kernels is deprecated.")
    print("I will initialize your x and xi grids.")
    print("The kernels will be initialized during the first retrieval.")
    set_x_xi_grids(nx, xi, grid_type, lagrange_order=5)
    return

def initialize_evolution_matrices(Q2=matrix_dict['Q2'], nlo=False):
    # TODO: rework docstring , or try to deprecate
    ''' Initializes the evolution matrices.
    This **MUST** be called before any method to get the evolution matrices.
    Moreover, initialize_kernels (above) **MUST** be called first.
    nx and nxi should be consistent with those used to initialize the kernels.
    INPUT:
    - Q2 .... numpy.array, with *at least* two values
         first value is the Q2 value at which evolutuon starts
    - nlo ... boolean, False by default; whether to include NLO corrections.
    OUTPUT:
    - None
    NOTES:
    - If the user wishes to change the Q2 array, this must be called again.
    - If the user wishes to change the nx or xi, initialize_kernels
      must be called to do this.
    - If the user wishes to toggle NLO, this must be called again.
    - If initialize_kernels has been called, any initialized evolution matrices
      have been erased and this routine must be called again.
    '''
    nQ2 = Q2.shape[0]
    assert(nQ2 >= 2)
    # If the kernels are not initialized, we must initialize them...
    if(matrix_dict['kernel_init']==False):
        f90src.make_kernels_wrap()
        matrix_dict['kernel_init'] = True
    # Reiniitialize the Q2 grid if something new has been plugged in
    if(not np.array_equal(Q2, matrix_dict['Q2'])):
        set_Q2_grid(Q2)
    # Finally, make the evolution matrices
    f90src.make_matrices_wrap(Q2, nlo)
    return

def pixelspace(nx, xi=0.5, grid_type=1):
    # TODO: replace by returning cached x space
    ''' The x space used by tiktaalik.
    1. grid_type=1, linear grid (default)
       The x values are the central values of nx intervals evenly dividing
       the domain [-1,1]. For instance, if nx=4, then [-1,1] is divided into the
       intervals [-1,-0.5], [-0.5,0], [0,0.5] and [0.5,1]. The midpoints of these
       are -0.75, -0.25, 0.25 and 0.75. These four midpoints are used as x values.
       Independent of xi.
    2. grid_type=2, log-linear-log
       x is broken down into the two DGLAP regions (x < -xi and x > xi) and the
       ERBL region (-xi < x < xi). nx/4 points are placed in each of the DGLAP
       regions, and are geometrically spaced, more closely around the x=-xi or
       x=xi endpoint. The other nx/2 points are linearly spaced in the ERBL region.
       The placement of the points in each region is done using midpoints
       (geometric rather than arithmetic in the DGLAP region) of interals,
       similarly to the linear grid.
    '''
    x = f90src.pixelspace_wrap(nx, xi, grid_type)
    return x

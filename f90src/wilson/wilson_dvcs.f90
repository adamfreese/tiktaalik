! wilson_dvcs.f90
!
! by Adam Freese
! part of the package tiktaalik for GPD evolution
!
! created June 30, 2025.
!
! Wilson coefficients taken from
!   Braun, Ji and Schoenleber
!   Physical Review Letters 129 (2022) 172001
!   Braun:2022bpn
!   arxiv:2207.06818

module wilson_dvcs
  use constants,      only: pi, CF, TF
  use gridspace,      only: push_forward
  use integration,    only: integrate2
  use kernels_common, only: theta_step, abslog, log2
  use pixelation,     only: interpixel
  use specfun,        only: dilog

  implicit none
  public

  integer,  parameter, private :: dp = kind(1d0)

  contains

    ! TODO

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Leading order

    function re_Cq0_dvcs_cst(xi) result(w)
        real(dp), intent(in) :: xi
        real(dp) :: w
        !
        w = log((1.+xi)/(1.-xi))
    end function re_Cq0_dvcs_cst

    function re_Cq0_dvcs_sub(x,xi) result(w)
        real(dp), intent(in) :: x, xi
        real(dp) :: w
        !
        w = 1./(xi-x)
    end function re_Cq0_dvcs_sub

    function im_Cq0_dvcs_cst(xi) result(w)
        real(dp), intent(in) :: xi
        real(dp) :: w
        !
        w = pi
    end function im_Cq0_dvcs_cst

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Next-to-leading order, quark part

    function re_Cq1_dvcs_reg(x,xi) result(w)
        real(dp), intent(in) :: x, xi
        real(dp) :: w
        !
        w = 3.*CF * abslog(0.5*(x+xi)/xi) / (xi-x)
    end function re_Cq1_dvcs_reg

    function re_Cq1_dvcs_sub(x,xi) result(w)
        real(dp), intent(in) :: x, xi
        real(dp) :: w
        !
        w = CF * log2(0.5*(xi-x)/xi)/(xi-x)
        w = w - 9.*CF*re_Cq0_dvcs_sub(x,xi)
    end function re_Cq1_dvcs_sub

    function re_Cq1_dvcs_cst(xi) result(w)
        real(dp), intent(in) :: xi
        real(dp) :: w
        !
        w = CF*(log(0.5*(1.+xi)/xi)**3 - log(0.5*(1.-xi)/xi)**3)/3. 
        w = w + CF*pi**2*log(0.5*(1.-xi)/xi)
        w = w - 9.*CF*re_Cq0_dvcs_cst(xi)
    end function re_Cq1_dvcs_cst

    function im_Cq1_dvcs_reg(x,xi) result(w)
        real(dp), intent(in) :: x, xi
        real(dp) :: w
        !
        w = 3.*CF * pi/(x+xi) * theta_step(x-xi)
    end function im_Cq1_dvcs_reg

    function im_Cq1_dvcs_sub(x,xi) result(w)
        real(dp), intent(in) :: x, xi
        real(dp) :: w
        !
        w = -2.*pi*abslog(0.5*(xi-x)/xi)/(xi-x) * theta_step(x-xi)
    end function im_Cq1_dvcs_sub

    function im_Cq1_dvcs_cst(xi) result(w)
        real(dp), intent(in) :: xi
        real(dp) :: w
        !
        w = pi*(log2(0.5*(1.-xi)/xi) - pi**2/3.)
        w = w - 9.*CF*im_Cq0_dvcs_cst(xi)
    end function im_Cq1_dvcs_cst

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Next-to-leading order, gluon part

    function re_CG1_dvcs_reg(x,xi) result(w)
        ! Does not yet contain factor sum(eq**2)
        real(dp), intent(in) :: x, xi
        real(dp) :: w
        real(dp) :: z, zbar
        !
        zbar = 0.5*(x+xi)/xi
        z = 0.5*(xi-x)/xi
        !w = TF*( 4.*z*abslog(zbar) + 8.*abslog(zbar) - 2.*log2(zbar) ) / (xi-x)**2
        ! New approach to A3
        w = TF*( 4.*z*abslog(zbar) - 2.*log2(zbar) ) / (xi-x)**2
    end function re_CG1_dvcs_reg

    function re_CG1_dvcs_sub(x,xi) result(w)
        ! Does not yet contain factor sum(eq**2)
        real(dp), intent(in) :: x, xi
        real(dp) :: w
        !
        w = TF * 2./xi*abslog(0.5*(xi-x)/xi) / (xi-x)
        ! New approach to A3
        w = w + TF * 8.*abslog(0.5*(xi+x)/xi) / (xi-x)**2
    end function re_CG1_dvcs_sub

    function re_CG1_dvcs_cst(xi) result(w)
        ! Does not yet contain factor sum(eq**2)
        real(dp), intent(in) :: xi
        real(dp) :: w
        !
        w = TF/xi * ( log2(0.5*(1.+xi)/xi) - log2(0.5*(1.-xi)/xi) + pi**2 )
        ! New approach to A3
        w = w + TF/xi * 4. * ( log(0.5*(1.-xi)/xi)/(1.+xi) - log(0.5*(1.+xi)/xi)/(1.-xi) )
    end function re_CG1_dvcs_cst

    function im_CG1_dvcs_reg(x,xi) result(w)
        ! Does not yet contain factor sum(eq**2)
        real(dp), intent(in) :: x, xi
        real(dp) :: w
        real(dp) :: z, zbar
        !
        w = 0.
        if(x > xi) then
          z = 0.5*(xi-x)/xi
          zbar = 0.5*(xi+x)/xi
          w = - pi * TF*( 4.*zbar + 8. - 4.*abslog(z) ) / (xi+x)**2
        endif
    end function im_CG1_dvcs_reg

    function im_CG1_dvcs_sub(x,xi) result(w)
        ! Does not yet contain factor sum(eq**2)
        real(dp), intent(in) :: x, xi
        real(dp) :: w
        !
        w = - TF/xi * 2.*pi * theta_step(x-xi) / (xi-x)
    end function im_CG1_dvcs_sub

    function im_CG1_dvcs_cst(xi) result(w)
        ! Does not yet contain factor sum(eq**2)
        real(dp), intent(in) :: xi
        real(dp) :: w
        !
        w = 2.*pi * TF/xi * log(0.5*(1.-xi)/xi)
        ! New approach to A3 releaved a missing term...?
        w = w - TF/xi * 4.*pi
    end function im_CG1_dvcs_cst

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Matrix-building routines
    ! TODO: move to separate matrices file (maybe combine with ulities?)

    !subroutine make_wilson_matrix(nx, nxi, matrix, kernel


    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Pixelated shifts
    ! TODO: move to a separate utilities file

    function pixel_coef(Kreg, Ksub, Kcst, xi, N, i, grid_type) result(shift)
        ! A modified pixel convolution method for Wilson coefficients.
        ! There's only one index since x is integrated out.
        real(dp), external   :: Kreg, Kcst, Ksub
        real(dp), intent(in) :: xi
        integer,  intent(in) :: N, i, grid_type
        real(dp) :: shift
        !
        shift = 0.0_dp
        shift = shift + pixel_coef_reg(Kreg, xi, N, i, grid_type)
        shift = shift + pixel_coef_sub(Ksub, xi, N, i, grid_type)
        shift = shift + pixel_coef_cst(Kcst, xi, N, i, grid_type)
    end function pixel_coef

    function pixel_coef_reg(kernel, xi, n_pixels, i, grid_type) result(shift)
        real(dp), external   :: kernel
        real(dp), intent(in) :: xi
        integer,  intent(in) :: n_pixels, i, grid_type
        real(dp) :: shift
        !
        shift = integrate2(integrand, xi, n_pixels, grid_type)
        return
        contains
          function integrand(x) result(intd)
              real(dp), intent(in) :: x
              real(dp) :: intd
              !
              real(dp) :: fx
              fx = interpixel(n_pixels, i-1, x, xi, grid_type)
              intd = kernel(x,xi)*fx
          end function integrand
    end function pixel_coef_reg

    function pixel_coef_sub(kernel, xi, n_pixels, i, grid_type) result(shift)
        real(dp), external   :: kernel
        real(dp), intent(in) :: xi
        integer,  intent(in) :: n_pixels, i, grid_type
        real(dp) :: shift
        !
        shift = integrate2(integrand, xi, n_pixels, grid_type)
        return
        contains
          function integrand(x) result(intd)
              real(dp), intent(in) :: x
              real(dp) :: intd
              !
              real(dp) :: fx, fxi
              fx  = interpixel(n_pixels, i-1,  x, xi, grid_type)
              fxi = interpixel(n_pixels, i-1, xi, xi, grid_type)
              intd = kernel(x,xi)*(fx-fxi)
          end function integrand
    end function pixel_coef_sub

    function pixel_coef_cst(kernel, xi, n_pixels, i, grid_type) result(shift)
        real(dp), external   :: kernel
        real(dp), intent(in) :: xi
        integer,  intent(in) :: n_pixels, i, grid_type
        real(dp) :: shift
        !
        real(dp) :: fxi
        fxi = interpixel(n_pixels, i-1, xi, xi, grid_type)
        shift = kernel(xi) * fxi
    end function pixel_coef_cst

end module wilson_dvcs

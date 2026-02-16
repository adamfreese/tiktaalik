! continuum_wilson.f90
!
! by Adam Freese
! part of the package tiktaalik for GPD evolution
!
! Created on July 25, 2025.

module continuum_wilson
  use alpha_qcd,     only: get_alpha_QCD, get_neff
  use constants,     only: i_, pi
  use integration,   only: adaptive_integrate2
  use kernels_common
  use wilson_dvcs

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)

  public :: continuum_cff_q, continuum_cff_g

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Public methods for testing

    function continuum_cff_q(gpd, xi, Q2, nlo) result(cff)
        real(dp), external   :: gpd
        real(dp), intent(in) :: xi, Q2
        logical,  intent(in) :: nlo
        complex(dp) :: cff
        !
        real(dp) :: ReCFF0, ImCFF0, ReCFF1, ImCFF1
        ReCFF0 = convolve(gpd, zero_func,       re_Cq0_dvcs_sub, re_Cq0_dvcs_cst, xi)
        ImCFF0 = convolve(gpd, zero_func,       zero_func,       im_Cq0_dvcs_cst, xi)
        ReCFF1 = convolve(gpd, re_Cq1_dvcs_reg, re_Cq1_dvcs_sub, re_Cq1_dvcs_cst, xi)
        ImCFF1 = convolve(gpd, im_Cq1_dvcs_reg, im_Cq1_dvcs_sub, im_Cq1_dvcs_cst, xi)
        cff = ReCFF0 + i_*ImCFF0
        if(nlo) then
          cff = cff + get_alpha_QCD(Q2)/(4.*pi)*( ReCFF1 + i_*ImCFF1 )
        endif
    end function continuum_cff_q

    function continuum_cff_g(gpd, xi, Q2, nlo) result(cff)
        real(dp), external   :: gpd
        real(dp), intent(in) :: xi, Q2
        logical,  intent(in) :: nlo
        complex(dp) :: cff
        !
        real(dp) :: ReCFF1, ImCFF1
        ReCFF1 = convolve(gpd, re_Cg1_dvcs_reg, re_Cg1_dvcs_sub, re_Cg1_dvcs_cst, xi)
        ImCFF1 = convolve(gpd, im_Cg1_dvcs_reg, im_Cg1_dvcs_sub, im_Cg1_dvcs_cst, xi)
        cff = 0.0_dp
        if(nlo) then
          cff = cff + get_alpha_QCD(Q2)/(4.*pi)*( ReCFF1 + i_*ImCFF1 )
        endif
    end function continuum_cff_g

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Continuous convolution

    function convolve(func, Kreg, Ksub, Kcst, xi) result(shift)
        real(dp), external   :: func, Kreg, Ksub, Kcst
        real(dp), intent(in) :: xi
        real(dp) :: shift
        !
        shift = 0.0_dp
        shift = shift + conv_reg(func, Kreg, xi)
        shift = shift + conv_sub(func, Ksub, xi)
        shift = shift + conv_cst(func, Kcst, xi)
    end function convolve

    function conv_reg(func, kernel, xi) result(shift)
        real(dp), external   :: func, kernel
        real(dp), intent(in) :: xi
        real(dp) :: shift
        !
        shift = adaptive_integrate2(integrand, xi)
        return
        contains
          function integrand(x) result(intd)
              real(dp), intent(in) :: x
              real(dp) :: intd
              !
              intd = kernel(x,xi)*func(x)
          end function integrand
    end function conv_reg

    function conv_sub(func, kernel, xi) result(shift)
        real(dp), external   :: func, kernel
        real(dp), intent(in) :: xi
        real(dp) :: shift
        !
        shift = adaptive_integrate2(integrand, xi)
        return
        contains
          function integrand(x) result(intd)
              real(dp), intent(in) :: x
              real(dp) :: intd
              !
              intd = kernel(x,xi)*(func(x) - func(xi))
          end function integrand
    end function conv_sub

    function conv_cst(func, kernel, xi) result(shift)
        real(dp), external   :: func, kernel
        real(dp), intent(in) :: xi
        real(dp) :: shift
        !
        shift = kernel(xi) * func(xi)
    end function conv_cst

end module continuum_wilson

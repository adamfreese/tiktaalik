! continuum_evolution.f90
!
! by Adam Freese
! part of the package tiktaalik for GPD evolution
!
! Created on May 24, 2024.

module continuum_evolution
  use alpha_qcd,     only: get_alpha_QCD, get_neff
  use integration,   only: adaptive_integrate
  use kernels_common
  use kvlo4tests

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)

  public :: continuum_shift_QQ, continuum_shift_QG, continuum_shift_GQ, continuum_shift_GG

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Public methods for testing

    function continuum_shift_QQ(func, x, xi, Q2, nlo, i_ns_type) result(shift)
        real(dp), external   :: func
        real(dp), intent(in) :: x, xi, Q2
        logical,  intent(in) :: nlo
        integer,  intent(in) :: i_ns_type
        real(dp) :: shift
        !
        real(dp) :: S1, S2, al2pi
        al2pi = get_alpha_QCD(Q2) / (2.*pi)
        S1 = convolve(func, zero_func, KV0_QQ_pls, KV0_QQ_cst, x, xi)
        S2 = 0.0_dp
        if(nlo) then
          select case(i_ns_type)
          case(1)
            S2 = convolve(func, zero_func, tKV1_NSp_pls, tKV1_NSp_cst, x, xi)
          case(-1)
            S2 = convolve(func, zero_func, tKV1_NSm_pls, tKV1_NSm_cst, x, xi)
          case(0)
            S2 = convolve(func, tKV1_QQ_reg, tKV1_NSp_pls, tKV1_NSp_cst, x, xi)
          end select
        endif
        shift = al2pi*S1 + al2pi**2*S2
    end function continuum_shift_QQ

    function continuum_shift_QG(func, x, xi, Q2, nlo) result(shift)
        real(dp), external   :: func
        real(dp), intent(in) :: x, xi, Q2
        logical,  intent(in) :: nlo
        real(dp) :: shift
        !
        real(dp) :: S1, S2, al2pi
        al2pi = get_alpha_QCD(Q2) / (2.*pi)
        S1 = convolve(func, KV0_QG_reg, zero_func, zero_func, x, xi)
        S2 = 0.0_dp
        if(nlo) then
          S2 = convolve(func, tKV1_QG_reg, zero_func, zero_func, x, xi)
        endif
        shift = al2pi*S1 + al2pi**2*S2
    end function continuum_shift_QG

    function continuum_shift_GQ(func, x, xi, Q2, nlo) result(shift)
        real(dp), external :: func
        real(dp), intent(in) :: x, xi, Q2
        logical,  intent(in) :: nlo
        real(dp) :: shift
        !
        real(dp) :: S1, S2, al2pi
        al2pi = get_alpha_QCD(Q2) / (2.*pi)
        S1 = convolve(func, KV0_GQ_reg, zero_func, zero_func, x, xi)
        S2 = 0.0_dp
        if(nlo) then
          S2 = convolve(func, tKV1_GQ_reg, zero_func, zero_func, x, xi)
        endif
        shift = al2pi*S1 + al2pi**2*S2
    end function continuum_shift_GQ

    function continuum_shift_GG(func, x, xi, Q2, nlo) result(shift)
        real(dp), external :: func
        real(dp), intent(in) :: x, xi, Q2
        logical,  intent(in) :: nlo
        real(dp) :: shift
        !
        real(dp) :: S1, S2, al2pi
        al2pi = get_alpha_QCD(Q2) / (2.*pi)
        S1 = convolve(func, KV0_GG_reg, KV0_GG_pls, KV0_GG_cst, x, xi)
        S2 = 0.0_dp
        if(nlo) then
            S2 = convolve(func, zero_func, tKV1_GG_pls, tKV1_GG_cst, x, xi)
        endif
        shift = al2pi*S1 + al2pi**2*S2
    end function continuum_shift_GG

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Continuous shifts

    function convolve(func, Kreg, Kpls, Kcst, x, xi) result(shift)
        real(dp), external   :: func, Kreg, Kpls, Kcst
        real(dp), intent(in) :: x, xi
        real(dp) :: shift
        !
        shift = 0.0_dp
        shift = shift + conv_reg(func, Kreg, x, xi)
        shift = shift + conv_pls(func, Kpls, x, xi)
        shift = shift + conv_cst(func, Kcst, x, xi)
    end function convolve

    function conv_reg(func, kernel, x, xi) result(shift)
        real(dp), external   :: func, kernel
        real(dp), intent(in) :: x, xi
        real(dp) :: shift
        !
        shift = adaptive_integrate(integrand, x, xi)
        return
        contains
          function integrand(y) result(intd)
              real(dp), intent(in) :: y
              real(dp) :: intd
              !
              intd = kernel(x,y,xi)*func(y)
          end function integrand
    end function conv_reg

    function conv_pls(func, kernel, x, xi) result(shift)
        real(dp), external   :: func, kernel
        real(dp), intent(in) :: x, xi
        real(dp) :: shift
        !
        shift = adaptive_integrate(integrand, x, xi)
        return
        contains
          function integrand(y) result(intd)
              real(dp), intent(in) :: y
              real(dp) :: intd
              !
              intd = kernel(x,y,xi)*(func(y) - func(x))
          end function integrand
    end function conv_pls

    function conv_cst(func, kernel, x, xi) result(shift)
        real(dp), external   :: func, kernel
        real(dp), intent(in) :: x, xi
        real(dp) :: shift
        !
        shift = kernel(x,xi) * func(x)
    end function conv_cst

end module continuum_evolution

! gk_zero.f90
!
! part of the package tiktaalik for GPD evolution
!
! Code added July 25, 2025 by Adam Freese.
!
! This implements a model due to Goloshkokov and Kroll.
! The specific reference I consulted is was:
!   P. Kroll, H. Moutarde, F. Sabatie
!   European Physical Journal C (2013) 73:2278
!   arxiv:1210.6975
!   Kroll:2012sm
! However, the original references for the GK model are:
!   S.V. Goloskokov and P. Kroll,
!   European Physical Journal C (2005) 42:281-301
!   arxiv:hep-ph/0501242
!   Goloskokov:2005sd
! and
!   S.V. Goloskokov and P. Kroll,
!   European Physical Journal C (2008) 53:367-384
!   arxiv:0708.3569
!   Goloskokov:2007nt
!
! This just gives the GPDs when xi=0.

module gk_zero

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)
  real(dp), parameter, private :: pi = acos(-1.0_dp)

  public :: H_n1_zero, H_n2_zero

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Forms for zero skewness (exact)
    ! See Eqs. (7) and (16)

    function H_n1_zero(x, t, C, delta, ap) result(H)
        real(dp), intent(in) :: x, t, delta, ap, C(0:3)
        real(dp) :: H
        !
        integer :: j
        H = 0.0_dp
        if(x <= 0.0_dp) return
        do j=0, 3, 1
          H = H + C(j)*sqrt(x)**j
        end do
        H = H * (1.-x)**3 * x**(-delta-ap*t)
    end function H_n1_zero

    function H_n2_zero(x, t, C, delta, ap, sgn) result(H)
        real(dp), intent(in) :: x, t, sgn, delta, ap, C(0:3)
        real(dp) :: H
        !
        integer :: j
        H = 0.0_dp
        do j=0, 3, 1
          H = H + C(j)*sqrt(abs(x))**j
        end do
        H = H * (1.-abs(x))**5 * abs(x)**(-delta-ap*t)
        if(x <= 0.0_dp) H = H * sgn
    end function H_n2_zero

end module gk_zero

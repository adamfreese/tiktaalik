! gk_smolxi.f90
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
! This module uses an approximate form for small xi,
! to mitigate numerical instabilities in the analytic results
! of gk_analytic when xi << |x|.

module gk_smol
  use gk_zero,     only: H_n1_zero, H_n2_zero

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)
  real(dp), parameter, private :: pi = acos(-1.0_dp)

  public :: H_n1_smol, H_n2_smol

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Forms for small skewness

    function H_n1_smol(x, xi, t, C, delta, ap) result(H)
        real(dp), intent(in) :: x, xi, t, delta, ap, C(0:3)
        real(dp) :: H
        !
        integer :: j
        H = H_n1_zero(x, t, C, delta, ap) + xi*H_n1_term(x, t, C, delta, ap)
    end function H_n1_smol

    function H_n2_smol(x, xi, t, C, delta, ap, sgn) result(H)
        real(dp), intent(in) :: x, xi, t, sgn, delta, ap, C(0:3)
        real(dp) :: H
        !
        integer :: j
        H = xi*H_n2_term(x, t, C, delta, ap)
        if(x < 0.0_dp) H = H * sgn
        H = H + H_n2_zero(x, t, C, delta, ap, sgn)
    end function H_n2_smol

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Implementation details

    function H_n1_term(x, t, C, delta, ap) result(H)
        ! The order-xi correction
        real(dp), intent(in) :: x, t, delta, ap, C(0:3)
        real(dp) :: H
        !
        integer :: j, w
        H = 0.0_dp
        do j=0, 3, 1
          do w=0, 3, 1
            H = H - pw(j, w, t, delta, ap)*C(j)*fw(w,1)*abs(x)**(pw(j, w, t, delta, ap)-1.)*(1.-abs(x))
          end do
        end do
        H = H * sign(1.0_dp, x) * 0.25
        H = H * 0.75 ! 3/4
    end function H_n1_term

    function H_n2_term(x, t, C, delta, ap) result(H)
        ! The order-xi correction
        real(dp), intent(in) :: x, t, delta, ap, C(0:3)
        real(dp) :: H
        !
        integer :: j, w
        H = 0.0_dp
        do j=0, 3, 1
          do w=0, 5, 1
            H = H - pw(j, w, t, delta, ap)*C(j)*fw(w,2)*abs(x)**(pw(j, w, t, delta, ap)-1.)*(1.-abs(x))
          end do
        end do
        H = H * sign(1.0_dp, x) / 6.
        H = H * 0.9375 ! 15/16
    end function H_n2_term

    function fw(w,n) result(f)
        integer, intent(in) ::  w, n
        real(dp) :: f
        f = 0.0_dp
        select case(n)
        case(1)
          select case(w)
          case(0)
            f = 1.
          case(1)
            f = -3.
          case(2)
            f = 3.
          case(3)
            f = -1.
          end select
        case(2)
          select case(w)
          case(0)
            f = 1.
          case(1)
            f = -5.
          case(2)
            f = 10.
          case(3)
            f = -10.
          case(4)
            f = 5.
          case(5)
            f = -1.
          end select
        end select
    end function fw

    function pw(j, w, t, delta, ap) result(p)
        ! Copied from gk_analytic
        integer,  intent(in) :: j, w
        real(dp), intent(in) :: t, delta, ap
        real(dp) :: p
        !
        p = real(j)/2. + real(w) - delta - ap*t
    end function pw

end module gk_smol

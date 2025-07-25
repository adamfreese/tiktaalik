! gk_analytic.f90
!
! part of the package tiktaalik for GPD evolution
!
! Code added July 30, 2024 by Adam Freese.
! Split off on July 25, 2025.
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
! This module uses analytic results for evaluating the (alpha,beta) integrals.
! The results are numerically unstable when xi << |x|, however.

module gk_analytic

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)
  real(dp), parameter, private :: pi = acos(-1.0_dp)

  public :: H_n1_analytic, H_n2_analytic

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Forms for moderate-to-large skewness

    function H_n1_analytic(x, xi, t, C, delta, ap) result(H)
        real(dp), intent(in) :: x, xi, t, delta, ap, C(0:3)
        real(dp) :: H
        !
        real(dp) :: zi, zf
        H = 0.0_dp
        zi = (x-xi)/(1.-xi)
        zf = (x+xi)/(1.+xi)
        if(x > xi) then
          H = H_n1_term(x, xi, t, zf, C, delta, ap) - H_n1_term(x, xi, t, zi, C, delta, ap)
        elseif( x > -xi) then
          H = H_n1_term(x, xi, t, zf, C, delta, ap)
        endif
    end function H_n1_analytic

    function H_n2_analytic(x, xi, t, C, delta, ap, sgn) result(H)
        real(dp), intent(in) :: x, xi, t, sgn, delta, ap, C(0:3)
        real(dp) :: H
        !
        real(dp) :: zi, zf, zm, zj
        H = 0.0_dp
        zi = (x-xi)/(1.-xi)
        zf = (x+xi)/(1.+xi)
        zj = (-x-xi)/(1.-xi)
        zm = (-x+xi)/(1.+xi)
        if(x > xi) then
          H = H_n2_term(x, xi, t, zf, C, delta, ap) - H_n2_term(x, xi, t, zi, C, delta, ap)
        elseif( x > -xi) then
          H = H_n2_term(x, xi, t, zf, C, delta, ap) + sgn*H_n2_term(-x, xi, t, zm, C, delta, ap)
        else
          H = sgn*(H_n2_term(-x, xi, t, zm, C, delta, ap) - H_n2_term(-x, xi, t, zj, C, delta, ap))
        endif
    end function H_n2_analytic

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Implementation details

    function H_n1_term(x, xi, t, z, C, delta, ap) result(H)
        real(dp), intent(in) :: x, xi, t, z, delta, ap, C(0:3)
        real(dp) :: H
        !
        integer :: j, w
        H = 0.0_dp
        do j=0, 3, 1
          do w=1, 3, 1
            H = H + C(j)*kw(x, xi, w)*z**pw(j, w, t, delta, ap) / pw(j, w, t, delta, ap)
          end do
        end do
        H = H * 0.75/xi ! 3/4
    end function H_n1_term

    function H_n2_term(x, xi, t, z, C, delta, ap) result(H)
        real(dp), intent(in) :: x, xi, t, z, delta, ap, C(0:3)
        real(dp) :: H
        !
        integer :: j, w
        H = 0.0_dp
        do j=0, 3, 1
          do w=1, 5, 1
            H = H + C(j)*ew(x, xi, w)*z**pw(j, w, t, delta, ap) / pw(j, w, t, delta, ap)
          end do
        end do
        H = H * 0.9375/xi ! 15/16
    end function H_n2_term

    function kw(x, xi, w) result(kappa)
        real(dp), intent(in) :: x, xi
        integer,  intent(in) :: w
        real(dp) :: kappa
        !
        select case(w)
        case(1)
          kappa = 1. - x**2/xi**2
        case(2)
          kappa = 2.*(x/xi**2-1.)
        case(3)
          kappa = 1. - 1./xi**2
        case default
          kappa = 0.0_dp
        end select
    end function kw

    function ew(x, xi, w) result(eta)
        real(dp), intent(in) :: x, xi
        integer,  intent(in) :: w
        real(dp) :: eta
        !
        select case(w)
        case(1)
          eta = kw(x,xi,1)**2
        case(2)
          eta = 2.*kw(x,xi,1)*kw(x,xi,2)
        case(3)
          eta = kw(x,xi,2)**2 + 2.*kw(x,xi,1)*kw(x,xi,3)
        case(4)
          eta = 2.*kw(x,xi,2)*kw(x,xi,3)
        case(5)
          eta = kw(x,xi,3)**2
        case default
          eta = 0.0_dp
        end select
    end function ew

    function pw(j, w, t, delta, ap) result(p)
        integer,  intent(in) :: j, w
        real(dp), intent(in) :: t, delta, ap
        real(dp) :: p
        !
        p = real(j)/2. + real(w) - delta - ap*t
    end function pw

end module gk_analytic

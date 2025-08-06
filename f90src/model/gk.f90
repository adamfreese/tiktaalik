! gk.f90
!
! part of the package tiktaalik for GPD evolution
!
! Code added July 30, 2024 by Adam Freese.
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

module gk
  use gk_analytic, only: H_n1_analytic, H_n2_analytic
  use gk_smol,     only: H_n1_smol,     H_n2_smol
  use gk_zero,     only: H_n1_zero,     H_n2_zero

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)
  real(dp), parameter, private :: pi = acos(-1.0_dp)
  real(dp), parameter, private :: eps = 1e-12_dp

  ! Valence
  real(dp), parameter, private :: a0_val = 0.48_dp
  real(dp), parameter, private :: ap_val = 0.9_dp ! GeV**-2
  real(dp), parameter, private :: c0_uv = 1.52_dp
  real(dp), parameter, private :: c1_uv = 2.88_dp
  real(dp), parameter, private :: c2_uv = -0.095_dp
  real(dp), parameter, private :: c3_uv = 0.0_dp
  real(dp), parameter, private :: c0_dv = 0.76_dp
  real(dp), parameter, private :: c1_dv = 3.11_dp
  real(dp), parameter, private :: c2_dv = -3.99_dp
  real(dp), parameter, private :: c3_dv = 0.0_dp

  ! Sea &  gluons
  real(dp), parameter, private :: a0_sea = 1.1_dp
  real(dp), parameter, private :: ap_sea = 0.15_dp ! GeV**-2
  real(dp), parameter, private :: c0_g = 2.23_dp
  real(dp), parameter, private :: c1_g = 5.43_dp
  real(dp), parameter, private :: c2_g = -34.0_dp
  real(dp), parameter, private :: c3_g = 40.6_dp
  real(dp), parameter, private :: c0_s = 0.123_dp
  real(dp), parameter, private :: c1_s = -0.327_dp
  real(dp), parameter, private :: c2_s = 0.692_dp
  real(dp), parameter, private :: c3_s = -0.486_dp
  real(dp), parameter, private :: kappa_s = 1.68_dp

  ! C matrix
  real(dp), parameter, private :: C(0:3,0:3) = reshape( &
      [ c0_g, c0_uv, c0_dv, c0_s, &
      &  c1_g, c1_uv, c1_dv, c1_s, &
      &  c2_g, c2_uv, c2_dv, c2_s, &
      &  c3_g, c3_uv, c3_dv, c3_s ], &
      [4,4] )

  ! Arrays for regge trajectories
  real(dp), parameter, private :: delta(0:3) = [ &
      a0_sea - 1., a0_val, a0_val, a0_sea ]
  real(dp), parameter, private :: ap(0:3) = [ &
      ap_sea, ap_val, ap_val, ap_sea ]

  public :: Hu, Hd, Hs, Hg, Hu_val, Hd_val, Hu_sea, Hd_sea

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Direct interfaces

    function Hu(x,xi,t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = Hu_val(x,xi,t) + Hu_sea(x,xi,t)
    end function Hu

    function Hd(x,xi,t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = Hd_val(x,xi,t) + Hd_sea(x,xi,t)
    end function Hd

    function Hs(x, xi, t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = H_n2(x, xi, t, C(3,:), delta(3), ap(3), -1.0_dp)
    end function Hs

    function Hg(x, xi, t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = H_n2(x, xi, t, C(0,:), delta(0), ap(0), 1.0_dp)
    end function Hg

    function Hu_val(x, xi, t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = H_n1(x, xi, t, C(1,:), delta(1), ap(1))
    end function Hu_val

    function Hd_val(x, xi, t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = H_n1(x, xi, t, C(2,:), delta(2), ap(2))
    end function Hd_val

    function Hu_sea(x, xi, t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = kappa_s * H_n2(x, xi, t, C(3,:), delta(3), ap(3), -1.0_dp)
    end function Hu_sea

    function Hd_sea(x, xi, t) result(H)
        real(dp), intent(in) :: x, xi, t
        real(dp) :: H
        !
        H = kappa_s * H_n2(x, xi, t, C(3,:), delta(3), ap(3), -1.0_dp)
    end function Hd_sea

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Switchers

    function H_n1(x, xi, t, C, delta, ap) result(H)
        real(dp), intent(in) :: x, xi, t, delta, ap, C(0:3)
        real(dp) :: H
        !
        if(xi < 0.01*abs(x)) then
          H = H_n1_zero(x, t, C, delta, ap)
        elseif( abs(x - xi) < eps) then
          H = H_n1_analytic(xi+eps, xi, t, C, delta, ap)
        elseif( abs(x + xi) < eps) then
          H = H_n1_analytic(-xi-eps, xi, t, C, delta, ap)
        else
          H = H_n1_analytic(x, xi, t, C, delta, ap)
        endif
    end function H_n1

    function H_n2(x, xi, t, C, delta, ap, sgn) result(H)
        real(dp), intent(in) :: x, xi, t, delta, ap, C(0:3), sgn
        real(dp) :: H
        !
        if(xi < 0.01*abs(x)) then
          H = H_n2_zero(x, t, C, delta, ap, sgn)
        elseif( abs(x - xi) < eps) then
          H = H_n2_analytic(xi+eps, xi, t, C, delta, ap, sgn)
        elseif( abs(x + xi) < eps) then
          H = H_n2_analytic(-xi-eps, xi, t, C, delta, ap, sgn)
        else
          H = H_n2_analytic(x, xi, t, C, delta, ap, sgn)
        endif
    end function H_n2


end module gk

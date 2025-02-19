! gridspace.f90
!
! by Adam Freese
! part of the package tiktaalik for GPD evolution
!
! Created February 14, 2025.

module gridspace
  implicit none
  private

  integer,  parameter :: dp = kind(1d0)
  real(dp), parameter :: lambda = 6.0_dp

  public :: linspace, geomspace, tanhspace, &
      & push_forward, pull_back, push_jacob

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Conversion between linear eta space and x space

    function push_forward(eta, xi, nx, grid_type) result(x)
        ! Grid types:
        ! - 1 ... linear
        ! - 2 ... log / linear / log
        real(dp), intent(in) :: eta, xi
        integer,  intent(in) :: nx, grid_type
        real(dp) :: x
        !
        select case(grid_type)
        case(1)
          ! Case 1 : linear grid type, same as in original paper
          x = eta
        case(2)
          ! Case 2 : log in DGLAP regions, linear in ERBL region
          if(eta < -0.5) then
            x = -(xi**2)**(1.+eta)
          elseif(eta <= 0.5) then
            x = 2.*xi*eta
          else
            x = (xi**2)**(1.-eta)
          endif
        case default
          ! Default to linear spacing given invalid grid_type
          x = eta
        end select
    end function push_forward

    function pull_back(x, xi, nx, grid_type) result(eta)
        ! Grid types:
        ! - 1 ... linear
        ! - 2 ... log / linear / log
        real(dp), intent(in) :: x, xi
        integer,  intent(in) :: nx, grid_type
        real(dp) :: eta
        !
        select case(grid_type)
        case(1)
          ! Case 1 : linear grid type, same as in original paper
          eta = x
        case(2)
          ! Case 2 : log in DGLAP regions, linear in ERBL region
          if(x < -xi) then
            eta = 0.5*log(-x/xi**2)/log(xi)
          elseif(x <= xi) then
            eta = 0.5*x/xi
          else
            eta = 0.5*log(xi**2/x)/log(xi)
          endif
        case default
          ! Default to linear spacing given invalid grid_type
          eta = x
        end select
    end function pull_back

    function push_jacob(eta, xi, nx, grid_type) result(J)
        ! Grid types:
        ! - 1 ... linear
        ! - 2 ... log / linear / log
        real(dp), intent(in) :: eta, xi
        integer,  intent(in) :: nx, grid_type
        real(dp) :: J
        !
        select case(grid_type)
        case(1)
          ! Case 1 : linear grid type, same as in original paper
          J = 1.0_dp
        case(2)
          ! Case 2 : log in DGLAP regions, linear in ERBL region
          if(eta < -0.5) then
            J = -(xi**2)**(1.+eta) * 2.*log(xi)
          elseif(eta <= 0.5) then
            J = 2.*xi
          else
            J = -(xi**2)**(1.-eta) * 2.*log(xi)
          endif
        case default
          ! Default to linear spacing given invalid grid_type
          J = 1.0_dp
        end select
    end function push_jacob

    !function pixelspace(xi, n_pixels, grid_type) result(xx)
    !    real(dp), intent(in) :: xi
    !    integer, intent(in)  :: n_pixels, grid_type
    !    real(dp) :: xx(n_pixels)
    !    !
    !    real(dp) :: endpoint
    !    integer :: i
    !    endpoint = 1. - 0.5/real(n_pixels)
    !    xx = linspace(-endpoint, endpoint, n_pixels)
    !    do i=1, n_pixels, 1
    !      xx(i) = push_forward(xx(i), xi, n_pixels, grid_type)
    !    end do
    !end function pixelspace

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Array-building routines

    function linspace(a, b, N) result(xx)
        ! Linearly-spaced array
        real(dp), intent(in) :: a, b
        integer,  intent(in) :: N
        real(dp) :: xx(N)
        !
        integer :: i
        do i=1, N, 1
          xx(i) = (b-a)*real(i-1)/real(N-1) + a
        end do
    end function linspace

    function geomspace(a, b, N) result(xx)
        ! Geometrically-spaced array, more dense on the small end
        real(dp), intent(in) :: a, b
        integer,  intent(in) :: N
        real(dp) :: xx(N)
        !
        integer :: i
        do i=1, N, 1
          xx(i) = (b/a)**(real(i-1)/real(N-1)) * a
        end do
    end function geomspace

    function tanhspace(a, b, N) result(xx)
        ! Geometrically-spaced, with increased density on both ends
        real(dp), intent(in) :: a, b
        integer,  intent(in) :: N
        real(dp) :: xx(N)
        !
        xx = linspace(atanh(2.*a-1.), atanh(2.*b-1.), N)
        xx = 0.5*(1. + tanh(xx))
    end function tanhspace

end module gridspace

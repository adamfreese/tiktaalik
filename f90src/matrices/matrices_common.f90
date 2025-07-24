! matrices_common.f90
!
! by Adam Freese
! part of the package tiktaalik for GPD evolution
!
! Created July 24, 2025
! Part of a breakup of the previous matevo.f90 file into smaller chunks

module matrices_common
  use gridspace, only: pixelspace
  use pixelation, only: initialize_lagrange_weights

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)
  real(dp), parameter, private :: epsilon = 1e-6_dp

  ! Cached information
  ! gt_ is grid type
  integer, private :: gt_, nx_, nxi_, nQ2_

  ! Cached arrays
  real(dp), allocatable, private :: xx_(:,:) ! 2D to allow for xi-dependent spacing
  real(dp), allocatable, private :: xi_(:)
  real(dp), allocatable, private :: Q2_(:)

  public :: get_nx, get_nxi, get_nQ2, get_grid_type, get_x, get_xi, get_Q2, &
      & initialize_x_xi, initialize_Q2

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Basic matrices

    function identity_matrix(N) result(M)
        integer, intent(in) :: N
        real(dp) :: M(N,N)
        !
        integer :: i
        M = 0.0_dp
        do i=1, N, 1
          M(i,i) = 1.0_dp
        end do
    end function identity_matrix

    function zero_matrix(N) result(M)
        integer, intent(in) :: N
        real(dp) :: M(N,N)
        !
        M = 0.0_dp
    end function zero_matrix

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Get methods

    function get_nx() result(nx)
        integer :: nx
        nx = nx_
    end function get_nx

    function get_nxi() result(nxi)
        integer :: nxi
        nxi = nxi_
    end function get_nxi

    function get_nQ2() result(nQ2)
        integer :: nQ2
        nQ2 = nQ2_
    end function get_nQ2

    function get_grid_type() result(grid_type)
        integer :: grid_type
        grid_type = gt_
    end function get_grid_type

    function get_x(nx, nxi) result(xx)
        integer, intent(in) :: nx, nxi
        real(dp) :: xx(nx,nxi)
        xx = xx_
    end function get_x

    function get_xi(nxi) result(xi)
        integer, intent(in) :: nxi
        real(dp) :: xi(nxi)
        if(nxi==nxi_) then
          xi = xi_
        else
          print *, "ERROR: incorrect xi grid size."
          print *, "Got:", nxi
          print *, "Expected:", nxi_
          xi = 0.0_dp
        endif
    end function get_xi

    function get_Q2(nQ2) result(Q2)
        integer, intent(in) :: nQ2
        real(dp) :: Q2(nQ2)
        if(nQ2==nQ2_) then
          Q2 = Q2_
        else
          print *, "ERROR: incorrect Q2 grid size."
          print *, "Got:", nQ2
          print *, "Expected:", nQ2_
          Q2 = 0.0_dp
        endif
    end function get_Q2

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Set methods

    subroutine initialize_x_xi(nx, nxi, xi_grid, grid_type, lagrange_order)
        integer,  intent(in) :: nx, nxi, grid_type, lagrange_order
        real(dp), intent(in) :: xi_grid(nxi)
        ! Deallocate any allocated arrays
        if(allocated(xx_)) deallocate(xx_)
        if(allocated(xi_)) deallocate(xi_)
        ! Cache passed integers
        nx_  = nx
        nxi_ = nxi
        gt_  = grid_type
        ! Allocate the cached arrays
        allocate(xx_(nx,nxi))
        allocate(xi_(nxi))
        ! Set the cached arrays
        xx_ = pixelspace(nx, nxi, xi_grid, grid_type)
        xi_ = xi_grid
        ! Move xi off of 0 or 1 to avoid singular behavior
        if(xi_(1)   <        epsilon) xi_(1)   =        epsilon
        if(xi_(nxi) > 1.0_dp-epsilon) xi_(nxi) = 1.0_dp-epsilon
        ! Finally, initialize the Lagrange weights
        call initialize_lagrange_weights(nx, 1+lagrange_order)
    end subroutine initialize_x_xi

    subroutine initialize_Q2(nQ2, Q2_array)
        integer,  intent(in) :: nQ2
        real(dp), intent(in) :: Q2_array(nQ2)
        ! Deallocate any allocated arrays
        if(allocated(Q2_)) deallocate(Q2_)
        ! Cache passed integers
        nQ2_ = nQ2
        ! Allocate the cached arrays
        allocate(Q2_(nQ2))
        ! Set the cached arrays
        Q2_ = Q2_array
    end subroutine initialize_Q2

end module matrices_common

! pixelation.f90
!
! by Adam Freese
! part of the package tiktaalik for GPD evolution
!
! Created on August 6, 2025, as a replacement for a previous file of the same name.

module pixelation
  use gridspace, only: pull_back

  implicit none
  private

  integer, parameter :: dp = kind(1d0)

  ! Interpolation stuff
  integer :: N_order_cache
  integer :: Nx_cache
  real(dp), dimension(:), allocatable :: w_table

  public :: interpixel, initialize_lagrange_weights

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! The interpixel

    function interpixel(n_pixels, i, x, xi, grid_type) result(f)
        ! Gives a function from polynomial interpolation of a pixel.
        ! One *MUST* call initialize_lagrange_weights before this!
        integer,  intent(in) :: n_pixels, i, grid_type
        real(dp), intent(in) :: x, xi
        real(dp) :: f
        !
        real(dp) :: eta
        ! Interpolate in the pulled back eta space
        eta = pull_back(x, xi, n_pixels, grid_type)
        if(grid_type==1) then
          call interpolate_lagrange_pixel_legacy(i, eta, f)
        else
          call interpolate_lagrange_pixel(i, eta, f)
        endif
    end function interpixel

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Initialization routines

    subroutine initialize_lagrange_weights(nx, norder)
        integer, intent(in) :: nx, norder
        !
        real(dp) :: dx
        integer :: i, j
        ! Keep track of nx and norder values
        N_order_cache = norder
        Nx_cache = nx
        ! Assume a linear spacing
        dx = 2. / real(nx-1)
        ! Fill the weight array
        if(allocated(w_table)) deallocate(w_table)
        allocate(w_table(0:norder))
        do i=0, norder, 1
          w_table(i) = 1.0_dp
          do j=0, norder, 1
            if(i.ne.j) w_table(i) = w_table(i) / (dx*real(i-j))
          end do
        end do
    end subroutine initialize_lagrange_weights

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Lagrange interpolation
    ! optimized for linearly spaced x grid

    subroutine interpolate_lagrange_pixel(ix, x, f)
        ! Piecewise polynomial interpolation from connecting Newton polynomials.
        real(dp), intent(in)  :: x
        integer,  intent(in)  :: ix
        real(dp), intent(out) :: f
        !
        integer :: iLoc, iStart, iEnd
        real(dp) :: x0
        iLoc = locate(x, Nx_cache)
        iStart = min( max(iLoc-(N_order_cache-0)/2, 0), Nx_cache-1-N_order_cache )
        iEnd   = iStart + N_order_cache
        ! New hack for grid_type=2 ... avoid DGLAP-ERBL crossing
        ! TODO: clean up later
        if(iStart < (nx_cache-1)/4 .and. iEnd > (nx_cache-1)/4) then
          if(iLoc < (nx_cache-1)/4) then
            iEnd = (nx_cache-1)/4
            iStart = iEnd - n_order_cache
          else
            iStart = (nx_cache-1)/4
            iEnd = iStart + n_order_cache
          endif
        elseif(iStart < 3*(nx_cache-1)/4 .and. iEnd > 3*(nx_cache-1)/4) then
          if(iLoc < 3*(nx_cache-1)/4) then
            iEnd = 3*(nx_cache-1)/4
            iStart = iEnd - n_order_cache
          else
            iStart = 3*(nx_cache-1)/4
            iEnd = iStart + n_order_cache
          endif
        endif
        if(ix < iStart .or. ix > iEnd) then
          f = 0.0_dp
          return
        endif
        x0 = -1. + real(2*iStart)/real(Nx_cache-1)
        f = lagrange_basis(x, x0, ix-iStart)
    end subroutine interpolate_lagrange_pixel

    subroutine interpolate_lagrange_pixel_legacy(ix, x, f)
        ! Piecewise polynomial interpolation from connecting Newton polynomials.
        real(dp), intent(in)  :: x
        integer,  intent(in)  :: ix
        real(dp), intent(out) :: f
        !
        integer :: iLoc, iStart, iEnd
        real(dp) :: x0
        iLoc = locate(x, Nx_cache)
        iStart = min( max(iLoc-(N_order_cache-1)/2, 0), Nx_cache+1-N_order_cache )
        iEnd   = iStart + N_order_cache
        if(ix < iStart .or. ix > iEnd) then
          f = 0.0_dp
          return
        endif
        x0 = -1. + real(2*iStart)/real(Nx_cache-1)
        f = lagrange_basis(x, x0, ix-iStart)
    end subroutine interpolate_lagrange_pixel_legacy

    function lagrange_basis(x, x0, i) result(f)
        real(dp), intent(in) :: x, x0
        integer,  intent(in) :: i
        real(dp) :: f
        !
        real(dp) :: dx, lfact(0:n_order_cache)
        integer :: j
        dx = 2. / real(nx_cache-1)
        do j=0, n_order_cache, 1
          lfact(j) = x - x0 - dx*real(j)
        end do
        lfact(i) = 1.
        f = w_table(i) * product(lfact)
    end function lagrange_basis

    function locate(x, nx) result(iloc)
        real(dp), intent(in) :: x
        integer,  intent(in) :: nx
        integer :: iloc
        !
        iloc = floor( 0.5*(real(nx-1)*(x+1.)) )
        if(iloc < 0)    iloc = 0
        if(iloc > nx-1) iloc = nx-1
    end function locate

end module pixelation

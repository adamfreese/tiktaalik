! model.f90
!
! by Adma Freese
! part of the package tiktaalik
!
! wrappers for f2py to access

module dummy
  use gk

  implicit none
  public

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Methods assuming independent x and xi arrays

    subroutine Hu_wrap(nx, nxi, nt, x, xi, t, f)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nt
        real(dp), intent(in)  :: x(nx), xi(nxi), t(nt)
        real(dp), intent(out) :: f(nx,nxi,nt)
        !
        integer :: ix, ixi, it
        !$OMP PARALLEL DO
        do ix=1, nx, 1
        do ixi=1, nxi, 1
        do it=1, nt, 1
          f(ix,ixi,it) = Hu(x(ix), xi(ixi), t(it))
        end do
        end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Hu_wrap

    subroutine Hd_wrap(nx, nxi, nt, x, xi, t, f)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nt
        real(dp), intent(in)  :: x(nx), xi(nxi), t(nt)
        real(dp), intent(out) :: f(nx,nxi,nt)
        !
        integer :: ix, ixi, it
        !$OMP PARALLEL DO
        do ix=1, nx, 1
        do ixi=1, nxi, 1
        do it=1, nt, 1
          f(ix,ixi,it) = Hd(x(ix), xi(ixi), t(it))
        end do
        end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Hd_wrap

    subroutine Hs_wrap(nx, nxi, nt, x, xi, t, f)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nt
        real(dp), intent(in)  :: x(nx), xi(nxi), t(nt)
        real(dp), intent(out) :: f(nx,nxi,nt)
        !
        integer :: ix, ixi, it
        !$OMP PARALLEL DO
        do ix=1, nx, 1
        do ixi=1, nxi, 1
        do it=1, nt, 1
          f(ix,ixi,it) = Hs(x(ix), xi(ixi), t(it))
        end do
        end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Hs_wrap

    subroutine Hg_wrap(nx, nxi, nt, x, xi, t, f)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nt
        real(dp), intent(in)  :: x(nx), xi(nxi), t(nt)
        real(dp), intent(out) :: f(nx,nxi,nt)
        !
        integer :: ix, ixi, it
        !$OMP PARALLEL DO
        do ix=1, nx, 1
        do ixi=1, nxi, 1
        do it=1, nt, 1
          f(ix,ixi,it) = Hg(x(ix), xi(ixi), t(it))
        end do
        end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Hg_wrap

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Methods allowing a 2D, xi-dependent x grid

    subroutine Hu_wrap_2d(nx, nxi, nt, x, xi, t, f)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nt
        real(dp), intent(in)  :: x(nx,nxi), xi(nxi), t(nt)
        real(dp), intent(out) :: f(nx,nxi,nt)
        !
        integer :: ix, ixi, it
        !$OMP PARALLEL DO
        do ix=1, nx, 1
        do ixi=1, nxi, 1
        do it=1, nt, 1
          f(ix,ixi,it) = Hu(x(ix,ixi), xi(ixi), t(it))
        end do
        end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Hu_wrap_2d

    subroutine Hd_wrap_2d(nx, nxi, nt, x, xi, t, f)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nt
        real(dp), intent(in)  :: x(nx,nxi), xi(nxi), t(nt)
        real(dp), intent(out) :: f(nx,nxi,nt)
        !
        integer :: ix, ixi, it
        !$OMP PARALLEL DO
        do ix=1, nx, 1
        do ixi=1, nxi, 1
        do it=1, nt, 1
          f(ix,ixi,it) = Hd(x(ix,ixi), xi(ixi), t(it))
        end do
        end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Hd_wrap_2d

    subroutine Hs_wrap_2d(nx, nxi, nt, x, xi, t, f)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nt
        real(dp), intent(in)  :: x(nx,nxi), xi(nxi), t(nt)
        real(dp), intent(out) :: f(nx,nxi,nt)
        !
        integer :: ix, ixi, it
        !$OMP PARALLEL DO
        do ix=1, nx, 1
        do ixi=1, nxi, 1
        do it=1, nt, 1
          f(ix,ixi,it) = Hs(x(ix,ixi), xi(ixi), t(it))
        end do
        end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Hs_wrap_2d

    subroutine Hg_wrap_2d(nx, nxi, nt, x, xi, t, f)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nt
        real(dp), intent(in)  :: x(nx,nxi), xi(nxi), t(nt)
        real(dp), intent(out) :: f(nx,nxi,nt)
        !
        integer :: ix, ixi, it
        !$OMP PARALLEL DO
        do ix=1, nx, 1
        do ixi=1, nxi, 1
        do it=1, nt, 1
          f(ix,ixi,it) = Hg(x(ix,ixi), xi(ixi), t(it))
        end do
        end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Hg_wrap_2d

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Deprecated methods

    !subroutine Hu_val_wrap(nx, nxi, nt, x, xi, t, f)
    !    integer,  parameter   :: dp = kind(1d0)
    !    integer,  intent(in)  :: nx, nxi, nt
    !    real(dp), intent(in)  :: x(nx), xi(nxi), t(nt)
    !    real(dp), intent(out) :: f(nx,nxi,nt)
    !    !
    !    integer :: ix, ixi, it
    !    !$OMP PARALLEL DO
    !    do ix=1, nx, 1
    !    do ixi=1, nxi, 1
    !    do it=1, nt, 1
    !      f(ix,ixi,it) = Hu_val(x(ix), xi(ixi), t(it))
    !    end do
    !    end do
    !    end do
    !    !$OMP END PARALLEL DO
    !end subroutine Hu_val_wrap

    !subroutine Hd_val_wrap(nx, nxi, nt, x, xi, t, f)
    !    integer,  parameter   :: dp = kind(1d0)
    !    integer,  intent(in)  :: nx, nxi, nt
    !    real(dp), intent(in)  :: x(nx), xi(nxi), t(nt)
    !    real(dp), intent(out) :: f(nx,nxi,nt)
    !    !
    !    integer :: ix, ixi, it
    !    !$OMP PARALLEL DO
    !    do ix=1, nx, 1
    !    do ixi=1, nxi, 1
    !    do it=1, nt, 1
    !      f(ix,ixi,it) = Hd_val(x(ix), xi(ixi), t(it))
    !    end do
    !    end do
    !    end do
    !    !$OMP END PARALLEL DO
    !end subroutine Hd_val_wrap

    !subroutine Hu_sea_wrap(nx, nxi, nt, x, xi, t, f)
    !    integer,  parameter   :: dp = kind(1d0)
    !    integer,  intent(in)  :: nx, nxi, nt
    !    real(dp), intent(in)  :: x(nx), xi(nxi), t(nt)
    !    real(dp), intent(out) :: f(nx,nxi,nt)
    !    !
    !    integer :: ix, ixi, it
    !    !$OMP PARALLEL DO
    !    do ix=1, nx, 1
    !    do ixi=1, nxi, 1
    !    do it=1, nt, 1
    !      f(ix,ixi,it) = Hu_sea(x(ix), xi(ixi), t(it))
    !    end do
    !    end do
    !    end do
    !    !$OMP END PARALLEL DO
    !end subroutine Hu_sea_wrap

    !subroutine Hd_sea_wrap(nx, nxi, nt, x, xi, t, f)
    !    integer,  parameter   :: dp = kind(1d0)
    !    integer,  intent(in)  :: nx, nxi, nt
    !    real(dp), intent(in)  :: x(nx), xi(nxi), t(nt)
    !    real(dp), intent(out) :: f(nx,nxi,nt)
    !    !
    !    integer :: ix, ixi, it
    !    !$OMP PARALLEL DO
    !    do ix=1, nx, 1
    !    do ixi=1, nxi, 1
    !    do it=1, nt, 1
    !      f(ix,ixi,it) = Hd_sea(x(ix), xi(ixi), t(it))
    !    end do
    !    end do
    !    end do
    !    !$OMP END PARALLEL DO
    !end subroutine Hd_sea_wrap

end module dummy

! matrices_wilson.f90
!
! by Adam Freese
! part of the package tiktaalik for GPD evolution
!
! Created July 24, 2025
! Part of a breakup of the previous matevo.f90 file into smaller chunks

module matrices_wilson
  use alpha_qcd, only: get_alpha_QCD, get_neff, get_charges_squared
  use constants, only: i_, pi
  use kernels_common
  use matrices_common
  use wilson_dvcs

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)

  ! Wilson coefficient matrices
  ! Indices are for (xi,x)
  real(dp), allocatable, dimension(:,:) :: re_Cq0_dvcs, im_Cq0_dvcs, &
      & re_Cq1_dvcs, im_Cq1_dvcs, re_CG1_dvcs, im_CG1_dvcs

  public :: make_wilson_matrices, &
      & Cq_dvcs, CG_dvcs

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Public methods to build the Wilson coefficient matrices

    subroutine make_wilson_matrices(nQ2)
        ! Public method to initialize the Wilson coefficient matrices
        ! IMPORTANT: the following must be called first!:
        ! - initialize_x_xi [in matrices_common]
        ! - initialize_Q2   [in matrices_common]
        ! The Python interface should automatically do this.
        integer, intent(in) :: nQ2
        !
        real(dp) :: Q2(nQ2)
        integer  :: nx, nxi, grid_type
        ! Retrieve relevant cahced variables
        Q2  = get_Q2(nQ2)
        nx  = get_nx()
        nxi = get_nxi()
        grid_type = get_grid_type()
        ! Make Wilson coefficient matrices
        call make_dvcs_Cq0(nx, nxi, nQ2, Q2, grid_type)
        call make_dvcs_Cq1(nx, nxi, nQ2, Q2, grid_type)
        call make_dvcs_CG1(nx, nxi, nQ2, Q2, grid_type)
    end subroutine make_wilson_matrices

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Public methods to return Wilson coefficient matrices

    function Cq_dvcs(nxi, nx, nQ2, l_nlo) result(M)
        integer,  intent(in) :: nxi, nx, nQ2
        logical,  intent(in) :: l_nlo
        complex(dp), dimension(nxi, nx, nQ2) :: M
        complex(dp), dimension(nxi, nx) :: M0, M1
        real(dp) :: Q2(nQ2)
        integer  :: iq
        Q2 = get_Q2(nQ2)
        M0 = re_Cq0_dvcs + i_*im_Cq0_dvcs
        M1 = re_Cq1_dvcs + i_*im_Cq1_dvcs
        !$OMP PARALLEL DO
        do iq=1, nQ2, 1
          if(l_nlo) then
            M(:,:,iq) = M0 + get_alpha_QCD(Q2(iq))/(4.*pi)*M1
          else
            M(:,:,iq) = M0
          endif
        end do
        !$OMP END PARALLEL DO
    end function Cq_dvcs

    function CG_dvcs(nxi, nx, nQ2, l_nlo) result(M)
        integer,  intent(in) :: nxi, nx, nQ2
        logical,  intent(in) :: l_nlo
        complex(dp), dimension(nxi, nx, nQ2) :: M
        complex(dp), dimension(nxi, nx) :: M0, M1
        real(dp) :: Q2(nQ2)
        integer  :: iq
        Q2 = get_Q2(nQ2)
        ! For now, only LO
        M1 = re_CG1_dvcs + i_*im_CG1_dvcs
        if(l_nlo) then
          !$OMP PARALLEL DO
          do iq=1, nQ2, 1
            M(:,:,iq) = get_charges_squared(Q2(iq))*get_alpha_QCD(Q2(iq))/(4.*pi)*M1
          end do
          !$OMP END PARALLEL DO
        else
          M = 0.0_dp
        endif
    end function CG_dvcs

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Private methods to make Wilson coefficient matrices

    subroutine make_dvcs_Cq0(nx, nxi, nQ2, Q2_array, grid_type)
        integer,  intent(in) :: nx, nxi, nQ2, grid_type
        real(dp), intent(in) :: Q2_array(nQ2)
        real(dp) :: xi(nxi)
        integer  :: ix, iz, iq, ic
        xi = get_xi(nxi)
        if(allocated(re_Cq0_dvcs)) deallocate(re_Cq0_dvcs)
        if(allocated(im_Cq0_dvcs)) deallocate(im_Cq0_dvcs)
        allocate(re_Cq0_dvcs(nxi,nx))
        allocate(im_Cq0_dvcs(nxi,nx))
        do iz=1, nxi, 1
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          re_Cq0_dvcs(iz,ix) = pixel_coef(zero_func, re_Cq0_dvcs_sub, re_Cq0_dvcs_cst, xi(iz), nx, ix, grid_type)
          im_Cq0_dvcs(iz,ix) = pixel_coef(zero_func, zero_func,       im_Cq0_dvcs_cst, xi(iz), nx, ix, grid_type)
        end do
        !$OMP END PARALLEL DO
        end do
    end subroutine make_dvcs_Cq0

    subroutine make_dvcs_Cq1(nx, nxi, nQ2, Q2_array, grid_type)
        integer,  intent(in) :: nx, nxi, nQ2, grid_type
        real(dp), intent(in) :: Q2_array(nQ2)
        real(dp) :: xi(nxi)
        integer  :: ix, iz, iq, ic
        xi = get_xi(nxi)
        if(allocated(re_Cq1_dvcs)) deallocate(re_Cq1_dvcs)
        if(allocated(im_Cq1_dvcs)) deallocate(im_Cq1_dvcs)
        allocate(re_Cq1_dvcs(nxi,nx))
        allocate(im_Cq1_dvcs(nxi,nx))
        do iz=1, nxi, 1
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          re_Cq1_dvcs(iz,ix) = pixel_coef(re_Cq1_dvcs_reg, re_Cq1_dvcs_sub, re_Cq1_dvcs_cst, xi(iz), nx, ix, grid_type)
          im_Cq1_dvcs(iz,ix) = pixel_coef(im_Cq1_dvcs_reg, im_Cq1_dvcs_sub, im_Cq1_dvcs_cst, xi(iz), nx, ix, grid_type)
        end do
        !$OMP END PARALLEL DO
        end do
    end subroutine make_dvcs_Cq1

    subroutine make_dvcs_CG1(nx, nxi, nQ2, Q2_array, grid_type)
        integer,  intent(in) :: nx, nxi, nQ2, grid_type
        real(dp), intent(in) :: Q2_array(nQ2)
        real(dp) :: xi(nxi)
        integer  :: ix, iz, iq, ic
        xi = get_xi(nxi)
        if(allocated(re_CG1_dvcs)) deallocate(re_CG1_dvcs)
        if(allocated(im_CG1_dvcs)) deallocate(im_CG1_dvcs)
        allocate(re_CG1_dvcs(nxi,nx))
        allocate(im_CG1_dvcs(nxi,nx))
        do iz=1, nxi, 1
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          re_CG1_dvcs(iz,ix) = pixel_coef(re_CG1_dvcs_reg, re_CG1_dvcs_sub, re_CG1_dvcs_cst, xi(iz), nx, ix, grid_type)
          im_CG1_dvcs(iz,ix) = pixel_coef(im_CG1_dvcs_reg, im_CG1_dvcs_sub, im_CG1_dvcs_cst, xi(iz), nx, ix, grid_type)
        end do
        !$OMP END PARALLEL DO
        end do
    end subroutine make_dvcs_CG1

end module matrices_wilson

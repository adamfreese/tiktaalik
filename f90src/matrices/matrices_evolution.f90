! matrices_evolution.f90
!
! by Adam Freese
! part of the package tiktaalik for GPD evolution
!
! Created July 24, 2025
! Part of a breakup of the previous matevo.f90 file into smaller chunks

module matrices_evolution
  use alpha_qcd,  only: get_alpha_QCD, get_neff
  use constants,  only: pi
  use convolution
  use kernels_common
  use kernels_lo
  use kernels_nlo
  use matrices_common

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)

  ! Flavor limits
  integer, parameter, private :: nfl_min = 3
  integer, parameter, private :: nfl_max = 5

  ! Kernel matrices
  ! Indices are for (x,y,xi,nfl)
  ! For singlet/gluon, the first two indices go up to 2*nx,
  ! with the first nx values being singlet quark and the last nx being gluon.
  ! LO
  real(dp), allocatable, dimension(:,:,:,:) :: K_NS_0, KV_SG_0, KA_SG_0
  ! NLO
  real(dp), allocatable, dimension(:,:,:,:) :: KV_NS_1p, KV_NS_1m, KV_SG_1, KA_SG_1

  ! Evolution matirces
  ! Indices are for (x,y,xi,Q2)
  ! For singlet/gluon, the first two indices go up to 2*nx,
  ! with the first nx values being singlet quark and the last nx being gluon.
  real(dp), allocatable, dimension(:,:,:,:) :: M_NS_pls, M_NS_min, MV_SG, MA_SG

  public :: make_kernels, make_evolution_matrices, &
      & evomat_V_NS, evomat_V_SG, evomat_A_NS, evomat_A_SG, &
      & kernel_V_QQ, kernel_V_QG, kernel_V_GQ, kernel_V_GG, &
      & kernel_A_QQ, kernel_A_QG, kernel_A_GQ, kernel_A_GG

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Public methods to build the kernels and evolution matrices.

    subroutine make_kernels()
        ! Public method to initialize the kernel matrices
        ! IMPORTANT: the following must be called first!:
        ! - initialize_x_xi [in matrices_common]
        ! The Python interface should automatically do this.
        integer :: nx, nxi, grid_type
        ! Retrieve relevant information from matrices_common
        nx  = get_nx()
        nxi = get_nxi()
        grid_type = get_grid_type()
        ! Initialize the kernel matrices
        call make_kernels_NS_0(nx, nxi, grid_type)
        call make_kernels_NS_1(nx, nxi, grid_type)
        call make_kernels_SG_0(nx, nxi, grid_type)
        call make_kernels_SG_1(nx, nxi, grid_type)
    end subroutine make_kernels

    subroutine make_evolution_matrices(nQ2, l_nlo)
        ! Public method to initialize the evolution matrices
        ! IMPORTANT: the following must be called first!:
        ! - initialize_x_xi [in matrices_common]
        ! - initialize_Q2   [in matrices_common]
        ! - make_kernels
        ! The Python interface should automatically do this.
        integer, intent(in) :: nQ2
        logical, intent(in) :: l_nlo
        !
        real(dp) :: Q2(nQ2)
        integer :: nx, nxi
        ! Retrieve relevant cahced variables
        Q2  = get_Q2(nQ2)
        nx  = get_nx()
        nxi = get_nxi()
        ! Make evolution matrices
        print *, "Flag A"
        call make_evomat_NS(  nx, nxi, nQ2, Q2, l_nlo)
        print *, "Flag B"
        call make_evomat_V_SG(nx, nxi, nQ2, Q2, l_nlo)
        print *, "Flag C"
        call make_evomat_A_SG(nx, nxi, nQ2, Q2, l_nlo)
        print *, "Flag D"
    end subroutine make_evolution_matrices

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Public methods to return kenrel matrices

    function kernel_V_QQ(Q2, nx, nxi, nfl, l_nlo, i_ns_type) result(K)
        real(dp), intent(in) :: Q2
        integer,  intent(in) :: nx, nxi, nfl, i_ns_type
        logical,  intent(in) :: l_nlo
        real(dp), dimension(nx,nx,nxi) :: K
        real(dp) :: al2pi
        al2pi = get_alpha_QCD(Q2) / (2.*pi)
        K = al2pi * K_NS_0(:,:,:,nfl)
        if(l_nlo) then
          select case(i_ns_type)
          case(1)
            K = K + al2pi**2 * KV_NS_1p(:,:,:,nfl)
          case(-1)
            K = K + al2pi**2 * KV_NS_1m(:,:,:,nfl)
          case(0)
            K = K + al2pi**2 * KV_SG_1(1:nx,1:nx,:,nfl)
          end select
        endif
    end function kernel_V_QQ

    function kernel_V_QG(Q2, nx, nxi, nfl, l_nlo) result(K)
        real(dp), intent(in) :: Q2
        integer,  intent(in) :: nx, nxi, nfl
        logical,  intent(in) :: l_nlo
        real(dp), dimension(nx,nx,nxi) :: K
        real(dp) :: al2pi
        al2pi = get_alpha_QCD(Q2) / (2.*pi)
        K = al2pi * KV_SG_0(1:nx,nx+1:2*nx,:,nfl)
        if(l_nlo) then
          K = K + al2pi**2 * KV_SG_1(1:nx,nx+1:2*nx,:,nfl)
        endif
    end function kernel_V_QG

    function kernel_V_GQ(Q2, nx, nxi, nfl, l_nlo) result(K)
        real(dp), intent(in) :: Q2
        integer,  intent(in) :: nx, nxi, nfl
        logical,  intent(in) :: l_nlo
        real(dp), dimension(nx,nx,nxi) :: K
        real(dp) :: al2pi
        al2pi = get_alpha_QCD(Q2) / (2.*pi)
        K = al2pi * KV_SG_0(nx+1:2*nx,1:nx,:,nfl)
        if(l_nlo) then
          K = K + al2pi**2 * KV_SG_1(nx+1:2*nx,1:nx,:,nfl)
        endif
    end function kernel_V_GQ

    function kernel_V_GG(Q2, nx, nxi, nfl, l_nlo) result(K)
        real(dp), intent(in) :: Q2
        integer,  intent(in) :: nx, nxi, nfl
        logical,  intent(in) :: l_nlo
        real(dp), dimension(nx,nx,nxi) :: K
        real(dp) :: al2pi
        al2pi = get_alpha_QCD(Q2) / (2.*pi)
        K = al2pi * KV_SG_0(nx+1:2*nx,nx+1:2*nx,:,nfl)
        if(l_nlo) then
          K = K + al2pi**2 * KV_SG_1(nx+1:2*nx,nx+1:2*nx,:,nfl)
        endif
    end function kernel_V_GG

    function kernel_A_QQ(Q2, nx, nxi, nfl, l_nlo, i_ns_type) result(K)
        ! NS kernels are the same as for V, except + and - are swapped.
        real(dp), intent(in) :: Q2
        integer,  intent(in) :: nx, nxi, nfl, i_ns_type
        logical,  intent(in) :: l_nlo
        real(dp), dimension(nx,nx,nxi) :: K
        real(dp) :: al2pi
        al2pi = get_alpha_QCD(Q2) / (2.*pi)
        K = al2pi * K_NS_0(:,:,:,nfl)
        if(l_nlo) then
          select case(i_ns_type)
          case(1)
            K = K + al2pi**2 * KV_NS_1m(:,:,:,nfl)
          case(-1)
            K = K + al2pi**2 * KV_NS_1p(:,:,:,nfl)
          case(0)
            K = K + al2pi**2 * KA_SG_1(1:nx,1:nx,:,nfl)
          end select
        endif
    end function kernel_A_QQ

    function kernel_A_QG(Q2, nx, nxi, nfl, l_nlo) result(K)
        real(dp), intent(in) :: Q2
        integer,  intent(in) :: nx, nxi, nfl
        logical,  intent(in) :: l_nlo
        real(dp), dimension(nx,nx,nxi) :: K
        real(dp) :: al2pi
        al2pi = get_alpha_QCD(Q2) / (2.*pi)
        K = al2pi * KA_SG_0(1:nx,nx+1:2*nx,:,nfl)
        if(l_nlo) then
          K = K + al2pi**2 * KA_SG_1(1:nx,nx+1:2*nx,:,nfl)
        endif
    end function kernel_A_QG

    function kernel_A_GQ(Q2, nx, nxi, nfl, l_nlo) result(K)
        real(dp), intent(in) :: Q2
        integer,  intent(in) :: nx, nxi, nfl
        logical,  intent(in) :: l_nlo
        real(dp), dimension(nx,nx,nxi) :: K
        real(dp) :: al2pi
        al2pi = get_alpha_QCD(Q2) / (2.*pi)
        K = al2pi * KA_SG_0(nx+1:2*nx,1:nx,:,nfl)
        if(l_nlo) then
          K = K + al2pi**2 * KA_SG_1(nx+1:2*nx,1:nx,:,nfl)
        endif
    end function kernel_A_GQ

    function kernel_A_GG(Q2, nx, nxi, nfl, l_nlo) result(K)
        real(dp), intent(in) :: Q2
        integer,  intent(in) :: nx, nxi, nfl
        logical,  intent(in) :: l_nlo
        real(dp), dimension(nx,nx,nxi) :: K
        real(dp) :: al2pi
        al2pi = get_alpha_QCD(Q2) / (2.*pi)
        K = al2pi * KA_SG_0(nx+1:2*nx,nx+1:2*nx,:,nfl)
        if(l_nlo) then
          K = K + al2pi**2 * KA_SG_1(nx+1:2*nx,nx+1:2*nx,:,nfl)
        endif
    end function kernel_A_GG

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Public methods to return evolution matrices

    function evomat_V_NS(nx, nxi, nQ2, nstype) result(M)
        integer,  intent(in) :: nx, nxi, nQ2, nstype
        real(dp), dimension(nx, nx, nxi, nQ2) :: M
        select case(nstype)
        case(1)
          M = M_NS_pls
        case(-1)
          M = M_NS_min
        end select
    end function evomat_V_NS

    function evomat_V_SG(nx, nxi, nQ2) result(M)
        integer,  intent(in) :: nx, nxi, nQ2
        real(dp), dimension(2*nx, 2*nx, nxi, nQ2) :: M
        M = MV_SG
    end function evomat_V_SG

    function evomat_A_NS(nx, nxi, nQ2, nstype) result(M)
        integer,  intent(in) :: nx, nxi, nQ2, nstype
        real(dp), dimension(nx, nx, nxi, nQ2) :: M
        select case(nstype)
        case(1)
          M = M_NS_min
        case(-1)
          M = M_NS_pls
        end select
    end function evomat_A_NS

    function evomat_A_SG(nx, nxi, nQ2) result(M)
        integer,  intent(in) :: nx, nxi, nQ2
        real(dp), dimension(2*nx, 2*nx, nxi, nQ2) :: M
        M = MA_SG
    end function evomat_A_SG

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Private methods to make kernel matrices

    subroutine make_kernels_NS_0(nx, nxi, grid_type)
        integer, intent(in) :: nx, nxi, grid_type
        integer  :: ix, iy, iz!, nfl
        real(dp) :: xi(nxi)
        if(allocated(K_NS_0)) deallocate(K_NS_0)
        allocate(K_NS_0(nx,nx,nxi,nfl_min:nfl_max))
        xi = get_xi(nxi)
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          do iy=1, nx, 1
            do iz=1, nxi, 1
              ! At leading order, no nfl dependence in QQ kernel.
              K_NS_0(ix,iy,iz,:) = pixel_conv(zero_func, K0_qq_pls, K0_qq_cst, xi(iz), nx, ix, iy, grid_type)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
    end subroutine make_kernels_NS_0

    subroutine make_kernels_NS_1(nx, nxi, grid_type)
        integer, intent(in) :: nx, nxi, grid_type
        if(allocated(KV_NS_1p)) deallocate(KV_NS_1p)
        if(allocated(KV_NS_1m)) deallocate(KV_NS_1m)
        allocate(KV_NS_1p(nx,nx,nxi,nfl_min:nfl_max))
        allocate(KV_NS_1m(nx,nx,nxi,nfl_min:nfl_max))
        call make_one_nlo_ns_kernel(nx, nxi, grid_type, &
            & zero_func, KV1_NSp_pls, KV1_NSp_cst, zero_func, KV1_NSp_pls_nfl, KV1_NSp_cst_nfl, &
            & KV_NS_1p)
        call make_one_nlo_ns_kernel(nx, nxi, grid_type, &
            & zero_func, KV1_NSm_pls, KV1_NSm_cst, zero_func, KV1_NSm_pls_nfl, KV1_NSm_cst_nfl, &
            & KV_NS_1m)
    end subroutine make_kernels_NS_1

    subroutine make_one_nlo_ns_kernel(nx, nxi, grid_type, freg0, fpls0, fcst0, freg1, fpls1, fcst1, kernel)
        integer, intent(in)   :: nx, nxi, grid_type
        real(dp), intent(out) :: kernel(nx,nx,nxi,nfl_min:nfl_max)
        real(dp), external :: freg0, fpls0, fcst0, freg1, fpls1, fcst1
        integer  :: ix, iy, iz, nfl
        real(dp) :: xi(nxi)
        real(dp), allocatable, dimension(:,:,:,:) :: k0, k1
        allocate(k0(nx,nx,nxi,nfl_min:nfl_max))
        allocate(k1(nx,nx,nxi,nfl_min:nfl_max))
        xi = get_xi(nxi)
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          do iy=1, nx, 1
            do iz=1, nxi, 1
              k0(ix,iy,iz,:) = pixel_conv(freg0, fpls0, fcst0, xi(iz), nx, ix, iy, grid_type)
              k1(ix,iy,iz,:) = pixel_conv(freg1, fpls1, fcst1, xi(iz), nx, ix, iy, grid_type)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
        do nfl=nfl_min, nfl_max, 1
          k1(:,:,:,nfl) = k1(:,:,:,nfl) * real(nfl)
        end do
        kernel = k0 + k1
        deallocate(k0, k1)
    end subroutine make_one_nlo_ns_kernel

    subroutine make_kernels_SG_0(nx, nxi, grid_type)
        integer, intent(in) :: nx, nxi, grid_type
        integer  :: ix, iy, iz, nfl
        real(dp) :: xi(nxi)
        real(dp), dimension(:,:,:), allocatable :: qq_nfl_0, qG_nfl_1, Gq_nfl_0, GG_nfl_0, GG_nfl_1
        if(allocated(KA_SG_0)) deallocate(KA_SG_0)
        if(allocated(KV_SG_0)) deallocate(KV_SG_0)
        allocate(KA_SG_0(2*nx,2*nx,nxi,nfl_min:nfl_max))
        allocate(KV_SG_0(2*nx,2*nx,nxi,nfl_min:nfl_max))
        ! Temporary sub-matrix arrays
        allocate(qq_nfl_0(nx,nx,nxi))
        allocate(qG_nfl_1(nx,nx,nxi))
        allocate(Gq_nfl_0(nx,nx,nxi))
        allocate(GG_nfl_0(nx,nx,nxi))
        allocate(GG_nfl_1(nx,nx,nxi))
        KA_SG_0 = 0.0_dp
        KV_SG_0 = 0.0_dp
        xi = get_xi(nxi)
        ! First, build up the A-type kernel
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          do iy=1, nx, 1
            do iz=1, nxi, 1
              qG_nfl_1(ix,iy,iz) = pixel_conv(KA0_qG_reg, zero_func,  zero_func,  xi(iz), nx, ix, iy, grid_type)
              Gq_nfl_0(ix,iy,iz) = pixel_conv(KA0_Gq_reg, zero_func,  zero_func,  xi(iz), nx, ix, iy, grid_type)
              GG_nfl_0(ix,iy,iz) = pixel_conv(zero_func,  KA0_GG_pls, KA0_GG_cst, xi(iz), nx, ix, iy, grid_type)
              GG_nfl_1(ix,iy,iz) = pixel_conv(zero_func,  zero_func,  KA0_GG_nfl, xi(iz), nx, ix, iy, grid_type)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
        ! For the QQ block, just copy over the NS kernel, since it's the same at LO.
        ! The NS kernel matrix is built before this method is called so this is safe.
        KA_SG_0(1:nx,1:nx,:,:) = K_NS_0(:,:,:,:)
        do nfl=nfl_min, nfl_max, 1
          !!KA_SG_0(1:nx,     1:nx,     :,nfl) = qq_nfl_0(:,:,:)
          KA_SG_0(1:nx,     nx+1:2*nx,:,nfl) = real(nfl)*qG_nfl_1(:,:,:)
          KA_SG_0(nx+1:2*nx,1:nx,     :,nfl) = Gq_nfl_0(:,:,:)
          KA_SG_0(nx+1:2*nx,nx+1:2*nx,:,nfl) = GG_nfl_0(:,:,:) + real(nfl)*GG_nfl_1(:,:,:)
        end do
        ! Next, the extra pieces for the V-type kernel
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          do iy=1, nx, 1
            do iz=1, nxi, 1
              qG_nfl_1(ix,iy,iz) = pixel_conv(KVmA0_qG_reg, zero_func, zero_func, xi(iz), nx, ix, iy, grid_type)
              Gq_nfl_0(ix,iy,iz) = pixel_conv(KVmA0_Gq_reg, zero_func, zero_func, xi(iz), nx, ix, iy, grid_type)
              GG_nfl_0(ix,iy,iz) = pixel_conv(KVmA0_GG_reg, zero_func, zero_func, xi(iz), nx, ix, iy, grid_type)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
        KV_SG_0(1:nx,1:nx,:,:) = 0.0_dp
        do nfl=nfl_min, nfl_max, 1
          KV_SG_0(1:nx,     nx+1:2*nx,:,nfl) = real(nfl)*qG_nfl_1(:,:,:)
          KV_SG_0(nx+1:2*nx,1:nx,     :,nfl) = Gq_nfl_0(:,:,:)
          KV_SG_0(nx+1:2*nx,nx+1:2*nx,:,nfl) = GG_nfl_0(:,:,:)
        end do
        ! And add the extra pieces to the A-type kernel
        KV_SG_0 = KV_SG_0 + KA_SG_0
        ! And we're done
        deallocate(qq_nfl_0)
        deallocate(qG_nfl_1)
        deallocate(Gq_nfl_0)
        deallocate(GG_nfl_0)
        deallocate(GG_nfl_1)
    end subroutine make_kernels_SG_0

    subroutine make_kernels_SG_1(nx, nxi, grid_type)
        integer, intent(in) :: nx, nxi, grid_type
        integer  :: ix, iy, iz, nfl
        real(dp) :: xi(nxi)
        real(dp), dimension(:,:,:), allocatable :: qq_nfl_1, qG_nfl_1, Gq_nfl_0, Gq_nfl_1, GG_nfl_0, GG_nfl_1
        if(allocated(KA_SG_1)) deallocate(KA_SG_1)
        if(allocated(KV_SG_1)) deallocate(KV_SG_1)
        allocate(KA_SG_1(2*nx,2*nx,nxi,nfl_min:nfl_max))
        allocate(KV_SG_1(2*nx,2*nx,nxi,nfl_min:nfl_max))
        allocate(qq_nfl_1(nx,nx,nxi))
        allocate(qG_nfl_1(nx,nx,nxi))
        allocate(Gq_nfl_0(nx,nx,nxi))
        allocate(Gq_nfl_1(nx,nx,nxi))
        allocate(GG_nfl_0(nx,nx,nxi))
        allocate(GG_nfl_1(nx,nx,nxi))
        KA_SG_1 = 0.0_dp
        KV_SG_1 = 0.0_dp
        xi = get_xi(nxi)
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !  V-type kernel first
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          do iy=1, nx, 1
            do iz=1, nxi, 1
              qq_nfl_1(ix,iy,iz) = pixel_conv(KV1_qq_reg,     zero_func,      zero_func,      xi(iz), nx, ix, iy, grid_type)
              qG_nfl_1(ix,iy,iz) = pixel_conv(KV1_qG_reg,     zero_func,      zero_func,      xi(iz), nx, ix, iy, grid_type)
              Gq_nfl_0(ix,iy,iz) = pixel_conv(KV1_Gq_reg,     zero_func,      zero_func,      xi(iz), nx, ix, iy, grid_type)
              Gq_nfl_1(ix,iy,iz) = pixel_conv(KV1_Gq_reg_nfl, zero_func,      zero_func,      xi(iz), nx, ix, iy, grid_type)
              GG_nfl_0(ix,iy,iz) = pixel_conv(zero_func,      KV1_GG_pls,     KV1_GG_cst,     xi(iz), nx, ix, iy, grid_type)
              GG_nfl_1(ix,iy,iz) = pixel_conv(zero_func,      KV1_GG_pls_nfl, KV1_GG_cst_nfl, xi(iz), nx, ix, iy, grid_type)
            end do
          end do
        end do
        ! Add nfl-independent and nfl-linear terms as needed.
        do nfl=nfl_min, nfl_max, 1
          KV_SG_1(1:nx,     1:nx,:,     nfl) = real(nfl)*qq_nfl_1(:,:,:)
          KV_SG_1(1:nx,     nx+1:2*nx,:,nfl) = real(nfl)*qG_nfl_1(:,:,:)
          KV_SG_1(nx+1:2*nx,1:nx,     :,nfl) = Gq_nfl_0(:,:,:) + real(nfl)*Gq_nfl_1(:,:,:)
          KV_SG_1(nx+1:2*nx,nx+1:2*nx,:,nfl) = GG_nfl_0(:,:,:) + real(nfl)*GG_nfl_1(:,:,:)
        end do
        ! For the QQ block, the plus-type NS kernel is a contribution, so we need to add it.
        ! (See Eq. (182).)
        ! The NS kernel matrix is built before this method is called so this is safe.
        KV_SG_1(1:nx,1:nx,:,:) = KV_SG_1(1:nx,1:nx,:,:) + KV_NS_1p(:,:,:,:)
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !  A-type kernel next
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          do iy=1, nx, 1
            do iz=1, nxi, 1
              qq_nfl_1(ix,iy,iz) = pixel_conv(KA1_qq_reg,     zero_func,      zero_func,      xi(iz), nx, ix, iy, grid_type)
              qG_nfl_1(ix,iy,iz) = pixel_conv(KA1_qG_reg,     zero_func,      zero_func,      xi(iz), nx, ix, iy, grid_type)
              Gq_nfl_0(ix,iy,iz) = pixel_conv(KA1_Gq_reg,     zero_func,      zero_func,      xi(iz), nx, ix, iy, grid_type)
              Gq_nfl_1(ix,iy,iz) = pixel_conv(KA1_Gq_reg_nfl, zero_func,      zero_func,      xi(iz), nx, ix, iy, grid_type)
              GG_nfl_0(ix,iy,iz) = pixel_conv(zero_func,      KA1_GG_pls,     KA1_GG_cst,     xi(iz), nx, ix, iy, grid_type)
              GG_nfl_1(ix,iy,iz) = pixel_conv(zero_func,      KA1_GG_pls_nfl, KA1_GG_cst_nfl, xi(iz), nx, ix, iy, grid_type)
            end do
          end do
        end do
        ! Add nfl-independent and nfl-linear terms as needed.
        do nfl=nfl_min, nfl_max, 1
          KA_SG_1(1:nx,     1:nx,:,     nfl) = real(nfl)*qq_nfl_1(:,:,:)
          KA_SG_1(1:nx,     nx+1:2*nx,:,nfl) = real(nfl)*qG_nfl_1(:,:,:)
          KA_SG_1(nx+1:2*nx,1:nx,     :,nfl) = Gq_nfl_0(:,:,:) + real(nfl)*Gq_nfl_1(:,:,:)
          KA_SG_1(nx+1:2*nx,nx+1:2*nx,:,nfl) = GG_nfl_0(:,:,:) + real(nfl)*GG_nfl_1(:,:,:)
        end do
        ! For the QQ block, the minus-type NS kernel is a contribution, so we need to add it.
        ! (See Eq. (178).)
        ! The NS kernel matrix is built before this method is called so this is safe.
        KA_SG_1(1:nx,1:nx,:,:) = KA_SG_1(1:nx,1:nx,:,:) + KV_NS_1m(:,:,:,:)
        ! And we're done
        deallocate(qq_nfl_1)
        deallocate(qG_nfl_1)
        deallocate(Gq_nfl_0)
        deallocate(Gq_nfl_1)
        deallocate(GG_nfl_0)
        deallocate(GG_nfl_1)
    end subroutine make_kernels_SG_1

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Private to make evolution matrices

    subroutine make_evomat_NS(nx, nxi, nQ2, Q2, l_nlo)
        ! Just one non-singlet evolution matrix at leading order.
        ! I'll need to break this into multiple routines,
        ! for both V-type and A-type, but also plus-type and minus-type,
        ! after going to NLO. But I'll burn that bridge when I come to it.
        integer,  intent(in) :: nx, nxi, nQ2
        real(dp), intent(in) :: Q2(nQ2)
        logical,  intent(in) :: l_nlo
        real(dp), dimension(nx,nx) :: idnx
        integer :: ix, ixi, iQ2, nfl
        print *, "Flag A1"
        if(allocated(M_NS_pls)) deallocate(M_NS_pls)
        if(allocated(M_NS_min)) deallocate(M_NS_min)
        print *, "Flag A2"
        allocate(M_NS_pls(nx,nx,nxi,nQ2))
        allocate(M_NS_min(nx,nx,nxi,nQ2))
        print *, "Flag A3"
        ! Build an identity matrix  ... TODO REPLACE
        idnx = 0.0_dp
        do ix=1, nx, 1
          idnx(ix,ix) = 1.0_dp
        end do
        ! Identity for no evolution at initial scale
        M_NS_pls = 0.0_dp
        M_NS_min = 0.0_dp
        do ixi=1, nxi, 1
          M_NS_pls(:,:,ixi,1) = idnx
          M_NS_min(:,:,ixi,1) = idnx
        end do
        print *, "Flag A4"
        ! Build evolution matrices for all other scales
        do iQ2=2, nQ2, 1
          print *, "Flag AQ", iQ2, nQ2
          !$OMP PARALLEL DO
          do ixi=1, nxi, 1
            ! First, an evolution matrix for prior Q2 step to current Q2 step
            nfl = get_neff(Q2(iQ2-1))
            M_NS_pls(:,:,ixi,iQ2) = idnx + &
                & rk4_NS(nx, nxi, Q2(iQ2-1), Q2(iQ2), K_NS_0(:,:,ixi,nfl), KV_NS_1p(:,:,ixi,nfl), l_nlo)
            M_NS_min(:,:,ixi,iQ2) = idnx + &
                & rk4_NS(nx, nxi, Q2(iQ2-1), Q2(iQ2), K_NS_0(:,:,ixi,nfl), KV_NS_1m(:,:,ixi,nfl), l_nlo)
            ! Then matrix multiplication to turn into matrix from initial Q2 to current Q2
            M_NS_pls(:,:,ixi,iQ2) = matmul(M_NS_pls(:,:,ixi,iQ2), M_NS_pls(:,:,ixi,iQ2-1))
            M_NS_min(:,:,ixi,iQ2) = matmul(M_NS_min(:,:,ixi,iQ2), M_NS_min(:,:,ixi,iQ2-1))
          end do
          !$OMP END PARALLEL DO
        end do
        print *, "Flag A5"
    end subroutine make_evomat_NS

    subroutine make_evomat_V_SG(nx, nxi, nQ2, Q2, l_nlo)
        integer,  intent(in) :: nx, nxi, nQ2
        real(dp), intent(in) :: Q2(nQ2)
        logical,  intent(in) :: l_nlo
        real(dp), dimension(2*nx,2*nx) :: id2nx
        integer :: ix, ixi, iQ2, nfl
        if(allocated(MV_SG)) deallocate(MV_SG)
        allocate(MV_SG(2*nx,2*nx,nxi,nQ2))
        ! Build an identity matrix
        id2nx = 0.0_dp
        do ix=1, 2*nx, 1
          id2nx(ix,ix) = 1.0_dp
        end do
        ! Identity for no evolution at initial scale
        MV_SG = 0.0_dp
        do ixi=1, nxi, 1
          MV_SG(:,:,ixi,1) = id2nx
        end do
        ! Build evolution matrices for all other scales
        do iQ2=2, nQ2, 1
          !$OMP PARALLEL DO
          do ixi=1, nxi, 1
            ! First, an evolution matrix for prior Q2 step to current Q2 step
            nfl = get_neff(Q2(iQ2-1))
            MV_SG(:,:,ixi,iQ2) = id2nx + &
                & rk4_SG(nx, nxi, Q2(iQ2-1), Q2(iQ2), &
                & KV_SG_0(:,:,ixi,nfl), KV_SG_1(:,:,ixi,nfl), l_nlo)
            ! Then matrix multiplication to turn into matrix from initial Q2 to current Q2
            MV_SG(:,:,ixi,iQ2) = matmul(MV_SG(:,:,ixi,iQ2), MV_SG(:,:,ixi,iQ2-1))
          end do
          !$OMP END PARALLEL DO
        end do
    end subroutine make_evomat_V_SG

    subroutine make_evomat_A_SG(nx, nxi, nQ2, Q2, l_nlo)
        integer,  intent(in) :: nx, nxi, nQ2
        real(dp), intent(in) :: Q2(nQ2)
        logical,  intent(in) :: l_nlo
        real(dp), dimension(2*nx,2*nx) :: id2nx
        integer :: ix, ixi, iQ2, nfl
        if(allocated(MA_SG)) deallocate(MA_SG)
        allocate(MA_SG(2*nx,2*nx,nxi,nQ2))
        ! Build an identity matrix
        id2nx = 0.0_dp
        do ix=1, 2*nx, 1
          id2nx(ix,ix) = 1.0_dp
        end do
        ! Identity for no evolution at initial scale
        MA_SG = 0.0_dp
        do ixi=1, nxi, 1
          MA_SG(:,:,ixi,1) = id2nx
        end do
        ! Build evolution matrices for all other scales
        do iQ2=2, nQ2, 1
          !$OMP PARALLEL DO
          do ixi=1, nxi, 1
            ! First, an evolution matrix for prior Q2 step to current Q2 step
            nfl = get_neff(Q2(iQ2-1))
            MA_SG(:,:,ixi,iQ2) = id2nx + &
                & rk4_SG(nx, nxi, Q2(iQ2-1), Q2(iQ2), &
                & KA_SG_0(:,:,ixi,nfl), KA_SG_1(:,:,ixi,nfl), l_nlo)
            ! Then matrix multiplication to turn into matrix from initial Q2 to current Q2
            MA_SG(:,:,ixi,iQ2) = matmul(MA_SG(:,:,ixi,iQ2), MA_SG(:,:,ixi,iQ2-1))
          end do
          !$OMP END PARALLEL DO
        end do
    end subroutine make_evomat_A_SG

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! RK4

    function rk4_NS(nx, nxi, Q2i, Q2f, K0, K1, l_nlo) result(K)
    !subroutine rk4_NS(nx, nxi, Q2i, Q2f, K0, K1, K)
        integer,  intent(in) :: nx, nxi
        real(dp), intent(in) :: Q2i, Q2f, K0(nx,nx), K1(nx,nx)
        logical,  intent(in) :: l_nlo
        !real(dp), intent(out) :: K(nx,nx)
        real(dp) :: K(nx,nx)
        !
        !real(dp), dimension(nx,nx) :: W1, W2, W3, W4, M_ini, M_mid, M_end
        real(dp), dimension(:,:), allocatable :: W1, W2, W3, W4, M_ini, M_mid, M_end
        real(dp) :: h, Q2m, a_ini, a_mid, a_end
        allocate(W1(nx,nx), W2(nx,nx), W3(nx,nx), W4(nx,nx), M_ini(nx,nx), M_mid(nx,nx), M_end(nx,nx))
        h = log(Q2f/Q2i)
        Q2m = sqrt(Q2i*Q2f) ! geometric mean for midpoint
        !
        a_ini = get_alpha_QCD(Q2i) / (2.*pi)
        a_mid = get_alpha_QCD(Q2m) / (2.*pi)
        a_end = get_alpha_QCD(Q2f) / (2.*pi)
        ! The three matrices over each subinterval
        M_ini = a_ini*K0
        M_mid = a_mid*K0
        M_end = a_end*K0
        ! Possible NLO corrections
        if(l_nlo) then
          M_ini = M_ini + a_ini**2*K1
          M_mid = M_mid + a_mid**2*K1
          M_end = M_end + a_end**2*K1
        endif
        ! The RK4 k-values
        W1 = M_ini
        W2 = M_mid + 0.5*h*matmul(M_mid,W1)
        W3 = M_mid + 0.5*h*matmul(M_mid,W2)
        W4 = M_end + h*matmul(M_end,W3)
        ! And the final formula for the shift
        K = h/6.*(W1 + 2.*W2 + 2.*W3 + W4)
        deallocate(W1, W2, W3, W4, M_ini, M_mid, M_end)
    end function rk4_NS
    !end subroutine rk4_NS

    function rk4_SG(nx, nxi, Q2i, Q2f, K0, K1, l_nlo) result(K)
    !subroutine rk4_SG(nx, nxi, Q2i, Q2f, K0, K1, K)
        integer,  intent(in) :: nx, nxi
        real(dp), intent(in) :: Q2i, Q2f, K0(2*nx,2*nx), K1(2*nx,2*nx)
        logical,  intent(in) :: l_nlo
        !real(dp), intent(out) :: :: K(2*nx,2*nx)
        real(dp) :: K(2*nx,2*nx)
        !
        !real(dp), dimension(2*nx,2*nx) :: W1, W2, W3, W4, M_ini, M_mid, M_end
        real(dp), dimension(:,:), allocatable :: W1, W2, W3, W4, M_ini, M_mid, M_end
        real(dp) :: h, Q2m, a_ini, a_mid, a_end
        allocate(W1(2*nx,2*nx), W2(2*nx,2*nx), W3(2*nx,2*nx), W4(2*nx,2*nx), M_ini(2*nx,2*nx), M_mid(2*nx,2*nx), M_end(2*nx,2*nx))
        h = log(Q2f/Q2i)
        Q2m = sqrt(Q2i*Q2f) ! geometric mean for midpoint
        !
        a_ini = get_alpha_QCD(Q2i) / (2.*pi)
        a_mid = get_alpha_QCD(Q2m) / (2.*pi)
        a_end = get_alpha_QCD(Q2f) / (2.*pi)
        ! The three matrices over each subinterval
        M_ini = a_ini*K0
        M_mid = a_mid*K0
        M_end = a_end*K0
        ! Possible NLO corrections
        if(l_nlo) then
          M_ini = M_ini + a_ini**2*K1
          M_mid = M_mid + a_mid**2*K1
          M_end = M_end + a_end**2*K1
        endif
        ! The RK4 k-values
        W1 = M_ini
        W2 = M_mid + 0.5*h*matmul(M_mid,W1)
        W3 = M_mid + 0.5*h*matmul(M_mid,W2)
        W4 = M_end + h*matmul(M_end,W3)
        ! And the final formula for the shift
        K = h/6.*(W1 + 2.*W2 + 2.*W3 + W4)
        deallocate(W1, W2, W3, W4, M_ini, M_mid, M_end)
    end function rk4_SG
    !end subroutine rk4_SG

end module matrices_evolution

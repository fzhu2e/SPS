!=================================================
! The model module of SPS-dynamic.
!-------------------------------------------------
! Version: 0.11
! Author: Zhu F.
! Email: lyricorpse@gmail.com
! Date: 2013-05-04 13:59:46 
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!-------------------------------------------------
! This module contains some settings about this model.
!=================================================
MODULE sp_module_model
USE sp_module_constant
IMPLICIT NONE
!=================================================
! Model settings.
!-------------------------------------------------
INTEGER, PARAMETER :: TimeScheme = 2       ! 1. Forward-backward; 2. Runge-Kutta;
INTEGER, PARAMETER :: AdvectionScheme = 5  ! 2. 2-order; 3. 3-order; 4. 4-order; 5. 5-order; 6. 6-order;
INTEGER, PARAMETER :: LateralBoundary = 1  ! 1. No-flux; 2. Periodic; 3. Open;
INTEGER, PARAMETER :: UpperBoundary = 1    ! 1. No-flux;
INTEGER, PARAMETER :: VertCoords = 1       ! 1. Height;

!=================================================
!-------------------------------------------------
! 1. Density current.
!-------------------------------------------------
INTEGER, PARAMETER :: RunCase = 1         ! 1. Density Current; 2. Thermal Bubble; 3. Internal gravity waves

INTEGER, PARAMETER :: nx = 512                              ! grid number along x-axis
INTEGER, PARAMETER :: nz = 64                               ! grid number along z-axis
REAL(kd), PARAMETER :: dx = 100.                        ! delta x (m)
REAL(kd), PARAMETER :: dz = 100.                       ! (m)

!INTEGER, PARAMETER :: nx = 1024                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 128                               ! grid number along z-axis
!REAL(kd), PARAMETER :: dx = 50.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 50.                       ! (m)

!INTEGER, PARAMETER :: nx = 2048                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 256                               ! grid number along z-axis
!REAL(kd), PARAMETER :: dx = 25.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 25.                       ! (m)

REAL(kd), PARAMETER :: dt = 0.1                               ! delta t (s)
INTEGER :: nstep = 9000

REAL(kd) :: Km = 75.        ! (m s)
REAL(kd) :: Kh = 75.        ! (K s)

!-------------------------------------------------
! 2. Thermal bubble.
!-------------------------------------------------
!INTEGER, PARAMETER :: RunCase = 2         ! 1. Density Current; 2. Thermal Bubble; 3. Internal gravity waves
!!INTEGER, PARAMETER :: nx = 200                              ! grid number along x-axis
!!INTEGER, PARAMETER :: nz = 100                               ! grid number along z-axis
!!REAL(kd), PARAMETER :: dx = 100.                        ! delta x (m)
!!REAL(kd), PARAMETER :: dz = 100.                       ! (m)

!INTEGER, PARAMETER :: nx = 400                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 200                               ! grid number along z-axis
!REAL(kd), PARAMETER :: dx = 50.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 50.                       ! (m)

!REAL(kd), PARAMETER :: dt = 0.1                               ! delta t (s)
!INTEGER :: nstep = 10000

!REAL(kd) :: Km = 20.        !(m s)
!REAL(kd) :: Kh = 20.        !(K s)

!-------------------------------------------------
! 3. Inertia gravity waves.
!-------------------------------------------------
!INTEGER, PARAMETER :: RunCase = 3         ! 1. Density Current; 2. Thermal Bubble; 3. Internal gravity waves

!!INTEGER, PARAMETER :: nx = 1200                              ! grid number along x-axis
!!INTEGER, PARAMETER :: nz = 400                               ! grid number along z-axis
!!REAL(kd), PARAMETER :: dx = 250.                        ! delta x (m)
!!REAL(kd), PARAMETER :: dz = 25.                       ! (m)

!!INTEGER, PARAMETER :: nx = 600                              ! grid number along x-axis
!!INTEGER, PARAMETER :: nz = 200                               ! grid number along z-axis
!!REAL(kd), PARAMETER :: dx = 500.                        ! delta x (m)
!!REAL(kd), PARAMETER :: dz = 50.                       ! (m)

!INTEGER, PARAMETER :: nx = 300                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 100                               ! grid number along z-axis
!REAL(kd), PARAMETER :: dx = 1000.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 100.                       ! (m)

!!INTEGER, PARAMETER :: nx = 150                              ! grid number along x-axis
!!INTEGER, PARAMETER :: nz = 50                               ! grid number along z-axis
!!REAL(kd), PARAMETER :: dx = 2000.                        ! delta x (m)
!!REAL(kd), PARAMETER :: dz = 200.                       ! (m)

!REAL(kd), PARAMETER :: dt = 0.1                               ! delta t (s)
!INTEGER :: nstep = 30000
!!INTEGER :: nstep = 1000

!REAL :: Km, Kh
!-------------------------------------------------
! 4. Schar mountain
!-------------------------------------------------
!INTEGER, PARAMETER :: RunCase = 4
!INTEGER, PARAMETER :: nx = 200                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 100                               ! grid number along z-axis
!!INTEGER, PARAMETER :: nx = 400                              ! grid number along x-axis
!!INTEGER, PARAMETER :: nz = 200                               ! grid number along z-axis
!REAL(kd), PARAMETER :: dx = 250.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 210.                       ! (m)

!!REAL(kd), PARAMETER :: pi_top = (p_top/p_0)**(R_d/C_pd)

!REAL(kd), PARAMETER :: dt = 0.12                               ! delta t (s)
!!INTEGER :: nstep = 300000
!INTEGER :: nstep = 10000
!!INTEGER :: nstep = 1000
!!INTEGER :: nstep = 300

!REAL :: Km, Kh
!-------------------------------------------------
INTEGER, PARAMETER :: its = 1
INTEGER, PARAMETER :: ite = its + nx - 1
INTEGER, PARAMETER :: kts = 1
INTEGER, PARAMETER :: kte = kts + nz - 1

INTEGER, PARAMETER :: halo = 4

INTEGER, PARAMETER :: ims = its - halo
INTEGER, PARAMETER :: ime = ite + halo
INTEGER, PARAMETER :: kms = kts - halo
INTEGER, PARAMETER :: kme = kte + halo

!-------------------------------------------------
! Vertical coordinates

!REAL(kd), PARAMETER :: ztop = nz*dz                         ! (m)
!REAL(kd), DIMENSION(ims:ime) :: zs = 0.                     ! (m) on u-grid
!REAL(kd), DIMENSION(ims:ime) :: zs_pi = 0.                  ! (m) on pi-grid
!REAL(kd), DIMENSION(ims:ime) :: PzsPx = 0.                  ! (m) on u-grid
!REAL(kd), DIMENSION(ims:ime) :: PzsPx_pi = 0.               ! (m) on pi-grid
!REAL(kd), DIMENSION(kms:kme) :: b                           !     on w-grid
!REAL(kd), DIMENSION(kms:kme) :: b_pi                        !     on pi-grid
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: z_hat               ! (m) on w-grid : Real Height in z_hat coords.
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: z_hat_pi            ! (m) on pi-grid
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: z_hat_u             ! (m) on u-grid
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: z_hat_vir             ! (m) on v-grid
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: PbPzhat             ! (m) on w-grid
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: PbPzhat_pi          ! (m) on pi-grid
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: PbPzhat_u           ! (m) on u-grid
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: PbPzhat_vir           ! (m) on v-grid
!REAL(kd), DIMENSION(ims:ime) :: VertA_u = 0.                !  Vertical Variation A on u-grid
!REAL(kd), DIMENSION(ims:ime) :: VertA_pi = 0.               !  Vertical Variation A on pi-grid
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: VertB_u = 1.        !  Vertical Variation B on u-grid
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: VertB_w = 1.        !  Vertical Variation B on w-grid
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: VertB_pi = 1.       !  Vertical Variation B on pi-grid
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: VertB_vir = 1.        !  Vertical Variation B on v-grid
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: VertC_u = 0.        !  Vertical Variation C on u-grid
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: VertC_w = 0.        !  Vertical Variation C on w-grid
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: VertC_pi = 0.       !  Vertical Variation C on pi-gridi
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: VertC_vir = 0.        !  Vertical Variation C on v-gridi
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime) :: xx      ! distance on u-grid along x-axis (m)
REAL(kd), DIMENSION(ims:ime) :: xpi     ! distance on pi-grid along x-axis (m)
REAL(kd), DIMENSION(kms:kme) :: zz      ! height on w-grid along z-axis (m) : Height of Model Level
REAL(kd), DIMENSION(kms:kme) :: zpi     ! height on pi-grid along z-axis (m)
!-------------------------------------------------
INTEGER :: imin, imax, kmin, kmax
!=================================================

CONTAINS
!=================================================
! Set calculate area of u-grid, w-grid, pi-grid, and v-grid.
! Fix the index of boundaries.
!=================================================
SUBROUTINE set_calc_area_u
IMPLICIT NONE
! u-grid (its + 1:ite - 1, kts:kte)
imin = its + 1
imax = ite - 1
kmin = kts
kmax = kte
END SUBROUTINE set_calc_area_u

SUBROUTINE set_calc_area_w
IMPLICIT NONE
! w-grid (it + 1:ite, kts + 1:kte)
imin = its + 1
imax = ite
kmin = kts + 1
kmax = kte
END SUBROUTINE set_calc_area_w

SUBROUTINE set_calc_area_pi
IMPLICIT NONE
! pi-grid (its + 1:ite, kts:kte)
imin = its + 1
imax = ite
kmin = kts
kmax = kte
END SUBROUTINE set_calc_area_pi

SUBROUTINE set_calc_area_vir
IMPLICIT NONE
! v-grid (its + 1:ite - 1, kts + 1:kte)
imin = its + 1
imax = ite - 1
kmin = kts + 1
kmax = kte
END SUBROUTINE set_calc_area_vir
!=================================================

!=================================================
! Arakawa-C grid.
! Set calculate area of u-grid, w-grid, pi-grid, and v-grid.
!=================================================
SUBROUTINE set_area_u
IMPLICIT NONE
imin = its
imax = ite
kmin = kts
kmax = kte
END SUBROUTINE set_area_u

SUBROUTINE set_area_w
IMPLICIT NONE
imin = its + 1
imax = ite
kmin = kts
kmax = kte + 1
END SUBROUTINE set_area_w

SUBROUTINE set_area_pi
IMPLICIT NONE
imin = its + 1
imax = ite
kmin = kts
kmax = kte
END SUBROUTINE set_area_pi

SUBROUTINE set_area_vir
IMPLICIT NONE
imin = its
imax = ite
kmin = kts
kmax = kte + 1
END SUBROUTINE set_area_vir

SUBROUTINE set_area_expand(num)
IMPLICIT NONE
INTEGER :: num
imin = its - num
imax = ite + num
kmin = kts - num
kmax = kte + num
END SUBROUTINE set_area_expand
!=================================================


!=================================================
END MODULE sp_module_model
!=================================================

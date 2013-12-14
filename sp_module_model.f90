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
INTEGER, PARAMETER :: TimeScheme = 1       ! 1. Runge-Kutta;
INTEGER, PARAMETER :: AdvectionScheme = 5  ! 2. 2-order; 3. 3-order; 4. 4-order; 5. 5-order; 6. 6-order;
INTEGER, PARAMETER :: LateralBoundary = 2  ! 1. No-flux; 2. Periodic; 3. Open;
INTEGER, PARAMETER :: UpperBoundary = 1    ! 1. No-flux;
INTEGER, PARAMETER :: VertCoords = 1       ! 1. Height;

!=================================================
!-------------------------------------------------
! 1. Density current.
!-------------------------------------------------
!INTEGER, PARAMETER :: RunCase = 1         ! 1. Density Current; 2. Thermal Bubble; 3. Internal gravity waves

!INTEGER, PARAMETER :: nx = 512                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 64                               ! grid number along z-axis
!REAL(kd), PARAMETER :: dx = 100.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 100.                       ! (m)

!!INTEGER, PARAMETER :: nx = 1024                              ! grid number along x-axis
!!INTEGER, PARAMETER :: nz = 128                               ! grid number along z-axis
!!REAL(kd), PARAMETER :: dx = 50.                        ! delta x (m)
!!REAL(kd), PARAMETER :: dz = 50.                       ! (m)

!!INTEGER, PARAMETER :: nx = 2048                              ! grid number along x-axis
!!INTEGER, PARAMETER :: nz = 256                               ! grid number along z-axis
!!REAL(kd), PARAMETER :: dx = 25.                        ! delta x (m)
!!REAL(kd), PARAMETER :: dz = 25.                       ! (m)

!REAL(kd), PARAMETER :: dt = 0.1                               ! delta t (s)
!INTEGER :: nstep = 9000

!REAL(kd) :: Km = 75.        ! (m s)
!REAL(kd) :: Kh = 75.        ! (K s)

!-------------------------------------------------
! 2. Thermal bubble.
!-------------------------------------------------
!INTEGER, PARAMETER :: RunCase = 2         ! 1. Density Current; 2. Thermal Bubble; 3. Internal gravity waves
!INTEGER, PARAMETER :: nx = 200                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 100                               ! grid number along z-axis
!REAL(kd), PARAMETER :: dx = 100.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 100.                       ! (m)

!!INTEGER, PARAMETER :: nx = 400                              ! grid number along x-axis
!!INTEGER, PARAMETER :: nz = 200                               ! grid number along z-axis
!!REAL(kd), PARAMETER :: dx = 50.                        ! delta x (m)
!!REAL(kd), PARAMETER :: dz = 50.                       ! (m)

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

!REAL :: Km = 0., Kh = 0.
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

!REAL :: Km = 0., Kh = 0.

!-------------------------------------------------
! 5. Wet bubble.
!-------------------------------------------------
INTEGER, PARAMETER :: RunCase = 5
INTEGER, PARAMETER :: nx = 100                              ! grid number along x-axis
INTEGER, PARAMETER :: nz = 200                               ! grid number along z-axis
REAL(kd), PARAMETER :: dx = 100.                        ! delta x (m)
REAL(kd), PARAMETER :: dz = 100.                       ! (m)

REAL(kd), PARAMETER :: dt = 0.1                               ! delta t (s)
INTEGER :: nstep = 10000

REAL(kd) :: Km = 20.        !(m s)
REAL(kd) :: Kh = 20.        !(K s)

!=================================================
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
REAL(kd), PARAMETER :: ztop = nz*dz                         ! (m)
!-------------------------------------------------
INTEGER :: imin, imax, kmin, kmax
!=================================================

CONTAINS

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

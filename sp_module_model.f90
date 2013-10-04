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
INTEGER, PARAMETER :: RunCase = 1         ! 1. Density Current; 2. Thermal Bubble; 3. Internal gravity waves
INTEGER, PARAMETER :: TimeScheme = 3      ! 1. Forward-backward; 2. Leapfrog; 3. Runge-Kutta; 99. Debug
INTEGER, PARAMETER :: AdvectionScheme = 2 ! 2. 2-order; 3. 3-order; 4. 4-order; 5. 5-order; 6. 6-order
INTEGER, PARAMETER :: LateralBoundary = 2 ! 1. Wall; 2. Periodic
!-------------------------------------------------
! 1. Density current.
!-------------------------------------------------
!INTEGER, PARAMETER :: nx = 512                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 64                               ! grid number along z-axis
!REAL(preci), PARAMETER :: dx = 100.                        ! delta x (m)
!REAL(preci), PARAMETER :: dz = 100.                       ! (m)

INTEGER, PARAMETER :: nx = 1024                              ! grid number along x-axis
INTEGER, PARAMETER :: nz = 128                               ! grid number along z-axis
REAL(preci), PARAMETER :: dx = 50.                        ! delta x (m)
REAL(preci), PARAMETER :: dz = 50.                       ! (m)

REAL(preci), PARAMETER :: dt = 0.1                               ! delta t (s)

!-------------------------------------------------
! 2. Thermal bubble.
!-------------------------------------------------
!INTEGER, PARAMETER :: nx = 200                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 100                               ! grid number along z-axis
!REAL(preci), PARAMETER :: dx = 100.                        ! delta x (m)
!REAL(preci), PARAMETER :: dz = 100.                       ! (m)

!INTEGER, PARAMETER :: nx = 400                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 200                               ! grid number along z-axis
!REAL(preci), PARAMETER :: dx = 50.                        ! delta x (m)
!REAL(preci), PARAMETER :: dz = 50.                       ! (m)

!REAL(preci), PARAMETER :: dt = 0.1                               ! delta t (s)

!-------------------------------------------------
! 3. Inertia gravity waves.
!-------------------------------------------------
!INTEGER, PARAMETER :: nx = 600                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 200                               ! grid number along z-axis
!REAL(preci), PARAMETER :: dx = 500.                        ! delta x (m)
!REAL(preci), PARAMETER :: dz = 50.                       ! (m)

!REAL(preci), PARAMETER :: dt = 0.1                               ! delta t (s)
!-------------------------------------------------
!REAL(preci), PARAMETER :: dx = 1000.                        ! delta x (m)
!REAL(preci), PARAMETER :: ztop = 35000.                     ! (m)
!REAL(preci), PARAMETER :: dz = ztop/REAL((nz - 1))          ! (m)
!REAL(preci), PARAMETER :: dz = 1000.                       ! (m)
!REAL(preci), PARAMETER :: pi_top = (p_top/p_0)**(R_d/C_pd)
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
!=================================================
END MODULE sp_module_model
!=================================================

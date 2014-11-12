!=================================================
! Super-Parametertization System (SPS)
!-------------------------------------------------
! Version: 0.2
! Author: Feng Zhu
! Email: zhuf.atmos@gmail.com
! Date: 2014-06-12 18:18:45
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!=================================================
MODULE sp_module_model
USE sp_module_constant
IMPLICIT NONE
!=================================================
! Model settings.
!-------------------------------------------------
INTEGER, PARAMETER :: TimeScheme = 1       ! 1. Runge-Kutta;
INTEGER, PARAMETER :: HoriAdv = 5          ! 2. 2-order; 5. 5-order;
INTEGER, PARAMETER :: VertAdv = 5          ! 2. 2-order; 5. 5-order;
INTEGER, PARAMETER :: LateralBoundary = 1  ! 1. No-flux; 2. Periodic; 3. Open;
INTEGER, PARAMETER :: UpperBoundary = 1    ! 1. No-flux; 2. Open
INTEGER, PARAMETER :: DampTop = 0          ! 0. Off; 1. On
INTEGER, PARAMETER :: DampLateral = 0      ! 0. Off; 1. On
INTEGER, PARAMETER :: VertCoords = 1       ! 1. Height;
INTEGER, PARAMETER :: VertLev = 0          ! 0. Uniform; 1. Uneven; 2. GRAPES

!--------------------- Uniform vertical layers (5-order VertAdv)
!REAL(kd), PARAMETER :: s = 0.50e4 ! Thickness of Rayleigh Layer (m)
!REAL(kd), PARAMETER :: tau0 = 0.02
!--------------------- GRAPES vertical layers (2-order VertAdv)
REAL(kd), PARAMETER :: s = 1.50e4 ! Thickness of Rayleigh Layer (m)
REAL(kd), PARAMETER :: tau0 = 0.02
!--------------------- GRAPES vertical layers (2-order VertAdv)
!REAL(kd), PARAMETER :: s = 1.50e4 ! Thickness of Rayleigh Layer (m)
!REAL(kd), PARAMETER :: tau0 = 1.0
!=================================================
!-------------------------------------------------
! 1. Density current.
!-------------------------------------------------
INTEGER, PARAMETER :: RunCase = 1         ! 1. Density Current; 2. Thermal Bubble; 3. Internal gravity waves
INTEGER, PARAMETER :: Vapor = 0

!INTEGER, PARAMETER :: nx = 266                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 32                               ! grid number along z-axis
!REAL(kd), PARAMETER :: dx = 200.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 200.                       ! (m)
!REAL(kd), PARAMETER :: dt = 0.3                               ! delta t (s)
!INTEGER :: nstep = 3000

INTEGER, PARAMETER :: nx = 531                              ! grid number along x-axis
INTEGER, PARAMETER :: nz = 64                               ! grid number along z-axis
REAL(kd), PARAMETER :: dx = 100.                        ! delta x (m)
REAL(kd), PARAMETER :: dz = 100.                       ! (m)
REAL(kd), PARAMETER :: dt = 0.2                               ! delta t (s)
INTEGER :: nstep = 4500

!INTEGER, PARAMETER :: nx = 1061                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 128                               ! grid number along z-axis
!REAL(kd), PARAMETER :: dx = 50.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 50.                       ! (m)
!REAL(kd), PARAMETER :: dt = 0.1                               ! delta t (s)
!INTEGER :: nstep = 9000

REAL(kd) :: Km = 75.        ! (m s)
REAL(kd) :: Kh = 75.        ! (K s)
REAL(kd) :: ztop = nz*dz                         ! (m)

!-------------------------------------------------
! 2. Thermal bubble.
!-------------------------------------------------
!INTEGER, PARAMETER :: RunCase = 2         ! 1. Density Current; 2. Thermal Bubble; 3. Internal gravity waves
!INTEGER, PARAMETER :: Vapor = 0

!!INTEGER, PARAMETER :: nx = 101                              ! grid number along x-axis
!!INTEGER, PARAMETER :: nz = 50                               ! grid number along z-axis
!!REAL(kd), PARAMETER :: dx = 200.                        ! delta x (m)
!!REAL(kd), PARAMETER :: dz = 200.                       ! (m)
!!REAL(kd), PARAMETER :: dt = 0.4                               ! delta t (s)
!!INTEGER :: nstep = 2500

!INTEGER, PARAMETER :: nx = 201                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 100                               ! grid number along z-axis
!REAL(kd), PARAMETER :: dx = 100.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 100.                       ! (m)
!REAL(kd), PARAMETER :: dt = 0.2                               ! delta t (s)
!INTEGER :: nstep = 5000

!!INTEGER, PARAMETER :: nx = 401                              ! grid number along x-axis
!!INTEGER, PARAMETER :: nz = 200                               ! grid number along z-axis
!!REAL(kd), PARAMETER :: dx = 50.                        ! delta x (m)
!!REAL(kd), PARAMETER :: dz = 50.                       ! (m)
!!REAL(kd), PARAMETER :: dt = 0.1                               ! delta t (s)
!!INTEGER :: nstep = 10000

!REAL(kd) :: Km = 0.        !(m s)
!REAL(kd) :: Kh = 0.        !(K s)
!REAL(kd) :: ztop = nz*dz                         ! (m)

!-------------------------------------------------
! 3. Inertia gravity waves.
!-------------------------------------------------
!INTEGER, PARAMETER :: RunCase = 3
!INTEGER, PARAMETER :: Vapor = 0

!!INTEGER, PARAMETER :: nx = 601                              ! grid number along x-axis
!!INTEGER, PARAMETER :: nz = 200                               ! grid number along z-axis
!!REAL(kd), PARAMETER :: dx = 500.                        ! delta x (m)
!!REAL(kd), PARAMETER :: dz = 50.                       ! (m)
!!REAL(kd), PARAMETER :: dt = 0.1                               ! delta t (s)
!!INTEGER :: nstep = 30000

!INTEGER, PARAMETER :: nx = 301                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 100                               ! grid number along z-axis
!REAL(kd), PARAMETER :: dx = 1000.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 100.                       ! (m)
!REAL(kd), PARAMETER :: dt = 0.2                               ! delta t (s)
!INTEGER :: nstep = 15000

!!INTEGER, PARAMETER :: nx = 151                              ! grid number along x-axis
!!INTEGER, PARAMETER :: nz = 50                               ! grid number along z-axis
!!REAL(kd), PARAMETER :: dx = 2000.                        ! delta x (m)
!!REAL(kd), PARAMETER :: dz = 200.                       ! (m)
!!REAL(kd), PARAMETER :: dt = 0.4                               ! delta t (s)
!!INTEGER :: nstep = 7500

!REAL :: Km = 0., Kh = 0.
!REAL(kd) :: ztop = nz*dz                         ! (m)
!-------------------------------------------------
! 4. Schar mountain
!-------------------------------------------------
!INTEGER, PARAMETER :: RunCase = 4
!INTEGER, PARAMETER :: Vapor = 0

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
!REAL(kd) :: ztop = nz*dz                         ! (m)

!-------------------------------------------------
! 5. Wet bubble
!-------------------------------------------------
!INTEGER, PARAMETER :: RunCase = 5
!INTEGER, PARAMETER :: Vapor = 1

!INTEGER, PARAMETER :: nx = 100                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 200                               ! grid number along z-axis
!REAL(kd), PARAMETER :: dx = 5.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 5.                       ! (m)

!REAL(kd), PARAMETER :: dt = 0.01                               ! delta t (s)
!INTEGER :: nstep = 10000

!!REAL(kd) :: Km = 20.        !(m s)
!!REAL(kd) :: Kh = 20.        !(K s)
!REAL(kd) :: Km = 0.        !(m s)
!REAL(kd) :: Kh = 0.        !(K s)
!REAL(kd) :: ztop = nz*dz                         ! (m)

!-------------------------------------------------
! 6a. Thunderstorm - Standard
!-------------------------------------------------
!INTEGER, PARAMETER :: RunCase = 6
!INTEGER, PARAMETER :: Vapor = 1

!!REAL(kd), PARAMETER :: s = 0.50e4 ! Thickness of Rayleigh Layer (m)
!!REAL(kd), PARAMETER :: tau0 = 0.02

!INTEGER, PARAMETER :: nx = 201                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 80                               ! grid number along z-axis
!REAL(kd), PARAMETER :: dx = 250.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 250.                       ! (m)
!REAL(kd), PARAMETER :: dt = 0.5                               ! delta t (s)
!INTEGER :: nstep = 3600

!REAL(kd) :: Km = 0.        !(m s)
!REAL(kd) :: Kh = 0.        !(K s)
!REAL(kd) :: ztop = nz*dz                         ! (m)

!-------------------------------------------------
! 6b. Thunderstorm - Plus
!-------------------------------------------------
!INTEGER, PARAMETER :: RunCase = 6
!INTEGER, PARAMETER :: Vapor = 1

!!REAL(kd), PARAMETER :: s = 1.50e4 ! Thickness of Rayleigh Layer (m)
!!REAL(kd), PARAMETER :: tau0 = 0.10

!INTEGER, PARAMETER :: nx = 201                              ! grid number along x-axis
!INTEGER, PARAMETER :: nz = 31                               ! grid number along z-axis
!REAL(kd), PARAMETER :: dx = 250.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 250.                       ! (m)
!REAL(kd), PARAMETER :: dt = 0.1                               ! delta t (s)
!INTEGER :: nstep = 18000
!REAL(kd) :: Km = 0.        !(m s)
!REAL(kd) :: Kh = 0.        !(K s)

!REAL(kd) :: ztop = nz*dz                         ! (m)

!-------------------------------------------------
! 99. Real Case
!-------------------------------------------------
!INTEGER, PARAMETER :: RunCase = 99
!INTEGER, PARAMETER :: Vapor = 1

!!REAL(kd), PARAMETER :: s = 0.50e4 ! Thickness of Rayleigh Layer (m)
!!REAL(kd), PARAMETER :: tau0 = 0.02

!INTEGER, PARAMETER :: nx = 31                              ! grid number along x-axis
!!INTEGER, PARAMETER :: nz = 80                               ! grid number along z-axis
!!REAL(kd), PARAMETER :: dx = 1000.                        ! delta x (m)
!REAL(kd), PARAMETER :: dx = 500.                        ! delta x (m)
!REAL(kd), PARAMETER :: dz = 250.                       ! (m)
!!REAL(kd), PARAMETER :: dt = 0.5                               ! delta t (s)
!!INTEGER :: nstep = 3600
!REAL(kd), PARAMETER :: dt = 0.1                               ! delta t (s)
!!INTEGER :: nstep = 900
!INTEGER :: nstep = 900*960


!REAL(kd) :: Km = 0.        !(m s)
!REAL(kd) :: Kh = 0.        !(K s)
!!REAL(kd) :: ztop = nz*dz                         ! (m)
!!INTEGER, PARAMETER :: nz = 30                               ! grid number along z-axis
!REAL(kd) :: ztop = 2.e4
!INTEGER, PARAMETER :: nz = 31                               ! grid number along z-axis

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
REAL(kd), DIMENSION(kms:kme) :: dk = undef
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

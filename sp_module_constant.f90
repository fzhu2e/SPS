!=================================================
! The constant module of SPS-dynamic.
!-------------------------------------------------
! Version: 0.11
! Author: Zhu F.
! Email: lyricorpse@gmail.com
! Date: 2013-05-04 13:59:46 
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!-------------------------------------------------
! This module contains some numerical and physical constants.
!=================================================
MODULE sp_module_constant
IMPLICIT NONE
!=================================================
! Numerical settings.
!-------------------------------------------------
INTEGER, PARAMETER :: preci = 8
REAL(preci), PARAMETER :: undef = -9999.
!=================================================
! Physical settings.
!-------------------------------------------------
!REAL(preci), PARAMETER :: PI_math = 3.1415926
REAL(preci), PARAMETER :: PI_math = 4.*ATAN(1.)
REAL(preci), PARAMETER :: g = 9.8            ! gravity (m/s^2)
REAL(preci), PARAMETER :: omega = 7.292E-5
REAL(preci), PARAMETER :: latitude = 45*PI_math/180.
REAL(preci), PARAMETER :: f = 2*omega*SIN(latitude)
REAL(preci), PARAMETER :: cs = 340.          ! speed of sound (m/s)
REAL(preci), PARAMETER :: Ts = 300.         ! (K)
REAL(preci), PARAMETER :: Cp = 1005.4152   ! 0.2403 (cal/g/K) = 0.2403*4.184*1000 (J/kg/K) = 1005.4152 (J/kg/K)
REAL(preci), PARAMETER :: p0 = 100000.      ! (Pa)
REAL(preci), PARAMETER :: p_top = 100.       ! (Pa)
REAL(preci), PARAMETER :: Rd = 287.04       ! 2.8704E+6 (erg/g/K) = 287.04 (J/kg/K)

! Density current.
!REAL(preci), PARAMETER :: Km = 75.        !(m^2/s)
!REAL(preci), PARAMETER :: Kh = 75.        !(m^2/s)

! Thermal bubble.
REAL(preci), PARAMETER :: Km = 20.        !(m^2/s)
REAL(preci), PARAMETER :: Kh = 20.        !(m^2/s)

! Inertia gravity waves.
!REAL(preci), PARAMETER :: Km = 0.        !(m^2/s)
!REAL(preci), PARAMETER :: Kh = 0.        !(m^2/s)
!=================================================
END MODULE sp_module_constant
!=================================================

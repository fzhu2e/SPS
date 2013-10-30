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
!REAL(preci), PARAMETER :: cs = 340.          ! speed of sound (m/s)
REAL(preci), PARAMETER :: cs = 300.          ! speed of sound (m/s)
REAL(preci), PARAMETER :: Cp = 1005.4152    ! 0.2403 (cal/g/K) = 0.2403*4.184*1000 (J/kg/K) = 1005.4152 (J/kg/K)
REAL(preci), PARAMETER :: Cv = 717.556    ! 0.1715 (cal/g/K) = 0.1715*4.184*1000 (J/kg/K) = 717.556 (J/kg/K)
REAL(preci), PARAMETER :: p0 = 100000.      ! (Pa)
REAL(preci), PARAMETER :: p_top = 100.       ! (Pa)
REAL(preci), PARAMETER :: Rd = 287.04       ! 2.8704E+6 (erg/g/K) = 287.04 (J/kg/K)
!=================================================
END MODULE sp_module_constant
!=================================================

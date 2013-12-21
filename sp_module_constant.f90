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
INTEGER, PARAMETER :: kd = KIND(1d0)
!REAL(kd), PARAMETER :: undef = -9999.
REAL(kd), PARAMETER :: undef = -999999999.
INTEGER, PARAMETER :: expand = 1
!=================================================
! Physical settings.
!-------------------------------------------------
REAL(kd), PARAMETER :: PI_math = 4.*ATAN(1.)  ! = 3.1415926
REAL(kd), PARAMETER :: g = 9.8                ! gravity (m/s^2)
REAL(kd), PARAMETER :: omega = 7.292E-5
REAL(kd), PARAMETER :: latitude = 45*PI_math/180.
REAL(kd), PARAMETER :: f = 2*omega*SIN(latitude)
REAL(kd), PARAMETER :: cs = 340.        ! speed of sound (m/s)
REAL(kd), PARAMETER :: Cp = 1005.4152    ! 0.2403 (cal/g/K) = 0.2403*4.184*1000 (J/kg/K) = 1005.4152 (J/kg/K)
REAL(kd), PARAMETER :: Cv = 717.556      ! 0.1715 (cal/g/K) = 0.1715*4.184*1000 (J/kg/K) = 717.556 (J/kg/K)
REAL(kd), PARAMETER :: p0 = 100000.      ! (Pa)
REAL(kd), PARAMETER :: p_top = 100.       ! (Pa)
REAL(kd), PARAMETER :: Rd = 287.04       ! 2.8704E+6 (erg/g/K) = 287.04 (J/kg/K)
!=================================================
! for WSM6
!-------------------------------------------------
REAL(kd), PARAMETER :: rhoair0 = 1.28
REAL(kd), PARAMETER :: rhosnow = 100.
REAL(kd), PARAMETER :: rhowater = 1000.

REAL(kd), PARAMETER :: cliq = 4190.
REAL(kd), PARAMETER :: rv = 461.6
REAL(kd), PARAMETER :: cpv = 4.*rv

REAL(kd), PARAMETER :: cpd = 1004.6

REAL(kd), PARAMETER :: qmin = 1.e-15

REAL(kd), PARAMETER :: xlv0 = 3.15e6
REAL(kd), PARAMETER :: t0c = 273.15

REAL(kd), PARAMETER :: den0 = rhoair0
REAL(kd), PARAMETER :: denr = rhowater

REAL(kd), PARAMETER :: xls = 2.834e6
REAL(kd), PARAMETER :: cice = 2.1060e3
REAL(kd), PARAMETER :: psat = 6.1078e2
REAL(kd), PARAMETER :: ep2 = rd/rv

REAL(kd), PARAMETER :: xlv = 2.5e6
REAL(kd), PARAMETER :: xlf0 = xls - xlv
!=================================================
END MODULE sp_module_constant
!=================================================

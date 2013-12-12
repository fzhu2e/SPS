!=================================================
! The interpolate module of SPS-dynamic-integrate
!-------------------------------------------------
! Version: 0.11
! Author: Zhu F.
! Email: lyricorpse@gmail.com
! Date: 2013-05-04 14:01:18 
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!=================================================
MODULE sp_module_interpolate
USE sp_module_constant
USE sp_module_model
USE sp_module_gridvar
USE sp_module_debug
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================
!/////////////////////////////////////////////////////////////////////
!=================================================
! Basic interpolate.
!=================================================
SUBROUTINE basic_interpolate(Main,uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
!=================================================
TYPE(mainvar), INTENT(IN) :: Main
TYPE(grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL u2w(Main%u,wGrid%u)
CALL u2pi(Main%u,piGrid%u)
CALL u2vir(Main%u,virGrid%u)

CALL w2pi(Main%w,piGrid%w)
CALL w2vir(Main%w,virGrid%w)

CALL pi2vir(Main%pi_1,virGrid%pi_1)

CALL w2u(wGrid%theta_M,uGrid%theta_M)

!CALL w2pi(Main%theta,piGrid%theta)

!CALL w2pi(Main%qv,piGrid%qv)
!CALL w2pi(Main%qc,piGrid%qc)
!CALL w2pi(Main%qr,piGrid%qr)

!CALL w2u(Main%qv,uGrid%qv)
!CALL w2u(Main%qc,uGrid%qc)
!CALL w2u(Main%qr,uGrid%qr)

!CALL w2vir(Main%qv,virGrid%qv)
!CALL w2vir(Main%qc,virGrid%qc)
!CALL w2vir(Main%qr,virGrid%qr)
!=================================================
END SUBROUTINE basic_interpolate
!=================================================
!/////////////////////////////////////////////////////////////////////
!=================================================
END MODULE sp_module_interpolate
!=================================================

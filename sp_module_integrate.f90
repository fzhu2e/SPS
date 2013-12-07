!=================================================
! The integrate module of SPS-dynamic.
!-------------------------------------------------
! Version: 0.11
! Author: Zhu F.
! Email: lyricorpse@gmail.com
! Date: 2013-05-04 13:59:46 
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!=================================================
MODULE sp_module_integrate
USE sp_module_constant
USE sp_module_model
USE sp_module_gridvar
USE sp_module_timeschemes
USE sp_module_boundary
USE sp_module_debug
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================
SUBROUTINE integrate(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE(grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
TYPE(mainvar) :: new
!=================================================
SELECT CASE (TimeScheme)
CASE (1)
	! Runge-Kutta Scheme
	CALL runge_kutta(uGrid,wGrid,piGrid,virGrid,new)
CASE DEFAULT
	STOP "Wrong time differencing scheme!!!"
END SELECT

! update u, w, pi_1, theta, theta_1
uGrid%u = new%u
wGrid%w = new%w
piGrid%pi_1 = new%pi_1
wGrid%theta = new%theta
!=================================================
END SUBROUTINE integrate
!=================================================


!=================================================
END MODULE sp_module_integrate
!=================================================

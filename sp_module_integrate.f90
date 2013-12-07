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
USE sp_module_timeschemes
USE sp_module_boundary
USE sp_module_debug
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================
SUBROUTINE integrate(u,w,pi_1,theta,theta_1)
IMPLICIT NONE
!=================================================
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: u        ! wind speed along x-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: w        ! wind speed along z-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_1     ! pi'
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_1  ! theta'
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: ud_u                    ! updated u
REAL(kd), DIMENSION(ims:ime,kms:kme) :: ud_w                    ! updated w
REAL(kd), DIMENSION(ims:ime,kms:kme) :: ud_pi_1                 ! updated ud_pi_1
REAL(kd), DIMENSION(ims:ime,kms:kme) :: ud_theta                ! updated ud_theta
REAL(kd), DIMENSION(ims:ime,kms:kme) :: ud_theta_1              ! updated ud_theta_1
!=================================================
SELECT CASE (TimeScheme)
CASE (2)
	! Runge-Kutta Scheme
	CALL runge_kutta( u,    w,    pi_1,    theta,    theta_1,  &
	                  ud_u, ud_w, ud_pi_1, ud_theta, ud_theta_1   )
CASE DEFAULT
	STOP "Wrong time differencing scheme!!!"
END SELECT

! update u, w, pi_1, theta, theta_1
u = ud_u
w = ud_w
pi_1 = ud_pi_1
theta = ud_theta
theta_1 = ud_theta_1
!=================================================
END SUBROUTINE integrate
!=================================================


!=================================================
END MODULE sp_module_integrate
!=================================================

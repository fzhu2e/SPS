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
USE sp_module_debug
IMPLICIT NONE
!=================================================
! Basic interpolated variations in Arakawa-C and Charney-Phillips grids.
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: u_pi, u_w, u_v
REAL(preci), DIMENSION(ims:ime,kms:kme) :: w_pi, w_u, w_v
REAL(preci), DIMENSION(ims:ime,kms:kme) :: pi_1_u, pi_1_w, pi_1_v
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_pi, theta_u, theta_v
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_0_pi, theta_0_u, theta_0_v
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_1_pi, theta_1_u, theta_1_v
REAL(preci), DIMENSION(ims:ime,kms:kme) :: rho_0_u, rho_0_w, rho_0_v
!=================================================
CONTAINS
!=================================================
!/////////////////////////////////////////////////////////////////////
!=================================================
! Basic interpolate.
!=================================================
SUBROUTINE basic_interpolate(u,w,pi_1,theta,theta_0,theta_1,rho_0)
IMPLICIT NONE
!=================================================
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: u
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: w
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_1
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_0
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_1
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0
!-------------------------------------------------
INTEGER :: i, k 
!=================================================
CALL debug_undef_all( u_pi,u_w,u_v,                               &
                      w_pi,w_u,w_v,                               &
                      pi_1_u,pi_1_w,pi_1_v,                       &
                      theta_pi,theta_u,theta_v,                   &
                      theta_0_pi,theta_0_u,theta_0_v,             &
                      theta_1_pi,theta_1_u,theta_1_v,             &
                      rho_0_u,rho_0_w,rho_0_v                     )
!CALL debug_ascii_output(u_pi)
!CALL debug_test_boundary(u_pi)
!CALL debug_ascii_output(pi_1)
!CALL debug_ascii_output(pi_1_u)
!=================================================
!CALL debug_test_boundary(w)
! To u-grid (pi, rho, w, theta) - Mid-vars can be calculated on boundaries. 
DO i = its, ite
	DO k = kts - 1, kte + 1
		pi_1_u(i,k) = (pi_1(i,k) + pi_1(i+1,k))/2.
		rho_0_u(i,k) = (rho_0(i,k) + rho_0(i+1,k))/2.
		w_u(i,k) = (w(i,k) + w(i+1,k) + w(i,k+1) + w(i+1,k+1))/4.
		theta_u(i,k) = (theta(i,k) + theta(i+1,k) + theta(i,k+1) + theta(i+1,k+1))/4.
		theta_0_u(i,k) = (theta_0(i,k) + theta_0(i+1,k) + theta_0(i,k+1) + theta_0(i+1,k+1))/4.
		theta_1_u(i,k) = (theta_1(i,k) + theta_1(i+1,k) + theta_1(i,k+1) + theta_1(i+1,k+1))/4.

		IF(pi_1(i,k) == undef .OR. pi_1(i+1,k) == undef) STOP "pi_1 is WRONG!!!"
		IF(rho_0(i,k) == undef .OR. rho_0(i+1,k) == undef) STOP "rho_0 is WRONG!!!"
		IF(w(i,k) == undef .OR. w(i+1,k) == undef .OR. w(i,k+1) == undef .OR. w(i+1,k+1) == undef) STOP "w_u is WRONG!!!"
		IF(theta(i,k) == undef .OR. theta(i+1,k) == undef .OR. theta(i,k+1) == undef .OR. theta(i+1,k+1) == undef) THEN
			!CALL debug_ascii_output(theta)
			WRITE(*,*) i, k
			STOP "theta_u is WRONG!!!"
		END IF
		!IF(theta(i,k) == undef .OR. theta(i+1,k) == undef .OR. theta(i,k+1) == undef .OR. theta(i+1,k+1) == undef) STOP "theta_u is WRONG!!!"
		IF(theta_0(i,k) == undef .OR. theta_0(i+1,k) == undef .OR. theta_0(i,k+1) == undef .OR. theta_0(i+1,k+1) == undef) STOP "theta_0_u is WRONG!!!"
		IF(theta_1(i,k) == undef .OR. theta_1(i+1,k) == undef .OR. theta_1(i,k+1) == undef .OR. theta_1(i+1,k+1) == undef) STOP "theta_0_u is WRONG!!!"
	END DO
END DO
! To pi-grid (u, w, theta)
DO i = its, ite + 1
	DO k = kts - 1, kte + 1
		u_pi(i,k) = (u(i-1,k) + u(i,k))/2.
		w_pi(i,k) = (w(i,k+1) + w(i,k))/2.
		theta_pi(i,k) = (theta(i,k+1) + theta(i,k))/2.
		theta_0_pi(i,k) = (theta_0(i,k+1) + theta_0(i,k))/2.
		theta_1_pi(i,k) = (theta_1(i,k+1) + theta_1(i,k))/2.

		IF(u(i-1,k) == undef .OR. u(i,k) == undef) STOP "u_pi is WRONG!!!"
		IF(w(i,k+1) == undef .OR. w(i,k) == undef) STOP "w_pi is WRONG!!!"
		IF(theta(i,k+1) == undef .OR. theta(i,k) == undef) STOP "theta_pi is WRONG!!!"
		IF(theta_0(i,k+1) == undef .OR. theta_0(i,k) == undef) STOP "theta_0_pi is WRONG!!!"
		IF(theta_1(i,k+1) == undef .OR. theta_1(i,k) == undef) STOP "theta_1_pi is WRONG!!!"
	END DO
END DO
! To w-grid (pi, rho, u)
DO i = its, ite + 1
	DO k = kts, kte + 1
		pi_1_w(i,k) = (pi_1(i,k-1) + pi_1(i,k))/2.
		rho_0_w(i,k) = (rho_0(i,k-1) + rho_0(i,k))/2.
		u_w(i,k) = (u(i,k) + u(i-1,k) + u(i,k-1) + u(i-1,k-1))/4.

		IF(pi_1(i,k-1) == undef .OR. pi_1(i,k) == undef) STOP "pi_1_w is WRONG!!!"
		IF(rho_0(i,k-1) == undef .OR. rho_0(i,k) == undef) STOP "rho_0_w is WRONG!!!"
		IF(u(i,k) == undef .OR. u(i-1,k) == undef .OR. u(i,k-1) == undef .OR. u(i-1,k-1) == undef) STOP "u_w is WRONG!!!"
	END DO
END DO
! To v(virtual)-grid (pi, rho, u, w, theta)
DO i = its, ite
	DO k = kts, kte + 1
		pi_1_v(i,k) = (pi_1(i,k) + pi_1(i + 1,k) + pi_1(i,k - 1) + pi_1(i + 1,k - 1))/4.
		rho_0_v(i,k) = (rho_0(i,k) + rho_0(i + 1,k) + rho_0(i,k - 1) + rho_0(i + 1,k - 1))/4.
		u_v(i,k) = (u(i,k - 1) + u(i,k))/2.
		w_v(i,k) = (w(i,k) + w(i + 1,k))/2.
		theta_v(i,k) = (theta(i,k) + theta(i + 1,k))/2.
		theta_0_v(i,k) = (theta_0(i,k) + theta_0(i + 1,k))/2.
		theta_1_v(i,k) = (theta_1(i,k) + theta_1(i + 1,k))/2.

		IF(pi_1(i,k) == undef .OR. pi_1(i+1,k) == undef .OR. pi_1(i,k-1) == undef .OR. pi_1(i+1,k-1) == undef) STOP "pi_1_v is WRONG!!!"
		IF(rho_0(i,k) == undef .OR. rho_0(i+1,k) == undef .OR. rho_0(i,k-1) == undef .OR. rho_0(i+1,k-1) == undef) STOP "rho_0_v is WRONG!!!"
		IF(u(i,k-1) == undef .OR. u(i,k) == undef) STOP "u_v is WRONG!!!"
		IF(w(i,k) == undef .OR. w(i+1,k) == undef) STOP "w_v is WRONG!!!"
		IF(theta(i,k) == undef .OR. theta(i+1,k) == undef) STOP "theta_v is WRONG!!!"
		IF(theta_0(i,k) == undef .OR. theta_0(i+1,k) == undef) STOP "theta_0_v is WRONG!!!"
		IF(theta_1(i,k) == undef .OR. theta_1(i+1,k) == undef) STOP "theta_1_v is WRONG!!!"
	END DO
END DO
!CALL debug_ascii_output(w)
!=================================================
END SUBROUTINE basic_interpolate
!=================================================
!/////////////////////////////////////////////////////////////////////

!=================================================
END MODULE sp_module_interpolate
!=================================================

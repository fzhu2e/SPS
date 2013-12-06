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
REAL(kd), DIMENSION(ims:ime,kms:kme) :: u_pi, u_w, u_vir
REAL(kd), DIMENSION(ims:ime,kms:kme) :: w_pi, w_u, w_vir
REAL(kd), DIMENSION(ims:ime,kms:kme) :: w_hat, w_hat_pi, w_hat_u, w_hat_vir ! Calculated, not interpolated.
REAL(kd), DIMENSION(ims:ime,kms:kme) :: pi_1_u, pi_1_w, pi_1_vir
REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_pi, theta_u, theta_vir
REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_0_pi, theta_0_u, theta_0_vir
REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_1_pi, theta_1_u, theta_1_vir
REAL(kd), DIMENSION(ims:ime,kms:kme) :: rho_0_u, rho_0_w, rho_0_vir
!=================================================
CONTAINS
!=================================================
!/////////////////////////////////////////////////////////////////////
!=================================================
! Basic interpolate.
!=================================================
SUBROUTINE basic_interpolate(u,w,pi_1,theta,theta_1)
IMPLICIT NONE
!=================================================
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_1
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_1
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL debug_undef_all( u_pi,u_w,u_vir,                               &
                      w_pi,w_u,w_vir,                               &
                      !w_hat_pi,w_hat_u,w_hat_vir,                   &
                      pi_1_u,pi_1_w,pi_1_vir,                       &
                      theta_pi,theta_u,theta_vir,                   &
                      theta_1_pi,theta_1_u,theta_1_vir              )
!CALL debug_ascii_output(u_pi)
!CALL debug_test_boundary(u_pi)
!CALL debug_ascii_output(pi_1)
!CALL debug_ascii_output(pi_1_u)
!=================================================
!CALL u2w(u,u_w)
!CALL u2pi(u,u_pi)
!CALL u2vir(u,u_vir)
!CALL w2pi(w,w_pi)
!CALL w2vir(w,w_vir)
!CALL pi2u(pi_1,pi_1_u)
!CALL pi2vir(pi_1,pi_1_vir)

!CALL w2pi(theta,theta_pi)
!CALL w2pi(wGrid%qv,piGrid%qv)
!CALL w2pi(wGrid%qc,piGrid%qc)
!CALL w2pi(wGrid%qr,piGrid%qr)

!CALL w2pi(wGrid%theta_v,piGrid%theta_v)
!CALL pi2w(piGrid%pi,wGrid%pi)
!=================================================

CALL set_area_u
CALL set_area_expand(expand)

!CALL debug_ascii_output(pi_1)
IF (ANY(pi_1(imin:imax+1,kmin:kmax) == undef)) STOP "pi_1_u is WRONG!!!"
IF (ANY(w(imin:imax+1,kmin:kmax) == undef)) STOP "w_u is WRONG!!!"
IF (ANY(theta(imin:imax+1,kmin:kmax+1) == undef)) STOP "theta_u is WRONG!!!"
IF (ANY(theta_1(imin:imax+1,kmin:kmax+1) == undef)) STOP "theta_0_u is WRONG!!!"

!IF (ANY(b_pi(kmin:kmax) == undef) .OR. ANY(VertA_u(imin:imax) == undef) .OR. ANY(VertB_u(imin:imax,kmin:kmax) == undef)) STOP "w_hat_u is WRONG!!!"

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		pi_1_u(i,k) = (pi_1(i,k) + pi_1(i+1,k))/2.
		w_u(i,k) = (w(i,k) + w(i+1,k) + w(i,k+1) + w(i+1,k+1))/4.
		theta_u(i,k) = (theta(i,k) + theta(i+1,k) + theta(i,k+1) + theta(i+1,k+1))/4.
		theta_1_u(i,k) = (theta_1(i,k) + theta_1(i+1,k) + theta_1(i,k+1) + theta_1(i+1,k+1))/4.

		!w_hat_u(i,k) = (w_u(i,k) - u(i,k)*b_pi(k)*VertA_u(i))/VertB_u(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL set_area_pi
CALL set_area_expand(expand)

IF (ANY(u(imin-1:imax,kmin:kmax) == undef)) STOP "u_pi is WRONG!!!"
IF (ANY(w(imin:imax,kmin:kmax+1) == undef)) STOP "w_pi is WRONG!!!"
IF (ANY(theta(imin:imax,kmin:kmax+1) == undef)) STOP "theta_pi is WRONG!!!"
IF (ANY(theta_1(imin:imax,kmin:kmax+1) == undef)) STOP "theta_1_pi is WRONG!!!"

!IF (ANY(b_pi(kmin:kmax) == undef) .OR. ANY(VertA_pi(imin:imax) == undef) .OR. ANY(VertB_pi(imin:imax,kmin:kmax) == undef)) STOP "w_hat_u is WRONG!!!"

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		u_pi(i,k) = (u(i-1,k) + u(i,k))/2.
		w_pi(i,k) = (w(i,k+1) + w(i,k))/2.
		theta_pi(i,k) = (theta(i,k+1) + theta(i,k))/2.
		theta_1_pi(i,k) = (theta_1(i,k+1) + theta_1(i,k))/2.

		!w_hat_pi(i,k) = (w_pi(i,k) - u_pi(i,k)*b_pi(k)*VertA_pi(i))/VertB_pi(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL set_area_w
CALL set_area_expand(expand)

IF (ANY(pi_1(imin:imax,kmin-1:kmax) == undef)) STOP "pi_1_w is WRONG!!!"
IF (ANY(u(imin-1:imax,kmin-1:kmax) == undef)) STOP "u_w is WRONG!!!"

!IF (ANY(b(kmin:kmax) == undef) .OR. ANY(VertA_pi(imin:imax) == undef) .OR. ANY(VertB_w(imin:imax,kmin:kmax) == undef)) STOP "w_hat_u is WRONG!!!"

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		pi_1_w(i,k) = (pi_1(i,k-1) + pi_1(i,k))/2.
		u_w(i,k) = (u(i,k) + u(i-1,k) + u(i,k-1) + u(i-1,k-1))/4.

		!w_hat(i,k) = (w(i,k) - u_w(i,k)*b(k)*VertA_pi(i))/VertB_w(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL set_area_vir
CALL set_area_expand(expand)

IF (ANY(pi_1(imin:imax+1,kmin-1:kmax) == undef)) STOP "pi_1_vir is WRONG!!!"
IF (ANY(u(imin:imax,kmin-1:kmax) == undef)) STOP "u_vir is WRONG!!!"
IF (ANY(w(imin:imax+1,kmin:kmax) == undef)) STOP "w_vir is WRONG!!!"
IF (ANY(theta(imin:imax+1,kmin:kmax) == undef)) STOP "theta_vir is WRONG!!!"
IF (ANY(theta_1(imin:imax+1,kmin:kmax) == undef)) STOP "theta_1_vir is WRONG!!!"

!IF (ANY(b(kmin:kmax) == undef) .OR. ANY(VertA_u(imin:imax) == undef) .OR. ANY(VertB_vir(imin:imax,kmin:kmax) == undef)) STOP "w_hat_u is WRONG!!!"
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		pi_1_vir(i,k) = (pi_1(i,k) + pi_1(i + 1,k) + pi_1(i,k - 1) + pi_1(i + 1,k - 1))/4.
		u_vir(i,k) = (u(i,k - 1) + u(i,k))/2.
		w_vir(i,k) = (w(i,k) + w(i + 1,k))/2.
		theta_vir(i,k) = (theta(i,k) + theta(i + 1,k))/2.
		theta_1_vir(i,k) = (theta_1(i,k) + theta_1(i + 1,k))/2.

		!w_hat_vir(i,k) = (w_vir(i,k) - u_vir(i,k)*b(k)*VertA_u(i))/VertB_vir(i,k)
	END DO
END DO
!OMP END PARALLEL DO

!CALL debug_ascii_output(w)
!=================================================
END SUBROUTINE basic_interpolate
!=================================================
!/////////////////////////////////////////////////////////////////////

!=================================================
SUBROUTINE u2pi(var_u,var_pi)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_pi
INTEGER :: i, k
CALL set_area_pi
CALL set_area_expand(expand)
DO k = kmin, kmax
	DO i = imin, imax
		var_pi(i,k) = (var_u(i-1,k) + var_u(i,k))/2.
	END DO
END DO
END SUBROUTINE u2pi
!=================================================

!=================================================
SUBROUTINE u2vir(var_u,var_vir)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_vir
INTEGER :: i, k
CALL set_area_vir
CALL set_area_expand(expand)
DO k = kmin, kmax
	DO i = imin, imax
		var_vir(i,k) = (var_u(i,k) + var_u(i,k-1))/2.
	END DO
END DO
END SUBROUTINE u2vir
!=================================================

!=================================================
SUBROUTINE u2w(var_u,var_w)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_u 
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_w 
INTEGER :: i, k
CALL set_area_w
CALL set_area_expand(expand)
DO k = kmin, kmax
	DO i = imin, imax
		var_w(i,k) = (var_u(i,k) + var_u(i,k-1) + var_u(i-1,k) + var_u(i-1,k-1))/4.
	END DO
END DO
END SUBROUTINE u2w
!=================================================

!=================================================
SUBROUTINE pi2u(var_pi,var_u)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_pi
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_u
INTEGER :: i, k
CALL set_area_u
CALL set_area_expand(expand)
DO k = kmin, kmax
	DO i = imin, imax
		var_u(i,k) = (var_pi(i,k) + var_pi(i+1,k))/2.
	END DO
END DO
END SUBROUTINE pi2u
!=================================================

!=================================================
SUBROUTINE pi2w(var_pi,var_w)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_pi
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_w
INTEGER :: i, k
CALL set_area_w
CALL set_area_expand(expand)
DO k = kmin, kmax
	DO i = imin, imax
		var_w(i,k) = (var_pi(i,k-1) + var_pi(i,k))/2.
	END DO
END DO
END SUBROUTINE pi2w
!=================================================

!=================================================
SUBROUTINE pi2vir(var_pi,var_vir)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_pi
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_vir
INTEGER :: i, k
CALL set_area_vir
CALL set_area_expand(expand)
DO k = kmin, kmax
	DO i = imin, imax
		var_vir(i,k) = (var_pi(i,k-1) + var_pi(i,k) + var_pi(i+1,k) + var_pi(i+1,k-1))/4.
	END DO
END DO
END SUBROUTINE pi2vir
!=================================================

!=================================================
SUBROUTINE w2vir(var_w,var_vir)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_vir
INTEGER :: i, k
CALL set_area_vir
CALL set_area_expand(expand)
DO k = kmin, kmax
	DO i = imin, imax
		var_vir(i,k) = (var_w(i,k) + var_w(i+1,k))/2.
	END DO
END DO
END SUBROUTINE w2vir
!=================================================

!=================================================
SUBROUTINE w2pi(var_w,var_pi)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_pi
INTEGER :: i, k
CALL set_area_pi
CALL set_area_expand(expand)
DO k = kmin, kmax
	DO i = imin, imax
		var_pi(i,k) = (var_w(i,k) + var_w(i,k+1))/2.
	END DO
END DO
END SUBROUTINE w2pi
!=================================================

!=================================================
SUBROUTINE w2u(var_w,var_u)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_u
INTEGER :: i, k
CALL set_area_u
CALL set_area_expand(expand)
DO k = kmin, kmax
	DO i = imin, imax
		var_u(i,k) = (var_w(i,k) + var_w(i,k+1) + var_w(i+1,k) + var_w(i+1,k+1))/4.
	END DO
END DO
END SUBROUTINE w2u
!=================================================

!=================================================
END MODULE sp_module_interpolate
!=================================================

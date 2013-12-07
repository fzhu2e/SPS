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
! Basic interpolated variations in Arakawa-C and Charney-Phillips grids.
!-------------------------------------------------
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: u_pi, u_w, u_vir
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: w_pi, w_u, w_vir
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: pi_1_u, pi_1_w, pi_1_vir
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_pi, theta_u, theta_vir
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_1_pi, theta_1_u, theta_1_vir
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
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: u
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: w
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_1
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_1
TYPE(mainvar), INTENT(IN) :: Main
TYPE(grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!-------------------------------------------------
INTEGER :: i, k
!=================================================
!CALL debug_undef_all( u_pi,u_w,u_vir,                               &
                      !w_pi,w_u,w_vir,                               &
                      !pi_1_u,pi_1_w,pi_1_vir,                       &
                      !theta_pi,theta_u,theta_vir,                   &
                      !theta_1_pi,theta_1_u,theta_1_vir              )
!=================================================
CALL u2w(Main%u,wGrid%u)
CALL u2pi(Main%u,piGrid%u)
CALL u2vir(Main%u,virGrid%u)

CALL w2pi(Main%w,piGrid%w)
CALL w2vir(Main%w,virGrid%w)

CALL pi2vir(Main%pi_1,virGrid%pi_1)

CALL w2pi(Main%theta,piGrid%theta)
CALL w2pi(Main%theta_1,piGrid%theta_1)
!=================================================
!CALL u2w(u,u_w)
!CALL u2pi(u,u_pi)
!CALL u2vir(u,u_vir)

!CALL w2pi(w,w_pi)
!CALL w2vir(w,w_vir)
!CALL w2u(w,w_u)

!CALL pi2u(pi_1,pi_1_u)
!CALL pi2w(pi_1,pi_1_w)
!CALL pi2vir(pi_1,pi_1_vir)

!CALL w2pi(theta,theta_pi)
!CALL w2u(theta,theta_u)
!CALL w2vir(theta,theta_vir)

!CALL w2pi(theta_1,theta_1_pi)
!CALL w2u(theta_1,theta_1_u)
!CALL w2vir(theta_1,theta_1_vir)
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
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_pi(i,k) = (var_u(i-1,k) + var_u(i,k))/2.
	END DO
END DO
!OMP END PARALLEL DO
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
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_vir(i,k) = (var_u(i,k) + var_u(i,k-1))/2.
	END DO
END DO
!OMP END PARALLEL DO
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
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_w(i,k) = (var_u(i,k) + var_u(i,k-1) + var_u(i-1,k) + var_u(i-1,k-1))/4.
	END DO
END DO
!OMP END PARALLEL DO
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
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_u(i,k) = (var_pi(i,k) + var_pi(i+1,k))/2.
	END DO
END DO
!OMP END PARALLEL DO
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
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_w(i,k) = (var_pi(i,k-1) + var_pi(i,k))/2.
	END DO
END DO
!OMP END PARALLEL DO
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
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_vir(i,k) = (var_pi(i,k-1) + var_pi(i,k) + var_pi(i+1,k) + var_pi(i+1,k-1))/4.
	END DO
END DO
!OMP END PARALLEL DO
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
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_vir(i,k) = (var_w(i,k) + var_w(i+1,k))/2.
	END DO
END DO
!OMP END PARALLEL DO
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
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_pi(i,k) = (var_w(i,k) + var_w(i,k+1))/2.
	END DO
END DO
!OMP END PARALLEL DO
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
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_u(i,k) = (var_w(i,k) + var_w(i,k+1) + var_w(i+1,k) + var_w(i+1,k+1))/4.
	END DO
END DO
!OMP END PARALLEL DO
END SUBROUTINE w2u
!=================================================

!=================================================
END MODULE sp_module_interpolate
!=================================================

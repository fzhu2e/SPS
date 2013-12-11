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


!=================================================
! Calculate virtual theta
!=================================================
SUBROUTINE calc_virTheta(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE (grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!=================================================
! theta_v, theta_M_0, theta_M_1
!-------------------------------------------------
wGrid%theta_v = wGrid%theta*(1. + 0.61*wGrid%qv)
wGrid%theta_M = wGrid%theta_v*(1. - wGrid%qc)
wGrid%theta_M_0 = wGrid%theta_0*(1. + 0.61*wGrid%qv)*(1. - wGrid%qc)
wGrid%theta_M_1 = wGrid%theta_M - wGrid%theta_M_0

!uGrid%theta_M_0 = uGrid%theta_0*(1. + 0.61*uGrid%qv)*(1. - uGrid%qc)
!piGrid%theta_M_0 = piGrid%theta_0*(1. + 0.61*piGrid%qv)*(1. - piGrid%qc)
!virGrid%theta_M_0 = virGrid%theta_0*(1. + 0.61*virGrid%qv)*(1. - virGrid%qc)
CALL w2u(wGrid%theta_M_0,uGrid%theta_M_0)
CALL w2pi(wGrid%theta_M_0,piGrid%theta_M_0)
!CALL w2vir(wGrid%theta_M_0,virGrid%theta_M_0)
!=================================================
END SUBROUTINE calc_virTheta
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

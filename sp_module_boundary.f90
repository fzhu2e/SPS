!=================================================
! Super-Parametertization System (SPS)
!-------------------------------------------------
! Version: 0.2
! Author: Feng Zhu
! Email: zhuf.atmos@gmail.com
! Date: 2014-06-12 18:18:45
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!=================================================
MODULE sp_module_boundary
USE sp_module_constant
USE sp_module_model
USE sp_module_gridvar
USE sp_module_debug
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================
! Initiate.
!=================================================
SUBROUTINE update_boundary(u, w, wGrid, pi_1, theta,                       &
                           qv, qc, qr, qi, qs, qg,                         &
                           rho_0_pi, rho_0_u, rho_0_w, rho_0_vir,          &
                           theta_0_w, pi_0_pi)
IMPLICIT NONE
TYPE(grid), INTENT(IN), OPTIONAL :: wGrid
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: u        ! wind speed along x-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: w        ! wind speed along z-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: pi_1     ! pi'
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: rho_0_pi, rho_0_u, rho_0_w, rho_0_vir, theta_0_w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: qv, qc, qr, qi, qs, qg
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: pi_0_pi
!-------------------------------------------------
INTEGER :: i, k
!=================================================

IF (PRESENT(u)) THEN
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_vector_lateral_u(u)
	CASE (2)
		CALL periodic_lateral_u(u)
	CASE (3)
		CALL open_lateral_u(u)
		!CALL open_lateral_u_vector(u)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_pi(u)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_pi(u)
	CASE (2)
		CALL open_top_pi(u)
	CASE DEFAULT
		STOP "Wrong vertical boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(w)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(w)
	CASE (2)
		CALL periodic_lateral_pi(w)
	CASE (3)
		CALL open_lateral_pi(w)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_vector_bottom_w(w,wGrid)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_vector_top_w(w)
	CASE (2)
		CALL open_top_w(w)
		!CALL open_top_w_vector(w)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

! pi-grid
IF (PRESENT(pi_1)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(pi_1)
	CASE (2)
		CALL periodic_lateral_pi(pi_1)
	CASE (3)
		CALL open_lateral_pi(pi_1)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_pi(pi_1)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_pi(pi_1)
	CASE (2)
		CALL open_top_pi(pi_1)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(theta)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(theta)
	CASE (2)
		CALL periodic_lateral_pi(theta)
	CASE (3)
		CALL open_lateral_pi(theta)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_w(theta)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(theta)
	CASE (2)
		CALL open_top_w(theta)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(qv)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(qv)
	CASE (2)
		CALL periodic_lateral_pi(qv)
	CASE (3)
		CALL open_lateral_pi(qv)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_w(qv)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(qv)
	CASE (2)
		CALL open_top_w(qv)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(qc)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(qc)
	CASE (2)
		CALL periodic_lateral_pi(qc)
	CASE (3)
		CALL open_lateral_pi(qc)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_w(qc)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(qc)
	CASE (2)
		CALL open_top_w(qc)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(qr)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(qr)
	CASE (2)
		CALL periodic_lateral_pi(qr)
	CASE (3)
		CALL open_lateral_pi(qr)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_w(qr)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(qr)
	CASE (2)
		CALL open_top_w(qr)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(qi)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(qi)
	CASE (2)
		CALL periodic_lateral_pi(qi)
	CASE (3)
		CALL open_lateral_pi(qi)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_w(qi)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(qi)
	CASE (2)
		CALL open_top_w(qi)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(qs)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(qs)
	CASE (2)
		CALL periodic_lateral_pi(qs)
	CASE (3)
		CALL open_lateral_pi(qs)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_w(qs)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(qs)
	CASE (2)
		CALL open_top_w(theta)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(qg)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(qg)
	CASE (2)
		CALL periodic_lateral_pi(qg)
	CASE (3)
		CALL open_lateral_pi(qg)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_w(qg)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(qg)
	CASE (2)
		CALL open_top_w(qg)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(rho_0_pi)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(rho_0_pi)
	CASE (2)
		CALL periodic_lateral_pi(rho_0_pi)
	CASE (3)
		CALL open_lateral_pi(rho_0_pi)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_pi(rho_0_pi)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_pi(rho_0_pi)
	CASE (2)
		CALL open_top_pi(rho_0_pi)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(rho_0_u)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_u(rho_0_u)
	CASE (2)
		CALL periodic_lateral_u(rho_0_u)
	CASE (3)
		CALL open_lateral_u(rho_0_u)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_pi(rho_0_u)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_pi(rho_0_u)
	CASE (2)
		CALL open_top_pi(rho_0_u)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(rho_0_w)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(rho_0_w)
	CASE (2)
		CALL periodic_lateral_pi(rho_0_w)
	CASE (3)
		CALL open_lateral_pi(rho_0_w)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_w(rho_0_w)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(rho_0_w)
	CASE (2)
		CALL open_top_w(rho_0_w)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(rho_0_vir)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_u(rho_0_vir)
	CASE (2)
		CALL periodic_lateral_u(rho_0_vir)
	CASE (3)
		CALL open_lateral_u(rho_0_vir)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_w(rho_0_vir)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(rho_0_vir)
	CASE (2)
		CALL open_top_w(rho_0_vir)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(theta_0_w)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(theta_0_w)
	CASE (2)
		CALL periodic_lateral_pi(theta_0_w)
	CASE (3)
		CALL open_lateral_pi(theta_0_w)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_w(theta_0_w)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(theta_0_w)
	CASE (2)
		CALL open_top_w(theta_0_w)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(pi_0_pi)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(pi_0_pi)
	CASE (2)
		CALL periodic_lateral_pi(pi_0_pi)
	CASE (3)
		CALL open_lateral_pi(pi_0_pi)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_pi(pi_0_pi)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_pi(pi_0_pi)
	CASE (2)
		CALL open_top_pi(pi_0_pi)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF
!=================================================
END SUBROUTINE update_boundary
!=================================================

!=================================================
! No Flux - Scalar - Bottom [w, pi]
!=================================================
SUBROUTINE no_flux_scalar_bottom_w(scalar)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: scalar
!-------------------------------------------------
CALL set_area_w
scalar(:,kms:kmin-1) = scalar(:,2*kmin-kms:kmin+1:-1)
END SUBROUTINE no_flux_scalar_bottom_w
!=================================================

!=================================================
SUBROUTINE no_flux_scalar_bottom_pi(scalar)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: scalar
!-------------------------------------------------
CALL set_area_pi
scalar(:,kms:kmin-1) = scalar(:,2*kmin-kms-1:kmin:-1)
END SUBROUTINE no_flux_scalar_bottom_pi
!=================================================

!=================================================
! No Flux - Scalar - Top [w, pi]
!=================================================
SUBROUTINE no_flux_scalar_top_w(scalar)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: scalar
!-------------------------------------------------
CALL set_area_w
scalar(:,kme:kmax+1:-1) = scalar(:,2*kmax-kme:kmax-1)
END SUBROUTINE no_flux_scalar_top_w
!=================================================

!=================================================
SUBROUTINE no_flux_scalar_top_pi(scalar)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: scalar
!-------------------------------------------------
CALL set_area_pi
scalar(:,kme:kmax+1:-1) = scalar(:,2*kmax-kme+1:kmax)
END SUBROUTINE no_flux_scalar_top_pi
!=================================================

!=================================================
! No Flux - Vector - Top [w]
!=================================================
SUBROUTINE no_flux_vector_top_w(vector)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: vector
!-------------------------------------------------
CALL set_area_w
vector(:,kme:kmax+1:-1) = vector(:,2*kmax-kme:kmax-1)
vector(:,kmin) = 0
vector(:,kmax) = 0
END SUBROUTINE no_flux_vector_top_w
!=================================================

!=================================================
! No Flux - Vector - Bottom [w]
!=================================================
SUBROUTINE no_flux_vector_bottom_w(vector,wGrid)
IMPLICIT NONE
TYPE(grid), INTENT(IN) :: wGrid
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: vector
!-------------------------------------------------
CALL set_area_w
vector(:,kms:kmin-1) = - vector(:,2*kmin-kms:kmin+1:-1)
vector(:,kmin) = wGrid%u(:,kmin)*wGrid%PzsPx(:)
END SUBROUTINE no_flux_vector_bottom_w
!=================================================

!=================================================
! No Flux - Scalar - Lateral [pi, u]
!=================================================
SUBROUTINE no_flux_scalar_lateral_pi(scalar)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: scalar
!-------------------------------------------------
CALL set_area_pi
scalar(ims:imin-1,:) = scalar(2*imin-ims-1:imin:-1,:)
scalar(imax+1:ime,:) = scalar(imax:2*imax-ime+1:-1,:)
END SUBROUTINE no_flux_scalar_lateral_pi
!=================================================

!=================================================
SUBROUTINE no_flux_scalar_lateral_u(scalar)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: scalar
!-------------------------------------------------
CALL set_area_u
scalar(ims:imin-1,:) = scalar(2*imin-ims:imin+1:-1,:)
scalar(imax+1:ime,:) = scalar(imax-1:2*imax-ime:-1,:)
END SUBROUTINE no_flux_scalar_lateral_u
!=================================================

!=================================================
! No Flux - Vector - Lateral [u]
!=================================================
SUBROUTINE no_flux_vector_lateral_u(vector)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: vector
!-------------------------------------------------
CALL set_area_u
vector(ims:imin-1,:) = - vector(2*imin-ims:imin+1:-1,:)
vector(imax+1:ime,:) = - vector(imax-1:2*imax-ime:-1,:)
vector(imin,:) = 0.
vector(imax,:) = 0.
END SUBROUTINE no_flux_vector_lateral_u
!=================================================

!=================================================
! Periodic - Lateral [pi, u]
!=================================================
SUBROUTINE periodic_lateral_pi(var)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: var
!-------------------------------------------------
CALL set_area_pi
var(ims:imin-1,:) = var(imax-(imin-1-ims):imax,:)
var(imax+1:ime,:) = var(imin:imin+ime-(imax+1),:)
END SUBROUTINE periodic_lateral_pi
!=================================================

!=================================================
SUBROUTINE periodic_lateral_u(var)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: var
!-------------------------------------------------
CALL set_area_u
var(ims:imin-1,:) = var(imax-(imin-1-ims):imax,:)
var(imax+1:ime,:) = var(imin:imin-(imax+1-ime),:)
END SUBROUTINE periodic_lateral_u
!=================================================

!=================================================
! Open - Lateral [pi, u]
!=================================================
SUBROUTINE open_lateral_pi(var)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: var
INTEGER :: i
!-------------------------------------------------
CALL set_area_pi
DO i = ims, imin-1
	var(i,:) = var(imin,:)
END DO
DO i = imax+1, ime
	var(i,:) = var(imax,:)
END DO
END SUBROUTINE open_lateral_pi
!=================================================

!=================================================
SUBROUTINE open_lateral_u(var)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: var
INTEGER :: i
!-------------------------------------------------
CALL set_area_u
DO i = ims, imin-1
	var(i,:) = var(imin,:)
END DO
DO i = imax+1, ime
	var(i,:) = var(imax,:)
END DO
END SUBROUTINE open_lateral_u
!=================================================

!=================================================
SUBROUTINE open_lateral_u_vector(var)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: var
INTEGER :: i
!-------------------------------------------------
CALL set_area_u
DO i = ims, imin
	var(i,:) = 0.
END DO
DO i = imax, ime
	var(i,:) = 0.
END DO
END SUBROUTINE open_lateral_u_vector
!=================================================

!=================================================
SUBROUTINE open_top_w(var)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: var
INTEGER :: k
!-------------------------------------------------
CALL set_area_w
DO k = kmax+1, kme
	var(:,k) = var(:,kmax)
END DO
END SUBROUTINE open_top_w
!=================================================

!=================================================
SUBROUTINE open_top_w_vector(var)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: var
INTEGER :: k
!-------------------------------------------------
CALL set_area_w
DO k = kmax, kme
	var(:,k) = 0.
END DO
END SUBROUTINE open_top_w_vector
!=================================================

!=================================================
SUBROUTINE open_top_pi(var)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: var
INTEGER :: k
!-------------------------------------------------
CALL set_area_pi
DO k = kmax+1, kme
	var(:,k) = var(:,kmax)
END DO
END SUBROUTINE open_top_pi
!=================================================

!=================================================
SUBROUTINE calc_tau_top(uGrid,wGrid)
IMPLICIT NONE
TYPE(grid), INTENT(INOUT) :: uGrid, wGrid
INTEGER :: i, k
!-------------------------------------------------
CALL set_area_u
DO k = kmin, kmax
	DO i = imin, imax
		IF (uGrid%zz(i,k) < ztop - s ) THEN
			uGrid%tau(i,k) = 0.
		ELSE
			IF (uGrid%tau(i,k) /= undef) THEN
				uGrid%tau(i,k) = MAX(tau0*((uGrid%zz(i,k) - (ztop - s))/s)**4,uGrid%tau(i,k))
			ELSE
				uGrid%tau(i,k) = tau0*((uGrid%zz(i,k) - (ztop - s))/s)**4
			END IF
		END IF
	END DO
END DO

CALL set_area_w
DO k = kmin, kmax
	DO i = imin, imax
		IF (wGrid%zz(i,k) < ztop - s ) THEN
			wGrid%tau(i,k) = 0.
		ELSE
			IF (wGrid%tau(i,k) /= undef) THEN
				wGrid%tau(i,k) = MAX(tau0*((wGrid%zz(i,k) - (ztop - s))/s)**4,wGrid%tau(i,k))
			ELSE
				wGrid%tau(i,k) = tau0*((wGrid%zz(i,k) - (ztop - s))/s)**4
			END IF
		END IF
	END DO
END DO

END SUBROUTINE calc_tau_top
!=================================================

!=================================================
SUBROUTINE calc_tau_lateral(uGrid,wGrid)
IMPLICIT NONE
TYPE(grid), INTENT(INOUT) :: uGrid, wGrid
INTEGER :: i, k
!-------------------------------------------------
CALL set_area_u
DO k = kmin, kmax
	DO i = imin, imax
		IF (uGrid%xx(i) <= s ) THEN
			IF (uGrid%tau(i,k) /= undef) THEN
				uGrid%tau(i,k) = MAX(tau0*((s - uGrid%xx(i))/s)**4,uGrid%tau(i,k))
			ELSE
				uGrid%tau(i,k) = tau0*((s - uGrid%xx(i))/s)**4
			END IF
		ELSE IF (uGrid%xx(i) >= nx*dx - s ) THEN
			IF (uGrid%tau(i,k) /= undef) THEN
				uGrid%tau(i,k) = MAX(tau0*((uGrid%xx(i) - (nx*dx - s))/s)**4,uGrid%tau(i,k))
			ELSE
				uGrid%tau(i,k) = tau0*((uGrid%xx(i) - (nx*dx - s))/s)**4
			END IF
		ELSE
			uGrid%tau(i,k) = 0.
		END IF
	END DO
END DO

CALL set_area_w
DO k = kmin, kmax
	DO i = imin, imax
		IF (wGrid%xx(i) <= s ) THEN
			IF (wGrid%tau(i,k) /= undef) THEN
				wGrid%tau(i,k) = MAX(tau0*((s - wGrid%xx(i))/s)**4,wGrid%tau(i,k))
			ELSE
				wGrid%tau(i,k) = tau0*((s - wGrid%xx(i))/s)**4
			END IF
		ELSE IF (wGrid%xx(i) >= nx*dx - s ) THEN
			IF (wGrid%tau(i,k) /= undef) THEN
				wGrid%tau(i,k) = MAX(tau0*((wGrid%xx(i) - (nx*dx - s))/s)**4,wGrid%tau(i,k))
			ELSE
				wGrid%tau(i,k) = tau0*((wGrid%xx(i) - (nx*dx - s))/s)**4
			END IF
		ELSE
			wGrid%tau(i,k) = 0.
		END IF
	END DO
END DO

END SUBROUTINE calc_tau_lateral
!=================================================

!=================================================
END MODULE sp_module_boundary
!=================================================

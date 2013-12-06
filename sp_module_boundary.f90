!=================================================
! The boundary module of SPS-dynamic.
!-------------------------------------------------
! Version: 0.01
! Author: Zhu F.
! Email: lyricorpse@gmail.com
! Date: 2013-04-20 12:20:45 
!=================================================
MODULE sp_module_boundary
USE sp_module_constant
USE sp_module_model
USE sp_module_debug
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================
! Initiate.
!=================================================
SUBROUTINE update_boundary(u,w,pi_1,theta,theta_1,                 &
                           theta_0,theta_0_pi,theta_0_u,theta_0_vir, &
                           rho_0,rho_0_w,rho_0_u,rho_0_vir           )
IMPLICIT NONE
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: u        ! wind speed along x-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: w        ! wind speed along z-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: pi_1     ! pi'
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta_1  ! theta'
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta_0  ! theta0
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta_0_pi
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta_0_u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta_0_vir
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: rho_0  ! rho0
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: rho_0_u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: rho_0_w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: rho_0_vir
!-------------------------------------------------
INTEGER :: i, k
!=================================================

IF (PRESENT(u)) THEN
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_virector_lateral_u(u)
	CASE (2)
		CALL periodic_lateral_u(u)
	CASE (3)
		CALL open_lateral_u(u)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_pi(u)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_pi(u)
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
	
	CALL no_flux_virector_bottom_w(w)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_virector_top_w(w)
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
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(theta_1)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(theta_1)
	CASE (2)
		CALL periodic_lateral_pi(theta_1)
	CASE (3)
		CALL open_lateral_pi(theta_1)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_w(theta_1)
	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(theta_1)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF


IF (PRESENT(theta_0)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(theta_0)
	CASE (2)
		CALL periodic_lateral_pi(theta_0)
	CASE (3)
		CALL open_lateral_pi(theta_0)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_w(theta_0)
	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(theta_0)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(theta_0_vir)) THEN
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_u(theta_0_vir)
	CASE (2)
		CALL periodic_lateral_u(theta_0_vir)
	CASE (3)
		CALL open_lateral_u(theta_0_vir)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT

	CALL no_flux_scalar_bottom_w(theta_0_vir)
	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(theta_0_vir)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
END IF

IF (PRESENT(theta_0_pi)) THEN
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(theta_0_pi)
	CASE (2)
		CALL periodic_lateral_pi(theta_0_pi)
	CASE (3)
		CALL open_lateral_pi(theta_0_pi)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT

	CALL no_flux_scalar_bottom_pi(theta_0_pi)
	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_pi(theta_0_pi)
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
END IF

IF (PRESENT(theta_0_u)) THEN
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_u(theta_0_u)
	CASE (2)
		CALL periodic_lateral_u(theta_0_u)
	CASE (3)
		CALL open_lateral_u(theta_0_u)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT

	CALL no_flux_scalar_bottom_pi(theta_0_u)
	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_pi(theta_0_u)
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
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(rho_0)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_pi(rho_0)
	CASE (2)
		CALL periodic_lateral_pi(rho_0)
	CASE (3)
		CALL open_lateral_pi(rho_0)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	CALL no_flux_scalar_bottom_pi(rho_0)

	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_pi(rho_0)
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
SUBROUTINE no_flux_virector_top_w(vector)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: vector
!-------------------------------------------------
CALL set_area_w
vector(:,kme:kmax+1:-1) = vector(:,2*kmax-kme:kmax-1)
vector(:,kmin) = 0
vector(:,kmax) = 0
END SUBROUTINE no_flux_virector_top_w
!=================================================

!=================================================
! No Flux - Vector - Bottom [w]
!=================================================
SUBROUTINE no_flux_virector_bottom_w(vector)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: vector
!-------------------------------------------------
CALL set_area_w
vector(:,kms:kmin-1) = - vector(:,2*kmin-kms:kmin+1:-1)
vector(:,kmin) = 0.
END SUBROUTINE no_flux_virector_bottom_w
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
SUBROUTINE no_flux_virector_lateral_u(vector)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: vector
!-------------------------------------------------
CALL set_area_u
vector(ims:imin-1,:) = - vector(2*imin-ims:imin+1:-1,:)
vector(imax+1:ime,:) = - vector(imax-1:2*imax-ime:-1,:)
vector(imin,:) = 0.
vector(imax,:) = 0.
END SUBROUTINE no_flux_virector_lateral_u
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
END MODULE sp_module_boundary
!=================================================

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
                           theta_0,theta_0_pi,theta_0_u,theta_0_v, &
                           rho_0,rho_0_w,rho_0_u,rho_0_v           )
IMPLICIT NONE
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: u        ! wind speed along x-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: w        ! wind speed along z-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: pi_1     ! pi'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta_1  ! theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta_0  ! theta0
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta_0_pi
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta_0_u
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta_0_v
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: rho_0  ! rho0
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: rho_0_u
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: rho_0_w
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: rho_0_v
!-------------------------------------------------
INTEGER :: i, k
!=================================================

IF (PRESENT(u)) THEN
	SELECT CASE (LateralBoundary)
	CASE (1)
		!u(ims:imin,:) = 0
		!u(imax:ime,:) = 0
		CALL no_flux_vector_lateral_u(u)
	CASE (2)
		!u(ims:imin,:) = u(ite+ims-its-1:ite-1,:)
		!u(imax:ime,:) = u(its+1:its+ime-ite+1,:)
		CALL periodic_lateral_u(u)
	CASE (4)
		DO i = ims, its
			u(i,:) = u(its+1,:)
		END DO
		DO i = ite, ime
			u(i,:) = u(ite-1,:)
		END DO
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	!u(:,kms:kts-1) = u(:,2*kts-kms-1:kts:-1) ! PuPz = 0
	CALL no_flux_scalar_bottom_pi(u)

	SELECT CASE (UpperBoundary)
	CASE (1)
		!u(:,kte+1:kme) = u(:,kte:2*kte-kme+1:-1) ! PuPz = 0
		CALL no_flux_scalar_top_pi(u)
	CASE (4)
		DO k = kte+1, kme
			u(:,k) = u(:,kte) ! PuPz = 0
		END DO
	CASE DEFAULT
		STOP "Wrong vertical boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(w)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		!w(its:ims:-1,:) = w(its+1:2*its-ims+1,:) ! PwPx = 0.
		!w(ite+1:ime,:) = w(ite:2*ite-ime+1:-1,:) ! PwPx = 0.
		CALL no_flux_scalar_lateral_pi(w)
	CASE (2)
		!w(ims:its,:) = w(ite+ims-its:ite,:)
		!w(ite+1:ime,:) = w(its+1:its+ime-ite,:)
		CALL periodic_lateral_pi(w)
	CASE (4)
		DO i = ims, its
			w(i,:) = w(its+1,:)
		END DO
		DO i = ite, ime
			w(i,:) = w(ite,:)
		END DO
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	!w(:,kms:kts) = 0  ! Lower boundary is always a wall.
	CALL no_flux_vector_bottom_w(w)

	SELECT CASE (UpperBoundary)
	CASE (1)
		!w(:,kte+1:kme) = 0
		CALL no_flux_vector_top_w(w)
	CASE (4)
		DO k = kte+1, kme
			w(:,k) = w(:,kte)
		END DO
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(theta)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		!theta(its:ims:-1,:) = theta(its+1:2*its-ims+1,:) ! PthetaPx = 0.
		!theta(ite+1:ime,:) = theta(ite:2*ite-ime+1:-1,:) ! PthetaPx = 0.
		CALL no_flux_scalar_lateral_pi(theta)
	CASE (2)
		!theta(ims:its,:) = theta(ite+ims-its:ite,:)
		!theta(ite+1:ime,:) = theta(its+1:its+ime-ite,:)
		CALL periodic_lateral_pi(theta)
	CASE (4)
		DO i = ims, its
			theta(i,:) = theta(its+1,:)
		END DO
		DO i = ite+1, ime
			theta(i,:) = theta(ite,:)
		END DO
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	!theta(:,kts) = theta(:,kts+1) ! PthetaPz = 0.
	!theta(:,kms:kts-1) = theta(:,2*kts-kms:kts+1:-1)
	CALL no_flux_scalar_bottom_w(theta)

	SELECT CASE (UpperBoundary)
	CASE (1)
		!theta(:,kte+1) = theta(:,kte) ! PthetaPz = 0.
		!theta(:,kte+2:kme) = theta(:,kte:2*kte-kme+2:-1)
		CALL no_flux_scalar_top_w(theta)
	CASE (4)
		DO k = kte+2, kme
			theta(:,k) = theta(:,kte+1)
		END DO
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(theta_1)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		!theta_1(its:ims:-1,:) = theta_1(its+1:2*its-ims+1,:) ! Ptheta_1Px = 0.
		!theta_1(ite+1:ime,:) = theta_1(ite:2*ite-ime+1:-1,:) ! Ptheta_1Px = 0.
		CALL no_flux_scalar_lateral_pi(theta_1)
	CASE (2)
		!theta_1(ims:its,:) = theta_1(ite+ims-its:ite,:)
		!theta_1(ite+1:ime,:) = theta_1(its+1:its+ime-ite,:)
		CALL periodic_lateral_pi(theta_1)
	CASE (4)
		DO i = ims, its
			theta_1(i,:) = theta_1(its+1,:)
		END DO
		DO i = ite+1, ime
			theta_1(i,:) = theta_1(ite,:)
		END DO
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	!theta_1(:,kts) = theta_1(:,kts+1) ! Ptheta_1Pz = 0.
	!theta_1(:,kms:kts-1) = theta_1(:,2*kts-kms:kts+1:-1)
	!theta_1(:,kte+1) = theta_1(:,kte) ! Ptheta_1Pz = 0.
	CALL no_flux_scalar_bottom_w(theta_1)
	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(theta_1)
		!theta_1(:,kte+2:kme) = theta_1(:,kte:2*kte-kme+2:-1)
	CASE (4)
		DO k = kte+2, kme
			theta_1(:,k) = theta_1(:,kte+1)
		END DO
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF


IF (PRESENT(theta_0)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		!theta_0(its:ims:-1,:) = theta_0(its+1:2*its-ims+1,:) ! Ptheta_0Px = 0.
		!theta_0(ite+1:ime,:) = theta_0(ite:2*ite-ime+1:-1,:) ! Ptheta_0Px = 0.
		CALL no_flux_scalar_lateral_pi(theta_0)
	CASE (2)
		!theta_0(ims:its,:) = theta_0(ite+ims-its:ite,:)
		!theta_0(ite+1:ime,:) = theta_0(its+1:its+ime-ite,:)
		CALL periodic_lateral_pi(theta_0)
	CASE (4)
		DO i = ims, its
			theta_0(i,:) = theta_0(its+1,:)
		END DO
		DO i = ite+1, ime
			theta_0(i,:) = theta_0(ite,:)
		END DO
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	!theta_0(:,kts) = theta_0(:,kts+1)
	!theta_0(:,kms:kts-1) = theta_0(:,2*kts-kms:kts+1:-1)
	!theta_0(:,kte+1) = theta_0(:,kte)
	CALL no_flux_scalar_bottom_w(theta_0)
	SELECT CASE (UpperBoundary)
	CASE (1)
		!theta_0(:,kte+2:kme) = theta_0(:,kte:2*kte-kme+2:-1)
		CALL no_flux_scalar_top_w(theta_0)
	CASE (4)
		DO k = kte+2, kme
			theta_0(:,k) = theta_0(:,kte+1)
		END DO
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(theta_0_v)) THEN
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_u(theta_0_v)
	CASE (2)
		CALL periodic_lateral_u(theta_0_v)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT

	CALL no_flux_scalar_bottom_w(theta_0_v)
	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(theta_0_v)
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

IF (PRESENT(rho_0_v)) THEN
	SELECT CASE (LateralBoundary)
	CASE (1)
		CALL no_flux_scalar_lateral_u(rho_0_v)
	CASE (2)
		CALL periodic_lateral_u(rho_0_v)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT

	CALL no_flux_scalar_bottom_w(rho_0_v)
	SELECT CASE (UpperBoundary)
	CASE (1)
		CALL no_flux_scalar_top_w(rho_0_v)
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
		!pi_1(its:ims:-1,:) = pi_1(its+1:2*its-ims+1,:) ! Ppi_1Px = 0.
		!pi_1(ite+1:ime,:) = pi_1(ite:2*ite-ime+1:-1,:) ! Ppi_1Px = 0.
		CALL no_flux_scalar_lateral_pi(pi_1)
	CASE (2)
		!pi_1(ims:its,:) = pi_1(ite+ims-its:ite,:)
		!pi_1(ite+1:ime,:) = pi_1(its+1:its+ime-ite,:)
		CALL periodic_lateral_pi(pi_1)
	CASE (4)
		DO i = ims, its
			pi_1(i,:) = pi_1(its+1,:)
		END DO
		DO i = ite+1, ime
			pi_1(i,:) = pi_1(ite,:)
		END DO
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	!pi_1(:,kms:kts-1) = pi_1(:,2*kts-kms-1:kts:-1) ! PuPz = 0
	CALL no_flux_scalar_bottom_pi(pi_1)

	SELECT CASE (UpperBoundary)
	CASE (1)
		!pi_1(:,kte+1:kme) = pi_1(:,kte:2*kte-kme+1:-1) ! PuPz = 0
		CALL no_flux_scalar_top_pi(pi_1)
	CASE (4)
		DO k = kte+1, kme
			pi_1(:,k) = pi_1(:,kte)
		END DO
	CASE DEFAULT
		STOP "Wrong upper boundary scheme!!!"
	END SELECT
	
END IF

IF (PRESENT(rho_0)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		!rho_0(its:ims:-1,:) = rho_0(its+1:2*its-ims+1,:) ! Ppi_1Px = 0.
		!rho_0(ite+1:ime,:) = rho_0(ite:2*ite-ime+1:-1,:) ! Ppi_1Px = 0.
		CALL no_flux_scalar_lateral_pi(rho_0)
	CASE (2)
		!rho_0(ims:its,:) = rho_0(ite+ims-its:ite,:)
		!rho_0(ite+1:ime,:) = rho_0(its+1:its+ime-ite,:)
		CALL periodic_lateral_pi(rho_0)
	CASE (4)
		DO i = ims, its
			rho_0(i,:) = rho_0(its+1,:)
		END DO
		DO i = ite+1, ime
			rho_0(i,:) = rho_0(ite,:)
		END DO
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	!rho_0(:,kms:kts-1) = rho_0(:,2*kts-kms-1:kts:-1) ! PuPz = 0
	CALL no_flux_scalar_bottom_pi(rho_0)

	SELECT CASE (UpperBoundary)
	CASE (1)
		!rho_0(:,kte+1:kme) = rho_0(:,kte:2*kte-kme+1:-1) ! PuPz = 0
		CALL no_flux_scalar_top_pi(rho_0)
	CASE (4)
		DO k = kte+1, kme
			rho_0(:,k) = rho_0(:,kte)
		END DO
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
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: scalar
!-------------------------------------------------
CALL set_area_w
scalar(:,kms:kmin-1) = scalar(:,2*kmin-kms:kmin+1:-1)
END SUBROUTINE no_flux_scalar_bottom_w
!=================================================

!=================================================
SUBROUTINE no_flux_scalar_bottom_pi(scalar)
IMPLICIT NONE
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: scalar
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
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: scalar
!-------------------------------------------------
CALL set_area_w
scalar(:,kme:kmax+1:-1) = scalar(:,2*kmax-kme:kmax-1)
END SUBROUTINE no_flux_scalar_top_w
!=================================================

!=================================================
SUBROUTINE no_flux_scalar_top_pi(scalar)
IMPLICIT NONE
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: scalar
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
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: vector
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
SUBROUTINE no_flux_vector_bottom_w(vector)
IMPLICIT NONE
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: vector
!-------------------------------------------------
CALL set_area_w
vector(:,kms:kmin-1) = - vector(:,2*kmin-kms:kmin+1:-1)
vector(:,kmin) = 0.
END SUBROUTINE no_flux_vector_bottom_w
!=================================================

!=================================================
! No Flux - Scalar - Lateral [pi, u]
!=================================================
SUBROUTINE no_flux_scalar_lateral_pi(scalar)
IMPLICIT NONE
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: scalar
!-------------------------------------------------
CALL set_area_pi
scalar(ims:imin-1,:) = scalar(2*imin-ims-1:imin:-1,:)
scalar(imax+1:ime,:) = scalar(imax:2*imax-ime+1:-1,:)
END SUBROUTINE no_flux_scalar_lateral_pi
!=================================================

!=================================================
SUBROUTINE no_flux_scalar_lateral_u(scalar)
IMPLICIT NONE
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: scalar
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
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: vector
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
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: var
!-------------------------------------------------
CALL set_area_pi
var(ims:imin-1,:) = var(imax-(imin-1-ims):imax,:)
var(imax+1:ime,:) = var(imin:imin+ime-(imax+1),:)
END SUBROUTINE periodic_lateral_pi
!=================================================

!=================================================
SUBROUTINE periodic_lateral_u(var)
IMPLICIT NONE
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: var
!-------------------------------------------------
CALL set_area_u
var(ims:imin-1,:) = var(imax-(imin-1-ims):imax,:)
var(imax+1:ime,:) = var(imin:imin-(imax+1-ime),:)
END SUBROUTINE periodic_lateral_u
!=================================================

!=================================================
END MODULE sp_module_boundary
!=================================================

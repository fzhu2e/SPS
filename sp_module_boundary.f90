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
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================
! Initiate.
!=================================================
SUBROUTINE update_boundary(u,w,pi_1,theta,theta_1,theta_0,rho_0)
IMPLICIT NONE
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: u        ! wind speed along x-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: w        ! wind speed along z-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: pi_1     ! pi'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta_1  ! theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: theta_0  ! theta0
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: rho_0  ! rho0
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! u-grid
!DO i = ims, its
!DO i = ite, ime
	!WRITE(*,*) i
!END DO
!WRITE(*,*) "----------------------------"
!DO i = ite + ims - its - 1, ite - 1
!DO i = its + 1, its + ime - ite + 1
	!WRITE(*,*) i
!END DO

IF (PRESENT(u)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		u(ims:its,:) = 0
		u(ite:ime,:) = 0
	CASE (2)
		u(ims:its,:) = u(ite+ims-its-1:ite-1,:)
		u(ite:ime,:) = u(its+1:its+ime-ite+1,:)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	u(:,kms:kts-1) = u(:,2*kts-kms-1:kts:-1) ! PuPz = 0
	u(:,kte+1:kme) = u(:,kte:2*kte-kme+1:-1) ! PuPz = 0
	
END IF

!DO k = kms, kts - 1
!DO k = kte+1, kme
	!WRITE(*,*) k
!END DO
!WRITE(*,*) "----------------------------"
!DO k = 2*kts - kms - 1, kts, - 1
!DO k = kte, 2*kte - kme + 1, - 1
	!WRITE(*,*) k
!END DO

! w-grid

!DO i = ims, its
!DO i = ite + 1, ime
	!WRITE(*,*) i
!END DO
!WRITE(*,*) "----------------------------"
!DO i = ite + ims - its, ite
!DO i = its + 1, its + ime - ite
	!WRITE(*,*) i
!END DO

IF (PRESENT(w)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		w(its:ims:-1,:) = w(its+1:2*its-ims+1,:) ! PwPx = 0.
		w(ite+1:ime,:) = w(ite:2*ite-ime+1:-1,:) ! PwPx = 0.
	CASE (2)
		w(ims:its,:) = w(ite+ims-its:ite,:)
		w(ite+1:ime,:) = w(its+1:its+ime-ite,:)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	w(:,kms:kts) = 0
	w(:,kte+1:kme) = 0
	
END IF

!DO k = kms, kts - 1
!DO k = kte + 2, kme
	!WRITE(*,*) k
!END DO
!WRITE(*,*) "----------------------------"
!DO k = 2*kts - kms, kts + 1, - 1
!DO k = kte, 2*kte - kme + 2, - 1
	!WRITE(*,*) k
!END DO

IF (PRESENT(theta)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		theta(its:ims:-1,:) = theta(its+1:2*its-ims+1,:) ! PthetaPx = 0.
		theta(ite+1:ime,:) = theta(ite:2*ite-ime+1:-1,:) ! PthetaPx = 0.
	CASE (2)
		theta(ims:its,:) = theta(ite+ims-its:ite,:)
		theta(ite+1:ime,:) = theta(its+1:its+ime-ite,:)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	!theta(:,kts) = theta(:,kts+1) ! PthetaPz = 0.
	theta(:,kms:kts-1) = theta(:,2*kts-kms:kts+1:-1)
	!theta(:,kte+1) = theta(:,kte) ! PthetaPz = 0.
	theta(:,kte+2:kme) = theta(:,kte:2*kte-kme+2:-1)
	
END IF

IF (PRESENT(theta_1)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		theta_1(its:ims:-1,:) = theta_1(its+1:2*its-ims+1,:) ! Ptheta_1Px = 0.
		theta_1(ite+1:ime,:) = theta_1(ite:2*ite-ime+1:-1,:) ! Ptheta_1Px = 0.
	CASE (2)
		theta_1(ims:its,:) = theta_1(ite+ims-its:ite,:)
		theta_1(ite+1:ime,:) = theta_1(its+1:its+ime-ite,:)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	!theta_1(:,kts) = theta_1(:,kts+1) ! Ptheta_1Pz = 0.
	theta_1(:,kms:kts-1) = theta_1(:,2*kts-kms:kts+1:-1)
	!theta_1(:,kte+1) = theta_1(:,kte) ! Ptheta_1Pz = 0.
	theta_1(:,kte+2:kme) = theta_1(:,kte:2*kte-kme+2:-1)
	
END IF


IF (PRESENT(theta_0)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		theta_0(its:ims:-1,:) = theta_0(its+1:2*its-ims+1,:) ! Ptheta_0Px = 0.
		theta_0(ite+1:ime,:) = theta_0(ite:2*ite-ime+1:-1,:) ! Ptheta_0Px = 0.
	CASE (2)
		theta_0(ims:its,:) = theta_0(ite+ims-its:ite,:)
		theta_0(ite+1:ime,:) = theta_0(its+1:its+ime-ite,:)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	!theta_0(:,kts) = theta_0(:,kts+1)
	theta_0(:,kms:kts-1) = theta_0(:,2*kts-kms:kts+1:-1)
	!theta_0(:,kte+1) = theta_0(:,kte)
	theta_0(:,kte+2:kme) = theta_0(:,kte:2*kte-kme+2:-1)
	
END IF

! pi-grid
IF (PRESENT(pi_1)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		pi_1(its:ims:-1,:) = pi_1(its+1:2*its-ims+1,:) ! Ppi_1Px = 0.
		pi_1(ite+1:ime,:) = pi_1(ite:2*ite-ime+1:-1,:) ! Ppi_1Px = 0.
	CASE (2)
		pi_1(ims:its,:) = pi_1(ite+ims-its:ite,:)
		pi_1(ite+1:ime,:) = pi_1(its+1:its+ime-ite,:)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	pi_1(:,kms:kts-1) = pi_1(:,2*kts-kms-1:kts:-1) ! PuPz = 0
	pi_1(:,kte+1:kme) = pi_1(:,kte:2*kte-kme+1:-1) ! PuPz = 0
	
END IF

IF (PRESENT(rho_0)) THEN
	
	SELECT CASE (LateralBoundary)
	CASE (1)
		rho_0(its:ims:-1,:) = rho_0(its+1:2*its-ims+1,:) ! Ppi_1Px = 0.
		rho_0(ite+1:ime,:) = rho_0(ite:2*ite-ime+1:-1,:) ! Ppi_1Px = 0.
	CASE (2)
		rho_0(ims:its,:) = rho_0(ite+ims-its:ite,:)
		rho_0(ite+1:ime,:) = rho_0(its+1:its+ime-ite,:)
	CASE DEFAULT
		STOP "Wrong lateral boundary scheme!!!"
	END SELECT
	
	rho_0(:,kms:kts-1) = rho_0(:,2*kts-kms-1:kts:-1) ! PuPz = 0
	rho_0(:,kte+1:kme) = rho_0(:,kte:2*kte-kme+1:-1) ! PuPz = 0
	
END IF
!=================================================
END SUBROUTINE update_boundary
!=================================================

!=================================================
END MODULE sp_module_boundary
!=================================================

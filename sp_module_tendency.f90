!=================================================
! The flux module of SPS-dynamic-integrate
!-------------------------------------------------
! Version: 0.15
! Author: Zhu F.
! Email: lyricorpse@gmail.com
! Date: 2013-05-04 13:59:46 
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!=================================================
MODULE sp_module_tendency
USE sp_module_constant
USE sp_module_model
USE sp_module_gridvar
USE sp_module_advection
USE sp_module_debug
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================
SUBROUTINE tendency_u(Main,tend_u,uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE(mainvar), INTENT(IN) :: Main
TYPE(grid), INTENT(IN) :: uGrid, wGrid, piGrid, virGrid
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: tend_u
!=================================================
REAL(kd), DIMENSION(ims:ime,kms:kme) :: Ppi_1Px_u = undef, Ppi_1Pzeta_u = undef

REAL(kd), DIMENSION(ims:ime,kms:kme) :: A_u = undef, P_u = undef, D_u = undef, S_u = undef
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL calc_advection_u(Main%u,A_u,uGrid,wGrid,piGrid,virGrid)

CALL ppx_u(Main%pi_1,Ppi_1Px_u)
CALL ppzeta_u(Main%pi_1,Ppi_1Pzeta_u)

IF (OpenTop == 2 .OR. OpenLateral == 2) THEN
	S_u = - uGrid%tau*Main%u
ELSE
	S_u = 0.
END IF

D_u = uGrid%Du

CALL set_area_u
!$OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		!P_u(i,k) = - Cp*uGrid%theta_M_0(i,k)*(Ppi_1Px_u(i,k) + uGrid%G(i,k)*Ppi_1Pzeta_u(i,k))
		P_u(i,k) = - Cp*uGrid%theta_M(i,k)*(Ppi_1Px_u(i,k) + uGrid%G(i,k)*Ppi_1Pzeta_u(i,k))

		tend_u(i,k) = A_u(i,k) + P_u(i,k) + D_u(i,k) + S_u(i,k)
	END DO
END DO
!$OMP END PARALLEL DO
!-------------------------------------------------
IF (ANY(ISNAN(tend_u(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITH tend_u!!!"
!=================================================
END SUBROUTINE tendency_u
!=================================================

!=================================================
SUBROUTINE tendency_w(Main,tend_w, uGrid, wGrid, piGrid, virGrid)
IMPLICIT NONE
TYPE(mainvar), INTENT(IN) :: Main
TYPE(grid), INTENT(IN) :: uGrid, wGrid, piGrid, virGrid
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: tend_w
!=================================================
REAL(kd), DIMENSION(ims:ime,kms:kme) :: A_w = undef, B_w = undef, P_w = undef, D_w = undef, S_w = undef
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: Ppi_1Pzeta_w = undef
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL calc_advection_w(Main%w,A_w,uGrid,wGrid,piGrid,virGrid)
CALL ppzeta_w(Main%pi_1,Ppi_1Pzeta_w)

IF (OpenTop == 2 .OR. OpenLateral == 2) THEN
	S_w = - wGrid%tau*Main%w
ELSE
	S_w = 0.
END IF

D_w = wGrid%Dw

CALL set_area_w
!$OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		B_w(i,k) = g*wGrid%theta_M_1(i,k)/wGrid%theta_M_0(i,k)

		!P_w(i,k) = - Cp*wGrid%theta_M_0(i,k)*wGrid%H(i)*Ppi_1Pzeta_w(i,k)
		P_w(i,k) = - Cp*wGrid%theta_M(i,k)*wGrid%H(i)*Ppi_1Pzeta_w(i,k)

		tend_w(i,k) = A_w(i,k) + B_w(i,k) + P_w(i,k) + D_w(i,k) + S_w(i,k)
	END DO
END DO
!$OMP END PARALLEL DO
!-------------------------------------------------
IF (ANY(ISNAN(tend_w(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITHT tend_w!!!"
!=================================================
END SUBROUTINE tendency_w
!=================================================

!=================================================
SUBROUTINE tendency_pi(Main,tend_pi_1,uGrid,wGrid,piGrid,virGrid )
IMPLICIT NONE
TYPE(mainvar), INTENT(IN) :: Main
TYPE(grid), INTENT(IN) :: uGrid, wGrid, piGrid, virGrid
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: tend_pi_1
!=================================================
REAL(kd), DIMENSION(ims:ime,kms:kme) :: A_pi = undef
REAL(kd), DIMENSION(ims:ime,kms:kme) :: Div_pi = undef
REAL(kd), DIMENSION(ims:ime,kms:kme) :: PuPx_pi = undef, PuPzeta_pi = undef, PwPzeta_pi = undef
REAL(kd), DIMENSION(ims:ime,kms:kme) :: pi = undef
REAL(kd), DIMENSION(ims:ime,kms:kme) :: Th_pi = undef
!-------------------------------------------------
REAL(kd) :: temp
INTEGER :: i, k
!=================================================
!-------------------------------------------------
CALL ppx_pi(Main%u,PuPx_pi)
CALL ppzeta_pi(wGrid%u,PuPzeta_pi)
CALL ppzeta_pi(Main%w,PwPzeta_pi)

pi = Main%pi_1 + piGrid%pi_0
CALL calc_advection_pi(pi,A_pi,uGrid,wGrid,piGrid,virGrid)
!CALL calc_advection_pi(Main%pi_1,A_pi,uGrid,wGrid,piGrid,virGrid)

CALL set_area_pi
!$OMP PARALLEL DO PRIVATE(temp)
DO k = kmin, kmax
	DO i = imin, imax
		temp = PuPx_pi(i,k) + piGrid%G(i,k)*PuPzeta_pi(i,k) + piGrid%H(i)*PwPzeta_pi(i,k)
		!Div_pi(i,k) = - cs*cs/Cp/piGrid%theta_M_0(i,k)*temp
		Div_pi(i,k) = - Rd*pi(i,k)/Cv*temp

		temp = (piGrid%Mtheta(i,k) + piGrid%Dtheta(i,k))*(1. + 0.61*piGrid%qv(i,k)) + (0.61*piGrid%theta(i,k)*(piGrid%Mqv(i,k) + piGrid%Dqv(i,k)))
		Th_pi(i,k) = - Rd*pi(i,k)/Cv/piGrid%theta_v(i,k)*temp

		tend_pi_1(i,k) = A_pi(i,k) + Div_pi(i,k) + Th_pi(i,k)
	END DO
END DO
!$OMP END PARALLEL DO
!-------------------------------------------------
IF (ANY(ISNAN(tend_pi_1(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITH tend_pi_1!!!"
!=================================================
END SUBROUTINE tendency_pi
!=================================================

!=================================================
SUBROUTINE tendency_q(flag,q,tend_q,uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
INTEGER, INTENT(IN) :: flag
TYPE(grid), INTENT(IN) :: uGrid, wGrid, piGrid, virGrid
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: q
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: tend_q
!=================================================
REAL(kd), DIMENSION(ims:ime,kms:kme) :: A_q = undef, D_q = undef, M_q = undef, S_q = undef
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL calc_advection_w(q,A_q,uGrid,wGrid,piGrid,virGrid)

CALL set_area_w
SELECT CASE (flag)
CASE (0)
	D_q = wGrid%Dtheta
	M_q = wGrid%Mtheta
	IF (OpenTop == 2 .OR. OpenLateral == 2) THEN
		S_q = - wGrid%tau*(wGrid%theta - wGrid%theta_0)
	ELSE
		S_q = 0.
	END IF
CASE (1)
	D_q = wGrid%Dqv
	M_q = wGrid%Mqv
	IF (OpenTop == 2 .OR. OpenLateral == 2) THEN
		S_q = - wGrid%tau*wGrid%qv
	ELSE
		S_q = 0.
	END IF
CASE (2)
	D_q = wGrid%Dqc
	M_q = wGrid%Mqc
	IF (OpenTop == 2 .OR. OpenLateral == 2) THEN
		S_q = - wGrid%tau*wGrid%qc
	ELSE
		S_q = 0.
	END IF
CASE (3)
	D_q = wGrid%Dqr
	M_q = wGrid%Mqr
	IF (OpenTop == 2 .OR. OpenLateral == 2) THEN
		S_q = - wGrid%tau*wGrid%qr
	ELSE
		S_q = 0.
	END IF
CASE (4)
	D_q = wGrid%Dqi
	M_q = wGrid%Mqi
	IF (OpenTop == 2 .OR. OpenLateral == 2) THEN
		S_q = - wGrid%tau*wGrid%qi
	ELSE
		S_q = 0.
	END IF
CASE (5)
	D_q = wGrid%Dqs
	M_q = wGrid%Mqs
	IF (OpenTop == 2 .OR. OpenLateral == 2) THEN
		S_q = - wGrid%tau*wGrid%qs
	ELSE
		S_q = 0.
	END IF
CASE (6)
	D_q = wGrid%Dqg
	M_q = wGrid%Mqg
	IF (OpenTop == 2 .OR. OpenLateral == 2) THEN
		S_q = - wGrid%tau*wGrid%qg
	ELSE
		S_q = 0.
	END IF
CASE DEFAULT
	STOP "WRONG flag of calc_advection_w"
END SELECT

!$OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		tend_q(i,k) = A_q(i,k) + D_q(i,k) + M_q(i,k) + S_q(i,k)
	END DO
END DO
!$OMP END PARALLEL DO
!-------------------------------------------------
IF (ANY(ISNAN(tend_q(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITH tend_theta!!!"
!=================================================
END SUBROUTINE tendency_q
!=================================================

!=================================================
END MODULE sp_module_tendency
!=================================================

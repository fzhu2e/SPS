!=================================================
! The flux module of SPS-dynamic-integrate
!-------------------------------------------------
! Version: 0.11
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
REAL(kd), DIMENSION(ims:ime,kms:kme) :: P2uPx2_u = undef, P2uPz2_u = undef

REAL(kd), DIMENSION(ims:ime,kms:kme) :: A_u = undef, P_u = undef, D_u = undef
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL calc_advection_u(Main%u,A_u,uGrid,wGrid,piGrid,virGrid)

CALL ppx_u(Main%pi_1,Ppi_1Px_u)
CALL ppzeta_u(Main%pi_1,Ppi_1Pzeta_u)

CALL set_area_u
IF (ANY(uGrid%theta_M_0(imin:imax,kmin:kmax) == undef)) STOP "tend_u"
IF (ANY(Ppi_1Px_u(imin:imax,kmin:kmax) == undef)) STOP "tend_u"
IF (ANY(uGrid%G(imin:imax,kmin:kmax) == undef)) STOP "tend_u"
IF (ANY(Ppi_1Pzeta_u(imin:imax,kmin:kmax) == undef)) STOP "tend_u"
IF (ANY(Main%u(imin-1:imax+1,kmin-1:kmax+1) == undef)) STOP "tend_u"
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		P_u(i,k) = - Cp*uGrid%theta_M_0(i,k)*(Ppi_1Px_u(i,k) + uGrid%G(i,k)*Ppi_1Pzeta_u(i,k))
		!P_u(i,k) = - Cp*uGrid%theta_M(i,k)*(Ppi_1Px_u(i,k) + uGrid%G(i,k)*Ppi_1Pzeta_u(i,k))

		P2uPx2_u(i,k) = (Main%u(i+1,k) + Main%u(i-1,k) - 2*Main%u(i,k))/dx/dx
		P2uPz2_u(i,k) = (Main%u(i,k+1) + Main%u(i,k-1) - 2*Main%u(i,k))/dz/dz
		D_u(i,k) = Km*(P2uPx2_u(i,k) + P2uPz2_u(i,k))

		tend_u(i,k) = A_u(i,k) + P_u(i,k) + D_u(i,k)
	END DO
END DO
!OMP END PARALLEL DO
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
REAL(kd), DIMENSION(ims:ime,kms:kme) :: A_w = undef, B_w = undef, P_w = undef, D_w = undef
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: fa, fb, fc, fd, fe, ff
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: Ppi_1Pzeta_w = undef

REAL(kd), DIMENSION(ims:ime,kms:kme) :: P2wPx2_w = undef, P2wPz2_w = undef
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL calc_advection_w(Main%w,A_w,uGrid,wGrid,piGrid,virGrid)
CALL ppzeta_w(Main%pi_1,Ppi_1Pzeta_w)

CALL set_area_w
IF (ANY(wGrid%theta_M_1(imin:imax,kmin:kmax) == undef)) STOP "tend_w"
IF (ANY(wGrid%theta_M_0(imin:imax,kmin:kmax) == undef)) STOP "tend_w"
IF (ANY(wGrid%H(imin:imax) == undef)) STOP "tend_w"
IF (ANY(Ppi_1Pzeta_w(imin:imax,kmin:kmax) == undef)) STOP "tend_w"
IF (ANY(Main%w(imin-1:imax+1,kmin-1:kmax+1) == undef)) STOP "tend_w"
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		B_w(i,k) = g*wGrid%theta_M_1(i,k)/wGrid%theta_M_0(i,k)

		P_w(i,k) = - Cp*wGrid%theta_M_0(i,k)*wGrid%H(i)*Ppi_1Pzeta_w(i,k)
		!P_w(i,k) = - Cp*wGrid%theta_M(i,k)*wGrid%H(i)*Ppi_1Pz_w(i,k)

		P2wPx2_w(i,k) = (Main%w(i+1,k) + Main%w(i-1,k) - 2*Main%w(i,k))/dx/dx
		P2wPz2_w(i,k) = (Main%w(i,k+1) + Main%w(i,k-1) - 2*Main%w(i,k))/dz/dz
		D_w(i,k) = Km*(P2wPx2_w(i,k) + P2wPz2_w(i,k))

		tend_w(i,k) = A_w(i,k) + B_w(i,k) + P_w(i,k) + D_w(i,k)
	END DO
END DO
!OMP END PARALLEL DO
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
REAL(kd), DIMENSION(ims:ime,kms:kme) :: A_pi_1 = undef
REAL(kd), DIMENSION(ims:ime,kms:kme) :: Div_pi_1 = undef
!-------------------------------------------------
REAL(kd) :: temp
REAL(kd), DIMENSION(ims:ime,kms:kme) :: PuPx_pi = undef, PuPzeta_pi = undef, PwPzeta_pi = undef
INTEGER :: i, k
!=================================================
!-------------------------------------------------
CALL ppx_pi(Main%u,PuPx_pi)
CALL ppzeta_pi(wGrid%u,PuPzeta_pi)
CALL ppzeta_pi(Main%w,PwPzeta_pi)

CALL calc_advection_pi(Main%pi_1,A_pi_1,uGrid,wGrid,piGrid,virGrid)

CALL set_area_pi
IF (ANY(PuPx_pi(imin:imax,kmin:kmax) == undef)) STOP "tend_pi"
IF (ANY(piGrid%G(imin:imax,kmin:kmax) == undef)) STOP "tend_pi"
IF (ANY(PuPzeta_pi(imin:imax,kmin:kmax) == undef)) STOP "tend_pi"
IF (ANY(piGrid%H(imin:imax) == undef)) STOP "tend_pi"
IF (ANY(PwPzeta_pi(imin:imax,kmin:kmax) == undef)) STOP "tend_pi"
IF (ANY(piGrid%theta_M_0(imin:imax,kmin:kmax) == undef)) STOP "tend_pi"
!OMP PARALLEL DO PRIVATE(temp)
DO k = kmin, kmax
	DO i = imin, imax
		temp = PuPx_pi(i,k) + piGrid%G(i,k)*PuPzeta_pi(i,k) + piGrid%H(i)*PwPzeta_pi(i,k)
		Div_pi_1(i,k) = - cs*cs/Cp/piGrid%theta_M_0(i,k)*temp
		!Div_pi_1(i,k) = - Rd*piGrid%pi_0(i,k)/Cv*temp
		tend_pi_1(i,k) = A_pi_1(i,k) + Div_pi_1(i,k)
	END DO
END DO
!OMP END PARALLEL DO
!-------------------------------------------------
IF (ANY(ISNAN(tend_pi_1(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITH tend_pi_1!!!"
!=================================================
END SUBROUTINE tendency_pi
!=================================================

!=================================================
SUBROUTINE tendency_theta(flag,Main,tend_theta,uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
INTEGER, INTENT(IN) :: flag
TYPE(mainvar), INTENT(IN) :: Main
TYPE(grid), INTENT(IN) :: uGrid, wGrid, piGrid, virGrid
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: tend_theta
!=================================================
REAL(kd), DIMENSION(ims:ime,kms:kme) :: A_theta = undef, D_theta = undef, M_theta = undef
REAL(kd), DIMENSION(ims:ime,kms:kme) :: A_qv = undef, D_qv = undef, M_qv = undef
REAL(kd), DIMENSION(ims:ime,kms:kme) :: A_qc = undef, D_qc = undef, M_qc = undef
REAL(kd), DIMENSION(ims:ime,kms:kme) :: A_qr = undef, D_qr = undef, M_qr = undef
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: P2thetaPx2_w = undef, P2thetaPz2_w = undef
INTEGER :: i, k
!=================================================
SELECT CASE (flag)
CASE (0)
	CALL calc_advection_w(Main%theta,A_theta,uGrid,wGrid,piGrid,virGrid)
	CALL set_area_w
	IF (ANY(Main%theta(imin-1:imax+1,kmin-1:kmax+1) == undef)) STOP "tend_theta"
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			P2thetaPx2_w(i,k) = (Main%theta(i+1,k) + Main%theta(i-1,k) - 2*Main%theta(i,k))/dx/dx
			P2thetaPz2_w(i,k) = (Main%theta(i,k+1) + Main%theta(i,k-1) - 2*Main%theta(i,k))/dz/dz
			D_theta(i,k) = Kh*(P2thetaPx2_w(i,k) + P2thetaPz2_w(i,k))
	
			M_theta(i,k) = 0.
			tend_theta(i,k) = A_theta(i,k) + D_theta(i,k) + M_theta(i,k)
		END DO
	END DO
	!OMP END PARALLEL DO

CASE (1)
	CALL calc_advection_w(Main%qv,A_qv,uGrid,wGrid,piGrid,virGrid)
	CALL set_area_w
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			D_qv(i,k) = 0.
			M_qv(i,k) = 0.
			tend_theta(i,k) = A_qv(i,k) + D_qv(i,k) + M_qv(i,k)
		END DO
	END DO
	!OMP END PARALLEL DO

CASE (2)
	CALL calc_advection_w(Main%qc,A_qc,uGrid,wGrid,piGrid,virGrid)
	CALL set_area_w
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			D_qc(i,k) = 0.
			M_qc(i,k) = 0.
			tend_theta(i,k) = A_qc(i,k) + D_qc(i,k) + M_qc(i,k)
		END DO
	END DO
	!OMP END PARALLEL DO

CASE (3)
	CALL calc_advection_w(Main%qr,A_qr,uGrid,wGrid,piGrid,virGrid)
	CALL set_area_w
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			D_qr(i,k) = 0.
			M_qr(i,k) = 0.
			tend_theta(i,k) = A_qr(i,k) + D_qr(i,k) + M_qr(i,k)
		END DO
	END DO
	!OMP END PARALLEL DO

CASE DEFAULT
	WRITE(*,*) "Wrong flag of tendency_theta"
END SELECT
!-------------------------------------------------
IF (ANY(ISNAN(tend_theta(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITH tend_theta!!!"
!=================================================
END SUBROUTINE tendency_theta
!=================================================

!=================================================
END MODULE sp_module_tendency
!=================================================

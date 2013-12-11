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
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		P_u(i,k) = - Cp*uGrid%theta_M_0(i,k)*(Ppi_1Px_u(i,k) + uGrid%G(i,k)*Ppi_1Pzeta_u(i,k))

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
REAL(kd), DIMENSION(ims:ime,kms:kme) :: Ppi_1Pz_w = undef

REAL(kd), DIMENSION(ims:ime,kms:kme) :: P2wPx2_w = undef, P2wPz2_w = undef
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL calc_advection_w(Main%w,A_w,uGrid,wGrid,piGrid,virGrid)
CALL ppzeta_w(Main%pi_1,Ppi_1Pz_w)

CALL set_area_w
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		B_w(i,k) = g*wGrid%theta_M_1(i,k)/wGrid%theta_M_0(i,k)

		P_w(i,k) = - Cp*wGrid%theta_M_0(i,k)*wGrid%H(i)*Ppi_1Pz_w(i,k)

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
REAL(kd), DIMENSION(ims:ime,kms:kme) :: F_pi
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: urhotheta_u = undef, wrhotheta_w = undef
REAL(kd), DIMENSION(ims:ime,kms:kme) :: PurhothetaPx_pi = undef, PwrhothetaPz_pi = undef
INTEGER :: i, k
!=================================================
! 5.1 F_pi = - c^2/(rho_0*theta_0^2)*(PurhothetaPx + PwrhothetaPz)
!-------------------------------------------------
CALL set_area_u
CALL set_area_expand(expand)
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		urhotheta_u(i,k) = Main%u(i,k)*uGrid%rho_0(i,k)*uGrid%theta_M_0(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL set_area_w
CALL set_area_expand(expand)
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		wrhotheta_w(i,k) = Main%w(i,k)*wGrid%rho_0(i,k)*wGrid%theta_M_0(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL ppx_pi(urhotheta_u,PurhothetaPx_pi)
CALL ppzeta_pi(wrhotheta_w,PwrhothetaPz_pi)

CALL set_area_pi
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		F_pi(i,k) = - cs*cs/Cp/piGrid%rho_0(i,k)/piGrid%theta_M_0(i,k)/piGrid%theta_M_0(i,k)*(PurhothetaPx_pi(i,k) + PwrhothetaPz_pi(i,k))
		tend_pi_1(i,k) = F_pi(i,k)
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
REAL(kd), DIMENSION(ims:ime,kms:kme) :: A_theta = undef, D_theta = undef
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: P2thetaPx2_w = undef, P2thetaPz2_w = undef
INTEGER :: i, k
!=================================================
SELECT CASE (flag)
CASE (0)
	CALL calc_advection_w(Main%theta,A_theta,uGrid,wGrid,piGrid,virGrid)
CASE (1)
	CALL calc_advection_w(Main%qv,A_theta,uGrid,wGrid,piGrid,virGrid)
CASE (2)
	CALL calc_advection_w(Main%qc,A_theta,uGrid,wGrid,piGrid,virGrid)
CASE (3)
	CALL calc_advection_w(Main%qr,A_theta,uGrid,wGrid,piGrid,virGrid)
CASE DEFAULT
	WRITE(*,*) "Wrong flag of tendency_theta"
END SELECT

CALL set_area_w
IF (flag == 0) THEN
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			P2thetaPx2_w(i,k) = (Main%theta(i+1,k) + Main%theta(i-1,k) - 2*Main%theta(i,k))/dx/dx
			P2thetaPz2_w(i,k) = (Main%theta(i,k+1) + Main%theta(i,k-1) - 2*Main%theta(i,k))/dz/dz
			D_theta(i,k) = Kh*(P2thetaPx2_w(i,k) + P2thetaPz2_w(i,k))
	
			tend_theta(i,k) = A_theta(i,k) + D_theta(i,k)
		END DO
	END DO
	!OMP END PARALLEL DO
ELSE
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			tend_theta(i,k) = A_theta(i,k)
		END DO
	END DO
	!OMP END PARALLEL DO
END IF

!-------------------------------------------------
IF (ANY(ISNAN(tend_theta(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITH tend_theta!!!"
!=================================================
END SUBROUTINE tendency_theta
!=================================================

!=================================================
END MODULE sp_module_tendency
!=================================================

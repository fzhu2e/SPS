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
USE sp_module_interpolate
USE sp_module_debug
IMPLICIT NONE
!=================================================
! Tendency term
REAL(kd), DIMENSION(ims:ime,kms:kme) ::    F_u,    F_w,    F_theta,    F_pi
REAL(kd), DIMENSION(ims:ime,kms:kme) :: tend_u, tend_w, tend_theta, tend_pi
!-------------------------------------------------
! Diffusion term
REAL(kd), DIMENSION(ims:ime,kms:kme) :: P2uPx2_u, P2uPz2_u
REAL(kd), DIMENSION(ims:ime,kms:kme) :: P2wPx2_w, P2wPz2_w
REAL(kd), DIMENSION(ims:ime,kms:kme) :: P2thetaPx2_w, P2thetaPz2_w
!-------------------------------------------------
! Conponents
!----------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: rhou_pi, rhouu_pi
REAL(kd), DIMENSION(ims:ime,kms:kme) :: rhow_v , rhowu_v
REAL(kd), DIMENSION(ims:ime,kms:kme) :: rhou_v , rhouw_v
REAL(kd), DIMENSION(ims:ime,kms:kme) :: rhow_pi, rhoww_pi
REAL(kd), DIMENSION(ims:ime,kms:kme) :: rhoutheta_v, rhowtheta_pi
REAL(kd), DIMENSION(ims:ime,kms:kme) :: urhotheta_u, wrhotheta_w
!----------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: PrhouPx_u, PrhouuPx_u, PrhowPz_u, PrhowuPz_u, Ppi_1Px_u
REAL(kd), DIMENSION(ims:ime,kms:kme) :: PrhouPx_w, PrhouwPx_w, PrhowPz_w, PrhowwPz_w, Ppi_1Pz_w
REAL(kd), DIMENSION(ims:ime,kms:kme) :: PrhouthetaPx_w, PrhowthetaPz_w
REAL(kd), DIMENSION(ims:ime,kms:kme) :: PurhothetaPx_pi, PwrhothetaPz_pi
!-------------------------------------------------
INTEGER, PARAMETER :: expand = 1
!=================================================
CONTAINS
!=================================================
SUBROUTINE tendency_u(u,rho_0,pi_1,F_u,tend_u)
IMPLICIT NONE
!=================================================
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0    ! density
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_1
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: F_u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: tend_u
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: uPuPx_u, wPuPz_u
REAL(kd), DIMENSION(ims:ime,kms:kme) :: fa, fb, fc, fd, fe, ff
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! 1. F_u = - u p.u/p.x - w p.u/p.z + fv   (fv = 0.)
! 1.1. - u p.u/p.x = - 1/rho (p.rhouu/p.x - u p.rhou/p.x)
! 1.2. - w p.u/p.z = - 1/rho (p.rhouw/p.z - u p.rhow/p.z)
!-------------------------------------------------
! pi-grid - Middle-vars can be calculated on boundaries.
!-------------------------------------------------
CALL set_area_pi
CALL set_area_expand(expand)

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		rhou_pi(i,k) = rho_0(i,k)*u_pi(i,k)
	END DO
END DO
!OMP END PARALLEL DO

SELECT CASE (AdvectionScheme)
CASE (5)
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			fa(i,k) = u(i-1,k) + u(i,k)
			fb(i,k) = u(i-2,k) + u(i+1,k)
			fc(i,k) = u(i-3,k) + u(i+2,k)
			rhouu_pi(i,k) = rho_0(i,k)*u_pi(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
			fd(i,k) = u(i,k) - u(i-1,k)
			fe(i,k) = u(i+1,k) - u(i-2,k)
			ff(i,k) = u(i+2,k) - u(i-3,k)
			rhouu_pi(i,k) = rhouu_pi(i,k) - ABS(u_pi(i,k))/60.*(10*fd(i,k) - 5*fe(i,k) + ff(i,k))
		END DO
	END DO
	!OMP END PARALLEL DO
CASE DEFAULT
	STOP "Wrong advection scheme!!!"
END SELECT
!-------------------------------------------------
! v(virtual)-grid
!-------------------------------------------------
CALL set_area_v
CALL set_area_expand(expand)

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		rhow_v(i,k) = rho_0_v(i,k)*w_v(i,k)
	END DO
END DO
!OMP END PARALLEL DO

SELECT CASE (AdvectionScheme)
CASE (5)
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			fa(i,k) = u(i,k) + u(i,k-1)
			fb(i,k) = u(i,k+1) + u(i,k-2)
			fc(i,k) = u(i,k+2) + u(i,k-3)
			rhowu_v(i,k) = rho_0_v(i,k)*w_v(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
			fd(i,k) = u(i,k) - u(i,k-1)
			fe(i,k) = u(i,k+1) - u(i,k-2)
			ff(i,k) = u(i,k+2) - u(i,k-3)
			rhowu_v(i,k) = rhowu_v(i,k) - ABS(w_v(i,k))/60.*(10*fd(i,k) - 5*fe(i,k) + ff(i,k))
		END DO
	END DO
	!OMP END PARALLEL DO
CASE DEFAULT
	STOP "Wrong advection scheme!!!"
END SELECT
!-------------------------------------------------
! u-grid
!-------------------------------------------------
CALL set_area_u
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		PrhouPx_u(i,k) = (rhou_pi(i+1,k) - rhou_pi(i,k))/dx
		PrhouuPx_u(i,k) = (rhouu_pi(i+1,k) - rhouu_pi(i,k))/dx
		PrhowPz_u(i,k) = (rhow_v(i,k+1) - rhow_v(i,k))/dz
		PrhowuPz_u(i,k) = (rhowu_v(i,k+1) - rhowu_v(i,k))/dz
		
		uPuPx_u(i,k) = 1./rho_0_u(i,k)*(PrhouuPx_u(i,k) - u(i,k)*PrhouPx_u(i,k))
		wPuPz_u(i,k) = 1./rho_0_u(i,k)*(PrhowuPz_u(i,k) - u(i,k)*PrhowPz_u(i,k))
	
		Ppi_1Px_u(i,k) = (pi_1(i + 1,k) - pi_1(i,k))/dx
	
		F_u(i,k) = - uPuPx_u(i,k) - wPuPz_u(i,k)
		tend_u(i,k) = F_u(i,k) - Cp*theta_0_u(i,k)*Ppi_1Px_u(i,k)
	END DO
END DO
!OMP END PARALLEL DO

IF (RunCase == 1 .OR. RunCase == 2) THEN
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			P2uPx2_u(i,k) = (u(i+1,k) + u(i-1,k) - 2*u(i,k))/dx/dx
			P2uPz2_u(i,k) = (u(i,k+1) + u(i,k-1) - 2*u(i,k))/dz/dz
			
			tend_u(i,k) = tend_u(i,k) + Km*(P2uPx2_u(i,k) + P2uPz2_u(i,k)) ! Add diffusion term.
		END DO
	END DO
	!OMP END PARALLEL DO
END IF
	
!-------------------------------------------------
IF (ANY(ISNAN(F_u(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITH F_u!!!"
!=================================================
END SUBROUTINE tendency_u
!=================================================

!=================================================
SUBROUTINE tendency_w(w,rho_0,theta_1,theta_0,pi_1,F_w,tend_w)
IMPLICIT NONE
!=================================================
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0    ! density
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_1
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_0
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_1
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: F_w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: tend_w
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: uPwPx_w, wPwPz_w
REAL(kd), DIMENSION(ims:ime,kms:kme) :: fa, fb, fc, fd, fe, ff
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! 2. F_w = - u p.w/p.x - w p.w/p.z + g(theta_1/theta_0)
! 2.1. - u p.w/p.x = - 1/rho (p.rhouw/p.x - w p.rhou/p.x)
! 2.2. - w p.w/p.z = - 1/rho (p.rhoww/p.z - w p.rhow/p.z)
!-------------------------------------------------
! v(virtual)-grid
!-------------------------------------------------
CALL set_area_v
CALL set_area_expand(expand)

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		rhou_v(i,k) = rho_0_v(i,k)*u_v(i,k)
	END DO
END DO
!OMP END PARALLEL DO

SELECT CASE (AdvectionScheme)
CASE (5)
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			fa(i,k) = w(i,k) + w(i+1,k)
			fb(i,k) = w(i-1,k) + w(i+2,k)
			fc(i,k) = w(i-2,k) + w(i+3,k)
			rhouw_v(i,k) = rho_0_v(i,k)*u_v(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
			fd(i,k) = w(i+1,k) - w(i,k)
			fe(i,k) = w(i+2,k) - w(i-1,k)
			ff(i,k) = w(i+3,k) - w(i-2,k)
			rhouw_v(i,k) = rhouw_v(i,k) - ABS(u_v(i,k))/60.*(10*fd(i,k) - 5*fe(i,k) + ff(i,k))
		END DO
	END DO
	!OMP END PARALLEL DO
CASE DEFAULT
	STOP "Wrong advection scheme!!!"
END SELECT
!-------------------------------------------------
! pi-grid
!-------------------------------------------------
CALL set_area_pi
CALL set_area_expand(expand)

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		rhow_pi(i,k) = rho_0(i,k)*w_pi(i,k)
	END DO
END DO
!OMP END PARALLEL DO

SELECT CASE (AdvectionScheme)
CASE (5)
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			fa(i,k) = w(i,k) + w(i,k+1)
			fb(i,k) = w(i,k-1) + w(i,k+2)
			fc(i,k) = w(i,k-2) + w(i,k+3)
			rhoww_pi(i,k) = rho_0(i,k)*w_pi(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
	
			fd(i,k) = w(i,k+1) - w(i,k)
			fe(i,k) = w(i,k+2) - w(i,k-1)
			ff(i,k) = w(i,k+3) - w(i,k-2)
			rhoww_pi(i,k) = rhoww_pi(i,k) - ABS(w_pi(i,k))/60.*(10*fd(i,k) - 5*fe(i,k) + ff(i,k))
		END DO
	END DO
	!OMP END PARALLEL DO
CASE DEFAULT
	STOP "Wrong advection scheme!!!"
END SELECT
!-------------------------------------------------
! w-grid 
!-------------------------------------------------
CALL set_area_w
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		PrhouPx_w(i,k) = (rhou_v(i,k) - rhou_v(i-1,k))/dx
		PrhouwPx_w(i,k) = (rhouw_v(i,k) - rhouw_v(i-1,k))/dx
		PrhowPz_w(i,k) = (rhow_pi(i,k) - rhow_pi(i,k-1))/dz
		PrhowwPz_w(i,k) = (rhoww_pi(i,k) - rhoww_pi(i,k-1))/dz

		uPwPx_w(i,k) = 1./rho_0_w(i,k)*(PrhouwPx_w(i,k) - w(i,k)*PrhouPx_w(i,k))
		wPwPz_w(i,k) = 1./rho_0_w(i,k)*(PrhowwPz_w(i,k) - w(i,k)*PrhowPz_w(i,k))

		Ppi_1Pz_w(i,k) = (pi_1(i,k) - pi_1(i,k - 1))/dz

		F_w(i,k) = - uPwPx_w(i,k) - wPwPz_w(i,k) + g*theta_1(i,k)/theta_0(i,k)
		tend_w(i,k) = F_w(i,k) - Cp*theta_0(i,k)*Ppi_1Pz_w(i,k)
	END DO
END DO
!OMP END PARALLEL DO

IF (RunCase == 1 .OR. RunCase == 2) THEN
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			P2wPx2_w(i,k) = (w(i+1,k) + w(i-1,k) - 2*w(i,k))/dx/dx
			P2wPz2_w(i,k) = (w(i,k+1) + w(i,k-1) - 2*w(i,k))/dz/dz
			
			tend_w(i,k) = tend_w(i,k) + Km*(P2wPx2_w(i,k) + P2wPz2_w(i,k)) ! Add diffusion term.
		END DO
	END DO
	!OMP END PARALLEL DO
END IF
!-------------------------------------------------
IF (ANY(ISNAN(F_w(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITHT F_w!!!"
!=================================================
END SUBROUTINE tendency_w
!=================================================

!=================================================
SUBROUTINE tendency_theta(u,w,rho_0,theta,F_theta,tend_theta)
IMPLICIT NONE
!=================================================
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0    ! density
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: F_theta
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: tend_theta
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: uPthetaPx_w, wPthetaPz_w
REAL(kd), DIMENSION(ims:ime,kms:kme) :: fa, fb, fc, fd, fe, ff
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! 3. F_theta = - u p.theta/p.x - w p.theta/p.z
! 3.1. - u p.theta/p.x = - 1/rho (p.rhoutheta/p.x - theta p.rhou/p.x)
! 3.2. - w p.w/p.z = - 1/rho (p.rhothetaw/p.z - theta p.rhow/p.z)
!-------------------------------------------------
! v(virtual)-grid
!-------------------------------------------------
CALL set_area_v
CALL set_area_expand(expand)

SELECT CASE (AdvectionScheme)
CASE (5)
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			fa(i,k) = theta(i,k) + theta(i+1,k)
			fb(i,k) = theta(i-1,k) + theta(i+2,k)
			fc(i,k) = theta(i-2,k) + theta(i+3,k)
			rhoutheta_v(i,k) = rho_0_v(i,k)*u_v(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
			fd(i,k) = theta(i+1,k) - theta(i,k)
			fe(i,k) = theta(i+2,k) - theta(i-1,k)
			ff(i,k) = theta(i+3,k) - theta(i-2,k)
			rhoutheta_v(i,k) = rhoutheta_v(i,k) - ABS(u_v(i,k))/60.*(10*fd(i,k) - 5*fe(i,k) + ff(i,k))
		END DO
	END DO
	!OMP END PARALLEL DO
CASE DEFAULT
	STOP "Wrong advection scheme!!!"
END SELECT
!-------------------------------------------------
! pi-grid
!-------------------------------------------------
CALL set_area_pi
CALL set_area_expand(expand)

SELECT CASE (AdvectionScheme)
CASE (5)
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			fa(i,k) = theta(i,k+1) + theta(i,k)
			fb(i,k) = theta(i,k+2) + theta(i,k-1)
			fc(i,k) = theta(i,k+3) + theta(i,k-2)
			rhowtheta_pi(i,k) = rho_0(i,k)*w_pi(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
			fd(i,k) = theta(i,k+1) - theta(i,k)
			fe(i,k) = theta(i,k+2) - theta(i,k-1)
			ff(i,k) = theta(i,k+3) - theta(i,k-2)
			rhowtheta_pi(i,k) = rhowtheta_pi(i,k) - ABS(w_hat_pi(i,k))/60.*(10*fd(i,k) - 5*fe(i,k) + ff(i,k))
		END DO
	END DO
	!OMP END PARALLEL DO
CASE DEFAULT
	STOP "Wrong advection scheme!!!"
END SELECT
!-------------------------------------------------
! w-grid - Theta on kts and kte+1 should also be updated.
!-------------------------------------------------
CALL set_area_w

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		PrhouthetaPx_w(i,k) = (rhoutheta_v(i,k) - rhoutheta_v(i-1,k))/dx
		PrhowthetaPz_w(i,k) = (rhowtheta_pi(i,k) - rhowtheta_pi(i,k-1))/dz
	END DO
END DO
!OMP END PARALLEL DO

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		uPthetaPx_w(i,k) = 1./rho_0_w(i,k)*(PrhouthetaPx_w(i,k) - theta(i,k)*PrhouPx_w(i,k))
		wPthetaPz_w(i,k) = 1./rho_0_w(i,k)*(PrhowthetaPz_w(i,k) - theta(i,k)*PrhowPz_w(i,k))
	
		F_theta(i,k) = - uPthetaPx_w(i,k) - wPthetaPz_w(i,k)
		tend_theta(i,k) = F_theta(i,k)
	END DO
END DO
!OMP END PARALLEL DO
	
IF (RunCase == 1 .OR. RunCase == 2) THEN
	!OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			P2thetaPx2_w(i,k) = (theta(i+1,k) + theta(i-1,k) - 2*theta(i,k))/dx/dx
			P2thetaPz2_w(i,k) = (theta(i,k+1) + theta(i,k-1) - 2*theta(i,k))/dz/dz
			
			tend_theta(i,k) = F_theta(i,k) + Kh*(P2thetaPx2_w(i,k) + P2thetaPz2_w(i,k)) ! Add diffusion term.
		END DO
	END DO
	!OMP END PARALLEL DO
END IF
!-------------------------------------------------
IF (ANY(ISNAN(F_theta(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITH F_theta!!!"
!=================================================
END SUBROUTINE tendency_theta
!=================================================

!=================================================
SUBROUTINE tendency_pi(u,w,pi_0,rho_0,theta_0,F_pi,tend_pi)
IMPLICIT NONE
!=================================================
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_0
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0    ! density
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_0
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: F_pi
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: tend_pi
INTEGER :: i, k
!=================================================
! 5.1 F_pi = - c^2/(rho_0*theta_0^2)*(PurhothetaPx + PwrhothetaPz)
!-------------------------------------------------
CALL set_area_u
CALL set_area_expand(expand)

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		urhotheta_u(i,k) = u(i,k)*rho_0_u(i,k)*theta_0_u(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL set_area_w
CALL set_area_expand(expand)

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		wrhotheta_w(i,k) = w(i,k)*rho_0_w(i,k)*theta_0(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL set_area_pi
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		PurhothetaPx_pi(i,k) = (urhotheta_u(i,k) - urhotheta_u(i - 1,k))/dx
		PwrhothetaPz_pi(i,k) = (wrhotheta_w(i,k + 1) - wrhotheta_w(i,k))/dz
		
		F_pi(i,k) = - cs*cs/Cp/rho_0(i,k)/theta_0_pi(i,k)/theta_0_pi(i,k)*(PurhothetaPx_pi(i,k) + PwrhothetaPz_pi(i,k))
		tend_pi(i,k) = F_pi(i,k)
	END DO
END DO
!OMP END PARALLEL DO
!-------------------------------------------------
IF (ANY(ISNAN(F_pi(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITH F_theta!!!"
!=================================================
END SUBROUTINE tendency_pi
!=================================================

!=================================================
END MODULE sp_module_tendency
!=================================================

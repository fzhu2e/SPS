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
REAL(preci), DIMENSION(ims:ime,kms:kme) :: F_u
REAL(preci), DIMENSION(ims:ime,kms:kme) :: tend_u

REAL(preci), DIMENSION(ims:ime,kms:kme) :: rhou_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: rhouu_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: PrhouPx_u
REAL(preci), DIMENSION(ims:ime,kms:kme) :: PrhouuPx_u

REAL(preci), DIMENSION(ims:ime,kms:kme) :: rhow_v
REAL(preci), DIMENSION(ims:ime,kms:kme) :: rhouw_v
REAL(preci), DIMENSION(ims:ime,kms:kme) :: PrhowPz_u
REAL(preci), DIMENSION(ims:ime,kms:kme) :: PrhouwPz_u

REAL(preci), DIMENSION(ims:ime,kms:kme) :: P2uPx2_u
REAL(preci), DIMENSION(ims:ime,kms:kme) :: P2uPz2_u
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: F_w
REAL(preci), DIMENSION(ims:ime,kms:kme) :: tend_w

REAL(preci), DIMENSION(ims:ime,kms:kme) :: rhou_v
!REAL(preci), DIMENSION(ims:ime,kms:kme) :: rhouw_v
REAL(preci), DIMENSION(ims:ime,kms:kme) :: PrhouPx_w
REAL(preci), DIMENSION(ims:ime,kms:kme) :: PrhouwPx_w

REAL(preci), DIMENSION(ims:ime,kms:kme) :: rhow_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: rhoww_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: PrhowPz_w
REAL(preci), DIMENSION(ims:ime,kms:kme) :: PrhowwPz_w

REAL(preci), DIMENSION(ims:ime,kms:kme) :: P2wPx2_w
REAL(preci), DIMENSION(ims:ime,kms:kme) :: P2wPz2_w
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: F_theta
REAL(preci), DIMENSION(ims:ime,kms:kme) :: tend_theta

!REAL(preci), DIMENSION(ims:ime,kms:kme) :: rhou_v
REAL(preci), DIMENSION(ims:ime,kms:kme) :: rhoutheta_v
!REAL(preci), DIMENSION(ims:ime,kms:kme) :: PrhouPx_w
REAL(preci), DIMENSION(ims:ime,kms:kme) :: PrhouthetaPx_w

!REAL(preci), DIMENSION(ims:ime,kms:kme) :: rhow_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: rhowtheta_pi
!REAL(preci), DIMENSION(ims:ime,kms:kme) :: PrhowPz_w
REAL(preci), DIMENSION(ims:ime,kms:kme) :: PrhowthetaPz_w

REAL(preci), DIMENSION(ims:ime,kms:kme) :: P2thetaPx2_w
REAL(preci), DIMENSION(ims:ime,kms:kme) :: P2thetaPz2_w
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: Ppi_1Px_u
REAL(preci), DIMENSION(ims:ime,kms:kme) :: Ppi_1Pz_w
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: F_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: tend_pi

REAL(preci), DIMENSION(ims:ime,kms:kme) :: urhotheta_u
REAL(preci), DIMENSION(ims:ime,kms:kme) :: wrhotheta_w
REAL(preci), DIMENSION(ims:ime,kms:kme) :: PurhothetaPx_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: PwrhothetaPz_pi
!=================================================
CONTAINS
!=================================================
SUBROUTINE tendency_u(u,rho_0,pi_1,F_u,tend_u)
IMPLICIT NONE
!=================================================
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: u
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0    ! density
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_1
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: F_u
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: tend_u
!-------------------------------------------------
REAL(preci) :: fa, fb, fc, fd, fe, ff
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! 1. F_u = - u p.u/p.x - w p.u/p.z + fv   (fv = 0.)
! 1.1. - u p.u/p.x = - 1/rho (p.rhouu/p.x - u p.rhou/p.x)
! 1.2. - w p.u/p.z = - 1/rho (p.rhouw/p.z - u p.rhow/p.z)
!-------------------------------------------------
! pi-grid - Middle-vars can be calculated on boundaries.
!-------------------------------------------------
DO i = its, ite + 1
	DO k = kts - 1, kte + 1
		rhou_pi(i,k) = rho_0(i,k)*u_pi(i,k)
		IF (rho_0(i,k) == undef .OR. u_pi(i,k) == undef) STOP "rhou_pi is WRONG!!!"
		SELECT CASE (AdvectionScheme)
		CASE (2)
			rhouu_pi(i,k) = rhou_pi(i,k)*u_pi(i,k)
			IF (rhou_pi(i,k) == undef .OR. u_pi(i,k) == undef) STOP "rhouu_pi is WRONG!!!"
		CASE (3)
			fa = u(i-1,k) + u(i,k)
			fb = u(i-2,k) + u(i+1,k)
			fc = u(i,k) - u(i-1,k)
			fd = u(i+1,k) - u(i-2,k)
			rhouu_pi(i,k) = rho_0(i,k)*(u_pi(i,k)/12.*(7*fa - fb) - ABS(u_pi(i,k))/12.*(3*fc - fd))
			IF(u(i-2,k) == undef .OR. u(i-1,k) == undef .OR. u(i,k) == undef .OR. u(i+1,k) == undef) STOP "rhouu_pi is WRONG!!!"
		CASE (4)
			fa = u(i-1,k) + u(i,k)
			fb = u(i-2,k) + u(i+1,k)
			rhouu_pi(i,k) = rho_0(i,k)*(u_pi(i,k)/12.*(7*fa - fb))
			IF(u(i-2,k) == undef .OR. u(i-1,k) == undef .OR. u(i,k) == undef .OR. u(i+1,k) == undef) STOP "rhouu_pi is WRONG!!!"
		CASE (5)
			fa = u(i-1,k) + u(i,k)
			fb = u(i-2,k) + u(i+1,k)
			fc = u(i-3,k) + u(i+2,k)
			rhouu_pi(i,k) = rho_0(i,k)*u_pi(i,k)/60.*(37*fa - 8*fb + fc)
			fd = u(i,k) - u(i-1,k)
			fe = u(i+1,k) - u(i-2,k)
			ff = u(i+2,k) - u(i-3,k)
			rhouu_pi(i,k) = rhouu_pi(i,k) - ABS(u_pi(i,k))/60.*(10*fd - 5*fe + ff)
			IF(u(i-2,k) == undef .OR. u(i-1,k) == undef .OR. u(i,k) == undef .OR. u(i+1,k) == undef .OR. u(i-3,k) == undef .OR. u(i+2,k) == undef) STOP "rhouu_pi is WRONG!!!"
		CASE (6)
			fa = u(i-1,k) + u(i,k)
			fb = u(i-2,k) + u(i+1,k)
			fc = u(i-3,k) + u(i+2,k)
			rhouu_pi(i,k) = rho_0(i,k)*u_pi(i,k)/60.*(37*fa - 8*fb + fc)
			IF(u(i-2,k) == undef .OR. u(i-1,k) == undef .OR. u(i,k) == undef .OR. u(i+1,k) == undef .OR. u(i-3,k) == undef .OR. u(i+2,k) == undef) STOP "rhouu_pi is WRONG!!!"
		CASE DEFAULT
			STOP "Wrong advection scheme!!!"
		END SELECT
	END DO
END DO
!-------------------------------------------------
! v(virtual)-grid
!-------------------------------------------------
!CALL debug_ascii_output(u)
!CALL debug_ascii_output(rho_0)
!CALL debug_ascii_output(rho_0_v)
DO i = its, ite
	DO k = kts, kte + 1
		rhow_v(i,k) = rho_0_v(i,k)*w_v(i,k)
		IF (rho_0_v(i,k) == undef .OR. w_v(i,k) == undef .OR. u_v(i,k) == undef) STOP "rhow_v or rhouw_v is WRONG!!!"
		SELECT CASE (AdvectionScheme)
		CASE (2)
			rhouw_v(i,k) = rhow_v(i,k)*u_v(i,k)
			IF (rhow_v(i,k) == undef .OR. u_v(i,k)== undef) STOP "rhouw_v is WRONG!!!"
		CASE (3)
			fa = u(i,k) + u(i,k-1)
			fb = u(i,k+1) + u(i,k-2)
			fc = u(i,k) - u(i,k-1)
			fd = u(i,k+1) - u(i,k-2)
			rhouw_v(i,k) = rho_0_v(i,k)*(w_v(i,k)/12.*(7*fa - fb) - ABS(w_v(i,k))/12.*(3*fc - fd))
			IF(u(i,k+1) == undef .OR. u(i,k) == undef .OR. u(i,k-1) == undef .OR. u(i,k-2) == undef) STOP "rhouw_v is WRONG!!!"
		CASE (4)
			fa = u(i,k) + u(i,k-1)
			fb = u(i,k+1) + u(i,k-2)
			rhouw_v(i,k) = rho_0_v(i,k)*(w_v(i,k)/12.*(7*fa - fb))
			IF(u(i,k+1) == undef .OR. u(i,k) == undef .OR. u(i,k-1) == undef .OR. u(i,k-2) == undef) STOP "rhouw_v is WRONG!!!"
		CASE (5)
			fa = u(i,k) + u(i,k-1)
			fb = u(i,k+1) + u(i,k-2)
			fc = u(i,k+2) + u(i,k-3)
			rhouw_v(i,k) = rho_0_v(i,k)*w_v(i,k)/60.*(37*fa - 8*fb + fc)
			fd = u(i,k) - u(i,k-1)
			fe = u(i,k+1) - u(i,k-2)
			ff = u(i,k+2) - u(i,k-3)
			rhouw_v(i,k) = rhouw_v(i,k) - ABS(w_v(i,k))/60.*(10*fd - 5*fe + ff)
			IF(u(i,k+1) == undef .OR. u(i,k) == undef .OR. u(i,k-1) == undef .OR. u(i,k-2) == undef .OR. u(i,k+2) == undef .OR. u(i,k-3) == undef) STOP "rhouw_v is WRONG!!!"
		CASE (6)
			fa = u(i,k) + u(i,k-1)
			fb = u(i,k+1) + u(i,k-2)
			fc = u(i,k+2) + u(i,k-3)
			rhouw_v(i,k) = rho_0_v(i,k)*w_v(i,k)/60.*(37*fa - 8*fb + fc)
			IF(u(i,k+1) == undef .OR. u(i,k) == undef .OR. u(i,k-1) == undef .OR. u(i,k-2) == undef .OR. u(i,k+2) == undef .OR. u(i,k-3) == undef) STOP "rhouw_v is WRONG!!!"
		CASE DEFAULT
			STOP "Wrong advection scheme!!!"
		END SELECT
	END DO
END DO
!rhow_v(its:ite,kts) = 0.      ! vertical boundary condition
!rhow_v(its:ite,kte+1) = 0.  ! vertical boundary condition
!rhouw_v(its:ite,kts) = 0.      ! vertical boundary condition
!rhouw_v(its:ite,kte+1) = 0.  ! vertical boundary condition
!rhow_v(its,kts:kte+1) = rhow_v(its+1,kts:kte+1)      ! lateral boundary condition
!rhow_v(ite,kts:kte+1) = rhow_v(ite-1,kts:kte+1)      ! lateral boundary condition
!rhouw_v(its,kts:kte+1) = rhouw_v(its+1,kts:kte+1)      ! lateral boundary condition
!rhouw_v(ite,kts:kte+1) = rhouw_v(ite-1,kts:kte+1)      ! lateral boundary condition
!CALL debug_ascii_output(w_v)
!CALL debug_test_boundary(rhow_v)
!CALL debug_test_boundary(rhouw_v)
!-------------------------------------------------
! u-grid
!-------------------------------------------------
DO i = its + 1, ite - 1
	DO k = kts, kte
		PrhouPx_u(i,k) = (rhou_pi(i+1,k) - rhou_pi(i,k))/dx
		PrhouuPx_u(i,k) = (rhouu_pi(i+1,k) - rhouu_pi(i,k))/dx
		PrhowPz_u(i,k) = (rhow_v(i,k+1) - rhow_v(i,k))/dz
		PrhouwPz_u(i,k) = (rhouw_v(i,k+1) - rhouw_v(i,k))/dz
		F_u(i,k) = - 1./rho_0_u(i,k)*(PrhouuPx_u(i,k) - u(i,k)*PrhouPx_u(i,k) + PrhouwPz_u(i,k) - u(i,k)*PrhowPz_u(i,k))
		!IF (ISNAN(F_u(i,k))) THEN
			!WRITE(*,*) u(i,k)
			!WRITE(*,*) rho_0_u(i,k)
			!WRITE(*,*) PrhouPx_u(i,k)
			!WRITE(*,*) PrhouuPx_u(i,k)
			!WRITE(*,*) PrhowPz_u(i,k)
			!WRITE(*,*) rhow_v(i,k+1)
			!WRITE(*,*) PrhouwPz_u(i,k)
			!STOP "WRONG F_u"
		!END IF
		Ppi_1Px_u(i,k) = (pi_1(i + 1,k) - pi_1(i,k))/dx
		tend_u(i,k) = F_u(i,k) - Cp*theta_0_u(i,k)*Ppi_1Px_u(i,k)
		
		P2uPx2_u(i,k) = (u(i+1,k) + u(i-1,k) - 2*u(i,k))/dx/dx
		P2uPz2_u(i,k) = (u(i,k+1) + u(i,k-1) - 2*u(i,k))/dz/dz
		IF (u(i+1,k) == undef .OR. u(i-1,k) == undef) STOP "P2uPx2_u is WRONG!!!"
		IF (u(i,k+1) == undef .OR. u(i,k-1) == undef) STOP "P2uPx2_u is WRONG!!!"
		
		tend_u(i,k) = tend_u(i,k) + Km*(P2uPx2_u(i,k) + P2uPz2_u(i,k)) ! Add diffusion term.
		
		IF (pi_1(i+1,k) == undef .OR. pi_1(i,k) == undef) STOP "Ppi_1Px_u is WRONG!!!"
		IF (rhou_pi(i+1,k) == undef .OR. rhou_pi(i,k) == undef) STOP "PrhouPx_u is WRONG!!!"
		IF (rhouu_pi(i+1,k) == undef .OR. rhouu_pi(i,k) == undef) STOP "PrhouuPx_u is WRONG!!!"
		IF (rhow_v(i,k+1) == undef .OR. rhow_v(i,k) == undef) STOP "PrhowPz_u is WRONG!!!"
		IF (rhouw_v(i,k+1) == undef .OR. rhouw_v(i,k) == undef) STOP "PrhouwPz_u is WRONG!!!"
	END DO
END DO
!-------------------------------------------------
IF (ANY(ISNAN(F_u(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITH F_u!!!"
!=================================================
END SUBROUTINE tendency_u
!=================================================

!=================================================
SUBROUTINE tendency_w(w,rho_0,theta_1,theta_0,pi_1,F_w,tend_w)
IMPLICIT NONE
!=================================================
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: w
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0    ! density
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_1
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_0
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_1
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: F_w
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: tend_w
!-------------------------------------------------
REAL(preci) :: fa, fb, fc, fd, fe, ff
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! 2. F_w = - u p.w/p.x - w p.w/p.z + g(theta_1/theta_0)
! 2.1. - u p.w/p.x = - 1/rho (p.rhouw/p.x - w p.rhou/p.x)
! 2.2. - w p.w/p.z = - 1/rho (p.rhoww/p.z - w p.rhow/p.z)
!-------------------------------------------------
! v(virtual)-grid
!-------------------------------------------------
DO i = its, ite
	DO k = kts, kte + 1
		rhou_v(i,k) = rho_0_v(i,k)*u_v(i,k)
		IF (rho_0_v(i,k) == undef .OR. u_v(i,k) == undef) STOP "rhou_v is WRONG!!!"
		SELECT CASE (AdvectionScheme)
		CASE (2)
			! rhouw_v have been calculated
		CASE (3)
			fa = w(i,k) + w(i+1,k)
			fb = w(i-1,k) + w(i+2,k)
			fc = w(i+1,k) - w(i,k)
			fd = w(i+2,k) - w(i-1,k)
			rhouw_v(i,k) = rho_0_v(i,k)*(u_v(i,k)/12.*(7*fa - fb) - ABS(u_v(i,k))/12.*(3*fc - fd))
			IF(w(i-1,k) == undef .OR. w(i,k) == undef .OR. w(i+1,k) == undef .OR. w(i+2,k) == undef) STOP "rhouw_v is WRONG!!!"
		CASE (4)
			fa = w(i,k) + w(i+1,k)
			fb = w(i-1,k) + w(i+2,k)
			rhouw_v(i,k) = rho_0_v(i,k)*(u_v(i,k)/12.*(7*fa - fb))
			IF(w(i-1,k) == undef .OR. w(i,k) == undef .OR. w(i+1,k) == undef .OR. w(i+2,k) == undef) STOP "rhouw_v is WRONG!!!"
		CASE (5)
			fa = w(i,k) + w(i+1,k)
			fb = w(i-1,k) + w(i+2,k)
			fc = w(i-2,k) + w(i+3,k)
			rhouw_v(i,k) = rho_0_v(i,k)*u_v(i,k)/60.*(37*fa - 8*fb + fc)
			fd = w(i+1,k) - w(i,k)
			fe = w(i+2,k) - w(i-1,k)
			ff = w(i+3,k) - w(i-2,k)
			rhouw_v(i,k) = rhouw_v(i,k) - ABS(u_v(i,k))/60.*(10*fd - 5*fe + ff)
			IF(w(i-1,k) == undef .OR. w(i,k) == undef .OR. w(i+1,k) == undef .OR. w(i+2,k) == undef .OR. w(i-2,k) == undef .OR. w(i+3,k) == undef) STOP "rhouw_v is WRONG!!!"
		CASE (6)
			fa = w(i,k) + w(i+1,k)
			fb = w(i-1,k) + w(i+2,k)
			fc = w(i-2,k) + w(i+3,k)
			rhouw_v(i,k) = rho_0_v(i,k)*u_v(i,k)/60.*(37*fa - 8*fb + fc)
			IF(w(i-1,k) == undef .OR. w(i,k) == undef .OR. w(i+1,k) == undef .OR. w(i+2,k) == undef .OR. w(i-2,k) == undef .OR. w(i+3,k) == undef) STOP "rhouw_v is WRONG!!!"
			
		CASE DEFAULT
			STOP "Wrong advection scheme!!!"
		END SELECT
	END DO
END DO
!rhou_v(its:ite,kts) = rhou_v(its:ite,kts+1)      ! vertical boundary condition
!rhou_v(its:ite,kte+1) = rhou_v(its:ite,kte)      ! vertical boundary condition
!rhou_v(its,kts:kte+1) = rhou_v(its+1,kts:kte+1)  ! lateral boundary condition
!rhou_v(ite,kts:kte+1) = rhou_v(ite-1,kts:kte+1)  ! lateral boundary condition
!-------------------------------------------------
! pi-grid
!-------------------------------------------------
DO i = its, ite + 1
	DO k = kts - 1, kte + 1
		rhow_pi(i,k) = rho_0(i,k)*w_pi(i,k)
		IF (rho_0(i,k) == undef .OR. w_pi(i,k) == undef) STOP "rhow_pi is WRONG!!!"
		SELECT CASE (AdvectionScheme)
		CASE (2)
			rhoww_pi(i,k) = rhow_pi(i,k)*w_pi(i,k)
			IF (rhow_pi(i,k) == undef .OR. w_pi(i,k) == undef) STOP "rhoww_pi is WRONG!!!"
		CASE (3)
			fa = w(i,k) + w(i,k+1)
			fb = w(i,k-1) + w(i,k+2)
			fc = w(i,k+1) - w(i,k)
			fd = w(i,k+2) - w(i,k-1)
			rhoww_pi(i,k) = rho_0(i,k)*(w_pi(i,k)/12.*(7*fa - fb) - ABS(w_pi(i,k))/12.*(3*fc - fd))
			IF(w(i,k+2) == undef .OR. w(i,k+1) == undef .OR. w(i,k) == undef .OR. w(i,k-1) == undef) STOP "rhoww_pi is WRONG!!!"
		CASE (4)
			fa = w(i,k) + w(i,k+1)
			fb = w(i,k-1) + w(i,k+2)
			rhoww_pi(i,k) = rho_0(i,k)*(w_pi(i,k)/12.*(7*fa - fb))
			IF(w(i,k+2) == undef .OR. w(i,k+1) == undef .OR. w(i,k) == undef .OR. w(i,k-1) == undef) STOP "rhoww_pi is WRONG!!!"
		CASE (5)
			fa = w(i,k) + w(i,k+1)
			fb = w(i,k-1) + w(i,k+2)
			fc = w(i,k-2) + w(i,k+3)
			rhoww_pi(i,k) = rho_0(i,k)*w_pi(i,k)/60.*(37*fa - 8*fb + fc)
			fd = w(i,k+1) - w(i,k)
			fe = w(i,k+2) - w(i,k-1)
			ff = w(i,k+3) - w(i,k-2)
			rhoww_pi(i,k) = rhoww_pi(i,k) - ABS(w_pi(i,k))/60.*(10*fd - 5*fe + ff)
			IF(w(i,k+2) == undef .OR. w(i,k+1) == undef .OR. w(i,k) == undef .OR. w(i,k-1) == undef .OR. w(i,k-2) == undef .OR. w(i,k+3) == undef) STOP "rhoww_pi is WRONG!!!"
		CASE (6)
			fa = w(i,k) + w(i,k+1)
			fb = w(i,k-1) + w(i,k+2)
			fc = w(i,k-2) + w(i,k+3)
			rhoww_pi(i,k) = rho_0(i,k)*w_pi(i,k)/60.*(37*fa - 8*fb + fc)
			IF(w(i,k+2) == undef .OR. w(i,k+1) == undef .OR. w(i,k) == undef .OR. w(i,k-1) == undef .OR. w(i,k-2) == undef .OR. w(i,k+3) == undef) STOP "rhoww_pi is WRONG!!!"
		CASE DEFAULT
			STOP "Wrong advection scheme!!!"
		END SELECT
	END DO
END DO
!-------------------------------------------------
! w-grid 
!-------------------------------------------------
DO i = its + 1, ite
	!DO k = kts + 1, kte
	DO k = kts, kte + 1
		PrhouPx_w(i,k) = (rhou_v(i,k) - rhou_v(i-1,k))/dx
		PrhouwPx_w(i,k) = (rhouw_v(i,k) - rhouw_v(i-1,k))/dx
		PrhowPz_w(i,k) = (rhow_pi(i,k) - rhow_pi(i,k-1))/dz
		PrhowwPz_w(i,k) = (rhoww_pi(i,k) - rhoww_pi(i,k-1))/dz
		F_w(i,k) = - 1./rho_0_w(i,k)*(PrhouwPx_w(i,k) - w(i,k)*PrhouPx_w(i,k) + PrhowwPz_w(i,k) - w(i,k)*PrhowPz_w(i,k)) + g*theta_1(i,k)/theta_0(i,k)
		IF (rhou_v(i,k) == undef .OR. rhou_v(i-1,k) == undef .OR. rhouw_v(i,k) == undef .OR. rhouw_v(i-1,k) == undef .OR. rhow_pi(i,k) == undef .OR. rhow_pi(i,k-1) == undef .OR. rhoww_pi(i,k) == undef .OR. rhoww_pi(i,k-1) == undef) STOP "F_w is WRONG!!!"
		Ppi_1Pz_w(i,k) = (pi_1(i,k) - pi_1(i,k - 1))/dz
		tend_w(i,k) = F_w(i,k) - Cp*theta_0(i,k)*Ppi_1Pz_w(i,k)
		
		P2wPx2_w(i,k) = (w(i+1,k) + w(i-1,k) - 2*w(i,k))/dx/dx
		P2wPz2_w(i,k) = (w(i,k+1) + w(i,k-1) - 2*w(i,k))/dz/dz
		IF (w(i+1,k) == undef .OR. w(i-1,k) == undef) STOP "P2uPx2_u is WRONG!!!"
		IF (w(i,k+1) == undef .OR. w(i,k-1) == undef) STOP "P2uPx2_u is WRONG!!!"
		
		tend_w(i,k) = tend_w(i,k) + Km*(P2wPx2_w(i,k) + P2wPz2_w(i,k)) ! Add diffusion term.
		
		IF (pi_1(i,k) == undef .OR. pi_1(i,k-1) == undef) STOP "Ppi_1Pz_w is WRONG!!!"
		IF (rhou_v(i,k) == undef .OR. rhou_v(i-1,k) == undef) STOP "PrhouPx_w is WRONG!!!"
		IF (rhouw_v(i,k) == undef .OR. rhouw_v(i-1,k) == undef) STOP "PrhouwPx_w is WRONG!!!"
		IF (rhow_pi(i,k) == undef .OR. rhow_pi(i,k-1) == undef) STOP "PrhowPz_w is WRONG!!!"
		IF (rhoww_pi(i,k) == undef .OR. rhoww_pi(i,k-1) == undef) STOP "PrhowwPz_w is WRONG!!!"
	END DO
END DO

!CALL debug_ascii_output(tend_w)
!CALL debug_test_boundary(F_w)
!-------------------------------------------------
IF (ANY(ISNAN(F_w(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITHT F_w!!!"
!=================================================
END SUBROUTINE tendency_w
!=================================================

!=================================================
SUBROUTINE tendency_theta(u,w,rho_0,theta,F_theta,tend_theta)
IMPLICIT NONE
!=================================================
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: u
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: w
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0    ! density
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: F_theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: tend_theta
!-------------------------------------------------
REAL(preci) :: fa, fb, fc, fd, fe, ff
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! 3. F_theta = - u p.theta/p.x - w p.theta/p.z
! 3.1. - u p.theta/p.x = - 1/rho (p.rhoutheta/p.x - theta p.rhou/p.x)
! 3.2. - w p.w/p.z = - 1/rho (p.rhothetaw/p.z - theta p.rhow/p.z)
!-------------------------------------------------
! v(virtual)-grid
!-------------------------------------------------
DO i = its, ite
	DO k = kts, kte + 1
		SELECT CASE (AdvectionScheme)
		CASE (2)
			rhoutheta_v(i,k) = rhou_v(i,k)*theta_v(i,k)
			IF (rhou_v(i,k) == undef .OR. theta_v(i,k) == undef) STOP "rhoutheta_v is WRONG!!!"
		CASE (3)
			fa = theta(i,k) + theta(i+1,k)
			fb = theta(i-1,k) + theta(i+2,k)
			fc = theta(i+1,k) - theta(i,k)
			fd = theta(i+2,k) - theta(i-1,k)
			rhoutheta_v(i,k) = rho_0_v(i,k)*(u_v(i,k)/12.*(7*fa - fb) - ABS(u_v(i,k))/12.*(3*fc - fd))
			IF (theta(i-1,k) == undef .OR. theta(i,k) == undef .OR. theta(i+1,k) == undef .OR. theta(i+2,k) == undef) STOP "rouutheta_v is WRONG!!!"
		CASE (4)
			fa = theta(i,k) + theta(i+1,k)
			fb = theta(i-1,k) + theta(i+2,k)
			rhoutheta_v(i,k) = rho_0_v(i,k)*(u_v(i,k)/12.*(7*fa - fb))
			IF (theta(i-1,k) == undef .OR. theta(i,k) == undef .OR. theta(i+1,k) == undef .OR. theta(i+2,k) == undef) STOP "rouutheta_v is WRONG!!!"
		CASE (5)
			fa = theta(i,k) + theta(i+1,k)
			fb = theta(i-1,k) + theta(i+2,k)
			fc = theta(i-2,k) + theta(i+3,k)
			rhoutheta_v(i,k) = rho_0_v(i,k)*u_v(i,k)/60.*(37*fa - 8*fb + fc)
			fd = theta(i+1,k) - theta(i,k)
			fe = theta(i+2,k) - theta(i-1,k)
			ff = theta(i+3,k) - theta(i-2,k)
			rhoutheta_v(i,k) = rhoutheta_v(i,k) - ABS(u_v(i,k))/60.*(10*fd - 5*fe + ff)
			IF (theta(i-1,k) == undef .OR. theta(i,k) == undef .OR. theta(i+1,k) == undef .OR. theta(i+2,k) == undef .OR. theta(i-2,k) == undef .OR. theta(i+3,k) == undef) STOP "rouutheta_v is WRONG!!!"
		CASE (6)
			fa = theta(i,k) + theta(i+1,k)
			fb = theta(i-1,k) + theta(i+2,k)
			fc = theta(i-2,k) + theta(i+3,k)
			rhoutheta_v(i,k) = rho_0_v(i,k)*u_v(i,k)/60.*(37*fa - 8*fb + fc)
			IF (theta(i-1,k) == undef .OR. theta(i,k) == undef .OR. theta(i+1,k) == undef .OR. theta(i+2,k) == undef .OR. theta(i-2,k) == undef .OR. theta(i+3,k) == undef) STOP "rouutheta_v is WRONG!!!"
		CASE DEFAULT
			STOP "Wrong advection scheme!!!"
		END SELECT
	END DO
END DO
!rhoutheta_v(its:ite,kts) = rhoutheta_v(its:ite,kts+1)      ! vertical boundary condition
!rhoutheta_v(its:ite,kte+1) = rhoutheta_v(its:ite,kte)      ! vertical boundary condition
!rhoutheta_v(its,kts:kte+1) = rhoutheta_v(its+1,kts:kte+1)  ! lateral boundary condition
!rhoutheta_v(ite,kts:kte+1) = rhoutheta_v(ite-1,kts:kte+1)  ! lateral boundary condition
!-------------------------------------------------
! pi-grid
!-------------------------------------------------
DO i = its, ite + 1
	DO k = kts - 1, kte + 1
		SELECT CASE (AdvectionScheme)
		CASE (2)
			rhowtheta_pi(i,k) = rhow_pi(i,k)*theta_pi(i,k)
			IF (rhow_pi(i,k) == undef .OR. theta_pi(i,k) == undef) STOP "rhowtheta_pi is WRONG!!!"
		CASE (3)
			fa = theta(i,k+1) + theta(i,k)
			fb = theta(i,k+2) + theta(i,k-1)
			fc = theta(i,k+1) - theta(i,k)
			fd = theta(i,k+2) - theta(i,k-1)
			rhowtheta_pi(i,k) = rho_0(i,k)*(w_pi(i,k)/12.*(7*fa - fb) - ABS(w_pi(i,k))/12.*(3*fc - fd))
			IF (theta(i,k-1) == undef .OR. theta(i,k) == undef .OR. theta(i,k+1) == undef .OR. theta(i,k+2) == undef) STOP "rhowtheta_pi is WRONG!!!"
		CASE (4)
			fa = theta(i,k+1) + theta(i,k)
			fb = theta(i,k+2) + theta(i,k-1)
			rhowtheta_pi(i,k) = rho_0(i,k)*(w_pi(i,k)/12.*(7*fa - fb))
			IF (theta(i,k-1) == undef .OR. theta(i,k) == undef .OR. theta(i,k+1) == undef .OR. theta(i,k+2) == undef) STOP "rhowtheta_pi is WRONG!!!"
		CASE (5)
			fa = theta(i,k+1) + theta(i,k)
			fb = theta(i,k+2) + theta(i,k-1)
			fc = theta(i,k+3) + theta(i,k-2)
			rhowtheta_pi(i,k) = rho_0(i,k)*w_pi(i,k)/60.*(37*fa - 8*fb + fc)
			fd = theta(i,k+1) - theta(i,k)
			fe = theta(i,k+2) - theta(i,k-1)
			ff = theta(i,k+3) - theta(i,k-2)
			rhowtheta_pi(i,k) = rhowtheta_pi(i,k) - ABS(w_pi(i,k))/60.*(10*fd - 5*fe + ff)
			IF (theta(i,k-1) == undef .OR. theta(i,k) == undef .OR. theta(i,k+1) == undef .OR. theta(i,k+2) == undef .OR. theta(i,k+3) == undef .OR. theta(i,k-2) == undef) STOP "rhowtheta_pi is WRONG!!!"
		CASE (6)
			fa = theta(i,k+1) + theta(i,k)
			fb = theta(i,k+2) + theta(i,k-1)
			fc = theta(i,k+3) + theta(i,k-2)
			rhowtheta_pi(i,k) = rho_0(i,k)*w_pi(i,k)/60.*(37*fa - 8*fb + fc)
			IF (theta(i,k-1) == undef .OR. theta(i,k) == undef .OR. theta(i,k+1) == undef .OR. theta(i,k+2) == undef .OR. theta(i,k+3) == undef .OR. theta(i,k-2) == undef) STOP "rhowtheta_pi is WRONG!!!"
		CASE DEFAULT
			STOP "Wrong advection scheme!!!"
		END SELECT
	END DO
END DO
!-------------------------------------------------
! w-grid - Theta on kts and kte+1 should also be updated.
!-------------------------------------------------
DO i = its + 1, ite
	!DO k = kts + 1, kte
	DO k = kts, kte + 1
		PrhouthetaPx_w(i,k) = (rhoutheta_v(i,k) - rhoutheta_v(i-1,k))/dx
		PrhowthetaPz_w(i,k) = (rhowtheta_pi(i,k) - rhowtheta_pi(i,k-1))/dz
		IF (rhoutheta_v(i,k) == undef .OR. rhoutheta_v(i-1,k) == undef) STOP "PrhouthetaPx_w is WRONG!!!"
		IF (rhowtheta_pi(i,k) == undef .OR. rhowtheta_pi(i,k-1) == undef) STOP "PrhowthetaPz_w is WRONG!!!"
		
		F_theta(i,k) = - 1./rho_0_w(i,k)*(PrhouthetaPx_w(i,k) - theta(i,k)*PrhouPx_w(i,k) + PrhowthetaPz_w(i,k) - theta(i,k)*PrhowPz_w(i,k))
		IF (rho_0_w(i,k) == undef .OR. PrhouthetaPx_w(i,k) == undef .OR. theta(i,k) == undef .OR. PrhouPx_w(i,k) == undef .OR. PrhowthetaPz_w(i,k) == undef .OR. PrhowPz_w(i,k) == undef) STOP "F_theta is WRONG!!!"

		tend_theta(i,k) = F_theta(i,k)
		
		P2thetaPx2_w(i,k) = (theta(i+1,k) + theta(i-1,k) - 2*theta(i,k))/dx/dx
		P2thetaPz2_w(i,k) = (theta(i,k+1) + theta(i,k-1) - 2*theta(i,k))/dz/dz
		IF (theta(i,k) == undef .OR. theta(i-1,k) == undef) STOP "P2thetaPx2_w is WRONG!!!"
		IF (theta(i,k+1) == undef .OR. theta(i,k-1) == undef) STOP "P2thetaPz2_w is WRONG!!!"
		
		tend_theta(i,k) = F_theta(i,k) + Kh*(P2thetaPx2_w(i,k) + P2thetaPz2_w(i,k)) ! Add diffusion term.
		
		IF (rhoutheta_v(i,k) == undef .OR. rhoutheta_v(i-1,k) == undef) STOP "PrhouthetaPx_w is WRONG!!!"
		IF (rhowtheta_pi(i,k) == undef .OR. rhowtheta_pi(i,k-1) == undef) STOP "PrhowthetaPz_w is WRONG!!!"
	END DO
END DO
!WRITE(*,*) PrhowPz_w(255,kts)
!WRITE(*,*) theta(255,kts)
!WRITE(*,*) PrhouthetaPx_w(255,kts)
!WRITE(*,*) rho_0_w(255,kts)
!WRITE(*,*) tend_theta(255,kts)
!-------------------------------------------------
IF (ANY(ISNAN(F_theta(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITH F_theta!!!"
!=================================================
END SUBROUTINE tendency_theta
!=================================================

!=================================================
SUBROUTINE tendency_pi(u,w,rho_0,theta_0,F_pi,tend_pi)
IMPLICIT NONE
!=================================================
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: u
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: w
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0    ! density
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_0
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: F_pi
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: tend_pi
INTEGER :: i, k
!=================================================
! 5.1 F_pi = - c^2/(rho_0*theta_0^2)*(PurhothetaPx + PwrhothetaPz)
!-------------------------------------------------
! u-grid
DO i = its, ite
	DO k = kts - 1, kte + 1
		urhotheta_u(i,k) = u(i,k)*rho_0_u(i,k)*theta_0_u(i,k)
		IF (u(i,k) == undef .OR. rho_0_u(i,k) == undef .OR. theta_0_u(i,k) == undef) STOP "urhotheta_u is WRONG!!!"
	END DO
END DO
!urhotheta_u(its,kts:kte) = 0.  ! lateral boundary condition
!urhotheta_u(ite,kts:kte) = 0.  ! lateral boundary condition

! w-grid
DO i = its, ite + 1
	DO k = kts, kte + 1
		wrhotheta_w(i,k) = w(i,k)*rho_0_w(i,k)*theta_0(i,k)
		IF (w(i,k) == undef .OR. rho_0_w(i,k) == undef .OR. theta_0(i,k) == undef) STOP "wrhotheta_w is WRONG!!!"
	END DO
END DO
!wrhotheta_w(its + 1:ite,kts) = 0.      ! vertical boundary condition
!wrhotheta_w(its + 1:ite,kte + 1) = 0.  ! vertical boundary condition

! To pi-grid
DO i = its + 1, ite
	DO k = kts, kte
		PurhothetaPx_pi(i,k) = (urhotheta_u(i,k) - urhotheta_u(i - 1,k))/dx
		PwrhothetaPz_pi(i,k) = (wrhotheta_w(i,k + 1) - wrhotheta_w(i,k))/dz
		F_pi(i,k) = - cs*cs/Cp/rho_0(i,k)/theta_0_pi(i,k)/theta_0_pi(i,k)*(PurhothetaPx_pi(i,k) + PwrhothetaPz_pi(i,k))
		tend_pi(i,k) = F_pi(i,k)
		IF (urhotheta_u(i,k) == undef .OR. urhotheta_u(i-1,k) == undef) STOP "PurhothetaPx_pi is WRONG!!!"
		IF (wrhotheta_w(i,k+1) == undef .OR. wrhotheta_w(i,k) == undef) STOP "PwrhothetaPz_pi is WRONG!!!"
	END DO
END DO
!-------------------------------------------------
IF (ANY(ISNAN(F_pi(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITH F_theta!!!"
!=================================================
END SUBROUTINE tendency_pi
!=================================================

!=================================================
END MODULE sp_module_tendency
!=================================================

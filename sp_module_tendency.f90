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
REAL(preci), DIMENSION(ims:ime,kms:kme) :: fa, fb, fc, fd, fe, ff
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! 1. F_u = - u p.u/p.x - w p.u/p.z + fv   (fv = 0.)
! 1.1. - u p.u/p.x = - 1/rho (p.rhouu/p.x - u p.rhou/p.x)
! 1.2. - w p.u/p.z = - 1/rho (p.rhouw/p.z - u p.rhow/p.z)
!-------------------------------------------------
! pi-grid - Middle-vars can be calculated on boundaries.
!-------------------------------------------------
CALL set_calc_area_pi
! ATTENTION: The calculated area includes the boundary layers.
imin = imin - 1
imax = imax + 1
kmin = kmin - 1
kmax = kmax + 1

IF (ANY(rho_0(imin:imax,kmin:kmax) == undef) .OR. ANY(u_pi(imin:imax,kmin:kmax) == undef)) STOP "rhou_pi is WRONG!!!"
FORALL(i = imin:imax, k = kmin:kmax)
	rhou_pi(i,k) = rho_0(i,k)*u_pi(i,k)
END FORALL

SELECT CASE (AdvectionScheme)

CASE (2)
	IF (ANY(rhou_pi(imin:imax,kmin:kmax) == undef) .OR. ANY(u_pi(imin:imax,kmin:kmax) == undef)) STOP "rhouu_pi is WRONG!!!"
	FORALL(i = imin:imax, k = kmin:kmax)
		rhouu_pi(i,k) = rhou_pi(i,k)*u_pi(i,k)
	END FORALL

CASE (3)
	IF (ANY(u(imin-2:imax+1,kmin:kmax) == undef)) STOP "rhouu_pi is WRONG!!!"
	FORALL(i = imin:imax, k = kmin:kmax)
		fa(i,k) = u(i-1,k) + u(i,k)
		fb(i,k) = u(i-2,k) + u(i+1,k)
		fc(i,k) = u(i,k) - u(i-1,k)
		fd(i,k) = u(i+1,k) - u(i-2,k)
		rhouu_pi(i,k) = rho_0(i,k)*(u_pi(i,k)/12.*(7*fa(i,k) - fb(i,k)) - ABS(u_pi(i,k))/12.*(3*fc(i,k) - fd(i,k)))
	END FORALL

CASE (4)
	IF (ANY(u(imin-2:imax+1,kmin:kmax) == undef)) STOP "rhouu_pi is WRONG!!!"
	FORALL(i = imin:imax, k = kmin:kmax)
		fa(i,k) = u(i-1,k) + u(i,k)
		fb(i,k) = u(i-2,k) + u(i+1,k)
		rhouu_pi(i,k) = rho_0(i,k)*(u_pi(i,k)/12.*(7*fa(i,k) - fb(i,k)))
	END FORALL

CASE (5)
	IF (ANY(u(imin-3:imax+2,kmin:kmax) == undef)) STOP "rhouu_pi is WRONG!!!"
	FORALL(i = imin:imax, k = kmin:kmax)
		fa(i,k) = u(i-1,k) + u(i,k)
		fb(i,k) = u(i-2,k) + u(i+1,k)
		fc(i,k) = u(i-3,k) + u(i+2,k)
		rhouu_pi(i,k) = rho_0(i,k)*u_pi(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
		fd(i,k) = u(i,k) - u(i-1,k)
		fe(i,k) = u(i+1,k) - u(i-2,k)
		ff(i,k) = u(i+2,k) - u(i-3,k)
		rhouu_pi(i,k) = rhouu_pi(i,k) - ABS(u_pi(i,k))/60.*(10*fd(i,k) - 5*fe(i,k) + ff(i,k))
	END FORALL

CASE (6)
	IF (ANY(u(imin-3:imax+2,kmin:kmax) == undef)) STOP "rhouu_pi is WRONG!!!"
	FORALL(i = imin:imax, k = kmin:kmax)
		fa(i,k) = u(i-1,k) + u(i,k)
		fb(i,k) = u(i-2,k) + u(i+1,k)
		fc(i,k) = u(i-3,k) + u(i+2,k)
		rhouu_pi(i,k) = rho_0(i,k)*u_pi(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
	END FORALL

CASE DEFAULT
	STOP "Wrong advection scheme!!!"
END SELECT

!-------------------------------------------------
! v(virtual)-grid
!-------------------------------------------------
CALL set_calc_area_v
! ATTENTION: The calculated area includes the boundary layers.
imin = imin - 1
imax = imax + 1
kmin = kmin - 1
kmax = kmax + 1

IF (ANY(rho_0_v(imin:imax,kmin:kmax) == undef) .OR. ANY(w_v(imin:imax,kmin:kmax) == undef) .OR. ANY(u_v(imin:imax,kmin:kmax) == undef)) STOP "rhow_v or rhouw_v is WRONG!!!"
FORALL (i = imin:imax, k = kmin:kmax)
	rhow_v(i,k) = rho_0_v(i,k)*w_v(i,k)
END FORALL

SELECT CASE (AdvectionScheme)
CASE (2)
	IF (ANY(rhow_v(imin:imax,kmin:kmax) == undef) .OR. ANY(u_v(imin:imax,kmin:kmax)== undef)) STOP "rhouw_v is WRONG!!!"
	FORALL (i = imin:imax, k = kmin:kmax)
		rhouw_v(i,k) = rhow_v(i,k)*u_v(i,k)
	END FORALL

CASE (3)
	IF (ANY(u(imin:imax,kmin-2:kmax+1) == undef)) STOP "rhouw_v is WRONG!!!"
	FORALL (i = imin:imax, k = kmin:kmax)
		fa(i,k) = u(i,k) + u(i,k-1)
		fb(i,k) = u(i,k+1) + u(i,k-2)
		fc(i,k) = u(i,k) - u(i,k-1)
		fd(i,k) = u(i,k+1) - u(i,k-2)
		rhouw_v(i,k) = rho_0_v(i,k)*(w_v(i,k)/12.*(7*fa(i,k) - fb(i,k)) - ABS(w_v(i,k))/12.*(3*fc(i,k) - fd(i,k)))
	END FORALL

CASE (4)
	IF (ANY(u(imin:imax,kmin-2:kmax+1) == undef)) STOP "rhouw_v is WRONG!!!"
	FORALL (i = imin:imax, k = kmin:kmax)
		fa(i,k) = u(i,k) + u(i,k-1)
		fb(i,k) = u(i,k+1) + u(i,k-2)
		rhouw_v(i,k) = rho_0_v(i,k)*(w_v(i,k)/12.*(7*fa(i,k) - fb(i,k)))
	END FORALL

CASE (5)
	IF (ANY(u(imin:imax,kmin-3:kmax+2) == undef)) STOP "rhouw_v is WRONG!!!"
	FORALL (i = imin:imax, k = kmin:kmax)
		fa(i,k) = u(i,k) + u(i,k-1)
		fb(i,k) = u(i,k+1) + u(i,k-2)
		fc(i,k) = u(i,k+2) + u(i,k-3)
		rhouw_v(i,k) = rho_0_v(i,k)*w_v(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
		fd(i,k) = u(i,k) - u(i,k-1)
		fe(i,k) = u(i,k+1) - u(i,k-2)
		ff(i,k) = u(i,k+2) - u(i,k-3)
		rhouw_v(i,k) = rhouw_v(i,k) - ABS(w_v(i,k))/60.*(10*fd(i,k) - 5*fe(i,k) + ff(i,k))
	END FORALL

CASE (6)
	IF (ANY(u(imin:imax,kmin-3:kmax+2) == undef)) STOP "rhouw_v is WRONG!!!"
	FORALL (i = imin:imax, k = kmin:kmax)
		fa(i,k) = u(i,k) + u(i,k-1)
		fb(i,k) = u(i,k+1) + u(i,k-2)
		fc(i,k) = u(i,k+2) + u(i,k-3)
		rhouw_v(i,k) = rho_0_v(i,k)*w_v(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
	END FORALL

CASE DEFAULT
	STOP "Wrong advection scheme!!!"
END SELECT
!-------------------------------------------------
! u-grid
!-------------------------------------------------
CALL set_calc_area_u

IF (ANY(u(imin-1:imax+1,k) == undef)) STOP "P2uPx2_u is WRONG!!!"
IF (ANY(u(imin:imax,kmin-1:kmax+1) == undef)) STOP "P2uPx2_u is WRONG!!!"

IF (ANY(pi_1(imin:imax+1,k) == undef)) STOP "Ppi_1Px_u is WRONG!!!"
IF (ANY(rhou_pi(imin:imax+1,k) == undef)) STOP "PrhouPx_u is WRONG!!!"
IF (ANY(rhouu_pi(imin:imax+1,k) == undef)) STOP "PrhouuPx_u is WRONG!!!"
IF (ANY(rhow_v(imin:imax,kmin:kmax+1) == undef)) STOP "PrhowPz_u is WRONG!!!"
IF (ANY(rhouw_v(imin:imax,kmin:kmax+1) == undef)) STOP "PrhouwPz_u is WRONG!!!"

FORALL (i = imin:imax, k = kmin:kmax)
	PrhouPx_u(i,k) = (rhou_pi(i+1,k) - rhou_pi(i,k))/dx
	PrhouuPx_u(i,k) = (rhouu_pi(i+1,k) - rhouu_pi(i,k))/dx
	PrhowPz_u(i,k) = (rhow_v(i,k+1) - rhow_v(i,k))/dz
	PrhouwPz_u(i,k) = (rhouw_v(i,k+1) - rhouw_v(i,k))/dz
	F_u(i,k) = - 1./rho_0_u(i,k)*(PrhouuPx_u(i,k) - u(i,k)*PrhouPx_u(i,k) + PrhouwPz_u(i,k) - u(i,k)*PrhowPz_u(i,k))
	Ppi_1Px_u(i,k) = (pi_1(i + 1,k) - pi_1(i,k))/dx
	tend_u(i,k) = F_u(i,k) - Cp*theta_0_u(i,k)*Ppi_1Px_u(i,k)
	
	P2uPx2_u(i,k) = (u(i+1,k) + u(i-1,k) - 2*u(i,k))/dx/dx
	P2uPz2_u(i,k) = (u(i,k+1) + u(i,k-1) - 2*u(i,k))/dz/dz
	
	tend_u(i,k) = tend_u(i,k) + Km*(P2uPx2_u(i,k) + P2uPz2_u(i,k)) ! Add diffusion term.
	
END FORALL
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
REAL(preci), DIMENSION(ims:ime,kms:kme) :: fa, fb, fc, fd, fe, ff
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! 2. F_w = - u p.w/p.x - w p.w/p.z + g(theta_1/theta_0)
! 2.1. - u p.w/p.x = - 1/rho (p.rhouw/p.x - w p.rhou/p.x)
! 2.2. - w p.w/p.z = - 1/rho (p.rhoww/p.z - w p.rhow/p.z)
!-------------------------------------------------
! v(virtual)-grid
!-------------------------------------------------
CALL set_calc_area_v
! ATTENTION: The calculated area includes the boundary layers.
imin = imin - 1
imax = imax + 1
kmin = kmin - 1
kmax = kmax + 1

IF (ANY(rho_0_v(imin:imax,kmin:kmax) == undef) .OR. ANY(u_v(imin:imax,kmin:kmax) == undef)) STOP "rhou_v is WRONG!!!"
FORALL (i = imin:imax, k = kmin:kmax)
	rhou_v(i,k) = rho_0_v(i,k)*u_v(i,k)
END FORALL

SELECT CASE (AdvectionScheme)
CASE (2)
	! rhouw_v have been calculated

CASE (3)
	IF (ANY(w(imin-1:imax+2,kmin:kmax) == undef)) STOP "rhouw_v is WRONG!!!"
	FORALL (i = imin:imax, k = kmin:kmax)
		fa(i,k) = w(i,k) + w(i+1,k)
		fb(i,k) = w(i-1,k) + w(i+2,k)
		fc(i,k) = w(i+1,k) - w(i,k)
		fd(i,k) = w(i+2,k) - w(i-1,k)
		rhouw_v(i,k) = rho_0_v(i,k)*(u_v(i,k)/12.*(7*fa(i,k) - fb(i,k)) - ABS(u_v(i,k))/12.*(3*fc(i,k) - fd(i,k)))
	END FORALL

CASE (4)
	IF (ANY(w(imin-1:imax+2,kmin:kmax) == undef)) STOP "rhouw_v is WRONG!!!"
	FORALL (i = imin:imax, k = kmin:kmax)
		fa(i,k) = w(i,k) + w(i+1,k)
		fb(i,k) = w(i-1,k) + w(i+2,k)
		rhouw_v(i,k) = rho_0_v(i,k)*(u_v(i,k)/12.*(7*fa(i,k) - fb(i,k)))
	END FORALL

CASE (5)
	IF (ANY(w(imin-2:imax+3,kmin:kmax) == undef)) STOP "rhouw_v is WRONG!!!"
	FORALL (i = imin:imax, k = kmin:kmax)
		fa(i,k) = w(i,k) + w(i+1,k)
		fb(i,k) = w(i-1,k) + w(i+2,k)
		fc(i,k) = w(i-2,k) + w(i+3,k)
		rhouw_v(i,k) = rho_0_v(i,k)*u_v(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
		fd(i,k) = w(i+1,k) - w(i,k)
		fe(i,k) = w(i+2,k) - w(i-1,k)
		ff(i,k) = w(i+3,k) - w(i-2,k)
		rhouw_v(i,k) = rhouw_v(i,k) - ABS(u_v(i,k))/60.*(10*fd(i,k) - 5*fe(i,k) + ff(i,k))
	END FORALL

CASE (6)
	IF (ANY(w(imin-2:imax+3,kmin:kmax) == undef)) STOP "rhouw_v is WRONG!!!"
	FORALL (i = imin:imax, k = kmin:kmax)
		fa(i,k) = w(i,k) + w(i+1,k)
		fb(i,k) = w(i-1,k) + w(i+2,k)
		fc(i,k) = w(i-2,k) + w(i+3,k)
		rhouw_v(i,k) = rho_0_v(i,k)*u_v(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
	END FORALL
	
CASE DEFAULT
	STOP "Wrong advection scheme!!!"
END SELECT
!-------------------------------------------------
! pi-grid
!-------------------------------------------------
CALL set_calc_area_pi
! ATTENTION: The calculated area includes the boundary layers.
imin = imin - 1
imax = imax + 1
kmin = kmin - 1
kmax = kmax + 1

IF (ANY(rho_0(imin:imax,kmin:kmax) == undef) .OR. ANY(w_pi(imin:imax,kmin:kmax) == undef)) STOP "rhow_pi is WRONG!!!"
FORALL (i = imin:imax, k = kmin:kmax)
	rhow_pi(i,k) = rho_0(i,k)*w_pi(i,k)
END FORALL

SELECT CASE (AdvectionScheme)

CASE (2)
	IF (ANY(rhow_pi(imin:imax,kmin:kmax) == undef) .OR. ANY(w_pi(imin:imax,kmin:kmax) == undef)) STOP "rhoww_pi is WRONG!!!"
	FORALL (i = imin:imax, k = kmin:kmax)
		rhoww_pi(i,k) = rhow_pi(i,k)*w_pi(i,k)
	END FORALL

CASE (3)
	IF (ANY(w(imin:imax,kmin-1:kmax+2) == undef)) STOP "rhoww_pi is WRONG!!!"
	FORALL (i = imin:imax, k = kmin:kmax)
		fa(i,k) = w(i,k) + w(i,k+1)
		fb(i,k) = w(i,k-1) + w(i,k+2)
		fc(i,k) = w(i,k+1) - w(i,k)
		fd(i,k) = w(i,k+2) - w(i,k-1)
		rhoww_pi(i,k) = rho_0(i,k)*(w_pi(i,k)/12.*(7*fa(i,k) - fb(i,k)) - ABS(w_pi(i,k))/12.*(3*fc(i,k) - fd(i,k)))
	END FORALL

CASE (4)
	IF (ANY(w(imin:imax,kmin-1:kmax+2) == undef)) STOP "rhoww_pi is WRONG!!!"
	FORALL (i = imin:imax, k = kmin:kmax)
		fa(i,k) = w(i,k) + w(i,k+1)
		fb(i,k) = w(i,k-1) + w(i,k+2)
		rhoww_pi(i,k) = rho_0(i,k)*(w_pi(i,k)/12.*(7*fa(i,k) - fb(i,k)))
	END FORALL

CASE (5)
	IF (ANY(w(imin:imax,kmin-2:kmax+3) == undef)) STOP "rhoww_pi is WRONG!!!"
	FORALL (i = imin:imax, k = kmin:kmax)
		fa(i,k) = w(i,k) + w(i,k+1)
		fb(i,k) = w(i,k-1) + w(i,k+2)
		fc(i,k) = w(i,k-2) + w(i,k+3)
		rhoww_pi(i,k) = rho_0(i,k)*w_pi(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
		fd(i,k) = w(i,k+1) - w(i,k)
		fe(i,k) = w(i,k+2) - w(i,k-1)
		ff(i,k) = w(i,k+3) - w(i,k-2)
		rhoww_pi(i,k) = rhoww_pi(i,k) - ABS(w_pi(i,k))/60.*(10*fd(i,k) - 5*fe(i,k) + ff(i,k))
	END FORALL

CASE (6)
	IF (ANY(w(imin:imax,kmin-2:kmax+3) == undef)) STOP "rhoww_pi is WRONG!!!"
	FORALL (i = imin:imax, k = kmin:kmax)
		fa(i,k) = w(i,k) + w(i,k+1)
		fb(i,k) = w(i,k-1) + w(i,k+2)
		fc(i,k) = w(i,k-2) + w(i,k+3)
		rhoww_pi(i,k) = rho_0(i,k)*w_pi(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
	END FORALL

CASE DEFAULT
	STOP "Wrong advection scheme!!!"
END SELECT
!-------------------------------------------------
! w-grid 
!-------------------------------------------------
CALL set_calc_area_w
! ATTENTION: The calculated area includes the boundary layers.
kmin = kmin - 1
kmax = kmax + 1

IF (ANY(rhou_v(imin-1:imax,kmin:kmax) == undef) .OR. ANY(rhouw_v(imin-1:imax,kmin:kmax) == undef) .OR. ANY(rhow_pi(imin:imax,kmin-1:kmax) == undef) .OR. ANY(rhoww_pi(imin:imax,kmin-1:kmax) == undef)) STOP "F_w is WRONG!!!"

IF (ANY(w(imin-1:imax+1,kmin:kmax) == undef)) STOP "P2wPx2_u is WRONG!!!"
IF (ANY(w(imin:imax,kmin-1:kmax+1) == undef)) STOP "P2wPz2_u is WRONG!!!"

IF (ANY(pi_1(imin:imax,kmin-1:kmax) == undef)) STOP "Ppi_1Pz_w is WRONG!!!"
IF (ANY(rhou_v(imin-1:imax,kmin:kmax) == undef)) STOP "PrhouPx_w is WRONG!!!"
IF (ANY(rhouw_v(imin-1:imax,kmin:kmax) == undef)) STOP "PrhouwPx_w is WRONG!!!"
IF (ANY(rhow_pi(imin:imax,kmin-1:kmax) == undef)) STOP "PrhowPz_w is WRONG!!!"
IF (ANY(rhoww_pi(imin:imax,kmin-1:kmax) == undef)) STOP "PrhowwPz_w is WRONG!!!"

FORALL( i = imin:imax, k = kmin:kmax)
	PrhouPx_w(i,k) = (rhou_v(i,k) - rhou_v(i-1,k))/dx
	PrhouwPx_w(i,k) = (rhouw_v(i,k) - rhouw_v(i-1,k))/dx
	PrhowPz_w(i,k) = (rhow_pi(i,k) - rhow_pi(i,k-1))/dz
	PrhowwPz_w(i,k) = (rhoww_pi(i,k) - rhoww_pi(i,k-1))/dz

	F_w(i,k) = - 1./rho_0_w(i,k)*(PrhouwPx_w(i,k) - w(i,k)*PrhouPx_w(i,k) + PrhowwPz_w(i,k) - w(i,k)*PrhowPz_w(i,k)) + g*theta_1(i,k)/theta_0(i,k)

	Ppi_1Pz_w(i,k) = (pi_1(i,k) - pi_1(i,k - 1))/dz
	tend_w(i,k) = F_w(i,k) - Cp*theta_0(i,k)*Ppi_1Pz_w(i,k)
	
	P2wPx2_w(i,k) = (w(i+1,k) + w(i-1,k) - 2*w(i,k))/dx/dx
	P2wPz2_w(i,k) = (w(i,k+1) + w(i,k-1) - 2*w(i,k))/dz/dz
	
	tend_w(i,k) = tend_w(i,k) + Km*(P2wPx2_w(i,k) + P2wPz2_w(i,k)) ! Add diffusion term.
	
END FORALL

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
REAL(preci), DIMENSION(ims:ime,kms:kme) :: fa, fb, fc, fd, fe, ff
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! 3. F_theta = - u p.theta/p.x - w p.theta/p.z
! 3.1. - u p.theta/p.x = - 1/rho (p.rhoutheta/p.x - theta p.rhou/p.x)
! 3.2. - w p.w/p.z = - 1/rho (p.rhothetaw/p.z - theta p.rhow/p.z)
!-------------------------------------------------
! v(virtual)-grid
!-------------------------------------------------
CALL set_calc_area_v
! ATTENTION: The calculated area includes the boundary layers.
imin = imin - 1
imax = imax + 1
kmin = kmin - 1
kmax = kmax + 1

SELECT CASE (AdvectionScheme)

CASE (2)
	IF (ANY(rhou_v(imin:imax,kmin:kmax) == undef) .OR. ANY(theta_v(imin:imax,kmin:kmax) == undef)) STOP "rhoutheta_v is WRONG!!!"
	FORALL( i = imin:imax, k = kmin:kmax)
		rhoutheta_v(i,k) = rhou_v(i,k)*theta_v(i,k)
	END FORALL

CASE (3)
	IF (ANY(theta(imin-1:imax+2,kmin:kmax) == undef)) STOP "rouutheta_v is WRONG!!!"
	FORALL( i = imin:imax, k = kmin:kmax)
		fa(i,k) = theta(i,k) + theta(i+1,k)
		fb(i,k) = theta(i-1,k) + theta(i+2,k)
		fc(i,k) = theta(i+1,k) - theta(i,k)
		fd(i,k) = theta(i+2,k) - theta(i-1,k)
		rhoutheta_v(i,k) = rho_0_v(i,k)*(u_v(i,k)/12.*(7*fa(i,k) - fb(i,k)) - ABS(u_v(i,k))/12.*(3*fc(i,k) - fd(i,k)))
	END FORALL

CASE (4)
	IF (ANY(theta(imin-1:imax+2,kmin:kmax) == undef)) STOP "rouutheta_v is WRONG!!!"
	FORALL( i = imin:imax, k = kmin:kmax)
		fa(i,k) = theta(i,k) + theta(i+1,k)
		fb(i,k) = theta(i-1,k) + theta(i+2,k)
		rhoutheta_v(i,k) = rho_0_v(i,k)*(u_v(i,k)/12.*(7*fa(i,k) - fb(i,k)))
	END FORALL

CASE (5)
	IF (ANY(theta(imin-1:imax+2,kmin:kmax) == undef)) STOP "rouutheta_v is WRONG!!!"
	FORALL( i = imin:imax, k = kmin:kmax)
		fa(i,k) = theta(i,k) + theta(i+1,k)
		fb(i,k) = theta(i-1,k) + theta(i+2,k)
		fc(i,k) = theta(i-2,k) + theta(i+3,k)
		rhoutheta_v(i,k) = rho_0_v(i,k)*u_v(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
		fd(i,k) = theta(i+1,k) - theta(i,k)
		fe(i,k) = theta(i+2,k) - theta(i-1,k)
		ff(i,k) = theta(i+3,k) - theta(i-2,k)
		rhoutheta_v(i,k) = rhoutheta_v(i,k) - ABS(u_v(i,k))/60.*(10*fd(i,k) - 5*fe(i,k) + ff(i,k))
	END FORALL

CASE (6)
	IF (ANY(theta(imin-1:imax+2,kmin:kmax) == undef)) STOP "rouutheta_v is WRONG!!!"
	FORALL( i = imin:imax, k = kmin:kmax)
		fa(i,k) = theta(i,k) + theta(i+1,k)
		fb(i,k) = theta(i-1,k) + theta(i+2,k)
		fc(i,k) = theta(i-2,k) + theta(i+3,k)
		rhoutheta_v(i,k) = rho_0_v(i,k)*u_v(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
	END FORALL

CASE DEFAULT
	STOP "Wrong advection scheme!!!"
END SELECT
!-------------------------------------------------
! pi-grid
!-------------------------------------------------
CALL set_calc_area_pi
! ATTENTION: The calculated area includes the boundary layers.
imin = imin - 1
imax = imax + 1
kmin = kmin - 1
kmax = kmax + 1

SELECT CASE (AdvectionScheme)
CASE (2)
IF (ANY(rhow_pi(imin:imax,kmin:kmax) == undef)) STOP "rhowtheta_pi is WRONG!!!"
FORALL( i = imin:imax, k = kmin:kmax)
	rhowtheta_pi(i,k) = rhow_pi(i,k)*theta_pi(i,k)
END FORALL

CASE (3)
IF (ANY(theta(imin:imax,kmin-1:kmax+2) == undef)) STOP "rhowtheta_pi is WRONG!!!"
FORALL( i = imin:imax, k = kmin:kmax)
	fa(i,k) = theta(i,k+1) + theta(i,k)
	fb(i,k) = theta(i,k+2) + theta(i,k-1)
	fc(i,k) = theta(i,k+1) - theta(i,k)
	fd(i,k) = theta(i,k+2) - theta(i,k-1)
	rhowtheta_pi(i,k) = rho_0(i,k)*(w_pi(i,k)/12.*(7*fa(i,k) - fb(i,k)) - ABS(w_pi(i,k))/12.*(3*fc(i,k) - fd(i,k)))
END FORALL

CASE (4)
IF (ANY(theta(imin:imax,kmin-1:kmax+2) == undef)) STOP "rhowtheta_pi is WRONG!!!"
FORALL( i = imin:imax, k = kmin:kmax)
	fa(i,k) = theta(i,k+1) + theta(i,k)
	fb(i,k) = theta(i,k+2) + theta(i,k-1)
	rhowtheta_pi(i,k) = rho_0(i,k)*(w_pi(i,k)/12.*(7*fa(i,k) - fb(i,k)))
END FORALL

CASE (5)
IF (ANY(theta(imin:imax,kmin-2:kmax+3) == undef)) STOP "rhowtheta_pi is WRONG!!!"
FORALL( i = imin:imax, k = kmin:kmax)
	fa(i,k) = theta(i,k+1) + theta(i,k)
	fb(i,k) = theta(i,k+2) + theta(i,k-1)
	fc(i,k) = theta(i,k+3) + theta(i,k-2)
	rhowtheta_pi(i,k) = rho_0(i,k)*w_pi(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
	fd(i,k) = theta(i,k+1) - theta(i,k)
	fe(i,k) = theta(i,k+2) - theta(i,k-1)
	ff(i,k) = theta(i,k+3) - theta(i,k-2)
	rhowtheta_pi(i,k) = rhowtheta_pi(i,k) - ABS(w_pi(i,k))/60.*(10*fd(i,k) - 5*fe(i,k) + ff(i,k))
END FORALL

CASE (6)
IF (ANY(theta(imin:imax,kmin-2:kmax+3) == undef)) STOP "rhowtheta_pi is WRONG!!!"
FORALL( i = imin:imax, k = kmin:kmax)
	fa(i,k) = theta(i,k+1) + theta(i,k)
	fb(i,k) = theta(i,k+2) + theta(i,k-1)
	fc(i,k) = theta(i,k+3) + theta(i,k-2)
	rhowtheta_pi(i,k) = rho_0(i,k)*w_pi(i,k)/60.*(37*fa(i,k) - 8*fb(i,k) + fc(i,k))
END FORALL

CASE DEFAULT
	STOP "Wrong advection scheme!!!"
END SELECT
!-------------------------------------------------
! w-grid - Theta on kts and kte+1 should also be updated.
!-------------------------------------------------
CALL set_calc_area_w
! ATTENTION: The calculated area includes the boundary layers.
kmin = kmin - 1
kmax = kmax + 1

IF (ANY(rhoutheta_v(imin-1:imax,kmin:kmax) == undef)) STOP "PrhouthetaPx_w is WRONG!!!"
IF (ANY(rhowtheta_pi(imin:imax,kmin-1:kmax) == undef)) STOP "PrhowthetaPz_w is WRONG!!!"

FORALL (i = imin:imax, k = kmin:kmax)
	PrhouthetaPx_w(i,k) = (rhoutheta_v(i,k) - rhoutheta_v(i-1,k))/dx
	PrhowthetaPz_w(i,k) = (rhowtheta_pi(i,k) - rhowtheta_pi(i,k-1))/dz
END FORALL

IF (ANY(rho_0_w(imin:imax,kmin:kmax) == undef) .OR. &
	ANY(PrhouthetaPx_w(imin:imax,kmin:kmax) == undef) .OR. &
	ANY(PrhouPx_w(imin:imax,kmin:kmax) == undef) .OR. &
	ANY(PrhowPz_w(imin:imax,kmin:kmax) == undef) .OR. &
	ANY(PrhowthetaPz_w(imin:imax,kmin:kmax) == undef)) STOP "F_theta is WRONG!!!"

IF (ANY(theta(imin-1:imax,kmin:kmax) == undef)) STOP "P2thetaPx2_w is WRONG!!!"
IF (ANY(theta(imin:imax,kmin-1:kmax+1) == undef)) STOP "P2thetaPz2_w is WRONG!!!"

FORALL (i = imin:imax, k = kmin:kmax)
	F_theta(i,k) = - 1./rho_0_w(i,k)*(PrhouthetaPx_w(i,k) - theta(i,k)*PrhouPx_w(i,k) + PrhowthetaPz_w(i,k) - theta(i,k)*PrhowPz_w(i,k))

	tend_theta(i,k) = F_theta(i,k)
	
	P2thetaPx2_w(i,k) = (theta(i+1,k) + theta(i-1,k) - 2*theta(i,k))/dx/dx
	P2thetaPz2_w(i,k) = (theta(i,k+1) + theta(i,k-1) - 2*theta(i,k))/dz/dz
	
	tend_theta(i,k) = F_theta(i,k) + Kh*(P2thetaPx2_w(i,k) + P2thetaPz2_w(i,k)) ! Add diffusion term.
END FORALL
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
CALL set_calc_area_u
! ATTENTION: The calculated area includes the boundary layers.
imin = imin - 1
imax = imax + 1
kmin = kmin - 1
kmax = kmax + 1

IF (ANY(u(imin:imax,kmin:kmax) == undef) .OR. ANY(rho_0_u(imin:imax,kmin:kmax) == undef) .OR. ANY(theta_0_u(imin:imax,kmin:kmax) == undef)) STOP "urhotheta_u is WRONG!!!"
FORALL (i = imin:imax, k = kmin:kmax)
	urhotheta_u(i,k) = u(i,k)*rho_0_u(i,k)*theta_0_u(i,k)
END FORALL

! w-grid
CALL set_calc_area_w
! ATTENTION: The calculated area includes the boundary layers.
imin = imin - 1
imax = imax + 1
kmin = kmin - 1
kmax = kmax + 1

IF (ANY(w(imin:imax,kmin:kmax) == undef) .OR. ANY(rho_0_w(imin:imax,kmin:kmax) == undef) .OR. ANY(theta_0(imin:imax,kmin:kmax) == undef)) STOP "wrhotheta_w is WRONG!!!"
FORALL (i = imin:imax, k = kmin:kmax)
	wrhotheta_w(i,k) = w(i,k)*rho_0_w(i,k)*theta_0(i,k)
END FORALL

! To pi-grid
CALL set_calc_area_pi

IF (ANY(urhotheta_u(imin-1:imax,kmin:kmax) == undef)) STOP "PurhothetaPx_pi is WRONG!!!"
IF (ANY(wrhotheta_w(imin:imax,kmin:kmax+1) == undef)) STOP "PwrhothetaPz_pi is WRONG!!!"
FORALL (i = imin:imax, k = kmin:kmax)
	PurhothetaPx_pi(i,k) = (urhotheta_u(i,k) - urhotheta_u(i - 1,k))/dx
	PwrhothetaPz_pi(i,k) = (wrhotheta_w(i,k + 1) - wrhotheta_w(i,k))/dz
	F_pi(i,k) = - cs*cs/Cp/rho_0(i,k)/theta_0_pi(i,k)/theta_0_pi(i,k)*(PurhothetaPx_pi(i,k) + PwrhothetaPz_pi(i,k))
	tend_pi(i,k) = F_pi(i,k)
END FORALL
!-------------------------------------------------
IF (ANY(ISNAN(F_pi(its:ite,kts:kte)))) STOP "SOMETHING IS WRONG WITH F_theta!!!"
!=================================================
END SUBROUTINE tendency_pi
!=================================================

!=================================================
END MODULE sp_module_tendency
!=================================================

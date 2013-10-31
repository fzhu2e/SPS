!=================================================
! The initial module of SPS-dynamic.
!-------------------------------------------------
! Version: 0.01
! Author: Zhu F.
! Email: lyricorpse@gmail.com
! Date: 2013-04-20 12:20:45 
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!=================================================
MODULE sp_module_initiate
USE sp_module_constant
USE sp_module_model
USE sp_module_interpolate
USE sp_module_debug
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================
! Initiate density current case.
!=================================================
SUBROUTINE initiate_dc(u,v,w,pi_1,pi_0,theta,theta_0,theta_1,rho_0)
IMPLICIT NONE
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: u        ! wind speed along x-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: v        ! wind speed along y-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: w        ! wind speed along z-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_1     ! pi'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_0 
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_0  ! theta = theta_0 + theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_1  ! theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: rho_0    ! density
!-------------------------------------------------
REAL(preci), PARAMETER :: x_c = 25.55*1000.  ! (m)
REAL(preci), PARAMETER :: z_c = 3.0*1000.     ! (m)
REAL(preci), PARAMETER :: r_x = 4*1000.  ! (m)
REAL(preci), PARAMETER :: r_z = 2*1000.  ! (m)
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_1_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_pi
!-------------------------------------------------
REAL(preci) :: L
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL debug_undef_all(pi,theta_1_pi,theta_pi)

CALL initiate_grid
!-------------------------------------------------
! theta_0, pi_0, rho_0
!-------------------------------------------------
CALL initiate_basic_state(theta_0,pi_0,rho_0)
!-------------------------------------------------
! u (on u-grid)
!-------------------------------------------------
CALL set_area_u
DO i = imin, imax
	DO k = kmin, kmax
		u(i,k) = 0.
	END DO
END DO

!-------------------------------------------------
! w (on w-grid)
!-------------------------------------------------
CALL set_area_w
DO i = imin, imax
	DO k = kmin, kmax
		w(i,k) = 0.
!-------------------------------------------------
! theta_1, theta (on w-grid)
!-------------------------------------------------
		L = SQRT((xpi(i) - x_c)*(xpi(i) - x_c)/r_x/r_x + (zz(k) - z_c)*(zz(k) - z_c)/r_z/r_z)
		IF (L <= 1.) THEN
			theta_1(i,k) = - 15./2.*(COS(PI_math*L) + 1.)  ! <= I want this.
		ELSE
			theta_1(i,k) = 0.
		END IF
		theta(i,k) = theta_0(i,k) + theta_1(i,k)   ! <= I want this.
	END DO
END DO

!-------------------------------------------------
! pi_1 (on pi-grid)
!-------------------------------------------------
CALL set_area_pi
DO i = imin, imax
	DO k = kmin, kmax
		!L = SQRT((xpi(i) - x_c)*(xpi(i) - x_c)/r_x/r_x + (zpi(k) - z_c)*(zpi(k) - z_c)/r_z/r_z)
		!IF (L <= 1.) THEN
			!theta_1_pi(i,k) = - 15./2.*(COS(PI_math*L) + 1.)  ! <= I want this.
		!ELSE
			!theta_1_pi(i,k) = 0.
		!END IF
		!theta_pi(i,k) = theta_0_pi(i,k) + theta_1_pi(i,k)
		pi_1(i,k) = 0.
	END DO
END DO
!=================================================
END SUBROUTINE initiate_dc
!=================================================

!=================================================
! Initiate thermal bubble case.
!=================================================
SUBROUTINE initiate_tb(u,v,w,pi_1,pi_0,theta,theta_0,theta_1,rho_0)
IMPLICIT NONE
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: u        ! wind speed along x-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: v        ! wind speed along y-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: w        ! wind speed along z-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_1     ! pi'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_0 
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_0  ! theta = theta_0 + theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_1  ! theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: rho_0    ! density
!-------------------------------------------------
REAL(preci), PARAMETER :: x_c = 10.0*1000. ! (m)
REAL(preci), PARAMETER :: z_c = 2.0*1000.  ! (m)
REAL(preci), PARAMETER :: R = 2.0*1000.    ! (m)
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_1_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_pi
!-------------------------------------------------
REAL(preci) :: L
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL debug_undef_all(pi,theta_1_pi,theta_pi)

CALL initiate_grid
!-------------------------------------------------
! theta_0, pi_0, rho_0
!-------------------------------------------------
CALL initiate_basic_state(theta_0,pi_0,rho_0)
!-------------------------------------------------
! u (on u-grid)
!-------------------------------------------------
CALL set_area_u
DO i = imin, imax
	DO k = kmin, kmax
		u(i,k) = 0.
	END DO
END DO

!-------------------------------------------------
! w (on w-grid)
!-------------------------------------------------
CALL set_area_w
DO i = imin, imax
	DO k = kmin, kmax
		w(i,k) = 0.
!-------------------------------------------------
! theta_1, theta (on w-grid)
!-------------------------------------------------
		L = SQRT((xpi(i) - x_c)*(xpi(i) - x_c) + (zz(k) - z_c)*(zz(k) - z_c))
		theta_1(i,k) = 2.*MAX(0.,1. - L/R)
		theta(i,k) = theta_0(i,k) + theta_1(i,k)
	END DO
END DO

!-------------------------------------------------
! pi_1 (on pi-grid)
!-------------------------------------------------
CALL set_area_pi
DO i = imin, imax
	DO k = kmin, kmax
		!L = SQRT((xpi(i) - x_c)*(xpi(i) - x_c) + (zpi(k) - z_c)*(zpi(k) - z_c))
		!theta_1_pi(i,k) = 2.*MAX(0.,1. - L/R)
		!theta_pi(i,k) = theta_0_pi(i,k) + theta_1_pi(i,k)
		pi_1(i,k) = 0.
	END DO
END DO
!=================================================
END SUBROUTINE initiate_tb
!=================================================

!=================================================
! Initiate inertia gravity waves.
!=================================================
SUBROUTINE initiate_igw(u,v,w,pi_1,pi_0,theta,theta_0,theta_1,rho_0)
IMPLICIT NONE
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: u        ! wind speed along x-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: v        ! wind speed along y-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: w        ! wind speed along z-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_1     ! pi'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_0 
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_0  ! theta = theta_0 + theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_1  ! theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: rho_0    ! density
!-------------------------------------------------
REAL(preci), PARAMETER :: x_c = 100.0*1000. ! (m)
REAL(preci), PARAMETER :: H = 10.0*1000.    ! (m)
REAL(preci), PARAMETER :: a = 5.0*1000.     ! (m)
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_1_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_pi
!-------------------------------------------------
REAL(preci) :: L
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL debug_undef_all(pi,theta_1_pi,theta_pi)

CALL initiate_grid
!-------------------------------------------------
! theta_0, pi_0, rho_0
!-------------------------------------------------
CALL initiate_basic_state(theta_0,pi_0,rho_0)
!-------------------------------------------------
! u (on u-grid)
!-------------------------------------------------
CALL set_area_u
DO i = imin, imax
	DO k = kmin, kmax
		u(i,k) = 20.
	END DO
END DO

!-------------------------------------------------
! w (on w-grid)
!-------------------------------------------------
CALL set_area_w
DO i = imin, imax
	DO k = kmin, kmax
		w(i,k) = 0.
!-------------------------------------------------
! theta_1, theta (on w-grid)
!-------------------------------------------------
		L = SIN(PI_math*zz(k)/H)/(1. + (xx(i) - x_c)*(xx(i) - x_c)/a/a)
		theta_1(i,k) = 0.01*L
		theta(i,k) = theta_0(i,k) + theta_1(i,k)
	END DO
END DO

!-------------------------------------------------
! pi_1 (on pi-grid)
!-------------------------------------------------
CALL set_area_pi
DO i = imin, imax
	DO k = kmin, kmax
		!L = SIN(PI_math*zpi(k)/H)/(1. + (xpi(i) - x_c)*(xpi(i) - x_c)/a/a)
		!theta_1_pi(i,k) = 0.01*L
		!theta_pi(i,k) = theta_0_pi(i,k) + theta_1_pi(i,k)
		pi_1(i,k) = 0.
	END DO
END DO
!=================================================
END SUBROUTINE initiate_igw
!=================================================


!=================================================
! Initiate Schar mountain case.
!=================================================
SUBROUTINE initiate_Sm(u,v,w,pi_1,pi_0,theta,theta_0,theta_1,rho_0)
IMPLICIT NONE
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: u        ! wind speed along x-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: v        ! wind speed along y-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: w        ! wind speed along z-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_1     ! pi'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_0 
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_0  ! theta = theta_0 + theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_1  ! theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: rho_0    ! density
!-------------------------------------------------
REAL(preci), PARAMETER :: h0 = 250.0            ! (m)
!REAL(preci), PARAMETER :: h0 = 50.0            ! (m)
REAL(preci), PARAMETER :: a0 = 5.0*1000.        ! (m)
REAL(preci), PARAMETER :: x_c = 25.0*1000.        ! (m)
REAL(preci), PARAMETER :: lambda0 = 4.0*1000.   ! (m)
REAL(preci), PARAMETER :: N0 = 0.01             ! (s-1)

REAL(preci), PARAMETER :: Ts = 280.             ! (K)
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_1_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_pi
!-------------------------------------------------
INTEGER :: i, k 
!=================================================
! To initiate:
! zs(i), zs_pi(i), PzsPx(i), PzsPx_pi(i)
! u(i,k), w(i,k)
! theta_0(i,k), pi_0(i,k), rho_0(i,k)
! theta(i,k), theta_1(i,k), pi_1(i,k)
!=================================================
CALL initiate_grid
!-------------------------------------------------
zs = undef
zs_pi = undef
PzsPx = undef
PzsPx_pi = undef

CALL set_area_u
DO i = imin, imax
	zs(i) = h0*EXP(-((xx(i) - x_c)/a0)**2)*COS(PI_math*(xx(i) - x_c)/lambda0)**2
	PzsPx(i) = - 2*h0*EXP(-((xx(i) - x_c)/a0)**2)*COS(PI_math*(xx(i) - x_c)/lambda0)*((xx(i) - x_c)/a0**2*COS(PI_math*(xx(i) - x_c)/lambda0) + PI_math/lambda0*SIN(PI_math*(xx(i) - x_c)/lambda0))
END DO

CALL set_area_pi
DO i = imin, imax
	zs_pi(i) = h0*EXP(-((xpi(i) - x_c)/a0)**2)*COS(PI_math*(xpi(i) - x_c)/lambda0)**2
	PzsPx_pi(i) = (zs(i) - zs(i-1))/dx
END DO
!-------------------------------------------------
CALL initiate_vertcoords
CALL initiate_basic_state(theta_0,pi_0,rho_0)
!-------------------------------------------------
CALL set_area_u
DO i = imin, imax
	DO k = kmin, kmax
		u(i,k) = 10.
	END DO
END DO

CALL set_area_w
DO i = imin, imax
	DO k = kmin, kmax
		w(i,k) = 0.
		theta_1(i,k) = 0.
		theta(i,k) = theta_0(i,k)
	END DO
END DO

CALL set_area_pi
DO i = imin, imax
	DO k = kmin, kmax
		pi_1(i,k) = 0.
	END DO
END DO
!=================================================
END SUBROUTINE initiate_Sm
!=================================================

!=================================================
! Initiate Grid Position
!=================================================
SUBROUTINE initiate_grid
IMPLICIT NONE
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! Undefine variations.
!-------------------------------------------------
xx = undef
xpi = undef
zz = undef
zpi = undef
!-------------------------------------------------

CALL set_area_u
DO i = imin, imax
	xx(i) = dx*(i - its)
END DO

CALL set_area_pi
DO i = imin, imax
	xpi(i) = (xx(i-1) + xx(i))/2.
END DO

CALL set_area_w
DO k = kmin, kmax
	zz(k) = dz*(k - kts)
END DO

CALL set_area_pi
DO k = kmin, kmax
	zpi(k) = (zz(k) + zz(k + 1))/2.
END DO

!-------------------------------------------------
! Unit Test
!-------------------------------------------------
!WRITE(*,*) xx
!WRITE(*,*) xpi
!WRITE(*,*) zz
!WRITE(*,*) zpi
!=================================================
END SUBROUTINE initiate_grid
!=================================================

!=================================================
! Initiate Basic State Atmosphere
! 
! To initiate:
! 
! theta_0, pi_0, rho_0
! theta_0_pi, theta_0_u, theta_0_v
! rho_0_u, rho_0_w, rho_0_v
!=================================================
SUBROUTINE initiate_basic_state(theta_0,pi_0,rho_0)
IMPLICIT NONE
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_0
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_0
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: rho_0
!-------------------------------------------------
REAL(preci), PARAMETER :: Ts = 300.         ! (K)
REAL(preci), PARAMETER :: N0 = 0.01         ! (s-1)
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: pi_0_w
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! Undefine variations.
!-------------------------------------------------
CALL debug_undef_all(theta_0,pi_0,rho_0)
CALL debug_undef_all(theta_0_pi,theta_0_u,theta_0_v)
CALL debug_undef_all(rho_0_u,rho_0_w,rho_0_v)

!-------------------------------------------------
! theta_0 (on w-grid)
!-------------------------------------------------
CALL set_area_w

DO i = imin, imax
	DO k = kmin, kmax
		IF (RunCase == 1 .OR. RunCase == 2) THEN
			theta_0(i,k) = Ts
			pi_0_w(i,k) = 1. - g*zz(k)/2./Cp/Ts
		ELSE IF (RunCase == 3) THEN
			theta_0(i,k) = Ts*EXP(N0*N0/g*zz(k))
			pi_0_w(i,k) = 1. + g*g/Cp/N0/N0/Ts*(EXP(-N0*N0*zz(k)/g) - 1.)
		ELSE IF (RunCase == 4) THEN
			theta_0(i,k) = (Ts - 20.)*EXP(N0*N0/g*z_hat(i,k))
			pi_0_w(i,k) = 1. + g*g/Cp/N0/N0/Ts*(EXP(-N0*N0*z_hat(i,k)/g) - EXP(-N0*N0*zs(i)/g))
		ELSE
			STOP "WRONG RunCase!!!"
		END IF
		rho_0_w(i,k) = p0/Rd/theta_0(i,k)*pi_0_w(i,k)*pi_0_w(i,k)**(Cp/Rd)
	END DO
END DO

CALL set_area_v
DO i = imin, imax
	DO k = kmin, kmax
		theta_0_v(i,k) = theta_0(imin+1,k)
		rho_0_v(i,k) = rho_0_w(imin+1,k)
	END DO
END DO

!-------------------------------------------------
! pi_0, rho_0, rho_0 (on pi-grid)
!-------------------------------------------------
CALL set_area_pi

DO i = imin, imax
	DO k = kmin, kmax
		IF (RunCase == 1 .OR. RunCase == 2) THEN
			theta_0_pi(i,k) = Ts
			pi_0(i,k) = 1. - g*zpi(k)/2./Cp/Ts
		ELSE IF (RunCase == 3) THEN
			theta_0_pi(i,k) = Ts*EXP(N0*N0/g*zpi(k))
			pi_0(i,k) = 1. + g*g/Cp/N0/N0/Ts*(EXP(-N0*N0*zpi(k)/g) - 1.)
		ELSE IF (RunCase == 4) THEN
			theta_0_pi(i,k) = (Ts - 20.)*EXP(N0*N0/g*z_hat_pi(i,k))
			pi_0(i,k) = 1. + g*g/Cp/N0/N0/Ts*(EXP(-N0*N0*z_hat_pi(i,k)/g) - EXP(-N0*N0*zs(i)/g))
		ELSE
			STOP "WRONG RunCase!!!"
		END IF
		rho_0(i,k) = p0/Rd/theta_0_pi(i,k)*pi_0(i,k)*pi_0(i,k)**(Cp/Rd)
	END DO
END DO

CALL set_area_u
DO i = imin, imax
	DO k = kmin, kmax
		theta_0_u(i,k) = theta_0_pi(imin+1,k)
		rho_0_u(i,k) = rho_0(imin+1,k)
	END DO
END DO

!=================================================
END SUBROUTINE initiate_basic_state
!=================================================

!=================================================
! Initiate Vertical Coordinates
! 
! To initiate:
! b(k), b_pi(k), VertA_u(i), VertA_pi(i)
! z_hat(i,k), z_hat_pi(i,k), z_hat_u(i,k), z_hat_v(i,k)  =  zz(k) + b(k)*zs(i)
! PbPzhat(i,k), PbPzhat_pi(i,k), PbPzhat_u(i,k), PbPhat_v(i,k)
! VertA_u(i), VertA_pi(i) = PzsPx(i)
! VertB_u(i,k), VertB_pi(i,k), VertB_u(i,k), VertB_v(i,k)
! VertC_u(i,k), VertC_pi(i,k), VertC_u(i,k), VertC_v(i,k)
!=================================================
SUBROUTINE initiate_vertcoords
IMPLICIT NONE
REAL(preci), PARAMETER :: sh = 3.0*1000.        ! (m)
!-------------------------------------------------
INTEGER :: i, k
!=================================================
VertA_u = undef
VertA_pi = undef

CALL set_area_u
DO i = imin, imax
	VertA_u(i) = PzsPx(i)
END DO

CALL set_area_pi
DO i = imin, imax
	VertA_pi(i) = PzsPx_pi(i)
END DO
!-------------------------------------------------
b = undef
b_pi = undef

CALL set_area_w
DO k = kmin, kmax
	b(k) = SINH((ztop - zz(k))/sh)/SINH(ztop/sh)
END DO

CALL set_area_pi
DO k = kmin, kmax
	b_pi(k) = SINH((ztop - zpi(k))/sh)/SINH(ztop/sh)
END DO
!-------------------------------------------------
CALL debug_undef_all(z_hat,z_hat_pi,z_hat_u,z_hat_v)
CALL debug_undef_all(PbPzhat,PbPzhat_pi,PbPzhat_u,PbPzhat_v)
CALL debug_undef_all(VertB_w,VertB_pi,VertB_u,VertB_v)
CALL debug_undef_all(VertC_w,VertC_pi,VertC_u,VertC_v)

CALL set_area_w
DO i = imin, imax
	DO k = kmin, kmax
		z_hat(i,k) = zz(k) + zs_pi(i)*b(k)
		PbPzhat(i,k) = - COSH((ztop - z_hat(i,k))/sh)/SINH(ztop/sh)/sh
		VertB_w(i,k) = 1 + zs_pi(i)*PbPzhat(i,k)
		VertC_w(i,k) = b(k)*VertA_pi(i)/VertB_w(i,k)
	END DO
END DO

CALL set_area_pi
DO i = imin, imax
	DO k = kmin, kmax
		z_hat_pi(i,k) = zpi(k) + zs_pi(i)*b_pi(k)
		PbPzhat_pi(i,k) = - COSH((ztop - z_hat_pi(i,k))/sh)/SINH(ztop/sh)/sh
		VertB_pi(i,k) = 1 + zs_pi(i)*PbPzhat_pi(i,k)
		VertC_pi(i,k) = b_pi(k)*VertA_pi(i)/VertB_pi(i,k)
	END DO
END DO

CALL set_area_u
DO i = imin, imax
	DO k = kmin, kmax
		z_hat_u(i,k) = zpi(k) + zs(i)*b_pi(k)
		PbPzhat_u(i,k) = - COSH((ztop - z_hat_u(i,k))/sh)/SINH(ztop/sh)/sh
		VertB_u(i,k) = 1 + zs(i)*PbPzhat_u(i,k)
		VertC_u(i,k) = b_pi(k)*VertA_u(i)/VertB_u(i,k)
	END DO
END DO

CALL set_area_v
DO i = imin, imax
	DO k = kmin, kmax
		z_hat_v(i,k) = zz(k) + zs(i)*b(k)
		PbPzhat_v(i,k) = - COSH((ztop - z_hat_v(i,k))/sh)/SINH(ztop/sh)/sh
		VertB_v(i,k) = 1 + zs(i)*PbPzhat_v(i,k)
		VertC_v(i,k) = b(k)*VertA_u(i)/VertB_v(i,k)
	END DO
END DO

OPEN(1, FILE="./output/modelvar_xpi.bin", FORM='binary', CONVERT='big_endian')
	WRITE(1) xpi(its:ite)
CLOSE(1)

OPEN(1, FILE="./output/modelvar_zhat.bin", FORM='binary', CONVERT='big_endian')
	WRITE(1) z_hat(its:ite,kts:kte)
CLOSE(1)

!=================================================
END SUBROUTINE initiate_vertcoords
!=================================================

!=================================================
END MODULE sp_module_initiate
!=================================================

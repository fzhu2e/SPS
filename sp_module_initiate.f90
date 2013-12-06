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
SUBROUTINE initiate_dc(u,w,pi_1,pi_0,theta,theta_0,theta_1,rho_0)
IMPLICIT NONE
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: u        ! wind speed along x-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: w        ! wind speed along z-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_1     ! pi'
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_0 
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_0  ! theta = theta_0 + theta'
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_1  ! theta'
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: rho_0    ! density
!-------------------------------------------------
REAL(kd), PARAMETER :: x_c = 25.55*1000.  ! (m)
REAL(kd), PARAMETER :: z_c = 3.0*1000.     ! (m)
REAL(kd), PARAMETER :: r_x = 4*1000.  ! (m)
REAL(kd), PARAMETER :: r_z = 2*1000.  ! (m)
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: pi
REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_1_pi
REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_pi
!-------------------------------------------------
REAL(kd) :: L
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
!SUBROUTINE initiate_tb(u,w,pi_1,pi_0,theta,theta_0,theta_1,rho_0)
!IMPLICIT NONE
!!-------------------------------------------------
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: u        ! wind speed along x-axis
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: w        ! wind speed along z-axis
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_1     ! pi'
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_0 
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_0  ! theta = theta_0 + theta'
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_1  ! theta'
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: rho_0    ! density
!!-------------------------------------------------
!REAL(kd), PARAMETER :: x_c = 10.0*1000. ! (m)
!REAL(kd), PARAMETER :: z_c = 2.0*1000.  ! (m)
!REAL(kd), PARAMETER :: R = 2.0*1000.    ! (m)
!!-------------------------------------------------
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: pi
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_1_pi
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_pi
!!-------------------------------------------------
!REAL(kd) :: L
!!-------------------------------------------------
!INTEGER :: i, k
!!=================================================
!CALL debug_undef_all(pi,theta_1_pi,theta_pi)

!CALL initiate_grid
!!-------------------------------------------------
!! theta_0, pi_0, rho_0
!!-------------------------------------------------
!CALL initiate_basic_state(theta_0,pi_0,rho_0)
!!-------------------------------------------------
!! u (on u-grid)
!!-------------------------------------------------
!CALL set_area_u
!DO i = imin, imax
	!DO k = kmin, kmax
		!u(i,k) = 0.
	!END DO
!END DO

!!-------------------------------------------------
!! w (on w-grid)
!!-------------------------------------------------
!CALL set_area_w
!DO i = imin, imax
	!DO k = kmin, kmax
		!w(i,k) = 0.
!!-------------------------------------------------
!! theta_1, theta (on w-grid)
!!-------------------------------------------------
		!L = SQRT((xpi(i) - x_c)*(xpi(i) - x_c) + (zz(k) - z_c)*(zz(k) - z_c))
		!theta_1(i,k) = 2.*MAX(0.,1. - L/R)
		!theta(i,k) = theta_0(i,k) + theta_1(i,k)
	!END DO
!END DO

!!-------------------------------------------------
!! pi_1 (on pi-grid)
!!-------------------------------------------------
!CALL set_area_pi
!DO i = imin, imax
	!DO k = kmin, kmax
		!!L = SQRT((xpi(i) - x_c)*(xpi(i) - x_c) + (zpi(k) - z_c)*(zpi(k) - z_c))
		!!theta_1_pi(i,k) = 2.*MAX(0.,1. - L/R)
		!!theta_pi(i,k) = theta_0_pi(i,k) + theta_1_pi(i,k)
		!pi_1(i,k) = 0.
	!END DO
!END DO
!!=================================================
!END SUBROUTINE initiate_tb
!=================================================

!=================================================
! Initiate inertia gravity waves.
!=================================================
!SUBROUTINE initiate_igw(u,w,pi_1,pi_0,theta,theta_0,theta_1,rho_0)
!IMPLICIT NONE
!!-------------------------------------------------
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: u        ! wind speed along x-axis
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: w        ! wind speed along z-axis
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_1     ! pi'
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_0 
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_0  ! theta = theta_0 + theta'
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_1  ! theta'
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: rho_0    ! density
!!-------------------------------------------------
!REAL(kd), PARAMETER :: x_c = 100.0*1000. ! (m)
!REAL(kd), PARAMETER :: H = 10.0*1000.    ! (m)
!REAL(kd), PARAMETER :: a = 5.0*1000.     ! (m)
!!-------------------------------------------------
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: pi
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_1_pi
!REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_pi
!!-------------------------------------------------
!REAL(kd) :: L
!!-------------------------------------------------
!INTEGER :: i, k
!!=================================================
!CALL debug_undef_all(pi,theta_1_pi,theta_pi)

!CALL initiate_grid
!!-------------------------------------------------
!! theta_0, pi_0, rho_0
!!-------------------------------------------------
!CALL initiate_basic_state(theta_0,pi_0,rho_0)
!!-------------------------------------------------
!! u (on u-grid)
!!-------------------------------------------------
!CALL set_area_u
!DO i = imin, imax
	!DO k = kmin, kmax
		!u(i,k) = 20.
	!END DO
!END DO

!!-------------------------------------------------
!! w (on w-grid)
!!-------------------------------------------------
!CALL set_area_w
!DO i = imin, imax
	!DO k = kmin, kmax
		!w(i,k) = 0.
!!-------------------------------------------------
!! theta_1, theta (on w-grid)
!!-------------------------------------------------
		!L = SIN(PI_math*zz(k)/H)/(1. + (xx(i) - x_c)*(xx(i) - x_c)/a/a)
		!theta_1(i,k) = 0.01*L
		!theta(i,k) = theta_0(i,k) + theta_1(i,k)
	!END DO
!END DO

!!-------------------------------------------------
!! pi_1 (on pi-grid)
!!-------------------------------------------------
!CALL set_area_pi
!DO i = imin, imax
	!DO k = kmin, kmax
		!!L = SIN(PI_math*zpi(k)/H)/(1. + (xpi(i) - x_c)*(xpi(i) - x_c)/a/a)
		!!theta_1_pi(i,k) = 0.01*L
		!!theta_pi(i,k) = theta_0_pi(i,k) + theta_1_pi(i,k)
		!pi_1(i,k) = 0.
	!END DO
!END DO
!!=================================================
!END SUBROUTINE initiate_igw
!=================================================


!=================================================
! Initiate Schar mountain case.
!=================================================
!SUBROUTINE initiate_Sm(u,w,pi_1,pi_0,theta,theta_0,theta_1,rho_0)
!IMPLICIT NONE
!!-------------------------------------------------
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: u        ! wind speed along x-axis
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: w        ! wind speed along z-axis
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_1     ! pi'
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_0 
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_0  ! theta = theta_0 + theta'
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_1  ! theta'
!REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: rho_0    ! density
!!-------------------------------------------------
!!REAL(kd), PARAMETER :: h0 = 250.0            ! (m)
!!!REAL(kd), PARAMETER :: h0 = 50.0            ! (m)
!!REAL(kd), PARAMETER :: a0 = 5.0*1000.        ! (m)
!!REAL(kd), PARAMETER :: x_c = 25.0*1000.        ! (m)
!!REAL(kd), PARAMETER :: lambda0 = 4.0*1000.   ! (m)
!!REAL(kd), PARAMETER :: N0 = 0.01             ! (s-1)

!!REAL(kd), PARAMETER :: Ts = 280.             ! (K)
!!!-------------------------------------------------
!!REAL(kd), DIMENSION(ims:ime,kms:kme) :: pi
!!REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_1_pi
!!REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_pi
!!!-------------------------------------------------
!!INTEGER :: i, k 
!!!=================================================
!!! To initiate:
!!! zs(i), zs_pi(i), PzsPx(i), PzsPx_pi(i)
!!! u(i,k), w(i,k)
!!! theta_0(i,k), pi_0(i,k), rho_0(i,k)
!!! theta(i,k), theta_1(i,k), pi_1(i,k)
!!!=================================================
!!CALL initiate_grid
!!!-------------------------------------------------
!!zs = undef
!!zs_pi = undef
!!PzsPx = undef
!!PzsPx_pi = undef

!!CALL set_area_u
!!DO i = imin, imax
	!!zs(i) = h0*EXP(-((xx(i) - x_c)/a0)**2)*COS(PI_math*(xx(i) - x_c)/lambda0)**2
	!!PzsPx(i) = - 2*h0*EXP(-((xx(i) - x_c)/a0)**2)*COS(PI_math*(xx(i) - x_c)/lambda0)*((xx(i) - x_c)/a0**2*COS(PI_math*(xx(i) - x_c)/lambda0) + PI_math/lambda0*SIN(PI_math*(xx(i) - x_c)/lambda0))
!!END DO

!!CALL set_area_pi
!!DO i = imin, imax
	!!zs_pi(i) = h0*EXP(-((xpi(i) - x_c)/a0)**2)*COS(PI_math*(xpi(i) - x_c)/lambda0)**2
	!!PzsPx_pi(i) = (zs(i) - zs(i-1))/dx
!!END DO
!!!-------------------------------------------------
!!CALL initiate_virertcoords
!!CALL initiate_basic_state(theta_0,pi_0,rho_0)
!!!-------------------------------------------------
!!CALL set_area_u
!!DO i = imin, imax
	!!DO k = kmin, kmax
		!!u(i,k) = 10.
	!!END DO
!!END DO

!!CALL set_area_w
!!DO i = imin, imax
	!!DO k = kmin, kmax
		!!w(i,k) = 0.
		!!theta_1(i,k) = 0.
		!!theta(i,k) = theta_0(i,k)
	!!END DO
!!END DO

!!CALL set_area_pi
!!DO i = imin, imax
	!!DO k = kmin, kmax
		!!pi_1(i,k) = 0.
	!!END DO
!!END DO
!!=================================================
!END SUBROUTINE initiate_Sm
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
! theta_0_pi, theta_0_u, theta_0_vir
! rho_0_u, rho_0_w, rho_0_vir
!=================================================
SUBROUTINE initiate_basic_state(theta_0,pi_0,rho_0)
IMPLICIT NONE
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_0
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_0
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: rho_0
!-------------------------------------------------
REAL(kd), PARAMETER :: Ts = 300.         ! (K)
REAL(kd), PARAMETER :: N0 = 0.01         ! (s-1)
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: pi_0_w
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! Undefine variations.
!-------------------------------------------------
CALL debug_undef_all(theta_0,pi_0,rho_0)
CALL debug_undef_all(theta_0_pi,theta_0_u,theta_0_vir)
CALL debug_undef_all(rho_0_u,rho_0_w,rho_0_vir)

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
		!ELSE IF (RunCase == 4) THEN
			!theta_0(i,k) = (Ts - 20.)*EXP(N0*N0/g*z_hat(i,k))
			!pi_0_w(i,k) = 1. + g*g/Cp/N0/N0/Ts*(EXP(-N0*N0*z_hat(i,k)/g) - EXP(-N0*N0*zs(i)/g))
		ELSE
			STOP "WRONG RunCase!!!"
		END IF
		rho_0_w(i,k) = p0/Rd/theta_0(i,k)*pi_0_w(i,k)*pi_0_w(i,k)**(Cp/Rd)
	END DO
END DO

CALL set_area_vir
DO i = imin, imax
	DO k = kmin, kmax
		theta_0_vir(i,k) = theta_0(imin+1,k)
		rho_0_vir(i,k) = rho_0_w(imin+1,k)
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
		!ELSE IF (RunCase == 4) THEN
			!theta_0_pi(i,k) = (Ts - 20.)*EXP(N0*N0/g*z_hat_pi(i,k))
			!pi_0(i,k) = 1. + g*g/Cp/N0/N0/Ts*(EXP(-N0*N0*z_hat_pi(i,k)/g) - EXP(-N0*N0*zs(i)/g))
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
END MODULE sp_module_initiate
!=================================================

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
SUBROUTINE initiate_dc(u,v,w,pi_1,theta,theta_0,theta_1,rho_0)
IMPLICIT NONE
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: u        ! wind speed along x-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: v        ! wind speed along y-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: w        ! wind speed along z-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_1     ! pi'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_0  ! theta = theta_0 + theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_1  ! theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: rho_0    ! density
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime) :: xx      ! distance on u-grid along x-axis (m)
REAL(preci), DIMENSION(ims:ime) :: xpi     ! distance on pi-grid along x-axis (m)
REAL(preci), DIMENSION(kms:kme) :: zz      ! height on w-grid along z-axis (m)
REAL(preci), DIMENSION(kms:kme) :: zpi     ! height on pi-grid along z-axis (m)
!-------------------------------------------------
REAL(preci), PARAMETER :: x_c = 25.55*1000.  ! (m)
REAL(preci), PARAMETER :: z_c = 3.0*1000.     ! (m)
REAL(preci), PARAMETER :: r_x = 4*1000.  ! (m)
REAL(preci), PARAMETER :: r_z = 2*1000.  ! (m)
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: pi_0 
REAL(preci), DIMENSION(ims:ime,kms:kme) :: pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_1_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_0_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_pi
!-------------------------------------------------
REAL(preci) :: L
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL debug_undef_all(pi_0,pi,theta_1_pi,theta_0,pi,theta_pi)
xx = undef
xpi = undef
zz = undef
zpi = undef
!=================================================
! 0. Calculate xx(its:ite), xpi(its+1:ite), zz(kts:kte+1), zpi(kts:kte).
!-------------------------------------------------
xx(its) = 0.
DO i = its + 1, ite
	xx(i) = dx*(i - its)
	xpi(i) = (xx(i-1) + xx(i))/2.
END DO
!WRITE(*,*) xpi
!WRITE(*,*) xx

DO k = kts, kte + 1
	zz(k) = dz*(k - kts)
END DO
DO k = kts, kte
	zpi(k) = (zz(k) + zz(k + 1))/2.
END DO
!WRITE(*,*) zz
!WRITE(*,*) zpi
!=================================================
! 1. u, w, pi_1, theta, theta_0, theta_1, rho_0
!-------------------------------------------------
! u-grid
u(its+1:ite-1,kts:kte) = 0.              ! <= I want this.

! w-grid
w(its+1:ite,kts+1:kte) = 0.              ! <= I want this.
!theta_0(its+1:ite,kts+1:kte) = Ts        ! <= I want this.
theta_0(its+1:ite,kts:kte+1) = Ts        ! <= I want this. Update boundary.

! Vertical boundary condition
!theta_0(:,kts-1) = theta_0(:,kts+1) ! Ptheta_0Pz = 0.
!theta_0(:,kts) = theta_0(:,kts+1) ! Ptheta_0Pz = 0.
!theta_0(:,kte+2) = theta_0(:,kte) ! Ptheta_0Pz = 0.
!theta_0(:,kte+1) = theta_0(:,kte) ! Ptheta_0Pz = 0.
! Lateral boundary condition
!theta_0(its,:) = theta_0(its+1,:) ! Ptheta_0Px = 0.
!theta_0(ite+1,:) = theta_0(ite,:) ! Ptheta_0Px = 0.

!DO k = kts + 1, kte
DO k = kts, kte + 1 ! Update boundary.
	DO i = its + 1, ite
		L = SQRT((xpi(i) - x_c)*(xpi(i) - x_c)/r_x/r_x + (zz(k) - z_c)*(zz(k) - z_c)/r_z/r_z)
		IF (L <= 1.) THEN
			theta_1(i,k) = - 15./2.*(COS(PI_math*L) + 1.)  ! <= I want this.
		ELSE
			theta_1(i,k) = 0.
		END IF
		theta(i,k) = theta_0(i,k) + theta_1(i,k)   ! <= I want this.
	END DO
END DO
!theta_1(its+1:ite,kte+1) = 0. ! Update boundary.

! pi-grid
pi_0(its+1:ite,kts) = 1.
theta_0_pi(its:ite+1,kts:kte+1) = Ts
rho_0(its+1:ite,kts) = p0/Rd/theta_0(its+1:ite,kts)*pi_0(its+1:ite,kts)*pi_0(its+1:ite,kts)**(Cp/Rd)     ! <= I want this.
DO k = kts + 1, kte
	DO i = its + 1, ite
		pi_0(i,k) = pi_0(i,k-1) - g/Cp/theta_0_pi(i,k)*dz                ! <= I want this.
		rho_0(i,k) = p0/Rd/theta_0_pi(i,k)*pi_0(i,k)*pi_0(i,k)**(Cp/Rd)     ! <= I want this.
	END DO
END DO
pi(its+1:ite,kte) = pi_0(its+1:ite,kte)
pi_1(its+1:ite,kte) = 0.
DO k = kte - 1, kts, - 1
	DO i = its + 1, ite
		L = SQRT((xpi(i) - x_c)*(xpi(i) - x_c)/r_x/r_x + (zpi(k) - z_c)*(zpi(k) - z_c)/r_z/r_z)
		IF (L <= 1.) THEN
			theta_1_pi(i,k) = - 15./2.*(COS(PI_math*L) + 1.)  ! <= I want this.
		ELSE
			theta_1_pi(i,k) = 0.
		END IF
		theta_pi(i,k) = theta_0_pi(i,k) + theta_1_pi(i,k)
		pi(i,k) = pi(i,k+1) + g/Cp/theta_pi(i,k)*dz                      ! <= I want this.
		pi_1(i,k) = pi(i,k) - pi_0(i,k)                                  ! <= I want this.
	END DO
END DO

!WRITE(*,*) pi(256,:)
!WRITE(*,*) pi_0(256,:)
!WRITE(*,*) rho_0(256,:)

! Vertical boundary condition
!rho_0(:,kts-1) = rho_0(:,kts) ! Prho_0Pz = 0.
!rho_0(:,kte+1) = rho_0(:,kte) ! Prho_0Pz = 0.
! Lateral boundary condition
!rho_0(its,:) = rho_0(its+1,:) ! Prho_0Px = 0.
!rho_0(ite+1,:) = rho_0(ite,:) ! Prho_0Px = 0.

!=================================================
!CALL debug_ascii_output(theta_1)
!CALL debug_ascii_output(theta_1_pi)
!CALL debug_ascii_output(pi_1)
!=================================================
END SUBROUTINE initiate_dc
!=================================================

!=================================================
! Initiate thermal bubble case.
!=================================================
SUBROUTINE initiate_tb(u,v,w,pi_1,theta,theta_0,theta_1,rho_0)
IMPLICIT NONE
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: u        ! wind speed along x-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: v        ! wind speed along y-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: w        ! wind speed along z-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_1     ! pi'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_0  ! theta = theta_0 + theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_1  ! theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: rho_0    ! density
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime) :: xx      ! distance on u-grid along x-axis (m)
REAL(preci), DIMENSION(ims:ime) :: xpi     ! distance on pi-grid along x-axis (m)
REAL(preci), DIMENSION(kms:kme) :: zz      ! height on w-grid along z-axis (m)
REAL(preci), DIMENSION(kms:kme) :: zpi     ! height on pi-grid along z-axis (m)
!-------------------------------------------------
REAL(preci), PARAMETER :: x_c = 10.0*1000. ! (m)
REAL(preci), PARAMETER :: z_c = 2.0*1000.  ! (m)
REAL(preci), PARAMETER :: R = 2.0*1000.    ! (m)

!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: pi_0 
REAL(preci), DIMENSION(ims:ime,kms:kme) :: pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_1_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_0_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_pi
!-------------------------------------------------
REAL(preci) :: L
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL debug_undef_all(pi_0,pi,theta_1_pi,theta_0,pi,theta_pi)
xx = undef
xpi = undef
zz = undef
zpi = undef
!=================================================
! 0. Calculate xx(its:ite), xpi(its+1:ite), zz(kts:kte+1), zpi(kts:kte).
!-------------------------------------------------
xx(its) = 0.
DO i = its + 1, ite
	xx(i) = dx*(i - its)
	xpi(i) = (xx(i-1) + xx(i))/2.
END DO

DO k = kts, kte + 1
	zz(k) = dz*(k - kts)
END DO
DO k = kts, kte
	zpi(k) = (zz(k) + zz(k + 1))/2.
END DO
!=================================================
! 1. u, w, pi_1, theta, theta_0, theta_1, rho_0
!-------------------------------------------------
! u-grid
u(its+1:ite-1,kts:kte) = 0.              ! <= I want this.

! w-grid
w(its+1:ite,kts+1:kte) = 0.              ! <= I want this.
!theta_0(its+1:ite,kts+1:kte) = Ts        ! <= I want this.
theta_0(its+1:ite,kts:kte+1) = Ts        ! <= I want this. Update boundary.

DO k = kts, kte + 1 ! Update boundary.
	DO i = its + 1, ite
		L = SQRT((xpi(i) - x_c)*(xpi(i) - x_c) + (zz(k) - z_c)*(zz(k) - z_c))
		theta_1(i,k) = 2.*MAX(0.,1. - L/R)
		theta(i,k) = theta_0(i,k) + theta_1(i,k)   ! <= I want this.
	END DO
END DO
!theta_1(its+1:ite,kte+1) = 0. ! Update boundary.

! pi-grid
pi_0(its+1:ite,kts) = 1.
theta_0_pi(its:ite+1,kts:kte+1) = Ts
rho_0(its+1:ite,kts) = p0/Rd/theta_0(its+1:ite,kts)*pi_0(its+1:ite,kts)*pi_0(its+1:ite,kts)**(Cp/Rd)     ! <= I want this.
DO k = kts + 1, kte
	DO i = its + 1, ite
		pi_0(i,k) = pi_0(i,k-1) - g/Cp/theta_0_pi(i,k)*dz                ! <= I want this.
		rho_0(i,k) = p0/Rd/theta_0_pi(i,k)*pi_0(i,k)*pi_0(i,k)**(Cp/Rd)     ! <= I want this.
	END DO
END DO
pi(its+1:ite,kte) = pi_0(its+1:ite,kte)
pi_1(its+1:ite,kte) = 0.
DO k = kte - 1, kts, - 1
	DO i = its + 1, ite
		L = SQRT((xpi(i) - x_c)*(xpi(i) - x_c) + (zpi(k) - z_c)*(zpi(k) - z_c))
		theta_1_pi(i,k) = 2.*MAX(0.,1. - L/R)
		theta_pi(i,k) = theta_0_pi(i,k) + theta_1_pi(i,k)
		pi(i,k) = pi(i,k+1) + g/Cp/theta_pi(i,k)*dz                      ! <= I want this.
		pi_1(i,k) = pi(i,k) - pi_0(i,k)                                  ! <= I want this.
	END DO
END DO
!=================================================
END SUBROUTINE initiate_tb
!=================================================

!=================================================
! Initiate inertia gravity waves.
!=================================================
SUBROUTINE initiate_igw(u,v,w,pi_1,theta,theta_0,theta_1,rho_0)
IMPLICIT NONE
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: u        ! wind speed along x-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: v        ! wind speed along y-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: w        ! wind speed along z-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: pi_1     ! pi'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_0  ! theta = theta_0 + theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: theta_1  ! theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT) :: rho_0    ! density
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime) :: xx      ! distance on u-grid along x-axis (m)
REAL(preci), DIMENSION(ims:ime) :: xpi     ! distance on pi-grid along x-axis (m)
REAL(preci), DIMENSION(kms:kme) :: zz      ! height on w-grid along z-axis (m)
REAL(preci), DIMENSION(kms:kme) :: zpi     ! height on pi-grid along z-axis (m)
!-------------------------------------------------
REAL(preci), PARAMETER :: x_c = 100.0*1000. ! (m)
REAL(preci), PARAMETER :: H = 10.0*1000.    ! (m)
REAL(preci), PARAMETER :: a = 5.0*1000.     ! (m)
REAL(preci), PARAMETER :: N0 = 0.01         ! (s-1)
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: pi_0 
REAL(preci), DIMENSION(ims:ime,kms:kme) :: pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_1_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_0_pi
REAL(preci), DIMENSION(ims:ime,kms:kme) :: theta_pi
!-------------------------------------------------
REAL(preci) :: L
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL debug_undef_all(pi_0,pi,theta_1_pi,theta_0_pi,theta_pi)
xx = undef
xpi = undef
zz = undef
zpi = undef
!=================================================
! 0. Calculate xx(its:ite), xpi(its+1:ite), zz(kts:kte+1), zpi(kts:kte).
!-------------------------------------------------
xx(its) = 0.
DO i = its + 1, ite
	xx(i) = dx*(i - its)
	xpi(i) = (xx(i-1) + xx(i))/2.
END DO

DO k = kts, kte + 1
	zz(k) = dz*(k - kts)
END DO
DO k = kts, kte
	zpi(k) = (zz(k) + zz(k + 1))/2.
END DO
!=================================================
! 1. u, w, pi_1, theta, theta_0, theta_1, rho_0
!-------------------------------------------------
! u-grid
u(its+1:ite-1,kts:kte) = 20.              ! <= I want this.

! w-grid
w(its+1:ite,kts+1:kte) = 0.              ! <= I want this.

DO k = kts, kte + 1 ! Update boundary.
	DO i = its + 1, ite
		theta_0(i,k) = Ts*EXP(N0*N0/g*zz(k))        ! <= I want this. Update boundary.
		L = SIN(PI_math*zz(k)/H)/(1. + (xx(i) - x_c)*(xx(i) - x_c)/a/a)
		theta_1(i,k) = 0.01*L
		theta(i,k) = theta_0(i,k) + theta_1(i,k)   ! <= I want this.
	END DO
END DO

! pi-grid
DO k = kts, kte
	DO i = its + 1, ite
		theta_0_pi(i,k) = Ts*EXP(N0*N0/g*zpi(k))        ! <= I want this. Update boundary.
		pi_0(i,k) = 1 + g*g/Cp/N0/N0/Ts*(EXP(- N0*N0*zpi(k)/g))
		rho_0(i,k) = p0/Rd/theta_0_pi(i,k)*pi_0(i,k)*pi_0(i,k)**(Cp/Rd)     ! <= I want this.
	END DO
END DO
CALL debug_ascii_output(pi_0)
pi(its+1:ite,kte) = pi_0(its+1:ite,kte)
pi_1(its+1:ite,kte) = 0.
DO k = kte - 1, kts, - 1
	DO i = its + 1, ite
		L = SIN(PI_math*zpi(k)/H)/(1. + (xpi(i) - x_c)*(xpi(i) - x_c)/a/a)
		theta_1_pi(i,k) = 0.01*L
		theta_pi(i,k) = theta_0_pi(i,k) + theta_1_pi(i,k)
		pi(i,k) = pi(i,k+1) + g/Cp/theta_pi(i,k)*dz                      ! <= I want this.
		pi_1(i,k) = pi(i,k) - pi_0(i,k)                                  ! <= I want this.
	END DO
END DO
!CALL debug_ascii_output(pi_0)
!CALL debug_ascii_output(pi_1)
!CALL debug_ascii_output(pi)
!WRITE(*,*) "----------------u---------------------"
!CALL debug_test_boundary(u)
!WRITE(*,*) "----------------w---------------------"
!CALL debug_test_boundary(w)
!WRITE(*,*) "----------------theta---------------------"
!CALL debug_test_boundary(theta)
!WRITE(*,*) "----------------theta_0---------------------"
!CALL debug_test_boundary(theta_0)
!WRITE(*,*) "----------------theta_1---------------------"
!CALL debug_test_boundary(theta_1)
!WRITE(*,*) "----------------pi_1---------------------"
!CALL debug_test_boundary(pi_1)
!WRITE(*,*) "----------------rho_0---------------------"
!CALL debug_test_boundary(rho_0)
!=================================================
END SUBROUTINE initiate_igw
!=================================================

!=================================================
END MODULE sp_module_initiate
!=================================================

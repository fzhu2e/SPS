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
USE sp_module_gridvar
USE sp_module_debug
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================
! Initiate density current case.
!=================================================
SUBROUTINE initiate_dc(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE (grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!-------------------------------------------------
REAL(kd), PARAMETER :: x_c = 25.55*1000.  ! (m)
REAL(kd), PARAMETER :: z_c = 3.0*1000.     ! (m)
REAL(kd), PARAMETER :: r_x = 4*1000.  ! (m)
REAL(kd), PARAMETER :: r_z = 2*1000.  ! (m)
REAL(kd) :: L
INTEGER :: i, k
!=================================================
! u, w, pi_1, theta, qc, qv, qr
!-------------------------------------------------
CALL set_area_u
DO k = kmin, kmax
	DO i = imin, imax
		uGrid%u(i,k) = 0
	END DO
END DO

CALL set_area_w
DO k = kmin, kmax
	DO i = imin, imax
		wGrid%w(i,k) = 0
		L = SQRT((wGrid%xx(i) - x_c)*(wGrid%xx(i) - x_c)/r_x/r_x + (wGrid%zz(i,k) - z_c)*(wGrid%zz(i,k) - z_c)/r_z/r_z)
		IF (L <= 1.) THEN
			wGrid%theta_1(i,k) = - 15./2.*(COS(PI_math*L) + 1.)  ! <= I want this.
		ELSE
			wGrid%theta_1(i,k) = 0.
		END IF
		wGrid%theta(i,k) = wGrid%theta_0(i,k) + wGrid%theta_1(i,k)   ! <= I want this.
	END DO
END DO

CALL set_area_pi
DO k = kmin, kmax
	DO i = imin, imax
		piGrid%pi_1(i,k) = 0.
		piGrid%pi(i,k) = piGrid%pi_0(i,k) + piGrid%pi_1(i,k)
	END DO
END DO
!=================================================
END SUBROUTINE initiate_dc
!=================================================

!=================================================
! Initiate thermal bubble case.
!=================================================
SUBROUTINE initiate_tb(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE (grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!-------------------------------------------------
REAL(kd), PARAMETER :: x_c = 10.0*1000. ! (m)
REAL(kd), PARAMETER :: z_c = 2.0*1000.  ! (m)
REAL(kd), PARAMETER :: R = 2.0*1000.    ! (m)
!-------------------------------------------------
REAL(kd) :: L
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL set_area_u
DO k = kmin, kmax
	DO i = imin, imax
		uGrid%u(i,k) = 0.
	END DO
END DO

CALL set_area_w
DO k = kmin, kmax
	DO i = imin, imax
		wGrid%w(i,k) = 0.
		L = SQRT((wGrid%xx(i) - x_c)*(wGrid%xx(i) - x_c) + (wGrid%zz(i,k) - z_c)*(wGrid%zz(i,k) - z_c))
		wGrid%theta_1(i,k) = 2.*MAX(0.,1. - L/R)
		wGrid%theta(i,k) = wGrid%theta_0(i,k) + wGrid%theta_1(i,k)
	END DO
END DO

CALL set_area_pi
DO k = kmin, kmax
	DO i = imin, imax
		piGrid%pi_1(i,k) = 0.
		piGrid%pi(i,k) = piGrid%pi_0(i,k) + piGrid%pi_1(i,k)
	END DO
END DO
!=================================================
END SUBROUTINE initiate_tb
!=================================================

!=================================================
! Initiate inertia gravity waves.
!=================================================
SUBROUTINE initiate_igw(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE (grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!-------------------------------------------------
REAL(kd), PARAMETER :: x_c = 100.0*1000. ! (m)
REAL(kd), PARAMETER :: H = 10.0*1000.    ! (m)
REAL(kd), PARAMETER :: a = 5.0*1000.     ! (m)
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: pi
REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_1_pi
REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_pi
!-------------------------------------------------
REAL(kd) :: L
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL set_area_u
DO k = kmin, kmax
	DO i = imin, imax
		uGrid%u(i,k) = 20.
	END DO
END DO

CALL set_area_w
DO k = kmin, kmax
	DO i = imin, imax
		wGrid%w(i,k) = 0.
		L = SIN(PI_math*wGrid%zz(i,k)/H)/(1. + (wGrid%xx(i) - x_c)*(wGrid%xx(i) - x_c)/a/a)
		wGrid%theta_1(i,k) = 0.01*L
		wGrid%theta(i,k) = wGrid%theta_0(i,k) + wGrid%theta_1(i,k)
	END DO
END DO

CALL set_area_pi
DO k = kmin, kmax
	DO i = imin, imax
		piGrid%pi_1(i,k) = 0.
		piGrid%pi(i,k) = piGrid%pi_0(i,k) + piGrid%pi_1(i,k)
	END DO
END DO
!=================================================
END SUBROUTINE initiate_igw
!=================================================

!=================================================
! Initiate Schar mountain case.
!=================================================
SUBROUTINE initiate_Sm(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE (grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!-------------------------------------------------
INTEGER :: i, k 
!=================================================
CALL set_area_u
DO k = kmin, kmax
	DO i = imin, imax
		uGrid%u(i,k) = 10.
	END DO
END DO

CALL set_area_w
DO k = kmin, kmax
	DO i = imin, imax
		wGrid%w(i,k) = 0.
		wGrid%theta_1(i,k) = 0.
		wGrid%theta(i,k) = wGrid%theta_0(i,k)
	END DO
END DO

CALL set_area_pi
DO k = kmin, kmax
	DO i = imin, imax
		piGrid%pi_1(i,k) = 0.
		piGrid%pi(i,k) = piGrid%pi_0(i,k) + piGrid%pi_1(i,k)
	END DO
END DO
!=================================================
END SUBROUTINE initiate_Sm
!=================================================

!=================================================
! Initiate wet bubble case.
!=================================================
SUBROUTINE initiate_wb(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE (grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!-------------------------------------------------
REAL(kd), PARAMETER :: x_c = 250. ! (m)
REAL(kd), PARAMETER :: z_c = 125.  ! (m)
REAL(kd), PARAMETER :: R = 75    ! (m)
!-------------------------------------------------
REAL(kd) :: L
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL set_area_u
DO k = kmin, kmax
	DO i = imin, imax
		uGrid%u(i,k) = 0.
	END DO
END DO

CALL set_area_w
DO k = kmin, kmax
	DO i = imin, imax
		wGrid%w(i,k) = 0.
		wGrid%theta_1(i,k) = 0.
		wGrid%theta(i,k) = wGrid%theta_0(i,k) + wGrid%theta_1(i,k)
		L = SQRT((wGrid%xx(i) - x_c)*(wGrid%xx(i) - x_c) + (wGrid%zz(i,k) - z_c)*(wGrid%zz(i,k) - z_c))
		!wGrid%qv(i,k) = 0.
		wGrid%qv(i,k) = 0.01*MAX(0.01,1. - L/R)
		wGrid%qc(i,k) = 0.
		wGrid%qr(i,k) = 0.
		wGrid%qi(i,k) = 0.
		wGrid%qs(i,k) = 0.
		wGrid%qg(i,k) = 0.
	END DO
END DO

CALL set_area_pi
DO k = kmin, kmax
	DO i = imin, imax
		piGrid%pi_1(i,k) = 0.
		piGrid%pi(i,k) = piGrid%pi_0(i,k) + piGrid%pi_1(i,k)
	END DO
END DO
!=================================================
END SUBROUTINE initiate_wb
!=================================================

!=================================================
! Initiate thunderstorm case.
!=================================================
SUBROUTINE initiate_th(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE (grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!-------------------------------------------------
REAL(kd), PARAMETER :: x_c = 25.0*1000. ! (m)
REAL(kd), PARAMETER :: z_c = 4.0*1000.  ! (m)
REAL(kd), PARAMETER :: R = 4.0*1000.    ! (m)
REAL(kd), PARAMETER :: S = 15.0*1000.    ! (m)
!-------------------------------------------------
REAL(kd) :: L
REAL(kd), DIMENSION(nz) :: qv
!-------------------------------------------------
INTEGER :: i, k
!=================================================
OPEN(1, FILE="./input/qv.input", STATUS='old')
DO k = 1, nz
	READ(1,*) qv(k)
END DO
CLOSE(1)

CALL set_area_pi
DO k = kmin, kmax
	DO i = ims, ime
		piGrid%qv(i,k) = qv(k)
	END DO
END DO

CALL set_area_w
wGrid%qv(:,kmin) = 1.400E-02
wGrid%qv(:,kmax) = 9.400E-05
DO k = kmin+1, kmax-1
	wGrid%qv(:,k) = (piGrid%qv(:,k-1) + piGrid%qv(:,k))/2.
END DO

!CALL set_area_w
!DO k = kmin, kmax
	!DO i = imin, imax
		!!wGrid%qv = 0.
		!IF (wGrid%zz(i,k) <= S) THEN
			!wGrid%qv(i,k) = 0.014*(1.-SIN(wGrid%zz(i,k)/S*PI_math/2.)**4)
		!ELSE
			!wGrid%qv(i,k) = 0.
		!END IF
	!END DO
!END DO
!WRITE(*,*) wGrid%qv(5,:)*1000
!CALL debug_SFSG

CALL set_area_u
DO k = kmin, kmax
	DO i = imin, imax
		!uGrid%u_0(i,k) = 0.
		uGrid%u_0(i,k) = - 1.2*MAX(0., REAL(11 - k))
		uGrid%u(i,k) = uGrid%u_0(i,k)
	END DO
END DO

CALL set_area_w
DO k = kmin, kmax
	DO i = imin, imax
		wGrid%w(i,k) = 0.
		L = SQRT((wGrid%xx(i) - x_c)*(wGrid%xx(i) - x_c) + (wGrid%zz(i,k) - z_c)*(wGrid%zz(i,k) - z_c))
		IF (L/R <= 1.) THEN
			wGrid%theta_1(i,k) = 3.*COS(.5*PI_math*L/R)**2  ! <= I want this.
		ELSE
			wGrid%theta_1(i,k) = 0.
		END IF
		wGrid%theta(i,k) = wGrid%theta_0(i,k) + wGrid%theta_1(i,k)
	END DO
END DO

CALL set_area_pi
DO k = kmin, kmax
	DO i = imin, imax
		piGrid%pi_1(i,k) = 0.
		piGrid%pi(i,k) = piGrid%pi_0(i,k) + piGrid%pi_1(i,k)
	END DO
END DO
!=================================================
END SUBROUTINE initiate_th
!=================================================

!=================================================
! Initiate Grid Position
!=================================================
SUBROUTINE initiate_grid(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE (grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! xx, zeta
!-------------------------------------------------
DO i = ims, ime
	uGrid%xx(i) = dx*(i - its)
	virGrid%xx(i) = uGrid%xx(i)
	piGrid%xx(i) = uGrid%xx(i) - dx/2.
	wGrid%xx(i) = piGrid%xx(i)
END DO

DO k = kms, kme
	wGrid%zeta(k) = dz*(k - kts)
	virGrid%zeta(k) = wGrid%zeta(k)
	piGrid%zeta(k) = wGrid%zeta(k) + dz/2.
	uGrid%zeta(k) = piGrid%zeta(k)
END DO
!=================================================
END SUBROUTINE initiate_grid
!=================================================

!=================================================
! Initiate Terrain
!=================================================
SUBROUTINE initiate_terrain(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE (grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!-------------------------------------------------
! for Sm
REAL(kd), PARAMETER :: h0 = 250.0            ! (m)
REAL(kd), PARAMETER :: a0 = 5.0*1000.        ! (m)
REAL(kd), PARAMETER :: x_c = 25.0*1000.        ! (m)
REAL(kd), PARAMETER :: lambda0 = 4.0*1000.   ! (m)
!-------------------------------------------------
INTEGER :: i, k
REAL(kd), PARAMETER :: s = 3000.  ! km
!=================================================
! zs, zz, G, H
!-------------------------------------------------
IF (RunCase /= 4) THEN
	DO i = ims, ime
		uGrid%zs(i) = 0.
		piGrid%zs(i) = 0.
	END DO

ELSE IF (RunCase == 4) THEN
	DO i = ims, ime
		uGrid%zs(i) = h0*EXP(-((uGrid%xx(i) - x_c)/a0)**2)*COS(PI_math*(uGrid%xx(i) - x_c)/lambda0)**2
		!uGrid%PzsPx(i) = - 2*h0*EXP(-((uGrid%xx(i) - x_c)/a0)**2)*COS(PI_math*(uGrid%xx(i) - x_c)/lambda0)*((uGrid%xx(i) - x_c)/a0**2*COS(PI_math*(uGrid%xx(i) - x_c)/lambda0) + PI_math/lambda0*SIN(PI_math*(uGrid%xx(i) - x_c)/lambda0))
	END DO
	
	DO i = ims, ime
		piGrid%zs(i) = h0*EXP(-((piGrid%xx(i) - x_c)/a0)**2)*COS(PI_math*(piGrid%xx(i) - x_c)/lambda0)**2
		!piGrid%PzsPx(i) = - 2*h0*EXP(-((piGrid%xx(i) - x_c)/a0)**2)*COS(PI_math*(piGrid%xx(i) - x_c)/lambda0)*((piGrid%xx(i) - x_c)/a0**2*COS(PI_math*(piGrid%xx(i) - x_c)/lambda0) + PI_math/lambda0*SIN(PI_math*(piGrid%xx(i) - x_c)/lambda0))
	END DO

ELSE
	STOP "WRONG RunCase!!!"
END IF

virGrid%zs = uGrid%zs
wGrid%zs = piGrid%zs

DO i = ims, ime - 1
	uGrid%PzsPx(i) = (piGrid%zs(i+1) - piGrid%zs(i))/dx
END DO

DO i = ims + 1, ime
	piGrid%PzsPx(i) = (piGrid%zs(i) - piGrid%zs(i-1))/dx
END DO

virGrid%PzsPx = uGrid%PzsPx
wGrid%PzsPx = piGrid%PzsPx

IF (VertCoords == 1) THEN
	DO k = kms, kme
		uGrid%b(k) = 1. - uGrid%zeta(k)/ztop
		wGrid%b(k) = 1. - wGrid%zeta(k)/ztop
		piGrid%b(k) = 1. - piGrid%zeta(k)/ztop
		virGrid%b(k) = 1. - virGrid%zeta(k)/ztop
	END DO

	DO i = ims, ime
		uGrid%H(i) = ztop/(ztop - uGrid%zs(i))
		wGrid%H(i) = ztop/(ztop - wGrid%zs(i))
		piGrid%H(i) = ztop/(ztop - piGrid%zs(i))
		virGrid%H(i) = ztop/(ztop - virGrid%zs(i))
	END DO

	DO k = kms, kme
		DO i = ims + 1, ime - 1
			uGrid%G(i,k) = (uGrid%zeta(k) - ztop)/(ztop - uGrid%zs(i))*uGrid%PzsPx(i)
			wGrid%G(i,k) = (wGrid%zeta(k) - ztop)/(ztop - uGrid%zs(i))*wGrid%PzsPx(i)
			piGrid%G(i,k) = (piGrid%zeta(k) - ztop)/(ztop - piGrid%zs(i))*piGrid%PzsPx(i)
			virGrid%G(i,k) = (virGrid%zeta(k) - ztop)/(ztop - virGrid%zs(i))*virGrid%PzsPx(i)
		END DO
	END DO

ELSE IF (VertCoords == 2) THEN
	DO k = kms, kme
		uGrid%b(k) = SINH((ztop - uGrid%zeta(k))/s)/SINH(ztop/s)
		wGrid%b(k) = SINH((ztop - wGrid%zeta(k))/s)/SINH(ztop/s)
		piGrid%b(k) = SINH((ztop - piGrid%zeta(k))/s)/SINH(ztop/s)
		virGrid%b(k) = SINH((ztop - virGrid%zeta(k))/s)/SINH(ztop/s)
	END DO

ELSE
	STOP "WRONG VertCoord!!!"
END IF

DO k = kms, kme
	DO i = ims, ime
		uGrid%zz(i,k) = uGrid%zeta(k) + uGrid%zs(i)*uGrid%b(k)
		wGrid%zz(i,k) = wGrid%zeta(k) + wGrid%zs(i)*wGrid%b(k)
		piGrid%zz(i,k) = piGrid%zeta(k) + piGrid%zs(i)*piGrid%b(k)
		virGrid%zz(i,k) = virGrid%zeta(k) + virGrid%zs(i)*virGrid%b(k)
	END DO
END DO
!=================================================
END SUBROUTINE initiate_terrain
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
SUBROUTINE initiate_basic_state(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
!-------------------------------------------------
TYPE (grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!-------------------------------------------------
REAL(kd) :: Ts = 300.
REAL(kd) :: N0 = 0.01         ! (s-1)

REAL(kd), DIMENSION(nz) :: theta,pressure ! sounding input in thunderstorm case
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! theta_0, pi_0, rho_0
!-------------------------------------------------

! theta_0, pi_0
!IF (RunCase == 1 .OR. RunCase == 2 .OR. RunCase == 6) THEN
IF (RunCase == 1 .OR. RunCase == 2) THEN
	CALL set_area_u
	DO k = kmin, kmax
		DO i = imin, imax
			uGrid%theta_0(i,k) = Ts
			uGrid%pi_0(i,k) = 1. - g*uGrid%zz(i,k)/2./Cp/Ts
		END DO
	END DO

	CALL set_area_w
	DO k = kmin, kmax
		DO i = imin, imax
			wGrid%theta_0(i,k) = Ts
			wGrid%pi_0(i,k) = 1. - g*wGrid%zz(i,k)/2./Cp/Ts
		END DO
	END DO

	CALL set_area_pi
	DO k = kmin, kmax
		DO i = imin, imax
			piGrid%theta_0(i,k) = Ts
			piGrid%pi_0(i,k) = 1. - g*piGrid%zz(i,k)/2./Cp/Ts
		END DO
	END DO

	CALL set_area_vir
	DO k = kmin, kmax
		DO i = imin, imax
			virGrid%theta_0(i,k) = Ts
			virGrid%pi_0(i,k) = 1. - g*virGrid%zz(i,k)/2./Cp/Ts
		END DO
	END DO


ELSE IF (RunCase /= 1 .AND. RunCase /= 2) THEN

	IF (RunCase == 5) THEN
		Ts = 270.
	END IF

	CALL set_area_u
	DO k = kmin, kmax
		DO i = imin, imax
			uGrid%theta_0(i,k) = Ts*EXP(N0*N0/g*uGrid%zz(i,k))
			uGrid%pi_0(i,k) = 1. + g*g/Cp/N0/N0/Ts*(EXP(-N0*N0*uGrid%zz(i,k)/g) - EXP(-N0*N0*uGrid%zs(i)/g))
		END DO
	END DO

	CALL set_area_w
	DO k = kmin, kmax
		DO i = imin, imax
			wGrid%theta_0(i,k) = Ts*EXP(N0*N0/g*wGrid%zz(i,k))
			wGrid%pi_0(i,k) = 1. + g*g/Cp/N0/N0/Ts*(EXP(-N0*N0*wGrid%zz(i,k)/g) - EXP(-N0*N0*wGrid%zs(i)/g))
		END DO
	END DO

	CALL set_area_pi
	DO k = kmin, kmax
		DO i = imin, imax
			piGrid%theta_0(i,k) = Ts*EXP(N0*N0/g*piGrid%zz(i,k))
			piGrid%pi_0(i,k) = 1. + g*g/Cp/N0/N0/Ts*(EXP(-N0*N0*piGrid%zz(i,k)/g) - EXP(-N0*N0*piGrid%zs(i)/g))
		END DO
	END DO

	CALL set_area_vir
	DO k = kmin, kmax
		DO i = imin, imax
			virGrid%theta_0(i,k) = Ts*EXP(N0*N0/g*virGrid%zz(i,k))
			virGrid%pi_0(i,k) = 1. + g*g/Cp/N0/N0/Ts*(EXP(-N0*N0*virGrid%zz(i,k)/g) - EXP(-N0*N0*virGrid%zs(i)/g))
		END DO
	END DO

	IF (RunCase == 6) THEN
		OPEN(1, FILE="./input/theta.input", STATUS='old')
		DO k = 1, nz
			READ(1,*) theta(k)
		END DO
		CLOSE(1)
		OPEN(1, FILE="./input/pressure.input", STATUS='old')
		DO k = 1, nz
			READ(1,*) pressure(k)
		END DO
		CLOSE(1)
		DO k = kmin, kmax
			DO i = ims, ime
				piGrid%theta_0(i,k) = theta(k)
				piGrid%pi_0(i,k) = (pressure(k)/p0)**(Rd/Cp)
			END DO
		END DO
		CALL set_area_w
		wGrid%theta_0(:,kmin) = 3.000E+02
		wGrid%theta_0(:,kmax) = 4.978E+02
		wGrid%pi_0(:,kmin) = 1
		wGrid%pi_0(:,kmax) = 0.4392
		DO k = kmin+1, kmax-1
			wGrid%theta_0(:,k) = (piGrid%theta_0(:,k-1) + piGrid%theta_0(:,k))/2.
			wGrid%pi_0(:,k) = (piGrid%pi_0(:,k-1) + piGrid%pi_0(:,k))/2.
		END DO
		virGrid%theta_0 = wGrid%theta_0
		uGrid%theta_0 = piGrid%theta_0
		virGrid%pi_0 = wGrid%pi_0
		uGrid%pi_0 = piGrid%pi_0
	END IF

ELSE
	STOP "WRONG RunCase!!!"
END IF

! rho_0
CALL set_area_u
DO k = kmin, kmax
	DO i = imin, imax
		uGrid%rho_0(i,k) = p0/Rd/uGrid%theta_0(i,k)*uGrid%pi_0(i,k)*uGrid%pi_0(i,k)**(Cp/Rd)
	END DO
END DO

CALL set_area_w
DO k = kmin, kmax
	DO i = imin, imax
		wGrid%rho_0(i,k) = p0/Rd/wGrid%theta_0(i,k)*wGrid%pi_0(i,k)*wGrid%pi_0(i,k)**(Cp/Rd)
	END DO
END DO

CALL set_area_pi
DO k = kmin, kmax
	DO i = imin, imax
		piGrid%rho_0(i,k) = p0/Rd/piGrid%theta_0(i,k)*piGrid%pi_0(i,k)*piGrid%pi_0(i,k)**(Cp/Rd)
	END DO
END DO

CALL set_area_vir
DO k = kmin, kmax
	DO i = imin, imax
		virGrid%rho_0(i,k) = p0/Rd/virGrid%theta_0(i,k)*virGrid%pi_0(i,k)*virGrid%pi_0(i,k)**(Cp/Rd)
	END DO
END DO
!=================================================
END SUBROUTINE initiate_basic_state
!=================================================

!=================================================
END MODULE sp_module_initiate
!=================================================

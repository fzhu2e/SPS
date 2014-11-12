!=================================================
! Super-Parametertization System (SPS)
!-------------------------------------------------
! Version: 0.2
! Author: Feng Zhu
! Email: zhuf.atmos@gmail.com
! Date: 2014-06-12 18:18:45
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
! u, w, pi_1, theta, qv
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
		wGrid%w(i,k) = 0
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
		wGrid%w(i,k) = 0
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
		wGrid%w(i,k) = 0
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
		wGrid%w(i,k) = 0
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
!REAL(kd), PARAMETER :: x_c = 37.5*1000. ! (m)
REAL(kd), PARAMETER :: z_c = 4.0*1000.  ! (m)
REAL(kd), PARAMETER :: R = 4.0*1000.    ! (m)
REAL(kd), PARAMETER :: DeltaTheta = 3.
REAL(kd), PARAMETER :: qv0 = 14.e-3
REAL(kd), PARAMETER :: S = 5.e3    ! (m)
!-------------------------------------------------
REAL(kd) :: L
REAL(kd), DIMENSION(nz) :: qv
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! qv
!-------------------------------------------------
!OPEN(1, FILE="./input/qv.input", STATUS='old')
!DO k = 1, nz
	!READ(1,*) qv(k)
!END DO
!CLOSE(1)

!CALL set_area_pi
!DO k = kmin, kmax
	!DO i = ims, ime
		!piGrid%qv(i,k) = qv(k)
	!END DO
!END DO

!CALL set_area_w
!wGrid%qv(:,kmin) = 1.400E-02
!wGrid%qv(:,kmax) = 9.400E-05
!DO k = kmin+1, kmax-1
	!wGrid%qv(:,k) = (piGrid%qv(:,k-1) + piGrid%qv(:,k))/2.
!END DO

!-----------------------------------------------------------

CALL set_area_w
DO k = kmin, kmax
	DO i = imin, imax
		IF (wGrid%zz(i,k) <= S) THEN
			wGrid%qv(i,k) = qv0*(1. - SIN(wGrid%zz(i,k)/S*PI_math/2.)**4)
			!wGrid%qv(i,k) = qv0*(1. - SIN(wGrid%zz(i,k)/S*PI_math/2.)**3)
			!wGrid%qv(i,k) = qv0*(1. - SIN(wGrid%zz(i,k)/S*PI_math/2.)**2)
			!wGrid%qv(i,k) = qv0*(1. - SIN(wGrid%zz(i,k)/S*PI_math/2.))
		ELSE
			wGrid%qv(i,k) = 0.
		END IF
	END DO
END DO
!=================================================

CALL set_area_u
DO k = kmin, kmax
	DO i = imin, imax
		!uGrid%u_0(i,k) = 0.
		uGrid%u_0(i,k) = - 12.*MAX(0.,1. - wGrid%zz(i,k)/2500.)
		uGrid%u(i,k) = uGrid%u_0(i,k)
	END DO
END DO

CALL set_area_w

DO k = kmin, kmax
	DO i = imin, imax
		wGrid%theta_M_0(i,k) = wGrid%theta_0(i,k)
		wGrid%theta_M(i,k) = wGrid%theta_M_0(i,k)
		wGrid%theta(i,k) = wGrid%theta_M(i,k)/(1. + 0.61*wGrid%qv(i,k))
	END DO
END DO

! Add bubble to theta
DO k = kmin, kmax
	DO i = imin, imax
		wGrid%w(i,k) = 0
		L = SQRT((wGrid%xx(i) - x_c)*(wGrid%xx(i) - x_c) + (wGrid%zz(i,k) - z_c)*(wGrid%zz(i,k) - z_c))
		IF (L/R <= 1.) THEN
			wGrid%theta_1(i,k) = DeltaTheta*COS(.5*PI_math*L/R)**2  ! <= I want this.
		ELSE
			wGrid%theta_1(i,k) = 0.
		END IF
		wGrid%theta(i,k) = wGrid%theta(i,k) + wGrid%theta_1(i,k)/(1. + 0.61*wGrid%qv(i,k))
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
REAL(kd) :: a, b
INTEGER, PARAMETER :: n_max = 34
REAL(kd), DIMENSION(1:n_max) :: aa = (/ &
     0.00000,     0.00000,     0.00000,    72.06358,   271.62646, &
   638.14893,  1195.88989,  1954.80518,  2911.56934,  4050.71777, &
  5345.91406,  6761.33984,  8253.21094,  9771.40625, 11261.23047, &
 12665.29297, 13925.51953, 14985.26953, 15791.59766, 16297.62109, &
 16465.00391, 16266.60938, 15689.20703, 14736.35547, 13431.39453, &
 11820.53906,  9976.13672,  8000.00000,  6000.00000,  4500.00000, &
  3000.00000,  2000.00000,  1500.00000,     0.00000 /)

REAL(kd), DIMENSION(1:n_max) :: bb = (/ &
 1.00000,     0.99228,     0.97299,     0.94421,     0.90772, &
 0.86496,     0.81717,     0.76541,     0.71059,     0.65356, &
 0.59508,     0.53589,     0.47669,     0.41820,     0.36115, &
 0.30624,     0.25419,     0.20570,     0.16144,     0.12200, &
 0.08789,     0.05949,     0.03697,     0.02032,     0.00919, &
 0.00292,     0.00039,     0.00000,     0.00000,     0.00000, &
 0.00000,     0.00000,     0.00000,     0.00000 /)
REAL(kd), DIMENSION(1:n_max) :: p, t0
REAL(kd), PARAMETER :: ps = 0.101325e6, ptop = 1000., t01=300., rr = 287., g = 9.81
REAL(kd) :: kk
!=================================================
! xx, zeta
!-------------------------------------------------
DO i = ims, ime
	uGrid%xx(i) = dx*(i - its)
	virGrid%xx(i) = uGrid%xx(i)
	piGrid%xx(i) = uGrid%xx(i) - dx/2.
	wGrid%xx(i) = piGrid%xx(i)
END DO

IF (VertLev == 0) THEN
	DO k = kms, kme
		dk(k) = dz
		wGrid%zeta(k) = dk(k)*(k - kts)
	END DO

ELSE IF (VertLev == 1) THEN
	CALL set_area_w
	wGrid%zeta(kmin) = 0.
	DO k = kmin+1, kmax-1
		a = REAL(k-kmin)/REAL(kmax-kmin)
		b = 0.75*a + 1.75*a*a*a - 1.5*a*a*a*a
		wGrid%zeta(kmax+1-k) = (1. - b)*ztop
	END DO
	wGrid%zeta(kmax) = ztop
	DO k = kmin, kmax-1
		dk(k) = wGrid%zeta(k+1) - wGrid%zeta(k)
	END DO
	DO k = kmin-1, kms, -1
		dk(k) = dk(kmin)
		wGrid%zeta(k) = wGrid%zeta(k+1) - dk(k)
	END DO
	DO k = kmax, kme
		dk(k) = dk(kmax-1)
		wGrid%zeta(k) = wGrid%zeta(k-1) + dk(k-1)
	END DO

ELSE IF (VertLev == 2) THEN
	CALL set_area_w
	IF (nz > n_max) THEN
		STOP "Too many vertical layers!!!"
	ELSE

		DO k = 1, n_max
			p(k) = aa(k) + bb(k)*ps
		END DO
		p(n_max) = ptop

		t0(kmin) = t01
		kk = - 6.5/1000.

		wGrid%zeta(kmin) = 0.
		DO k = kmin+1, kmax
			wGrid%zeta(k) = wGrid%zeta(k-1) + t0(k-1)/kk*(EXP(kk*rr/g*LOG(p(k-1)/p(k))) - 1.)
			IF (wGrid%zeta(k) < 13000.) THEN
				kk = - 6.0/1000.
			ELSE
				kk = 2./1000.
			END IF
			t0(k) = t0(k-1) + kk*(wGrid%zeta(k) - wGrid%zeta(k-1))
		END DO
		ztop = wGrid%zeta(kmax)
		WRITE(*,*) "Modify ztop:            ", ztop
		DO k = kmin-1, kms, -1
			wGrid%zeta(k) = - wGrid%zeta(kmin+kmin-k)
		END DO
		DO k = kmax+1, kme
			wGrid%zeta(k) = ztop + ztop - wGrid%zeta(kmax - (k - kmax))
		END DO
		DO k = kms, kme-1
			dk(k) = wGrid%zeta(k+1) - wGrid%zeta(k)
		END DO
		dk(kme) = dk(kme-1)
	END IF

ELSE IF (VertLev == 99) THEN
	CALL set_area_w
	DO k = kms, kme
		dk(k) = dz
	END DO
	wGrid%zeta(kmin) = 0.
	DO k = kmin-1, kms, -1
		wGrid%zeta(k) = wGrid%zeta(k+1) - dk(k)
	END DO
	DO k = kmin+1, kme
		wGrid%zeta(k) = wGrid%zeta(k-1) + dk(k-1)
	END DO
END IF

DO k = kms, kme
	virGrid%zeta(k) = wGrid%zeta(k)
	piGrid%zeta(k) = wGrid%zeta(k) + dk(k)/2.
	uGrid%zeta(k) = piGrid%zeta(k)
END DO

WRITE(*,*) wGrid%zeta
OPEN(1, FILE="./output/modelvar_zeta.bin", FORM='binary', CONVERT='big_endian')
	WRITE(1) wGrid%zeta
CLOSE(1)
!WRITE(*,*) piGrid%zeta
!CALL debug_SFSG
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

REAL(kd), DIMENSION(nz) :: theta, pressure ! sounding input in thunderstorm case
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! theta_0, pi_0, rho_0
!-------------------------------------------------

! theta_0, pi_0
!IF (RunCase == 1 .OR. RunCase == 2 .OR. RunCase == 6) THEN
!IF (RunCase == 1 .OR. RunCase == 2 .OR. RunCase == 99) THEN
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

	!IF (RunCase == 6) THEN
		!OPEN(1, FILE="./input/theta.input", STATUS='old')
		!DO k = 1, nz
			!READ(1,*) theta(k)
		!END DO
		!CLOSE(1)

		!OPEN(1, FILE="./input/pressure.input", STATUS='old')
		!DO k = 1, nz
			!READ(1,*) pressure(k)
		!END DO
		!CLOSE(1)

		!CALL set_area_pi
		!DO k = kmin, kmax
			!DO i = ims, ime
				!piGrid%theta_M(i,k) = theta(k)
				!piGrid%pi(i,k) = (pressure(k)/p0)**(Rd/Cp)
			!END DO
		!END DO

		!CALL set_area_w
		!wGrid%theta_M(:,kmin) = 3.000E+02
		!wGrid%theta_M(:,kmax) = 4.978E+02
		!wGrid%pi(:,kmin) = 1
		!wGrid%pi(:,kmax) = 0.4392

		!DO k = kmin+1, kmax-1
			!wGrid%theta_M(:,k) = (piGrid%theta_M(:,k-1) + piGrid%theta_M(:,k))/2.
			!wGrid%pi(:,k) = (piGrid%pi(:,k-1) + piGrid%pi(:,k))/2.
		!END DO

		!virGrid%theta_M = wGrid%theta_M
		!uGrid%theta_M = piGrid%theta_M
		!virGrid%pi = wGrid%pi
		!uGrid%pi = piGrid%pi
	!END IF

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
! Initiate Real Case.
!=================================================
SUBROUTINE initiate_real(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
!-------------------------------------------------
TYPE (grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!-------------------------------------------------
INTEGER :: i, k
REAL(kd), PARAMETER :: DeltaTime = 180.
INTEGER, PARAMETER :: n = 3

REAL(kd), DIMENSION(kts:kte) :: u_ls, v_ls, pi_ls
REAL(kd), DIMENSION(kts:kte+1) :: w_ls, th_ls, qv_ls, qc_ls, qr_ls, qi_ls, qs_ls, qg_ls
REAL(kd) :: deltaTheta
REAL(kd), DIMENSION(ims:ime) :: ran
!=================================================
OPEN(1, FILE="./input/real_case/v_ls.dat", STATUS='old')
DO k = kts, kte
	READ(1,*) v_ls(k)
	DO i = ims, ime
		uGrid%v(i,k) = v_ls(k)
	END DO
END DO
CLOSE(1)

OPEN(1, FILE="./input/real_case/u_ls.dat", STATUS='old')
DO k = kts, kte
	READ(1,*) u_ls(k)
	DO i = ims, ime
        !uGrid%u_0(i,k) = u_ls(k)
        uGrid%u_0(i,k) = 0.
		uGrid%u(i,k) = uGrid%u_0(i,k)
	END DO
END DO
CLOSE(1)
OPEN(1, FILE="./input/real_case/w_ls.dat", STATUS='old')
DO k = kts, kte + 1
	READ(1,*) w_ls(k)
	DO i = ims, ime
        !wGrid%w(i,k) = w_ls(k)
        wGrid%w(i,k) = 0.
	END DO
END DO
CLOSE(1)
OPEN(1, FILE="./input/real_case/pi_ls.dat", STATUS='old')
DO k = kts, kte
	READ(1,*) pi_ls(k)
	DO i = ims, ime
		piGrid%pi(i,k) = pi_ls(k)
	END DO
END DO
CLOSE(1)
OPEN(1, FILE="./input/real_case/th_ls.dat", STATUS='old')
DO k = kts, kte + 1
	READ(1,*) th_ls(k)
	DO i = ims, ime
		wGrid%theta(i,k) = th_ls(k)
	END DO
END DO
CLOSE(1)
OPEN(1, FILE="./input/real_case/qv_ls.dat", STATUS='old')
DO k = kts, kte + 1
	READ(1,*) qv_ls(k)
	DO i = ims, ime
		wGrid%qv(i,k) = qv_ls(k)
	END DO
END DO
CLOSE(1)
OPEN(1, FILE="./input/real_case/qc_ls.dat", STATUS='old')
DO k = kts, kte + 1
	READ(1,*) qc_ls(k)
	DO i = ims, ime
		wGrid%qc(i,k) = qc_ls(k)
	END DO
END DO
CLOSE(1)
OPEN(1, FILE="./input/real_case/qr_ls.dat", STATUS='old')
DO k = kts, kte + 1
	READ(1,*) qr_ls(k)
	DO i = ims, ime
		wGrid%qr(i,k) = qr_ls(k)
	END DO
END DO
CLOSE(1)
OPEN(1, FILE="./input/real_case/qi_ls.dat", STATUS='old')
DO k = kts, kte + 1
	READ(1,*) qi_ls(k)
	DO i = ims, ime
		wGrid%qi(i,k) = qi_ls(k)
	END DO
END DO
CLOSE(1)
OPEN(1, FILE="./input/real_case/qs_ls.dat", STATUS='old')
DO k = kts, kte + 1
	READ(1,*) qs_ls(k)
	DO i = ims, ime
		wGrid%qs(i,k) = qs_ls(k)
	END DO
END DO
CLOSE(1)
OPEN(1, FILE="./input/real_case/qg_ls.dat", STATUS='old')
DO k = kts, kte + 1
	READ(1,*) qg_ls(k)
	DO i = ims, ime
		wGrid%qg(i,k) = qg_ls(k)
	END DO
END DO
CLOSE(1)
!-------------------------------------------------
CALL set_area_w
CALL RANDOM_SEED()
CALL RANDOM_NUMBER(ran(imin:imax))
deltaTheta = 2.
ran(imin:imax) = ran(imin:imax) - .5
ran(imin:imax) = ran(imin:imax)*deltaTheta
DO i = imin, imax
    wGrid%theta(i,1) = wGrid%theta(i,1) + ran(i)
END DO
!-------------------------------------------------

CALL set_area_pi
DO k = kmin, kmax
	DO i = imin, imax
		piGrid%pi_1(i,k) = piGrid%pi(i,k) - piGrid%pi_0(i,k)
	END DO
END DO

CALL set_area_w
DO k = kmin, kmax
	DO i = imin, imax
		wGrid%theta_1(i,k) = piGrid%theta(i,k) - piGrid%theta_0(i,k)
	END DO
END DO
!-------------------------------------------------
wGrid%rain = 0.
wGrid%rainncv = 0.
wGrid%sr = 0.
wGrid%snow = 0.
wGrid%snowncv = 0.
wGrid%graupel = 0.
wGrid%graupelncv = 0.
wGrid%cldfra = 0.
wGrid%Mtheta = 0.
wGrid%Mqv = 0.
wGrid%Mqc = 0.
wGrid%Mqr = 0.
wGrid%Mqi = 0.
wGrid%Mqs = 0.
wGrid%Mqg = 0.
!-------------------------------------------------

!CALL set_area_u
!DO k = kmin, kmax
	!uGrid%forcing_u(k) = u_ls(k)/DeltaTime
!END DO
uGrid%forcing_u = 0.
wGrid%forcing_theta = 0.
wGrid%forcing_qv = 0.
wGrid%forcing_qc = 0.
wGrid%forcing_qr = 0.
wGrid%forcing_qi = 0.
wGrid%forcing_qs = 0.
wGrid%forcing_qg = 0.

!=================================================
!CALL debug_ascii_output(wGrid%theta,"theta")
!CALL debug_ascii_output(wGrid%theta_0,"theta_0")
!CALL debug_ascii_output(wGrid%theta_1,"theta_1")
!CALL debug_ascii_output(piGrid%pi_1,"pi_1")
!CALL debug_ascii_output(wGrid%pi_1,"rho_0_w")
!CALL debug_ascii_output(uGrid%pi_1,"rho_0_u")
!CALL debug_ascii_output(piGrid%pi_1,"rho_0_pi")
!CALL debug_ascii_output(virGrid%pi_1,"rho_0_vir")
!CALL debug_SFSG
END SUBROUTINE initiate_real
!=================================================

!=================================================
END MODULE sp_module_initiate
!=================================================

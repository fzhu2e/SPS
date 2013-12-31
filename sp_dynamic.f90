!=================================================
! The dynamic core of Super-Parametertization System (SPS).
!-------------------------------------------------
! Version: 0.11
! Author: Zhu F.
! Email: lyricorpse@gmail.com
! Date: 2013-05-04 12:47:26 
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!=================================================
PROGRAM sp_dynamic
USE sp_module_constant
USE sp_module_model
USE sp_module_gridvar
USE sp_module_initiate
USE sp_module_boundary
USE sp_module_integrate
USE sp_module_output
USE sp_module_debug
IMPLICIT NONE
!=================================================
TYPE(grid) :: uGrid, wGrid, piGrid, virGrid
!-------------------------------------------------
INTEGER :: i, k
INTEGER :: t_start, t_end, rate
REAL(kd) :: t_lapse, t_left, t_all
!=================================================
! Initial an ideal case.
!-------------------------------------------------
WRITE(*,*) "====================="
WRITE(*,*) " Initial model..."
WRITE(*,*) "====================="
WRITE(*,*) " Kind:            ", kd
WRITE(*,*) " RunCase:         ", RunCase
WRITE(*,*) " TimeScheme:      ", TimeScheme
WRITE(*,*) " AdvectionScheme: ", AdvectionScheme
WRITE(*,*) " LateralBoundary: ", LateralBoundary
WRITE(*,*) " UpperBoundary:   ", UpperBoundary
WRITE(*,*) "---------------------"
WRITE(*,*) " nstep: ", nstep
WRITE(*,*) " nx/nz: ", nx, nz
WRITE(*,"(1X,A9,2F9.2)") " dx/dz: ", dx, dz
WRITE(*,"(1X,A9,F9.2)") "    dt: ", dt
WRITE(*,*) "---------------------"
WRITE(*,"(1X,A9,2F9.2)") " Km/Kh: ", Km, Kh
WRITE(*,*) "====================="
WRITE(*,*)

CALL initiate_grid(uGrid,wGrid,piGrid,virGrid)
CALL initiate_terrain(uGrid,wGrid,piGrid,virGrid)
CALL initiate_basic_state(uGrid,wGrid,piGrid,virGrid)
!-------------------------------------------------
! Initiate.
!-------------------------------------------------
WRITE(*,*) "====================="
WRITE(*,*) " Initial case..."
WRITE(*,*) "====================="
WRITE(*,*)
SELECT CASE (RunCase)
CASE (1)
	CALL initiate_dc(uGrid,wGrid,piGrid,virGrid)
CASE (2)
	CALL initiate_tb(uGrid,wGrid,piGrid,virGrid)
CASE (3)
	CALL initiate_igw(uGrid,wGrid,piGrid,virGrid)
CASE (4)
	CALL initiate_Sm(uGrid,wGrid,piGrid,virGrid)
CASE (5)
	CALL initiate_wb(uGrid,wGrid,piGrid,virGrid)
CASE (6)
	CALL initiate_th(uGrid,wGrid,piGrid,virGrid)
CASE DEFAULT
	STOP "Wrong ideal case!!!"
END SELECT
IF (Vapor == 0) THEN
	wGrid%qv = 0.
END IF

CALL set_area_w
DO k = kmin, kmax
	DO i = imin, imax
		wGrid%qc = 0.
		wGrid%qr = 0.
		wGrid%qi = 0.
		wGrid%qs = 0.
		wGrid%qg = 0.
		wGrid%rain = 0.
		wGrid%rainncv = 0.
		wGrid%sr = 0.
		wGrid%snow = 0.
		wGrid%snowncv = 0.
		wGrid%graupel = 0.
		wGrid%graupelncv = 0.
		wGrid%Mtheta = 0.
		wGrid%Mqv = 0.
		wGrid%Mqc = 0.
		wGrid%Mqr = 0.
		wGrid%Mqi = 0.
		wGrid%Mqs = 0.
		wGrid%Mqg = 0.
	END DO
END DO

CALL update_boundary(uGrid%u,wGrid%w,wGrid,piGrid%pi_1,wGrid%theta,                   &
                     wGrid%qv,wGrid%qc,wGrid%qr,wGrid%qi,wGrid%qs,wGrid%qg,           &
                     piGrid%rho_0,uGrid%rho_0,wGrid%rho_0,virGrid%rho_0,              &
                     wGrid%theta_0, piGrid%pi_0)

wGrid%theta_M_0 = wGrid%theta_0
piGrid%theta_M_0 = piGrid%theta_0
uGrid%theta_M_0 = uGrid%theta_0

CALL calc_virTheta(uGrid,wGrid,piGrid,virGrid)

!=================================================
! Adjust pi_1
!-------------------------------------------------
!CALL set_area_w
!wGrid%pi(imin:imax,kmax) = wGrid%pi_0(imin:imax,kmax)
!DO k = kmax-1, kmin, -1
	!DO i = imin, imax
		!wGrid%pi(i,k) = wGrid%pi(i,k+1) + g*dz/Cp/wGrid%theta(i,k)
	!END DO
!END DO

!CALL set_area_pi
!DO k = kmin, kmax
	!DO i = imin, imax
		!piGrid%pi(i,k) = (wGrid%pi(i,k) + wGrid%pi(i,k+1))/2.
		!piGrid%pi_1(i,k) = piGrid%pi(i,k) - piGrid%pi_0(i,k)
	!END DO
!END DO
!CALL update_boundary(pi_1=piGrid%pi_1,wGrid=wGrid)
!=================================================

IF (Vapor == 0) THEN
	CALL output(0,uGrid%u,wGrid%w,piGrid%pi_1,wGrid%theta_M_1,wGrid%theta_M, wGrid%theta)
	
ELSE
	CALL output(0,uGrid%u,wGrid%w,piGrid%pi_1,wGrid%theta_M_1,wGrid%theta_M, wGrid%theta, &
	              wGrid%qv,wGrid%qc,wGrid%qr,wGrid%qi,wGrid%qs,wGrid%qg,                  &
	              wGrid%rain,wGrid%snow,wGrid%graupel     )
END IF
!=================================================
! Integrate.
!-------------------------------------------------
IF (OpenTop == 2) THEN
	CALL calc_tau_top(uGrid,wGrid)
END IF
IF (OpenLateral == 2) THEN
	CALL calc_tau_lateral(uGrid,wGrid)
END IF
!CALL debug_ascii_output(uGrid%tau,"tau_u")
!CALL debug_ascii_output(wGrid%tau,"tau_w")
!CALL debug_SFSG

CALL wsm6init(rhoair0,rhowater,rhosnow,cliq,cpv)
t_all = 0.
DO i = 1, nstep
	CALL SYSTEM_CLOCK(t_start,rate)
	CALL integrate(i,uGrid,wGrid,piGrid,virGrid) ! main integrate module
	CALL update_boundary(uGrid%u,wGrid%w,wGrid,piGrid%pi_1,wGrid%theta,       &
	                     wGrid%qv,wGrid%qc,wGrid%qr,wGrid%qi,wGrid%qs,wGrid%qg)
	CALL calc_virTheta(uGrid,wGrid,piGrid,virGrid)
	IF (MOD(i,100) == 0.) THEN
		IF (Vapor == 0) THEN
			CALL output(1,uGrid%u,wGrid%w,piGrid%pi_1,wGrid%theta_M_1,wGrid%theta_M, wGrid%theta)
		ELSE
			CALL output(1,uGrid%u,wGrid%w,piGrid%pi_1,wGrid%theta_M_1,wGrid%theta_M, wGrid%theta, &
			              wGrid%qv,wGrid%qc,wGrid%qr,wGrid%qi,wGrid%qs,wGrid%qg,                  &
	                      wGrid%rain,wGrid%snow,wGrid%graupel     )
		END IF
	END IF
	
	CALL SYSTEM_CLOCK(t_end)
	t_lapse = REAL(t_end - t_start)/REAL(rate)
	!t_left = t_lapse*(nstep - i)/60./60.  ! unit: hour
	t_left = t_lapse*(nstep - i)/60.  ! unit: minute
	t_all = t_all + t_lapse
	!WRITE(*,"('Step/nStep -- time lapse/left: ',2X,I6,'/ ',I6,' --',F12.6,' sec/',1X,F6.3,' hr')") , i, nstep, t_lapse, t_left
	WRITE(*,"('Step/nStep -- time lapse/left: ',2X,I6,'/ ',I6,' --',F12.6,' sec/',1X,F6.3,' min')") , i, nstep, t_lapse, t_left
END DO
!=================================================
! Finish.
!-------------------------------------------------
IF (Vapor == 0) THEN
	CALL output(99,uGrid%u,wGrid%w,piGrid%pi_1,wGrid%theta_M_1,wGrid%theta_M, wGrid%theta)
ELSE
	CALL output(99,uGrid%u,wGrid%w,piGrid%pi_1,wGrid%theta_M_1,wGrid%theta_M, wGrid%theta, &
	              wGrid%qv,wGrid%qc,wGrid%qr,wGrid%qi,wGrid%qs,wGrid%qg,                  &
                  wGrid%rain,wGrid%snow,wGrid%graupel     )
END IF
WRITE(*,*)
WRITE(*,*) "====================="
WRITE(*,*) " Finish!!!"
WRITE(*,*) "====================="
WRITE(*,*) " nstep: ", nstep
WRITE(*,*) " nx/nz: ", nx, nz
WRITE(*,*) " dx/dz: ", dx, dz
WRITE(*,*) "    dt: ", dt
WRITE(*,*) "---------------------"
WRITE(*,*) " RunCase:         ", RunCase
WRITE(*,*) " TimeScheme:      ", TimeScheme
WRITE(*,*) " AdvectionScheme: ", AdvectionScheme
WRITE(*,*) " LateralBoundary: ", LateralBoundary
WRITE(*,*) " UpperBoundary:   ", UpperBoundary
WRITE(*,*) "---------------------"
WRITE(*,*) " Km/Kh: ", Km, Kh
WRITE(*,*) "---------------------"
!WRITE(*,*) " TIME: ", t_all/60./60., "hr"
WRITE(*,*) " TIME: ", t_all/60., "min"
WRITE(*,*) "====================="
WRITE(*,*)
!-------------------------------------------------
!=================================================
END PROGRAM sp_dynamic
!=================================================

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
REAL(kd), DIMENSION(ims:ime,kms:kme) :: u        ! wind speed along x-axis
REAL(kd), DIMENSION(ims:ime,kms:kme) :: w        ! wind speed along z-axis
REAL(kd), DIMENSION(ims:ime,kms:kme) :: pi_1     ! pi'
REAL(kd), DIMENSION(ims:ime,kms:kme) :: pi_0 
REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta
REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_0  ! theta = theta_0 + theta'
REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_1  ! theta'
REAL(kd), DIMENSION(ims:ime,kms:kme) :: rho_0    ! density
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

CALL debug_undef_all(u,w,pi_1,pi_0,theta,theta_0,theta_1,rho_0)
!CALL debug_ascii_output(u)
!-------------------------------------------------
! Initiate.
!-------------------------------------------------
WRITE(*,*) "====================="
WRITE(*,*) " Initial case..."
WRITE(*,*) "====================="
WRITE(*,*)
SELECT CASE (RunCase)
CASE (1)
	CALL initiate_dc(u,w,pi_1,pi_0,theta,theta_0,theta_1,rho_0)  ! initiate the DC case
!CASE (2)
	!CALL initiate_tb(u,w,pi_1,pi_0,theta,theta_0,theta_1,rho_0)  ! initiate the TB case
CASE (3)
	CALL initiate_igw(u,w,pi_1,pi_0,theta,theta_0,theta_1,rho_0)  ! initiate the IGW case
!CASE (4)
	!CALL initiate_Sm(u,w,pi_1,pi_0,theta,theta_0,theta_1,rho_0)  ! initiate the IGW case
CASE DEFAULT
	STOP "Wrong ideal case!!!"
END SELECT

CALL update_boundary(u,w,pi_1,theta,theta_1,                  &
                     theta_0,theta_0_pi,theta_0_u,theta_0_vir,  &
                     rho_0,rho_0_w,rho_0_u,rho_0_vir            )

!CALL debug_ascii_output(u)

CALL output(0,u,w,pi_1,theta_1)                               ! output the initial fields
!-------------------------------------------------
!CALL debug_ascii_output(pi)
!CALL debug_ascii_output(pi_0)
!CALL debug_ascii_output(pi_1)
!CALL debug_ascii_output(theta)
!CALL debug_ascii_output(theta_1)
!CALL debug_ascii_output(rho_0)
!=================================================
! Integrate.
!-------------------------------------------------
t_all = 0.
DO i = 1, nstep
	CALL SYSTEM_CLOCK(t_start,rate)
	CALL integrate(u,w,pi_1,pi_0,theta,theta_0,theta_1,rho_0) ! main integrate module
	CALL update_boundary(u,w,pi_1,theta,theta_1)
	!IF (MOD(i,1000) == 0.) THEN
	!IF (MOD(i,200) == 0.) THEN
	IF (MOD(i,100) == 0.) THEN
		CALL output(1,u,w,pi_1,theta_1)               ! output the fields at each time step
	END IF
	
	CALL SYSTEM_CLOCK(t_end)
	t_lapse = REAL(t_end - t_start)/REAL(rate)
	t_left = t_lapse*(nstep - i)/60./60.  ! unit: hour
	t_all = t_all + t_lapse
	WRITE(*,"('Step/nStep -- time lapse/left: ',2X,I6,'/ ',I6,' --',F12.6,' sec/',1X,F6.3,' hr')") , i, nstep, t_lapse, t_left
END DO
!=================================================
! Finish.
!-------------------------------------------------
CALL output(99,u,w,pi_1,theta_1)                   ! finish
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
WRITE(*,*) " TIME: ", t_all/60./60., "hr"
WRITE(*,*) "====================="
WRITE(*,*)
!-------------------------------------------------
!=================================================
END PROGRAM sp_dynamic
!=================================================

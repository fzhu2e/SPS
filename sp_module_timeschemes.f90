!=================================================
! The timeschemes module of SPS-dynamic-integrate
!-------------------------------------------------
! Version: 0.11
! Author: Zhu F.
! Email: lyricorpse@gmail.com
! Date: 2013-05-04 13:59:46 
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!=================================================
MODULE sp_module_timeschemes
USE sp_module_constant
USE sp_module_model
USE sp_module_tendency
USE sp_module_interpolate
USE sp_module_boundary
USE sp_module_math
USE sp_module_debug
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================

!=================================================
! Forward-backward Scheme.
!-------------------------------------------------
! 1. Update u, w, theta with value at n-time point.
! 2. Update pi' with the updated u, w.
!=================================================
SUBROUTINE forward_backward( DeltaT,v,theta_0,rho_0,                     &
                             old_u,old_w,old_pi_1,old_theta,old_theta_1, &
                             new_u,new_w,new_pi_1,new_theta,new_theta_1  )
IMPLICIT NONE
!=================================================
REAL(preci), INTENT(IN) :: DeltaT
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: v            ! wind speed along y-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_0      ! theta = theta_0 + theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0        ! density
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_u        ! wind speed along x-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_w        ! wind speed along z-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_pi_1     ! pi'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_theta_1  ! theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_u
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_w
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_pi_1
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_theta_1
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL basic_interpolate(old_u,old_w,old_pi_1,old_theta,theta_0,old_theta_1,rho_0)

CALL debug_undef_all( new_u,new_w,new_pi_1,new_theta,new_theta_1, &
                      F_u,rhou_pi,rhouu_pi,                       &
                      PrhouPx_u,PrhouuPx_u,                       &
                      rhow_v,rhouw_v,                             &
                      PrhowPz_u,PrhouwPz_u,                       &
                      F_w,rhou_v,PrhouPx_w,PrhouwPx_w,            &
                      rhow_pi,rhoww_pi,PrhowPz_w,PrhowwPz_w,      &
                      F_theta,rhoutheta_v,PrhowthetaPz_w,         &
                      Ppi_1Px_u,Ppi_1Pz_w,                        &
                      F_pi,urhotheta_u,wrhotheta_w,               &
                      PurhothetaPx_pi,PwrhothetaPz_pi,            &
                      tend_u,tend_w,tend_theta,tend_pi,           &
                      P2uPx2_u,P2uPz2_u,                          &
                      P2wPx2_w,P2wPz2_w,                          &
                      P2thetaPx2_w,P2thetaPz2_w                   )

CALL tendency_u(old_u,rho_0,old_pi_1,F_u,tend_u)
CALL tendency_w(old_w,rho_0,old_theta_1,theta_0,old_pi_1,F_w,tend_w)
CALL tendency_theta(old_u,old_w,rho_0,old_theta,F_theta,tend_theta)
!=================================================
! 4. Update u, w, theta
!-------------------------------------------------
! u-grid
DO i = its + 1, ite -1
	DO k = kts, kte
		new_u(i,k) = old_u(i,k) + DeltaT*tend_u(i,k)
	END DO
END DO
! w-grid
DO i = its + 1, ite
	DO k = kts + 1, kte
		new_w(i,k) = old_w(i,k) + DeltaT*tend_w(i,k)
		new_theta(i,k) = old_theta(i,k) + DeltaT*tend_theta(i,k)
		new_theta_1(i,k) = new_theta(i,k) - theta_0(i,k)
	END DO
END DO
CALL update_boundary(new_u,new_w)

!DO k = kme, kms, -1
	!WRITE(*,"(2F15.9)") tend_theta(255,k), tend_w(255,k)
!END DO
!CALL debug_ascii_output(Ppi_1Pz_w)
!CALL debug_test_boundary(Ppi_1Pz_w)
!CALL debug_ascii_output(Cp*theta_0*Ppi_1Pz_w)
!CALL debug_ascii_output(theta_0*Ppi_1Pz_w)
!CALL debug_test_boundary(Cp*theta_0*Ppi_1Pz_w)
!CALL debug_ascii_output(F_w)
!CALL debug_ascii_output(F_w - Cp*theta_0*Ppi_1Pz_w)
!CALL debug_ascii_output(F_w - theta_0*Ppi_1Pz_w)

!=================================================
! 5.2 Update pi_1
!-------------------------------------------------
CALL tendency_pi(new_u,new_w,rho_0,theta_0,F_pi,tend_pi)
! pi-grid
DO i = its + 1, ite
	DO k = kts, kte
		new_pi_1(i,k) = old_pi_1(i,k) + DeltaT*tend_pi(i,k)
	END DO
END DO
!=================================================
END SUBROUTINE forward_backward
!=================================================


!=================================================
! Leapfrog scheme.
!-------------------------------------------------
! var(n+1) = var(n-1) + DeltaT*tend_var(n)
!=================================================
SUBROUTINE leapfrog( DeltaT,v,theta_0,rho_0,                        &
                     old_u,old_w,old_pi_1,old_theta,old_theta_1,    &
                     new_u,new_w,new_pi_1,new_theta,new_theta_1     )
IMPLICIT NONE
!=================================================
REAL(preci), INTENT(IN) :: DeltaT
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: v            ! wind speed along y-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_0      ! theta = theta_0 + theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0        ! density
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_u        ! wind speed along x-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_w        ! wind speed along z-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_pi_1     ! pi'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_theta_1  ! theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_u
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_w
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_pi_1
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_theta_1
!-------------------------------------------------
INTEGER :: i, k
!=================================================

!=================================================
END SUBROUTINE leapfrog
!=================================================

!=================================================
! Runge-Kutta 3 order scheme.
!-------------------------------------------------
! phi* = phi(n) + dt/3.*tend(phi(n))
! phi** = phi(n) + dt/2.*tend(phi*)
! phi(n+1) = phi(n) + dt*tend(phi**)
!=================================================
SUBROUTINE runge_kutta( DeltaT,v,theta_0,rho_0,                        &
                        old_u,old_w,old_pi_1,old_theta,old_theta_1,    &
                        new_u,new_w,new_pi_1,new_theta,new_theta_1     )
IMPLICIT NONE
!=================================================
REAL(preci), INTENT(IN) :: DeltaT
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: v
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_0
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_u
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_w
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_pi_1
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_theta_1
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_u
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_w
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_pi_1
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_theta_1
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid1_u
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid1_w
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid1_pi_1
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid1_theta
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid1_theta_1
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid2_u
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid2_w
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid2_pi_1
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid2_theta
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid2_theta_1
!-------------------------------------------------
INTEGER :: i, k, imin, imax, kmin, kmax
!=================================================
! Step 1. phi* = phi(n) + dt/3.*tend(phi(n))
!-------------------------------------------------
CALL basic_interpolate(old_u,old_w,old_pi_1,old_theta,theta_0,old_theta_1,rho_0)

CALL debug_undef_all( mid1_u,mid1_w,mid1_pi_1,mid1_theta,mid1_theta_1, &
                      F_u,rhou_pi,rhouu_pi,                       &
                      PrhouPx_u,PrhouuPx_u,                       &
                      rhow_v,rhouw_v,                             &
                      PrhowPz_u,PrhouwPz_u,                       &
                      F_w,rhou_v,PrhouPx_w,PrhouwPx_w,            &
                      rhow_pi,rhoww_pi,PrhowPz_w,PrhowwPz_w,      &
                      F_theta,rhoutheta_v,PrhowthetaPz_w,         &
                      Ppi_1Px_u,Ppi_1Pz_w,                        &
                      F_pi,urhotheta_u,wrhotheta_w,               &
                      PurhothetaPx_pi,PwrhothetaPz_pi,            &
                      tend_u,tend_w,tend_theta,tend_pi,           &
                      P2uPx2_u,P2uPz2_u,                          &
                      P2wPx2_w,P2wPz2_w,                          &
                      P2thetaPx2_w,P2thetaPz2_w                   )

CALL tendency_u(old_u,rho_0,old_pi_1,F_u,tend_u)
CALL tendency_w(old_w,rho_0,old_theta_1,theta_0,old_pi_1,F_w,tend_w)
CALL tendency_theta(old_u,old_w,rho_0,old_theta,F_theta,tend_theta)

! u-grid (its + 1:ite - 1, kts:kte)
imin = its + 1
imax = ite - 1
kmin = kts
kmax = kte

FORALL (i = imin:imax, k = kmin:kmax)
	mid1_u(i,k) = old_u(i,k) + DeltaT/3.*tend_u(i,k)
END FORALL

! w-grid (it + 1:ite, kts + 1:kte)
imin = its + 1
imax = ite
kmin = kts + 1
kmax = kte
! ATTENTION: The calculated area includes the boundary layers.
kmin = kmin - 1
kmax = kmax + 1

FORALL (i = imin:imax, k = kmin:kmax)
	mid1_w(i,k) = old_w(i,k) + DeltaT/3.*tend_w(i,k)
	mid1_theta(i,k) = old_theta(i,k) + DeltaT/3.*tend_theta(i,k)
	mid1_theta_1(i,k) = mid1_theta(i,k) - theta_0(i,k)
END FORALL

CALL update_boundary(mid1_u,mid1_w)

CALL tendency_pi(mid1_u,mid1_w,rho_0,theta_0,F_pi,tend_pi)

! pi-grid (its + 1:ite, kts:kte)
imin = its + 1
imax = ite
kmin = kts
kmax = kte

FORALL (i = imin:imax, k = kmin:kmax)
	mid1_pi_1(i,k) = old_pi_1(i,k) + DeltaT/3.*tend_pi(i,k)
END FORALL

CALL update_boundary(mid1_u,mid1_w,mid1_pi_1,mid1_theta,mid1_theta_1)

!=================================================
! Step 2. phi** = phi(n) + dt/2.*tend(phi*)
!-------------------------------------------------
CALL basic_interpolate(mid1_u,mid1_w,mid1_pi_1,mid1_theta,theta_0,mid1_theta_1,rho_0)

CALL debug_undef_all( mid2_u,mid2_w,mid2_pi_1,mid2_theta,mid2_theta_1, &
                      F_u,rhou_pi,rhouu_pi,                       &
                      PrhouPx_u,PrhouuPx_u,                       &
                      rhow_v,rhouw_v,                             &
                      PrhowPz_u,PrhouwPz_u,                       &
                      F_w,rhou_v,PrhouPx_w,PrhouwPx_w,            &
                      rhow_pi,rhoww_pi,PrhowPz_w,PrhowwPz_w,      &
                      F_theta,rhoutheta_v,PrhowthetaPz_w,         &
                      Ppi_1Px_u,Ppi_1Pz_w,                        &
                      F_pi,urhotheta_u,wrhotheta_w,               &
                      PurhothetaPx_pi,PwrhothetaPz_pi,            &
                      tend_u,tend_w,tend_theta,tend_pi,           &
                      P2uPx2_u,P2uPz2_u,                          &
                      P2wPx2_w,P2wPz2_w,                          &
                      P2thetaPx2_w,P2thetaPz2_w                   )

CALL tendency_u(mid1_u,rho_0,mid1_pi_1,F_u,tend_u)
CALL tendency_w(mid1_w,rho_0,mid1_theta_1,theta_0,mid1_pi_1,F_w,tend_w)
CALL tendency_theta(mid1_u,mid1_w,rho_0,mid1_theta,F_theta,tend_theta)

! u-grid (its + 1:ite - 1, kts:kte)
imin = its + 1
imax = ite - 1
kmin = kts
kmax = kte

FORALL (i = imin:imax, k = kmin:kmax)
	mid2_u(i,k) = old_u(i,k) + DeltaT/2.*tend_u(i,k)
END FORALL

! w-grid (it + 1:ite, kts + 1:kte)
imin = its + 1
imax = ite
kmin = kts + 1
kmax = kte
! ATTENTION: The calculated area includes the boundary laye.
kmin = kmin - 1
kmax = kmax + 1

FORALL (i = imin:imax, k = kmin:kmax)
	mid2_w(i,k) = old_w(i,k) + DeltaT/2.*tend_w(i,k)
	mid2_theta(i,k) = old_theta(i,k) + DeltaT/2.*tend_theta(i,k)
	mid2_theta_1(i,k) = mid1_theta(i,k) - theta_0(i,k)
END FORALL

CALL update_boundary(mid2_u,mid2_w)

CALL tendency_pi(mid2_u,mid2_w,rho_0,theta_0,F_pi,tend_pi)

! pi-grid (its + 1:ite, kts:kte)
imin = its + 1
imax = ite
kmin = kts
kmax = kte

FORALL (i = imin:imax, k = kmin:kmax)
	mid2_pi_1(i,k) = old_pi_1(i,k) + DeltaT/2.*tend_pi(i,k)
END FORALL

CALL update_boundary(mid2_u,mid2_w,mid2_pi_1,mid2_theta,mid2_theta_1)

!=================================================
! Step 3. phi(n+1) = phi(n) + dt*tend(phi**)
!-------------------------------------------------
CALL basic_interpolate(mid2_u,mid2_w,mid2_pi_1,mid2_theta,theta_0,mid2_theta_1,rho_0)

CALL debug_undef_all( new_u,new_w,new_pi_1,new_theta,new_theta_1, &
                      F_u,rhou_pi,rhouu_pi,                       &
                      PrhouPx_u,PrhouuPx_u,                       &
                      rhow_v,rhouw_v,                             &
                      PrhowPz_u,PrhouwPz_u,                       &
                      F_w,rhou_v,PrhouPx_w,PrhouwPx_w,            &
                      rhow_pi,rhoww_pi,PrhowPz_w,PrhowwPz_w,      &
                      F_theta,rhoutheta_v,PrhowthetaPz_w,         &
                      Ppi_1Px_u,Ppi_1Pz_w,                        &
                      F_pi,urhotheta_u,wrhotheta_w,               &
                      PurhothetaPx_pi,PwrhothetaPz_pi,            &
                      tend_u,tend_w,tend_theta,tend_pi,           &
                      P2uPx2_u,P2uPz2_u,                          &
                      P2wPx2_w,P2wPz2_w,                          &
                      P2thetaPx2_w,P2thetaPz2_w                   )

CALL tendency_u(mid2_u,rho_0,mid2_pi_1,F_u,tend_u)
CALL tendency_w(mid2_w,rho_0,mid2_theta_1,theta_0,mid2_pi_1,F_w,tend_w)
CALL tendency_theta(mid2_u,mid2_w,rho_0,mid2_theta,F_theta,tend_theta)

! u-grid (its + 1:ite - 1, kts:kte)
imin = its + 1
imax = ite - 1
kmin = kts
kmax = kte

FORALL (i = imin:imax, k = kmin:kmax)
	new_u(i,k) = old_u(i,k) + DeltaT*tend_u(i,k)
END FORALL

! w-grid (it + 1:ite, kts + 1:kte)
imin = its + 1
imax = ite
kmin = kts + 1
kmax = kte
! ATTENTION: The calculated area includes the boundary laye.
kmin = kmin - 1
kmax = kmax + 1

FORALL (i = imin:imax, k = kmin:kmax)
	new_w(i,k) = old_w(i,k) + DeltaT*tend_w(i,k)
	new_theta(i,k) = old_theta(i,k) + DeltaT*tend_theta(i,k)
	new_theta_1(i,k) = mid1_theta(i,k) - theta_0(i,k)
END FORALL

CALL update_boundary(new_u,new_w)

CALL tendency_pi(new_u,new_w,rho_0,theta_0,F_pi,tend_pi)

! pi-grid (its + 1:ite, kts:kte)
imin = its + 1
imax = ite
kmin = kts
kmax = kte

FORALL (i = imin:imax, k = kmin:kmax)
	new_pi_1(i,k) = old_pi_1(i,k) + DeltaT*tend_pi(i,k)
END FORALL
!=================================================
END SUBROUTINE runge_kutta
!=================================================

!=================================================
! Debug mode.
!=================================================
SUBROUTINE debug_integrate( MODE,DeltaT,v,theta_0,rho_0,                &
                            old_u,old_w,old_pi_1,old_theta,old_theta_1, &
                            new_u,new_w,new_pi_1,new_theta,new_theta_1  )
IMPLICIT NONE
!=================================================
INTEGER, INTENT(IN) :: MODE
REAL(preci), INTENT(IN) :: DeltaT
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: v            ! wind speed along y-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_0      ! theta = theta_0 + theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0        ! density
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_u        ! wind speed along x-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_w        ! wind speed along z-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_pi_1     ! pi'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_theta_1  ! theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_u
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_w
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_pi_1
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_theta
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_theta_1
!-------------------------------------------------
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid1_u
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid1_w
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid1_pi_1
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid1_theta
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid1_theta_1
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid2_u
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid2_w
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid2_pi_1
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid2_theta
REAL(preci), DIMENSION(ims:ime,kms:kme) :: mid2_theta_1
!-------------------------------------------------
INTEGER :: i, k
!=================================================
SELECT CASE (MODE)
CASE (1) ! Forward-backward Scheme
	CALL basic_interpolate(old_u,old_w,old_pi_1,old_theta,theta_0,old_theta_1,rho_0)
	
	CALL debug_undef_all( new_u,new_w,new_pi_1,new_theta,new_theta_1, &
	                      F_u,rhou_pi,rhouu_pi,                       &
	                      PrhouPx_u,PrhouuPx_u,                       &
	                      rhow_v,rhouw_v,                             &
	                      PrhowPz_u,PrhouwPz_u,                       &
	                      F_w,rhou_v,PrhouPx_w,PrhouwPx_w,            &
	                      rhow_pi,rhoww_pi,PrhowPz_w,PrhowwPz_w,      &
	                      F_theta,rhoutheta_v,PrhowthetaPz_w,         &
	                      Ppi_1Px_u,Ppi_1Pz_w,                        &
	                      F_pi,urhotheta_u,wrhotheta_w,               &
	                      PurhothetaPx_pi,PwrhothetaPz_pi,            &
	                      tend_u,tend_w,tend_theta,tend_pi,           &
	                      P2uPx2_u,P2uPz2_u,                          &
	                      P2wPx2_w,P2wPz2_w,                          &
                          P2thetaPx2_w,P2thetaPz2_w                   )
	
	CALL tendency_u(old_u,rho_0,old_pi_1,F_u,tend_u)
	CALL tendency_w(old_w,rho_0,old_theta_1,theta_0,old_pi_1,F_w,tend_w)
	CALL tendency_theta(old_u,old_w,rho_0,old_theta,F_theta,tend_theta)
	
	! w-grid
	DO i = its + 1, ite
		DO k = kts + 1, kte
			new_theta(i,k) = old_theta(i,k) + DeltaT*tend_theta(i,k)
			new_theta_1(i,k) = new_theta(i,k) - theta_0(i,k)
		END DO
	END DO
	new_u = old_u
	new_w = old_w
	new_pi_1 = old_pi_1
CASE (2) ! Forward-backward Scheme: only advtion of theta.
	CALL debug_undef_all( u_pi,u_w,u_v,                               &
	                      w_pi,w_u,w_v,                               &
	                      theta_pi,theta_u,theta_v,                   &
	                      theta_0_pi,theta_0_u,theta_0_v,             &
	                      theta_1_pi,theta_1_u,theta_1_v,             &
	                      rho_0_u,rho_0_w,rho_0_v                     )
	! v-grid
	DO i = its + 1, ite - 1
		DO k = kts + 1, kte
			rho_0_v(i,k) = (rho_0(i,k) + rho_0(i,k-1) + rho_0(i+1,k) + rho_0(i+1,k-1))/4.
			u_v(i,k) = (old_u(i,k) + old_u(i,k-1))/2.
			theta_v(i,k) = (old_theta(i,k) + old_theta(i+1,k))/2.
			IF (old_u(i,k) == undef .OR. old_u(i,k-1) == undef) STOP "u_v is WRONG!!!"
			IF (old_theta(i,k) == undef .OR. old_theta(i+1,k) == undef) STOP "theta_v is WRONG!!!"
		END DO
	END DO
	! w-grid
	DO i = its + 1, ite
		DO k = kts + 1, kte
			rho_0_w(i,k) = (rho_0(i,k) + rho_0(i,k-1))/2.
			IF (rho_0(i,k) == undef .OR. rho_0(i,k-1) == undef) STOP "rho_0_w is WRONG!!!"
		END DO
	END DO
	! pi-grid
	DO i = its + 1, ite
		DO k = kts, kte
			w_pi(i,k) = (old_w(i,k) + old_w(i,k+1))/2.
			theta_pi(i,k) = (old_theta(i,k) + old_theta(i,k+1))/2.
			IF (old_w(i,k) == undef .OR. old_w(i,k+1) == undef) STOP "w_pi is WRONG!!!"
			IF (old_theta(i,k) == undef .OR. old_theta(i,k+1) == undef) STOP "theta_pi is WRONG!!!"
		END DO
	END DO
	!CALL basic_interpolate(old_u,old_w,old_pi_1,old_theta,theta_0,old_theta_1,rho_0)
	CALL debug_undef_all( rhou_v,rhoutheta_v,rhow_pi,rhowtheta_pi,                &
                          PrhouPx_w,PrhouthetaPx_w,PrhowPz_w,PrhowthetaPz_w,      &
                          F_theta,tend_theta,new_theta,new_theta_1,new_u,new_w,new_pi_1)
	! v-grid
	DO i = its + 1, ite - 1
		DO k = kts + 1, kte
			rhou_v(i,k) = rho_0_v(i,k)*u_v(i,k)
			rhoutheta_v(i,k) = rho_0_v(i,k)*u_v(i,k)*theta_v(i,k)
			IF (rho_0_v(i,k) == undef .OR. u_v(i,k) == undef .OR. theta_v(i,k) == undef) STOP "rhou_v is WRONG!!!"
		END DO
	END DO
	rhou_v(its:ite,kts) = rhou_v(its:ite,kts+1)      ! vertical boundary condition
	rhou_v(its:ite,kte+1) = rhou_v(its:ite,kte)      ! vertical boundary condition
	rhou_v(its,kts:kte+1) = rhou_v(its+1,kts:kte+1)  ! lateral boundary condition
	rhou_v(ite,kts:kte+1) = rhou_v(ite-1,kts:kte+1)  ! lateral boundary condition
	! pi-grid
	DO i = its + 1, ite
		DO k = kts, kte
			rhow_pi(i,k) = rho_0(i,k)*w_pi(i,k)
			rhowtheta_pi(i,k) = rho_0(i,k)*w_pi(i,k)*theta_pi(i,k)
			IF (rho_0(i,k) == undef .OR. w_pi(i,k) == undef .OR. theta_pi(i,k) == undef) STOP "rhou_v is WRONG!!!"
		END DO
	END DO
	rhoutheta_v(its:ite,kts) = rhoutheta_v(its:ite,kts+1)      ! vertical boundary condition
	rhoutheta_v(its:ite,kte+1) = rhoutheta_v(its:ite,kte)      ! vertical boundary condition
	rhoutheta_v(its,kts:kte+1) = rhoutheta_v(its+1,kts:kte+1)  ! lateral boundary condition
	rhoutheta_v(ite,kts:kte+1) = rhoutheta_v(ite-1,kts:kte+1)  ! lateral boundary condition
	! w-grid
	DO i = its + 1, ite
		DO k = kts + 1, kte
			PrhouPx_w(i,k) = (rhou_v(i,k) - rhou_v(i-1,k))/dx
			PrhouthetaPx_w(i,k) = (rhoutheta_v(i,k) - rhoutheta_v(i-1,k))/dx
			PrhowPz_w(i,k) = (rhow_pi(i,k) - rhow_pi(i,k-1))/dz
			PrhowthetaPz_w(i,k) = (rhowtheta_pi(i,k) - rhowtheta_pi(i,k-1))/dz

			IF (rhou_v(i,k) == undef .OR. rhou_v(i-1,k) == undef) STOP "PrhouPx_w is WRONG!!!"
			IF (rhoutheta_v(i,k) == undef .OR. rhoutheta_v(i-1,k) == undef) STOP "PrhouthetaPx_w is WRONG!!!"
			IF (rhow_pi(i,k) == undef .OR. rhow_pi(i,k-1) == undef) STOP "PrhowPz_w is WRONG!!!"
			IF (rhowtheta_pi(i,k) == undef .OR. rhowtheta_pi(i,k-1) == undef) STOP "PrhowthetaPz_w is WRONG!!!"

			F_theta(i,k) = - 1./rho_0_w(i,k)*(PrhouthetaPx_w(i,k) - old_theta(i,k)*PrhouPx_w(i,k) + PrhowthetaPz_w(i,k) - old_theta(i,k)*PrhowPz_w(i,k))
			tend_theta(i,k) = F_theta(i,k)
			new_theta(i,k) = old_theta(i,k) + DeltaT*tend_theta(i,k)
			!new_theta(i,k) = (old_theta(i-1,k) + old_theta(i+1,k))/2. + DeltaT*tend_theta(i,k) ! Lax Scheme
			!new_theta(i,k) = old_theta(i,k) + DeltaT*tend_theta(i,k) + (old_u(i,k)*DeltaT/dx)**2/2.*(old_theta(i+1,k) - 2*old_theta(i,k) + old_theta(i-1,k))! Lax-Wendorf Scheme
			new_theta_1(i,k) = new_theta(i,k) - theta_0(i,k)
		END DO
	END DO
	new_u = old_u
	new_w = old_w
	new_pi_1 = old_pi_1
	
CASE (3) ! Runge-Kutta Scheme
	!=================================================
	! Step 1. phi* = phi(n) + dt/3.*tend(phi(n))
	!-------------------------------------------------
	CALL basic_interpolate(old_u,old_w,old_pi_1,old_theta,theta_0,old_theta_1,rho_0)
	
	CALL debug_undef_all( mid1_u,mid1_w,mid1_pi_1,mid1_theta,mid1_theta_1, &
	                      F_u,rhou_pi,rhouu_pi,                       &
	                      PrhouPx_u,PrhouuPx_u,                       &
	                      rhow_v,rhouw_v,                             &
	                      PrhowPz_u,PrhouwPz_u,                       &
	                      F_w,rhou_v,PrhouPx_w,PrhouwPx_w,            &
	                      rhow_pi,rhoww_pi,PrhowPz_w,PrhowwPz_w,      &
	                      F_theta,rhoutheta_v,PrhowthetaPz_w,         &
	                      Ppi_1Px_u,Ppi_1Pz_w,                        &
	                      F_pi,urhotheta_u,wrhotheta_w,               &
	                      PurhothetaPx_pi,PwrhothetaPz_pi,            &
	                      tend_u,tend_w,tend_theta,tend_pi            )
	
	CALL tendency_u(old_u,rho_0,old_pi_1,F_u,tend_u)
	CALL tendency_w(old_w,rho_0,old_theta_1,theta_0,old_pi_1,F_w,tend_w)
	CALL tendency_theta(old_u,old_w,rho_0,old_theta,F_theta,tend_theta)
	! w-grid
	DO i = its + 1, ite
		DO k = kts + 1, kte
			mid1_theta(i,k) = old_theta(i,k) + DeltaT*tend_theta(i,k)
			mid1_theta_1(i,k) = new_theta(i,k) - theta_0(i,k)
		END DO
	END DO
	mid1_u = old_u
	mid1_w = old_w
	mid1_pi_1 = old_pi_1
	CALL update_boundary(mid1_u,mid1_w,mid1_pi_1,mid1_theta,mid1_theta_1)
	!=================================================
	! Step 2. phi** = phi(n) + dt/2.*tend(phi*)
	!-------------------------------------------------
	CALL basic_interpolate(mid1_u,mid1_w,mid1_pi_1,mid1_theta,theta_0,mid1_theta_1,rho_0)
	
	CALL debug_undef_all( mid2_u,mid2_w,mid2_pi_1,mid2_theta,mid2_theta_1, &
	                      F_u,rhou_pi,rhouu_pi,                       &
	                      PrhouPx_u,PrhouuPx_u,                       &
	                      rhow_v,rhouw_v,                             &
	                      PrhowPz_u,PrhouwPz_u,                       &
	                      F_w,rhou_v,PrhouPx_w,PrhouwPx_w,            &
	                      rhow_pi,rhoww_pi,PrhowPz_w,PrhowwPz_w,      &
	                      F_theta,rhoutheta_v,PrhowthetaPz_w,         &
	                      Ppi_1Px_u,Ppi_1Pz_w,                        &
	                      F_pi,urhotheta_u,wrhotheta_w,               &
	                      PurhothetaPx_pi,PwrhothetaPz_pi,            &
	                      tend_u,tend_w,tend_theta,tend_pi            )
	
	CALL tendency_u(mid1_u,rho_0,mid1_pi_1,F_u,tend_u)
	CALL tendency_w(mid1_w,rho_0,mid1_theta_1,theta_0,mid1_pi_1,F_w,tend_w)
	CALL tendency_theta(mid1_u,mid1_w,rho_0,mid1_theta,F_theta,tend_theta)
	! w-grid
	DO i = its + 1, ite
		DO k = kts + 1, kte
			mid2_theta(i,k) = old_theta(i,k) + DeltaT*tend_theta(i,k)
			mid2_theta_1(i,k) = new_theta(i,k) - theta_0(i,k)
		END DO
	END DO
	mid2_u = mid1_u
	mid2_w = mid1_w
	mid2_pi_1 = mid1_pi_1
	CALL update_boundary(mid2_u,mid2_w,mid2_pi_1,mid2_theta,mid2_theta_1)
	!=================================================
	! Step 3. phi(n+1) = phi(n) + dt*tend(phi**)
	!-------------------------------------------------
	CALL basic_interpolate(mid2_u,mid2_w,mid2_pi_1,mid2_theta,theta_0,mid2_theta_1,rho_0)
	
	CALL debug_undef_all( new_u,new_w,new_pi_1,new_theta,new_theta_1, &
	                      F_u,rhou_pi,rhouu_pi,                       &
	                      PrhouPx_u,PrhouuPx_u,                       &
	                      rhow_v,rhouw_v,                             &
	                      PrhowPz_u,PrhouwPz_u,                       &
	                      F_w,rhou_v,PrhouPx_w,PrhouwPx_w,            &
	                      rhow_pi,rhoww_pi,PrhowPz_w,PrhowwPz_w,      &
	                      F_theta,rhoutheta_v,PrhowthetaPz_w,         &
	                      Ppi_1Px_u,Ppi_1Pz_w,                        &
	                      F_pi,urhotheta_u,wrhotheta_w,               &
	                      PurhothetaPx_pi,PwrhothetaPz_pi,            &
	                      tend_u,tend_w,tend_theta,tend_pi            )
	
	CALL tendency_u(mid2_u,rho_0,mid2_pi_1,F_u,tend_u)
	CALL tendency_w(mid2_w,rho_0,mid2_theta_1,theta_0,mid2_pi_1,F_w,tend_w)
	CALL tendency_theta(mid2_u,mid2_w,rho_0,mid2_theta,F_theta,tend_theta)
	! w-grid
	DO i = its + 1, ite
		DO k = kts + 1, kte
			new_theta(i,k) = old_theta(i,k) + DeltaT*tend_theta(i,k)
			new_theta_1(i,k) = new_theta(i,k) - theta_0(i,k)
		END DO
	END DO
	new_u = mid2_u
	new_w = mid2_w
	new_pi_1 = mid2_pi_1
	CALL update_boundary(new_u,new_w,new_pi_1,new_theta,new_theta_1)
	
CASE DEFAULT
	STOP "Wrong time differencing scheme!!!"
END SELECT
!=================================================
END SUBROUTINE debug_integrate
!=================================================

!=================================================
END MODULE sp_module_timeschemes
!=================================================

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
SUBROUTINE forward_backward( DeltaT,v,theta_0,pi_0,rho_0,                &
                             old_u,old_w,old_pi_1,old_theta,old_theta_1, &
                             new_u,new_w,new_pi_1,new_theta,new_theta_1  )
IMPLICIT NONE
!=================================================
REAL(preci), INTENT(IN) :: DeltaT
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: v            ! wind speed along y-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_0      ! theta = theta_0 + theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_0
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
CALL tendency_pi(new_u,new_w,pi_0,rho_0,theta_0,F_pi,tend_pi)
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
! Runge-Kutta 3 order scheme.
!-------------------------------------------------
! phi* = phi(n) + dt/3.*tend(phi(n))
! phi** = phi(n) + dt/2.*tend(phi*)
! phi(n+1) = phi(n) + dt*tend(phi**)
!=================================================
SUBROUTINE runge_kutta( DeltaT,v,theta_0,pi_0,rho_0,                   &
                        old_u,old_w,old_pi_1,old_theta,old_theta_1,    &
                        new_u,new_w,new_pi_1,new_theta,new_theta_1     )
IMPLICIT NONE
!=================================================
REAL(preci), INTENT(IN) :: DeltaT
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: v
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_0
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_0
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
INTEGER :: i, k
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
CALL set_calc_area_u

FORALL (i = imin:imax, k = kmin:kmax)
	mid1_u(i,k) = old_u(i,k) + DeltaT/3.*tend_u(i,k)
END FORALL

! w-grid (it + 1:ite, kts + 1:kte)
CALL set_calc_area_w
! ATTENTION: The calculated area includes the boundary layers.
kmin = kmin - 1
kmax = kmax + 1

FORALL (i = imin:imax, k = kmin:kmax)
	mid1_w(i,k) = old_w(i,k) + DeltaT/3.*tend_w(i,k)
	mid1_theta(i,k) = old_theta(i,k) + DeltaT/3.*tend_theta(i,k)
	mid1_theta_1(i,k) = mid1_theta(i,k) - theta_0(i,k)
END FORALL

CALL update_boundary(mid1_u,mid1_w)

CALL tendency_pi(mid1_u,mid1_w,pi_0,rho_0,theta_0,F_pi,tend_pi)

! pi-grid (its + 1:ite, kts:kte)
CALL set_calc_area_pi

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
CALL set_calc_area_u

FORALL (i = imin:imax, k = kmin:kmax)
	mid2_u(i,k) = old_u(i,k) + DeltaT/2.*tend_u(i,k)
END FORALL

! w-grid (it + 1:ite, kts + 1:kte)
CALL set_calc_area_w
! ATTENTION: The calculated area includes the boundary laye.
kmin = kmin - 1
kmax = kmax + 1

FORALL (i = imin:imax, k = kmin:kmax)
	mid2_w(i,k) = old_w(i,k) + DeltaT/2.*tend_w(i,k)
	mid2_theta(i,k) = old_theta(i,k) + DeltaT/2.*tend_theta(i,k)
	mid2_theta_1(i,k) = mid1_theta(i,k) - theta_0(i,k)
END FORALL

CALL update_boundary(mid2_u,mid2_w)

CALL tendency_pi(mid2_u,mid2_w,pi_0,rho_0,theta_0,F_pi,tend_pi)

! pi-grid (its + 1:ite, kts:kte)
CALL set_calc_area_pi

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
CALL set_calc_area_u

FORALL (i = imin:imax, k = kmin:kmax)
	new_u(i,k) = old_u(i,k) + DeltaT*tend_u(i,k)
END FORALL

! w-grid (it + 1:ite, kts + 1:kte)
CALL set_calc_area_w
! ATTENTION: The calculated area includes the boundary laye.
kmin = kmin - 1
kmax = kmax + 1

FORALL (i = imin:imax, k = kmin:kmax)
	new_w(i,k) = old_w(i,k) + DeltaT*tend_w(i,k)
	new_theta(i,k) = old_theta(i,k) + DeltaT*tend_theta(i,k)
	new_theta_1(i,k) = mid1_theta(i,k) - theta_0(i,k)
END FORALL

CALL update_boundary(new_u,new_w)

CALL tendency_pi(new_u,new_w,pi_0,rho_0,theta_0,F_pi,tend_pi)

! pi-grid (its + 1:ite, kts:kte)
CALL set_calc_area_pi

FORALL (i = imin:imax, k = kmin:kmax)
	new_pi_1(i,k) = old_pi_1(i,k) + DeltaT*tend_pi(i,k)
END FORALL
!=================================================
END SUBROUTINE runge_kutta
!=================================================

!=================================================
END MODULE sp_module_timeschemes
!=================================================

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
REAL(kd), INTENT(IN) :: DeltaT
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: v            ! wind speed along y-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_0      ! theta = theta_0 + theta'
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_0
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0        ! density
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_u        ! wind speed along x-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_w        ! wind speed along z-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_pi_1     ! pi'
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_theta
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_theta_1  ! theta'
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_pi_1
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_theta
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_theta_1
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL basic_interpolate(old_u,old_w,old_pi_1,old_theta,old_theta_1)

CALL debug_undef_all(new_u,new_w,new_pi_1,new_theta,new_theta_1)

CALL debug_undef_all(   F_u,   F_w,   F_theta,   F_pi, &
                     tend_u,tend_w,tend_theta,tend_pi  )
CALL debug_undef_all(    P2uPx2_u,    P2uPz2_u,   &
                         P2wPx2_w,    P2wPz2_w,   &
                     P2thetaPx2_w,P2thetaPz2_w    )
CALL debug_undef_all( rhou_pi, rhouu_pi,          &
                      rhow_vir , rhowu_vir,           &
                      rhou_vir , rhouw_vir,           &
                      rhow_pi, rhoww_pi,          &
                      rhoutheta_vir, rhowtheta_pi,  &
                      urhotheta_u, wrhotheta_w    )
CALL debug_undef_all(PrhouPx_u, PrhouuPx_u, PrhowPz_u, PrhowuPz_u, Ppi_1Px_u,         &
                     PrhouPx_w, PrhouwPx_w, PrhowPz_w, PrhowwPz_w, Ppi_1Pz_w,         &
                     PrhouthetaPx_w, PrhowthetaPz_w, PurhothetaPx_pi, PwrhothetaPz_pi )

CALL tendency_u(old_u,rho_0,old_pi_1,F_u,tend_u)
CALL tendency_w(old_w,rho_0,old_theta_1,theta_0,old_pi_1,F_w,tend_w)
CALL tendency_theta(old_u,old_w,rho_0,old_theta,F_theta,tend_theta)
!=================================================
! Update u, w, theta
!-------------------------------------------------
! u-grid
!OMP PARALLEL DO
DO k = kts, kte
	DO i = its + 1, ite -1
		new_u(i,k) = old_u(i,k) + DeltaT*tend_u(i,k)
	END DO
END DO
!OMP END PARALLEL DO

! w-grid
!OMP PARALLEL DO
DO k = kts + 1, kte
	DO i = its + 1, ite
		new_w(i,k) = old_w(i,k) + DeltaT*tend_w(i,k)
		new_theta(i,k) = old_theta(i,k) + DeltaT*tend_theta(i,k)
		new_theta_1(i,k) = new_theta(i,k) - theta_0(i,k)
	END DO
END DO
!OMP END PARALLEL DO
CALL update_boundary(new_u,new_w)
!=================================================
! Update pi_1
!-------------------------------------------------
CALL tendency_pi(new_u,new_w,pi_0,rho_0,theta_0,F_pi,tend_pi)
! pi-grid
!OMP PARALLEL DO
DO k = kts, kte
	DO i = its + 1, ite
		new_pi_1(i,k) = old_pi_1(i,k) + DeltaT*tend_pi(i,k)
	END DO
END DO
!OMP END PARALLEL DO
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
REAL(kd), INTENT(IN) :: DeltaT
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: v
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_0
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_0
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: rho_0
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_pi_1
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_theta
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: old_theta_1
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_pi_1
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_theta
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: new_theta_1
!-------------------------------------------------
REAL(kd), DIMENSION(ims:ime,kms:kme) :: mid1_u
REAL(kd), DIMENSION(ims:ime,kms:kme) :: mid1_w
REAL(kd), DIMENSION(ims:ime,kms:kme) :: mid1_pi_1
REAL(kd), DIMENSION(ims:ime,kms:kme) :: mid1_theta
REAL(kd), DIMENSION(ims:ime,kms:kme) :: mid1_theta_1
REAL(kd), DIMENSION(ims:ime,kms:kme) :: mid2_u
REAL(kd), DIMENSION(ims:ime,kms:kme) :: mid2_w
REAL(kd), DIMENSION(ims:ime,kms:kme) :: mid2_pi_1
REAL(kd), DIMENSION(ims:ime,kms:kme) :: mid2_theta
REAL(kd), DIMENSION(ims:ime,kms:kme) :: mid2_theta_1
!-------------------------------------------------
INTEGER :: i, k
!=================================================
! Step 1. phi* = phi(n) + dt/3.*tend(phi(n))
!-------------------------------------------------
CALL basic_interpolate(old_u,old_w,old_pi_1,old_theta,old_theta_1)

CALL debug_undef_all(mid1_u,mid1_w,mid1_pi_1,mid1_theta,mid1_theta_1)

CALL debug_undef_all(   F_u,   F_w,   F_theta,   F_pi, &
                     tend_u,tend_w,tend_theta,tend_pi  )
CALL debug_undef_all(    P2uPx2_u,    P2uPz2_u,   &
                         P2wPx2_w,    P2wPz2_w,   &
                     P2thetaPx2_w,P2thetaPz2_w    )
CALL debug_undef_all( rhou_pi, rhouu_pi,          &
                      rhow_vir , rhowu_vir,           &
                      rhou_vir , rhouw_vir,           &
                      rhow_pi, rhoww_pi,          &
                      rhoutheta_vir, rhowtheta_pi,  &
                      urhotheta_u, wrhotheta_w    )
CALL debug_undef_all(PrhouPx_u, PrhouuPx_u, PrhowPz_u, PrhowuPz_u, Ppi_1Px_u,         &
                     PrhouPx_w, PrhouwPx_w, PrhowPz_w, PrhowwPz_w, Ppi_1Pz_w,         &
                     PrhouthetaPx_w, PrhowthetaPz_w, PurhothetaPx_pi, PwrhothetaPz_pi )

CALL tendency_u(old_u,rho_0,old_pi_1,F_u,tend_u)
CALL tendency_w(old_w,rho_0,old_theta_1,theta_0,old_pi_1,F_w,tend_w)
CALL tendency_theta(old_u,old_w,rho_0,old_theta,F_theta,tend_theta)

! u-grid (its + 1:ite - 1, kts:kte)
CALL set_area_u
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		mid1_u(i,k) = old_u(i,k) + DeltaT/3.*tend_u(i,k)
	END DO
END DO
!OMP END PARALLEL DO

! w-grid (it + 1:ite, kts + 1:kte)
CALL set_area_w

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		mid1_w(i,k) = old_w(i,k) + DeltaT/3.*tend_w(i,k)
		mid1_theta(i,k) = old_theta(i,k) + DeltaT/3.*tend_theta(i,k)
		mid1_theta_1(i,k) = mid1_theta(i,k) - theta_0(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL update_boundary(mid1_u,mid1_w)
CALL basic_interpolate(mid1_u,mid1_w,old_pi_1,old_theta,old_theta_1) !!!

CALL tendency_pi(mid1_u,mid1_w,pi_0,rho_0,theta_0,F_pi,tend_pi)

! pi-grid (its + 1:ite, kts:kte)
CALL set_area_pi

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		mid1_pi_1(i,k) = old_pi_1(i,k) + DeltaT/3.*tend_pi(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL update_boundary(mid1_u,mid1_w,mid1_pi_1,mid1_theta,mid1_theta_1)

!=================================================
! Step 2. phi** = phi(n) + dt/2.*tend(phi*)
!-------------------------------------------------
CALL basic_interpolate(mid1_u,mid1_w,mid1_pi_1,mid1_theta,mid1_theta_1)

CALL debug_undef_all( mid2_u,mid2_w,mid2_pi_1,mid2_theta,mid2_theta_1)

CALL debug_undef_all(   F_u,   F_w,   F_theta,   F_pi, &
                     tend_u,tend_w,tend_theta,tend_pi  )
CALL debug_undef_all(    P2uPx2_u,    P2uPz2_u,   &
                         P2wPx2_w,    P2wPz2_w,   &
                     P2thetaPx2_w,P2thetaPz2_w    )
CALL debug_undef_all( rhou_pi, rhouu_pi,          &
                      rhow_vir , rhowu_vir,           &
                      rhou_vir , rhouw_vir,           &
                      rhow_pi, rhoww_pi,          &
                      rhoutheta_vir, rhowtheta_pi,  &
                      urhotheta_u, wrhotheta_w    )
CALL debug_undef_all(PrhouPx_u, PrhouuPx_u, PrhowPz_u, PrhowuPz_u, Ppi_1Px_u,         &
                     PrhouPx_w, PrhouwPx_w, PrhowPz_w, PrhowwPz_w, Ppi_1Pz_w,         &
                     PrhouthetaPx_w, PrhowthetaPz_w, PurhothetaPx_pi, PwrhothetaPz_pi )

CALL tendency_u(mid1_u,rho_0,mid1_pi_1,F_u,tend_u)
CALL tendency_w(mid1_w,rho_0,mid1_theta_1,theta_0,mid1_pi_1,F_w,tend_w)
CALL tendency_theta(mid1_u,mid1_w,rho_0,mid1_theta,F_theta,tend_theta)

! u-grid (its + 1:ite - 1, kts:kte)
CALL set_area_u

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		mid2_u(i,k) = old_u(i,k) + DeltaT/2.*tend_u(i,k)
	END DO
END DO
!OMP END PARALLEL DO

! w-grid (it + 1:ite, kts + 1:kte)
CALL set_area_w

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		mid2_w(i,k) = old_w(i,k) + DeltaT/2.*tend_w(i,k)
		mid2_theta(i,k) = old_theta(i,k) + DeltaT/2.*tend_theta(i,k)
		mid2_theta_1(i,k) = mid1_theta(i,k) - theta_0(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL update_boundary(mid2_u,mid2_w)
CALL basic_interpolate(mid2_u,mid2_w,mid1_pi_1,mid1_theta,mid1_theta_1)  !!!

CALL tendency_pi(mid2_u,mid2_w,pi_0,rho_0,theta_0,F_pi,tend_pi)

! pi-grid (its + 1:ite, kts:kte)
CALL set_area_pi

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		mid2_pi_1(i,k) = old_pi_1(i,k) + DeltaT/2.*tend_pi(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL update_boundary(mid2_u,mid2_w,mid2_pi_1,mid2_theta,mid2_theta_1)

!=================================================
! Step 3. phi(n+1) = phi(n) + dt*tend(phi**)
!-------------------------------------------------
CALL basic_interpolate(mid2_u,mid2_w,mid2_pi_1,mid2_theta,mid2_theta_1)

CALL debug_undef_all( new_u,new_w,new_pi_1,new_theta,new_theta_1)

CALL debug_undef_all(   F_u,   F_w,   F_theta,   F_pi, &
                     tend_u,tend_w,tend_theta,tend_pi  )
CALL debug_undef_all(    P2uPx2_u,    P2uPz2_u,   &
                         P2wPx2_w,    P2wPz2_w,   &
                     P2thetaPx2_w,P2thetaPz2_w    )
CALL debug_undef_all( rhou_pi, rhouu_pi,          &
                      rhow_vir , rhowu_vir,           &
                      rhou_vir , rhouw_vir,           &
                      rhow_pi, rhoww_pi,          &
                      rhoutheta_vir, rhowtheta_pi,  &
                      urhotheta_u, wrhotheta_w    )
CALL debug_undef_all(PrhouPx_u, PrhouuPx_u, PrhowPz_u, PrhowuPz_u, Ppi_1Px_u,         &
                     PrhouPx_w, PrhouwPx_w, PrhowPz_w, PrhowwPz_w, Ppi_1Pz_w,         &
                     PrhouthetaPx_w, PrhowthetaPz_w, PurhothetaPx_pi, PwrhothetaPz_pi )

CALL tendency_u(mid2_u,rho_0,mid2_pi_1,F_u,tend_u)
CALL tendency_w(mid2_w,rho_0,mid2_theta_1,theta_0,mid2_pi_1,F_w,tend_w)
CALL tendency_theta(mid2_u,mid2_w,rho_0,mid2_theta,F_theta,tend_theta)

! u-grid (its + 1:ite - 1, kts:kte)
CALL set_area_u

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		new_u(i,k) = old_u(i,k) + DeltaT*tend_u(i,k)
	END DO
END DO
!OMP END PARALLEL DO

! w-grid (it + 1:ite, kts + 1:kte)
CALL set_area_w

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		new_w(i,k) = old_w(i,k) + DeltaT*tend_w(i,k)
		new_theta(i,k) = old_theta(i,k) + DeltaT*tend_theta(i,k)
		new_theta_1(i,k) = mid1_theta(i,k) - theta_0(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL update_boundary(new_u,new_w)
CALL basic_interpolate(new_u,new_w,mid2_pi_1,mid2_theta,mid2_theta_1)  !!!

CALL tendency_pi(new_u,new_w,pi_0,rho_0,theta_0,F_pi,tend_pi)

! pi-grid (its + 1:ite, kts:kte)
CALL set_area_pi

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		new_pi_1(i,k) = old_pi_1(i,k) + DeltaT*tend_pi(i,k)
	END DO
END DO
!OMP END PARALLEL DO
!=================================================
END SUBROUTINE runge_kutta
!=================================================

!=================================================
END MODULE sp_module_timeschemes
!=================================================

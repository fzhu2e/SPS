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
USE sp_module_gridvar
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
SUBROUTINE runge_kutta( uGrid, wGrid, piGrid, virGrid, new )
IMPLICIT NONE
TYPE(grid), INTENT(IN) :: uGrid, wGrid, piGrid, virGrid
TYPE(mainvar), INTENT(OUT) :: new
TYPE(mainvar) :: old, mid1, mid2
!=================================================
INTEGER :: i, k
!=================================================
! Step 1. phi* = phi(n) + dt/3.*tend(phi(n))
!-------------------------------------------------
old%u = uGrid%u
old%w = wGrid%w
old%pi_1 = piGrid%pi_1
old%theta = wGrid%theta
old%theta_1 = wGrid%theta_1

CALL update(old%u,old%w,old%pi_1,old%theta,old%theta_1,     &
            old%u,old%w,old%pi_1,old%theta,old%theta_1,     &
            mid1%u,mid1%w,mid1%pi_1,mid1%theta,mid1%theta_1,&
            3 )

!=================================================
! Step 2. phi** = phi(n) + dt/2.*tend(phi*)
!-------------------------------------------------
CALL update(old%u,old%w,old%pi_1,old%theta,old%theta_1,     &
            mid1%u,mid1%w,mid1%pi_1,mid1%theta,mid1%theta_1,&
            mid2%u,mid2%w,mid2%pi_1,mid2%theta,mid2%theta_1,&
            2 )

!=================================================
! Step 3. phi(n+1) = phi(n) + dt*tend(phi**)
!-------------------------------------------------
CALL update(old%u,old%w,old%pi_1,old%theta,old%theta_1,     &
            mid2%u,mid2%w,mid2%pi_1,mid2%theta,mid2%theta_1,&
            new%u,new%w,new%pi_1,new%theta,new%theta_1,     &
            1 )

!=================================================
END SUBROUTINE runge_kutta
!=================================================

!=================================================
SUBROUTINE update(A_u,A_w,A_pi_1,A_theta,A_theta_1, &
                  B_u,B_w,B_pi_1,B_theta,B_theta_1, &
                  C_u,C_w,C_pi_1,C_theta,C_theta_1, &
                  deno )
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: A_u, A_w, A_pi_1, A_theta, A_theta_1
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: B_u, B_w, B_pi_1, B_theta, B_theta_1
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: C_u, C_w, C_pi_1, C_theta, C_theta_1
INTEGER, INTENT(IN) :: deno
INTEGER :: i, k
!=================================================
CALL basic_interpolate(A_u,A_w,A_pi_1,A_theta,A_theta_1)

CALL debug_undef_all(C_u,C_w,C_pi_1,C_theta,C_theta_1)

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

CALL tendency_u(B_u,B_pi_1,tend_u)
CALL tendency_w(B_w,B_theta_1,B_pi_1,tend_w)
CALL tendency_theta(B_u,B_w,B_theta,tend_theta)

! u-grid (its + 1:ite - 1, kts:kte)
CALL set_area_u
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		C_u(i,k) = A_u(i,k) + dt/REAL(deno)*tend_u(i,k)
	END DO
END DO
!OMP END PARALLEL DO

! w-grid (it + 1:ite, kts + 1:kte)
CALL set_area_w

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		C_w(i,k) = A_w(i,k) + dt/REAL(deno)*tend_w(i,k)
		C_theta(i,k) = A_theta(i,k) + dt/REAL(deno)*tend_theta(i,k)
		C_theta_1(i,k) = C_theta(i,k) - theta_0(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL update_boundary(C_u,C_w)
CALL basic_interpolate(C_u,C_w,A_pi_1,A_theta,A_theta_1) !!!

CALL tendency_pi(C_u,C_w,tend_pi)

! pi-grid (its + 1:ite, kts:kte)
CALL set_area_pi

!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		C_pi_1(i,k) = A_pi_1(i,k) + dt/REAL(deno)*tend_pi(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL update_boundary(C_u,C_w,C_pi_1,C_theta,C_theta_1)

!=================================================
END SUBROUTINE update
!=================================================

!=================================================
END MODULE sp_module_timeschemes
!=================================================

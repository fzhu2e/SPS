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
SUBROUTINE runge_kutta(uGrid,wGrid,piGrid,virGrid,new)
IMPLICIT NONE
TYPE(grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
TYPE(mainvar), INTENT(OUT) :: new
TYPE(mainvar) :: old, mid1, mid2
!=================================================
INTEGER :: i, k
!=================================================
old%u = uGrid%u
old%w = wGrid%w
old%pi_1 = piGrid%pi_1
old%theta = wGrid%theta
old%theta_1 = wGrid%theta_1
!=================================================
! Step 1. phi* = phi(n) + dt/3.*tend(phi(n))
!-------------------------------------------------
CALL update(old,old,mid1,3,uGrid,wGrid,piGrid,virGrid)

!=================================================
! Step 2. phi** = phi(n) + dt/2.*tend(phi*)
!-------------------------------------------------
CALL update(old,mid1,mid2,2,uGrid,wGrid,piGrid,virGrid)

!=================================================
! Step 3. phi(n+1) = phi(n) + dt*tend(phi**)
!-------------------------------------------------
CALL update(old,mid2,new,1,uGrid,wGrid,piGrid,virGrid)

!=================================================
END SUBROUTINE runge_kutta
!=================================================

!=================================================
SUBROUTINE update(A,B,C,deno,uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE(mainvar), INTENT(IN) :: A, B
INTEGER, INTENT(IN) :: deno
TYPE(grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
TYPE(mainvar), INTENT(OUT) :: C
INTEGER :: i, k
!=================================================
CALL basic_interpolate(A,uGrid,wGrid,piGrid,virGrid)

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

CALL tendency_u(B%u,B%pi_1,tend_u,uGrid,wGrid,piGrid,virGrid)
CALL tendency_w(B%w,B%theta_1,B%pi_1,tend_w,uGrid,wGrid,piGrid,virGrid)
CALL tendency_theta(B%u,B%w,B%theta,tend_theta,uGrid,wGrid,piGrid,virGrid)

CALL set_area_u
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		C%u(i,k) = A%u(i,k) + dt/REAL(deno)*tend_u(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL set_area_w
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		C%w(i,k) = A%w(i,k) + dt/REAL(deno)*tend_w(i,k)
		C%theta(i,k) = A%theta(i,k) + dt/REAL(deno)*tend_theta(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL update_boundary(C%u,C%w)
CALL basic_interpolate(C,uGrid,wGrid,piGrid,virGrid)

CALL tendency_pi(C%u,C%w,tend_pi,uGrid,wGrid,piGrid,virGrid)

CALL set_area_pi
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		C%pi_1(i,k) = A%pi_1(i,k) + dt/REAL(deno)*tend_pi(i,k)
	END DO
END DO
!OMP END PARALLEL DO

CALL update_boundary(C%u,C%w,C%pi_1,C%theta)
C%theta_1 = C%theta - wGrid%theta_0
!=================================================
END SUBROUTINE update
!=================================================

!=================================================
END MODULE sp_module_timeschemes
!=================================================

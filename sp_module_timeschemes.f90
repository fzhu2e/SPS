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
old%qv = wGrid%qv
old%qc = wGrid%qc
old%qr = wGrid%qr
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

REAL(kd), DIMENSION(ims:ime,kms:kme) :: tend_u = undef, tend_w = undef
REAL(kd), DIMENSION(ims:ime,kms:kme) :: tend_pi_1 = undef, tend_theta = undef
REAL(kd), DIMENSION(ims:ime,kms:kme) :: tend_qv = undef, tend_qc = undef, tend_qr = undef

INTEGER :: i, k
!=================================================
C = A
CALL basic_interpolate(A,uGrid,wGrid,piGrid,virGrid)

!-------------------------------------------------
CALL tendency_u(B,tend_u,uGrid,wGrid,piGrid,virGrid)
CALL tendency_w(B,tend_w,uGrid,wGrid,piGrid,virGrid)

CALL tendency_theta(0,B,tend_theta,uGrid,wGrid,piGrid,virGrid)

CALL tendency_theta(1,B,tend_qv,uGrid,wGrid,piGrid,virGrid)
CALL tendency_theta(2,B,tend_qc,uGrid,wGrid,piGrid,virGrid)
CALL tendency_theta(3,B,tend_qr,uGrid,wGrid,piGrid,virGrid)

!-------------------------------------------------
CALL set_area_u
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		C%u(i,k) = A%u(i,k) + dt/REAL(deno)*tend_u(i,k)
	END DO
END DO
!OMP END PARALLEL DO

!-------------------------------------------------
CALL set_area_w
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		C%w(i,k) = A%w(i,k) + dt/REAL(deno)*tend_w(i,k)

		C%theta(i,k) = A%theta(i,k) + dt/REAL(deno)*tend_theta(i,k)
		wGrid%theta_M(i,k) = C%theta(i,k)*(1. + 0.61*wGrid%qv(i,k))*(1. - wGrid%qc(i,k))
		wGrid%theta_M_1(i,k) = wGrid%theta_M(i,k) - wGrid%theta_M_0(i,k)

		C%qv(i,k) = A%qv(i,k) + dt/REAL(deno)*tend_qv(i,k)
		C%qc(i,k) = A%qc(i,k) + dt/REAL(deno)*tend_qc(i,k)
		C%qr(i,k) = A%qr(i,k) + dt/REAL(deno)*tend_qr(i,k)
	END DO
END DO
!OMP END PARALLEL DO

!-------------------------------------------------
CALL update_boundary(C%u,C%w,wGrid=wGrid)
CALL basic_interpolate(C,uGrid,wGrid,piGrid,virGrid)
!CALL debug_SFSG

!-------------------------------------------------
CALL tendency_pi(C,tend_pi_1,uGrid,wGrid,piGrid,virGrid)

CALL set_area_pi
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		C%pi_1(i,k) = A%pi_1(i,k) + dt/REAL(deno)*tend_pi_1(i,k)
	END DO
END DO
!OMP END PARALLEL DO

!-------------------------------------------------
CALL update_boundary(C%u,C%w,C%pi_1,C%theta,C%qv,C%qc,C%qr,wGrid)
!=================================================
END SUBROUTINE update
!=================================================

!=================================================
END MODULE sp_module_timeschemes
!=================================================

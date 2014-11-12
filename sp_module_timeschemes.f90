!=================================================
! Super-Parametertization System (SPS)
!-------------------------------------------------
! Version: 0.2
! Author: Feng Zhu
! Email: zhuf.atmos@gmail.com
! Date: 2014-06-12 18:18:45
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
old%qi = wGrid%qi
old%qs = wGrid%qs
old%qg = wGrid%qg
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
REAL(kd), DIMENSION(ims:ime,kms:kme) :: tend_qi = undef, tend_qs = undef, tend_qg = undef

INTEGER :: i, k
!=================================================
C = A
CALL basic_interpolate(A,uGrid,wGrid,piGrid,virGrid)

!-------------------------------------------------
CALL tendency_u(B,tend_u,uGrid,wGrid,piGrid,virGrid)
CALL set_area_u
!$OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		C%u(i,k) = A%u(i,k) + dt/REAL(deno)*tend_u(i,k)
	END DO
END DO
!$OMP END PARALLEL DO

CALL tendency_w(B,tend_w,uGrid,wGrid,piGrid,virGrid)
CALL tendency_q(0,B%theta,tend_theta,uGrid,wGrid,piGrid,virGrid)

CALL set_area_w
!$OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		C%w(i,k) = A%w(i,k) + dt/REAL(deno)*tend_w(i,k)

		C%theta(i,k) = A%theta(i,k) + dt/REAL(deno)*tend_theta(i,k)
		wGrid%theta_M(i,k) = C%theta(i,k)*(1. + 0.61*wGrid%qv(i,k))*(1. - wGrid%qc(i,k) - wGrid%qi(i,k))
		wGrid%theta_M_1(i,k) = wGrid%theta_M(i,k) - wGrid%theta_M_0(i,k)
	END DO
END DO
!$OMP END PARALLEL DO

IF (Vapor /= 0) THEN
	CALL tendency_q(1,B%qv,tend_qv,uGrid,wGrid,piGrid,virGrid)
	CALL tendency_q(2,B%qc,tend_qc,uGrid,wGrid,piGrid,virGrid)
	CALL tendency_q(3,B%qr,tend_qr,uGrid,wGrid,piGrid,virGrid)
	CALL tendency_q(4,B%qi,tend_qi,uGrid,wGrid,piGrid,virGrid)
	CALL tendency_q(5,B%qs,tend_qs,uGrid,wGrid,piGrid,virGrid)
	CALL tendency_q(6,B%qg,tend_qg,uGrid,wGrid,piGrid,virGrid)
END IF
IF (Vapor == 0) THEN
	C%qv = 0.
	C%qc = 0.
	C%qr = 0.
	C%qi = 0.
	C%qs = 0.
	C%qg = 0.
ELSE
	!$OMP PARALLEL DO
	DO k = kmin, kmax
		DO i = imin, imax
			C%qv(i,k) = A%qv(i,k) + dt/REAL(deno)*tend_qv(i,k)
			C%qc(i,k) = A%qc(i,k) + dt/REAL(deno)*tend_qc(i,k)
			C%qr(i,k) = A%qr(i,k) + dt/REAL(deno)*tend_qr(i,k)
			C%qi(i,k) = A%qi(i,k) + dt/REAL(deno)*tend_qi(i,k)
			C%qs(i,k) = A%qs(i,k) + dt/REAL(deno)*tend_qs(i,k)
			C%qg(i,k) = A%qg(i,k) + dt/REAL(deno)*tend_qg(i,k)
		END DO
	END DO
	!$OMP END PARALLEL DO
END IF

!-------------------------------------------------
CALL update_boundary(C%u,C%w,wGrid)
CALL basic_interpolate(C,uGrid,wGrid,piGrid,virGrid)

CALL tendency_pi(C,tend_pi_1,uGrid,wGrid,piGrid,virGrid)
CALL set_area_pi
!$OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		C%pi_1(i,k) = A%pi_1(i,k) + dt/REAL(deno)*tend_pi_1(i,k)
	END DO
END DO
!$OMP END PARALLEL DO

!-------------------------------------------------
CALL update_boundary(C%u,C%w,wGrid,C%pi_1,C%theta,&
                     C%qv,C%qc,C%qr,C%qi,C%qs,C%qg)
!=================================================
END SUBROUTINE update
!=================================================

!=================================================
END MODULE sp_module_timeschemes
!=================================================

!=================================================
! The physics module of SPS-dynamic-integrate.
!-------------------------------------------------
! Version: 0.11
! Author: Zhu F.
! Email: lyricorpse@gmail.com
! Date: 2013-05-16 19:08:57 
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!=================================================
MODULE sp_module_physics
USE sp_module_constant
USE sp_module_model
USE sp_module_gridvar
USE sp_module_debug
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================

!=================================================
! Microphysics - WSM6 
!=================================================
SUBROUTINE subgrid(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
!-------------------------------------------------
TYPE(grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid

REAL(kd) :: temp_a, temp_b
INTEGER :: i, k
!=================================================
CALL set_area_u
!OMP PARALLEL DO PRIVATE(temp_a,temp_b)
DO k = kmin, kmax
	DO i = imin, imax
		temp_a = (uGrid%u(i+1,k) + uGrid%u(i-1,k) - 2*uGrid%u(i,k))/dx/dx
		temp_b = (uGrid%u(i,k+1) + uGrid%u(i,k-1) - 2*uGrid%u(i,k))/dz/dz
		uGrid%Du(i,k) = Km*(temp_a + temp_b)
	END DO
END DO
!OMP END PARALLEL DO

CALL set_area_w
!OMP PARALLEL DO PRIVATE(temp_a,temp_b)
DO k = kmin, kmax
	DO i = imin, imax
		temp_a = (wGrid%w(i+1,k) + wGrid%w(i-1,k) - 2*wGrid%w(i,k))/dx/dx
		temp_b = (wGrid%w(i,k+1) + wGrid%w(i,k-1) - 2*wGrid%w(i,k))/dz/dz
		wGrid%Dw(i,k) = Km*(temp_a + temp_b)

		temp_a = (wGrid%theta(i+1,k) + wGrid%theta(i-1,k) - 2*wGrid%theta(i,k))/dx/dx
		temp_b = (wGrid%theta(i,k+1) + wGrid%theta(i,k-1) - 2*wGrid%theta(i,k))/dz/dz
		wGrid%Dtheta(i,k) = Kh*(temp_a + temp_b)

		temp_a = (wGrid%qv(i+1,k) + wGrid%qv(i-1,k) - 2*wGrid%qv(i,k))/dx/dx
		temp_b = (wGrid%qv(i,k+1) + wGrid%qv(i,k-1) - 2*wGrid%qv(i,k))/dz/dz
		wGrid%Dqv(i,k) = Kh*(temp_a + temp_b)

		temp_a = (wGrid%qc(i+1,k) + wGrid%qc(i-1,k) - 2*wGrid%qc(i,k))/dx/dx
		temp_b = (wGrid%qc(i,k+1) + wGrid%qc(i,k-1) - 2*wGrid%qc(i,k))/dz/dz
		wGrid%Dqc(i,k) = Kh*(temp_a + temp_b)

		temp_a = (wGrid%qr(i+1,k) + wGrid%qr(i-1,k) - 2*wGrid%qr(i,k))/dx/dx
		temp_b = (wGrid%qr(i,k+1) + wGrid%qr(i,k-1) - 2*wGrid%qr(i,k))/dz/dz
		wGrid%Dqr(i,k) = Kh*(temp_a + temp_b)

		temp_a = (wGrid%qi(i+1,k) + wGrid%qi(i-1,k) - 2*wGrid%qi(i,k))/dx/dx
		temp_b = (wGrid%qi(i,k+1) + wGrid%qi(i,k-1) - 2*wGrid%qi(i,k))/dz/dz
		wGrid%Dqi(i,k) = Kh*(temp_a + temp_b)

		temp_a = (wGrid%qs(i+1,k) + wGrid%qs(i-1,k) - 2*wGrid%qs(i,k))/dx/dx
		temp_b = (wGrid%qs(i,k+1) + wGrid%qs(i,k-1) - 2*wGrid%qs(i,k))/dz/dz
		wGrid%Dqs(i,k) = Kh*(temp_a + temp_b)

		temp_a = (wGrid%qg(i+1,k) + wGrid%qg(i-1,k) - 2*wGrid%qg(i,k))/dx/dx
		temp_b = (wGrid%qg(i,k+1) + wGrid%qg(i,k-1) - 2*wGrid%qg(i,k))/dz/dz
		wGrid%Dqg(i,k) = Kh*(temp_a + temp_b)
	END DO
END DO
!OMP END PARALLEL DO
!=================================================
END SUBROUTINE subgrid
!=================================================

!=================================================
! Microphysics - WSM6 
!=================================================
SUBROUTINE mp_wsm6(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
!-------------------------------------------------
TYPE(grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
INTEGER :: i, k
!=================================================
wGrid%Mtheta = 0.
wGrid%Mqv = 0.
wGrid%Mqc = 0.
wGrid%Mqr = 0.
wGrid%Mqi = 0.
wGrid%Mqs = 0.
wGrid%Mqg = 0.
!=================================================
END SUBROUTINE mp_wsm6
!=================================================

!=================================================
END MODULE sp_module_physics
!=================================================

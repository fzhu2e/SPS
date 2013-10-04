!=================================================
! The math module of SPS-dynamic-integrate.
!-------------------------------------------------
! Version: 0.11
! Author: Zhu F.
! Email: lyricorpse@gmail.com
! Date: 2013-05-04 20:05:59 
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!=================================================
MODULE sp_module_math
USE sp_module_constant
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================

!=================================================
! Tri-Diagonal Matrix (no measure with aberant matrix)
!=================================================
SUBROUTINE TDMatrix(a,b,c,d,x,x_start,x_end,z_start,z_end)
IMPLICIT NONE
INTEGER, INTENT(IN) :: x_start, x_end, z_start, z_end
REAL(preci), DIMENSION(x_start:x_end,z_start:z_end), INTENT(INOUT) :: a, b, c, d
REAL(preci), DIMENSION(x_start:x_end,z_start:z_end), INTENT(INOUT) :: x
!-------------------------------------------------
!REAL, DIMENSION(x_start:x_end,z_start:z_end) :: temp
REAL, DIMENSION(x_start:x_end,z_start+1:z_end) :: m
INTEGER :: k
!=================================================
!WRITE(*,"(6f11.5)") b(2,1),  c(2,1),  0     ,  0     ,  0     ,  d(2,1)
!WRITE(*,"(6f11.5)") a(2,2),  b(2,2),  c(2,2),  0     ,  0     ,  d(2,2)
!WRITE(*,"(6f11.5)") 0     ,  a(2,3),  b(2,3),  c(2,3),  0     ,  d(2,3)  
!WRITE(*,"(6f11.5)") 0     ,  0     ,  a(2,4),  b(2,4),  c(2,4),  d(2,4)
!WRITE(*,"(6f11.5)") 0     ,  0     ,  0     ,  a(2,4),  b(2,4),  d(2,4)
DO k = z_start, z_end - 1
	m(:,k+1) = a(:,k+1)/b(:,k)
	a(:,k+1) = 0.
	b(:,k+1) = b(:,k+1) - c(:,k)*m(:,k+1)
	d(:,k+1) = d(:,k+1) - d(:,k)*m(:,k+1)
END DO

x(:,z_end) = d(:,z_end)/b(:,z_end)

DO k = z_end-1, z_start, -1
	x(:,k) = (d(:,k)-c(:,k)*x(:,k+1))/b(:,k)
END DO
!WRITE(*,*) "=========================="
!WRITE(*,"(6f11.5)") b(2,1),  c(2,1),  0     ,  0     ,  0     ,  d(2,1)
!WRITE(*,"(6f11.5)") a(2,2),  b(2,2),  c(2,2),  0     ,  0     ,  d(2,2)
!WRITE(*,"(6f11.5)") 0     ,  a(2,3),  b(2,3),  c(2,3),  0     ,  d(2,3)  
!WRITE(*,"(6f11.5)") 0     ,  0     ,  a(2,4),  b(2,4),  c(2,4),  d(2,4)
!WRITE(*,"(6f11.5)") 0     ,  0     ,  0     ,  a(2,4),  b(2,4),  d(2,4)
!WRITE(*,*) "=========================="
!WRITE(*,"(5f11.5)") x(2,1),  x(2,2),  x(2,3), x(2,4), x(2,5)
!=================================================
END SUBROUTINE TDMatrix
!=================================================

!=================================================
END MODULE sp_module_math
!=================================================

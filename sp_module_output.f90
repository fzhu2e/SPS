!=================================================
! The output module of SPS-dynamic.
!-------------------------------------------------
! Version: 0.01
! Author: Zhu F.
! Email: lyricorpse@gmail.com
! Date: 2013-04-20 12:20:45 
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!=================================================
MODULE sp_module_output
USE sp_module_constant
USE sp_module_model
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================
! Output the fields.
!=================================================
SUBROUTINE output(flag,u,w,theta_1,pi_1)
IMPLICIT NONE
!-------------------------------------------------
INTEGER,INTENT(IN) :: flag
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: u        ! wind speed along x-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: w        ! wind speed along z-axis
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_1  ! theta'
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_1     ! pi'
!=================================================
SELECT CASE(flag)
CASE (0)
	OPEN(1, FILE="./output/modelvar_u.bin", FORM='binary', CONVERT='big_endian')
	OPEN(2, FILE="./output/modelvar_w.bin", FORM='binary', CONVERT='big_endian')
	OPEN(3, FILE="./output/modelvar_theta_1.bin", FORM='binary', CONVERT='big_endian')
	OPEN(4, FILE="./output/modelvar_pi_1.bin", FORM='binary', CONVERT='big_endian')
	WRITE(1) u(its:ite,kts:kte)
	WRITE(2) w(its:ite,kts:kte)
	WRITE(3) theta_1(its:ite,kts:kte)
	WRITE(4) pi_1(its:ite,kts:kte)
	!WRITE(1) u(ims:ime,kms:kme)
	!WRITE(2) w(ims:ime,kms:kme)
	!WRITE(3) theta_1(ims:ime,kms:kme)
	!WRITE(4) pi_1(its:ime,kms:kme)
CASE (1)
	WRITE(1) u(its:ite,kts:kte)
	WRITE(2) w(its:ite,kts:kte)
	WRITE(3) theta_1(its:ite,kts:kte)
	WRITE(4) pi_1(its:ite,kts:kte)
	!WRITE(1) u(ims:ime,kms:kme)
	!WRITE(2) w(ims:ime,kms:kme)
	!WRITE(3) theta_1(ims:ime,kms:kme)
	!WRITE(4) pi_1(its:ime,kms:kme)
CASE (99)
	CLOSE(1)
	CLOSE(2)
	CLOSE(3)
	CLOSE(4)
CASE DEFAULT
	WRITE(*,*) "==================================="
	WRITE(*,*) " WARNING: Illegal flag of output."
	WRITE(*,*) "==================================="
END SELECT
!=================================================
END SUBROUTINE output
!=================================================
END MODULE sp_module_output
!=================================================

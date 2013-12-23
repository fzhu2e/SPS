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
SUBROUTINE output(flag,u,w,pi_1,theta_M_1,theta_M,theta,qv,qc,qr,qi,qs,qg)
IMPLICIT NONE
!-------------------------------------------------
INTEGER,INTENT(IN) :: flag
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: u        ! wind speed along x-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: w        ! wind speed along z-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_1     ! pi'
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_M_1
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_M
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: qv,qc,qr,qi,qs,qg
!=================================================
SELECT CASE(flag)
CASE (0)
	OPEN(1, FILE="./output/modelvar_u.bin", FORM='binary', CONVERT='big_endian')
	WRITE(1) u

	OPEN(2, FILE="./output/modelvar_w.bin", FORM='binary', CONVERT='big_endian')
	WRITE(2) w

	OPEN(3, FILE="./output/modelvar_pi_1.bin", FORM='binary', CONVERT='big_endian')
	WRITE(3) pi_1

	OPEN(4, FILE="./output/modelvar_theta_M_1.bin", FORM='binary', CONVERT='big_endian')
	WRITE(4) theta_M_1

	OPEN(5, FILE="./output/modelvar_theta_M.bin", FORM='binary', CONVERT='big_endian')
	WRITE(5) theta_M

	OPEN(6, FILE="./output/modelvar_theta.bin", FORM='binary', CONVERT='big_endian')
	WRITE(6) theta

	OPEN(7, FILE="./output/modelvar_qv.bin", FORM='binary', CONVERT='big_endian')
	WRITE(7) qv

	OPEN(8, FILE="./output/modelvar_qc.bin", FORM='binary', CONVERT='big_endian')
	WRITE(8) qc

	OPEN(9, FILE="./output/modelvar_qr.bin", FORM='binary', CONVERT='big_endian')
	WRITE(9) qr

	OPEN(10, FILE="./output/modelvar_qi.bin", FORM='binary', CONVERT='big_endian')
	WRITE(10) qi

	OPEN(11, FILE="./output/modelvar_qs.bin", FORM='binary', CONVERT='big_endian')
	WRITE(11) qs

	OPEN(12, FILE="./output/modelvar_qg.bin", FORM='binary', CONVERT='big_endian')
	WRITE(12) qg
CASE (1)
	WRITE(1) u
	WRITE(2) w
	WRITE(3) pi_1
	WRITE(4) theta_M_1
	WRITE(5) theta_M
	WRITE(6) theta
	WRITE(7) qv
	WRITE(8) qc
	WRITE(9) qr
	WRITE(10) qi
	WRITE(11) qs
	WRITE(12) qg
CASE (99)
	CLOSE(1)
	CLOSE(2)
	CLOSE(3)
	CLOSE(4)
	CLOSE(5)
	CLOSE(6)
	CLOSE(7)
	CLOSE(8)
	CLOSE(9)
	CLOSE(10)
	CLOSE(11)
	CLOSE(12)
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

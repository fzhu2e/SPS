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
SUBROUTINE output(flag,u,w,pi_1,theta_1,qv,qc,qr)
IMPLICIT NONE
!-------------------------------------------------
INTEGER,INTENT(IN) :: flag
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN), OPTIONAL :: u        ! wind speed along x-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN), OPTIONAL :: w        ! wind speed along z-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN), OPTIONAL :: pi_1     ! pi'
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN), OPTIONAL :: theta_1  ! theta'
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN), OPTIONAL :: qv,qc,qr
!=================================================
SELECT CASE(flag)
CASE (0)
	IF (PRESENT(u)) THEN
		OPEN(1, FILE="./output/modelvar_u.bin", FORM='binary', CONVERT='big_endian')
		WRITE(1) u(its:ite,kts:kte)
	END IF
	IF (PRESENT(w)) THEN
		OPEN(2, FILE="./output/modelvar_w.bin", FORM='binary', CONVERT='big_endian')
		WRITE(2) w(its:ite,kts:kte)
	END IF
	IF (PRESENT(pi_1)) THEN
		OPEN(3, FILE="./output/modelvar_pi_1.bin", FORM='binary', CONVERT='big_endian')
		WRITE(3) pi_1(its:ite,kts:kte)
	END IF
	IF (PRESENT(theta_1)) THEN
		OPEN(4, FILE="./output/modelvar_theta_1.bin", FORM='binary', CONVERT='big_endian')
		WRITE(4) theta_1(its:ite,kts:kte)
	END IF
	IF (PRESENT(qv)) THEN
		OPEN(5, FILE="./output/modelvar_qv.bin", FORM='binary', CONVERT='big_endian')
		WRITE(5) qv(its:ite,kts:kte)
	END IF
	IF (PRESENT(qc)) THEN
		OPEN(6, FILE="./output/modelvar_qc.bin", FORM='binary', CONVERT='big_endian')
		WRITE(6) qc(its:ite,kts:kte)
	END IF
	IF (PRESENT(qr)) THEN
		OPEN(7, FILE="./output/modelvar_qr.bin", FORM='binary', CONVERT='big_endian')
		WRITE(7) qr(its:ite,kts:kte)
	END IF
CASE (1)
	IF (PRESENT(u)) THEN
		WRITE(1) u(its:ite,kts:kte)
	END IF
	IF (PRESENT(w)) THEN
		WRITE(2) w(its:ite,kts:kte)
	END IF
	IF (PRESENT(pi_1)) THEN
		WRITE(3) pi_1(its:ite,kts:kte)
	END IF
	IF (PRESENT(theta_1)) THEN
		WRITE(4) theta_1(its:ite,kts:kte)
	END IF
	IF (PRESENT(qv)) THEN
		WRITE(5) qv(its:ite,kts:kte)
	END IF
	IF (PRESENT(qc)) THEN
		WRITE(6) qc(its:ite,kts:kte)
	END IF
	IF (PRESENT(qr)) THEN
		WRITE(7) qr(its:ite,kts:kte)
	END IF
CASE (99)
	IF (PRESENT(u)) THEN
		CLOSE(1)
	END IF
	IF (PRESENT(w)) THEN
		CLOSE(2)
	END IF
	IF (PRESENT(pi_1)) THEN
		CLOSE(3)
	END IF
	IF (PRESENT(theta_1)) THEN
		CLOSE(4)
	END IF
	IF (PRESENT(qv)) THEN
		CLOSE(5)
	END IF
	IF (PRESENT(qc)) THEN
		CLOSE(6)
	END IF
	IF (PRESENT(qr)) THEN
		CLOSE(7)
	END IF
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

!=================================================
! Super-Parametertization System (SPS)
!-------------------------------------------------
! Version: 0.2
! Author: Feng Zhu
! Email: zhuf.atmos@gmail.com
! Date: 2014-06-12 18:18:45
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
SUBROUTINE output(flag,u,w,pi_1,theta_1,theta_M_1,theta_M,theta, &
                  qv,qc,qr,qi,qs,qg,rain,snow,graupel,cldfra)
IMPLICIT NONE
!-------------------------------------------------
INTEGER,INTENT(IN) :: flag
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: u        ! wind speed along x-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: w        ! wind speed along z-axis
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: pi_1     ! pi'
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_1
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_M_1
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta_M
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: theta
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN), OPTIONAL :: qv,qc,qr,qi,qs,qg
REAL(kd), DIMENSION(ims:ime), INTENT(IN), OPTIONAL :: rain, snow, graupel
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN), OPTIONAL :: cldfra
!=================================================
SELECT CASE(flag)
CASE (0)
	OPEN(1, FILE="./output/modelvar_u.bin", FORM='binary', CONVERT='big_endian')
	WRITE(1) u

	OPEN(2, FILE="./output/modelvar_w.bin", FORM='binary', CONVERT='big_endian')
	WRITE(2) w

	OPEN(3, FILE="./output/modelvar_pi_1.bin", FORM='binary', CONVERT='big_endian')
	WRITE(3) pi_1

	OPEN(4, FILE="./output/modelvar_theta_1.bin", FORM='binary', CONVERT='big_endian')
	WRITE(4) theta_1

	OPEN(5, FILE="./output/modelvar_theta_M_1.bin", FORM='binary', CONVERT='big_endian')
	WRITE(5) theta_M_1

	OPEN(6, FILE="./output/modelvar_theta_M.bin", FORM='binary', CONVERT='big_endian')
	WRITE(6) theta_M

	OPEN(7, FILE="./output/modelvar_theta.bin", FORM='binary', CONVERT='big_endian')
	WRITE(7) theta

	IF (PRESENT(qv)) THEN
		OPEN(8, FILE="./output/modelvar_qv.bin", FORM='binary', CONVERT='big_endian')
		WRITE(8) qv
	END IF

	IF (PRESENT(qc)) THEN
		OPEN(9, FILE="./output/modelvar_qc.bin", FORM='binary', CONVERT='big_endian')
		WRITE(9) qc
	END IF

	IF (PRESENT(qr)) THEN
		OPEN(10, FILE="./output/modelvar_qr.bin", FORM='binary', CONVERT='big_endian')
		WRITE(10) qr
	END IF

	IF (PRESENT(qi)) THEN
		OPEN(11, FILE="./output/modelvar_qi.bin", FORM='binary', CONVERT='big_endian')
		WRITE(11) qi
	END IF

	IF (PRESENT(qs)) THEN
		OPEN(12, FILE="./output/modelvar_qs.bin", FORM='binary', CONVERT='big_endian')
		WRITE(12) qs
	END IF

	IF (PRESENT(qg)) THEN
		OPEN(13, FILE="./output/modelvar_qg.bin", FORM='binary', CONVERT='big_endian')
		WRITE(13) qg
	END IF

	IF (PRESENT(rain)) THEN
		OPEN(14, FILE="./output/modelvar_rain.bin", FORM='binary', CONVERT='big_endian')
		WRITE(14) rain
	END IF

	IF (PRESENT(snow)) THEN
		OPEN(15, FILE="./output/modelvar_snow.bin", FORM='binary', CONVERT='big_endian')
		WRITE(15) snow
	END IF

	IF (PRESENT(graupel)) THEN
		OPEN(16, FILE="./output/modelvar_graupel.bin", FORM='binary', CONVERT='big_endian')
		WRITE(16) graupel
	END IF

	IF (PRESENT(cldfra)) THEN
		OPEN(17, FILE="./output/modelvar_cldfra.bin", FORM='binary', CONVERT='big_endian')
		WRITE(17) cldfra
	END IF

CASE (1)
	WRITE(1) u
	WRITE(2) w
	WRITE(3) pi_1
	WRITE(4) theta_1
	WRITE(5) theta_M_1
	WRITE(6) theta_M
	WRITE(7) theta
	IF (PRESENT(qv)) WRITE(8) qv
	IF (PRESENT(qc)) WRITE(9) qc
	IF (PRESENT(qr)) WRITE(10) qr
	IF (PRESENT(qi)) WRITE(11) qi
	IF (PRESENT(qs)) WRITE(12) qs
	IF (PRESENT(qg)) WRITE(13) qg
	IF (PRESENT(rain)) WRITE(14) rain
	IF (PRESENT(snow)) WRITE(15) snow
	IF (PRESENT(graupel)) WRITE(16) graupel
	IF (PRESENT(graupel)) WRITE(17) cldfra
CASE (99)
	CLOSE(1)
	CLOSE(2)
	CLOSE(3)
	CLOSE(4)
	CLOSE(5)
	CLOSE(6)
	CLOSE(7)
	IF (PRESENT(qv)) CLOSE(8)
	IF (PRESENT(qc)) CLOSE(9)
	IF (PRESENT(qr)) CLOSE(10)
	IF (PRESENT(qi)) CLOSE(11)
	IF (PRESENT(qs)) CLOSE(12)
	IF (PRESENT(qg)) CLOSE(13)
	IF (PRESENT(rain)) CLOSE(14)
	IF (PRESENT(snow)) CLOSE(15)
	IF (PRESENT(graupel)) CLOSE(16)
	IF (PRESENT(cldfra)) CLOSE(17)
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

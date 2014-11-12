!=================================================
! Super-Parametertization System (SPS)
!-------------------------------------------------
! Version: 0.2
! Author: Feng Zhu
! Email: zhuf.atmos@gmail.com
! Date: 2014-06-12 18:18:45
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!=================================================
MODULE sp_module_cldfra
USE sp_module_constant
USE sp_module_model
USE sp_module_gridvar
USE sp_module_debug
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================

!=================================================
SUBROUTINE calc_cldfra(wGrid)
IMPLICIT NONE
!-------------------------------------------------
TYPE(grid), INTENT(INOUT) :: wGrid
!-------------------------------------------------
REAL(kd), PARAMETER :: thresh = 1.0e-6
INTEGER :: i, k
!=================================================
CALL set_area_w
!$OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		IF (wGrid%qc(i,k) + wGrid%qi(i,k) > thresh) THEN
			wGrid%cldfra(i,k) = 1.
		ELSE
			wGrid%cldfra(i,k) = 0.
		END IF
	END DO
END DO
!$OMP END PARALLEL DO
!=================================================
END SUBROUTINE calc_cldfra
!=================================================

!=================================================
SUBROUTINE calc_cldfra2(wGrid)
IMPLICIT NONE
!-------------------------------------------------
TYPE(grid), INTENT(INOUT) :: wGrid
INTEGER :: i, k
!-------------------------------------------------
REAL(kd) :: RHUM, tc, esw, esi, weight, qvsw, qvsi, qvs_weight, QIMID, QWMID, QCLD, DENOM, ARG, SUBSAT

REAL(kd), PARAMETER :: ALPHA0=100., GAMMA=0.49, QCLDMIN=1.E-12, PEXP=0.25, RHGRID=1.0
REAL(kd), PARAMETER :: SVP1=0.61078
REAL(kd), PARAMETER :: SVP2=17.2693882
REAL(kd), PARAMETER :: SVPI2=21.8745584
REAL(kd), PARAMETER :: SVP3=35.86
REAL(kd), PARAMETER :: SVPI3=7.66
REAL(kd), PARAMETER :: SVPT0=273.15
REAL(kd), PARAMETER :: r_d = 287.
REAL(kd), PARAMETER :: r_v = 461.6
REAL(kd), PARAMETER :: ep_2=r_d/r_v
!-------------------------------------------------
REAL(kd) :: t, p

!=================================================
CALL set_area_w
DO k = kmin, kmax
	DO i = imin, imax
		t = wGrid%theta(i,k)*wGrid%pi(i,k)
		p = EXP(Cp/Rd*LOG(wGrid%pi(i,k)))*p0
		tc = t - SVPT0
		esw = 1000.0 * SVP1 * EXP( SVP2  * tc / ( t - SVP3  ) )
		esi = 1000.0 * SVP1 * EXP( SVPI2 * tc / ( t - SVPI3 ) )
		QVSW = EP_2 * esw / ( p - esw )
		QVSI = EP_2 * esi / ( p - esi )
		!QCLD = wGrid%qc(i,k) + wGrid%qi(i,k) + wGrid%qs(i,k)
		QCLD = wGrid%qc(i,k) + wGrid%qi(i,k) + wGrid%qs(i,k) + wGrid%qg(i,k)
		IF (QCLD < QCLDMIN) THEN
			weight = 0.
			wGrid%cldfra(i,k) = 0.
		ELSE
			!weight = (wGrid%qi(i,k) + wGrid%qs(i,k))/QCLD
			weight = (wGrid%qi(i,k) + wGrid%qs(i,k) + wGrid%qg(i,k))/QCLD
			qvs_weight = (1. - weight)*QVSW + weight*QVSI
			RHUM = wGrid%qv(i,k)/qvs_weight
			IF (RHUM >= RHGRID) THEN
				wGrid%cldfra(i,k) = 1.
			ELSE
				SUBSAT = MAX(1.e-10, RHGRID*qvs_weight - wGrid%qv(i,k))
				DENOM = (SUBSAT)**GAMMA
				ARG = MAX(-6.9, -ALPHA0*QCLD/DENOM)
				RHUM = MAX(1.e-10, RHUM)
				wGrid%cldfra(i,k) = (RHUM/RHGRID)**PEXP*(1.-EXP(ARG))
				IF (wGrid%cldfra(i,k) < .01) wGrid%cldfra = 0.
			END IF
		END IF
	END DO
END DO
!=================================================
END SUBROUTINE calc_cldfra2
!=================================================

!=================================================
END MODULE sp_module_cldfra
!=================================================

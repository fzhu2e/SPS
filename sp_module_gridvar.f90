!=================================================
MODULE sp_module_gridvar
USE sp_module_constant
USE sp_module_model
USE sp_module_debug
IMPLICIT NONE
!=================================================
TYPE :: mainvar
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: u = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: w = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: pi_1 = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: qc = undef, qv = undef, qr = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: qi = undef, qs = undef, qg = undef
END TYPE mainvar
!=================================================

!=================================================
TYPE, EXTENDS(mainvar) :: grid
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: pi = undef, pi_0 = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_0 = undef, theta_1 = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_v = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: theta_M = undef, theta_M_0 = undef, theta_M_1 = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: rho_0 = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: qt = undef  ! qt = qv + qc + qr + qi + qs + qg

	REAL(kd), DIMENSION(ims:ime,kms:kme) :: Du = undef, Dw = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: Dtheta = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: Dqv = undef, Dqc = undef, Dqr = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: Dqi = undef, Dqs = undef, Dqg = undef

	REAL(kd), DIMENSION(ims:ime,kms:kme) :: Mtheta = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: Mqv = undef, Mqc = undef, Mqr = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: Mqi = undef, Mqs = undef, Mqg = undef
	
	REAL(kd), DIMENSION(ims:ime) :: xx = undef
	REAL(kd), DIMENSION(kms:kme) :: zeta = undef               ! model height
	REAL(kd), DIMENSION(ims:ime) :: zs = undef                 ! terrain height
	REAL(kd), DIMENSION(ims:ime) :: PzsPx = undef
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: zz = undef         ! real height

	REAL(kd), DIMENSION(kms:kme) :: b = undef                  ! decay coeff of terrain
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: G = undef
	REAL(kd), DIMENSION(ims:ime) :: H = undef

!=================================================
! Sponge Layer
!-------------------------------------------------
	REAL(kd), DIMENSION(ims:ime,kms:kme) :: tau = undef

END TYPE grid
!=================================================

!=================================================
CONTAINS
!=================================================

!=================================================
SUBROUTINE u2pi(var_u,var_pi)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_pi
INTEGER :: i, k
CALL set_area_pi
CALL set_area_expand(expand)
!IF (ANY(var_u(imin-1:imax,kmin:kmax) == undef)) STOP "u2pi WRONG!!!"
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_pi(i,k) = (var_u(i-1,k) + var_u(i,k))/2.
	END DO
END DO
!OMP END PARALLEL DO
END SUBROUTINE u2pi
!=================================================

!=================================================
SUBROUTINE u2vir(var_u,var_vir)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_vir
INTEGER :: i, k
CALL set_area_vir
CALL set_area_expand(expand)
!IF (ANY(var_u(imin:imax,kmin-1:kmax) == undef)) STOP "u2vir WRONG!!!"
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_vir(i,k) = (var_u(i,k) + var_u(i,k-1))/2.
	END DO
END DO
!OMP END PARALLEL DO
END SUBROUTINE u2vir
!=================================================

!=================================================
SUBROUTINE u2w(var_u,var_w)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_u 
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_w 
INTEGER :: i, k
CALL set_area_w
CALL set_area_expand(expand)
!IF (ANY(var_u(imin-1:imax,kmin-1:kmax) == undef)) STOP "u2w WRONG!!!"
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_w(i,k) = (var_u(i,k) + var_u(i,k-1) + var_u(i-1,k) + var_u(i-1,k-1))/4.
	END DO
END DO
!OMP END PARALLEL DO
END SUBROUTINE u2w
!=================================================

!=================================================
SUBROUTINE pi2u(var_pi,var_u)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_pi
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_u
INTEGER :: i, k
CALL set_area_u
CALL set_area_expand(expand)
!IF (ANY(var_pi(imin:imax+1,kmin:kmax) == undef)) STOP "pi2u WRONG!!!"
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_u(i,k) = (var_pi(i,k) + var_pi(i+1,k))/2.
	END DO
END DO
!OMP END PARALLEL DO
END SUBROUTINE pi2u
!=================================================

!=================================================
SUBROUTINE pi2w(var_pi,var_w)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_pi
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_w
INTEGER :: i, k
CALL set_area_w
CALL set_area_expand(expand)
!IF (ANY(var_pi(imin:imax,kmin-1:kmax) == undef)) STOP "pi2w WRONG!!!"
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_w(i,k) = (var_pi(i,k-1) + var_pi(i,k))/2.
	END DO
END DO
!OMP END PARALLEL DO
END SUBROUTINE pi2w
!=================================================

!=================================================
SUBROUTINE pi2vir(var_pi,var_vir)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_pi
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_vir
INTEGER :: i, k
CALL set_area_vir
CALL set_area_expand(expand)
!IF (ANY(var_pi(imin:imax+1,kmin-1:kmax) == undef)) STOP "pi2vir WRONG!!!"
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_vir(i,k) = (var_pi(i,k-1) + var_pi(i,k) + var_pi(i+1,k) + var_pi(i+1,k-1))/4.
	END DO
END DO
!OMP END PARALLEL DO
END SUBROUTINE pi2vir
!=================================================

!=================================================
SUBROUTINE w2vir(var_w,var_vir)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_vir
INTEGER :: i, k
CALL set_area_vir
CALL set_area_expand(expand)
!IF (ANY(var_w(imin:imax+1,kmin:kmax) == undef)) STOP "w2vir WRONG!!!"
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_vir(i,k) = (var_w(i,k) + var_w(i+1,k))/2.
	END DO
END DO
!OMP END PARALLEL DO
END SUBROUTINE w2vir
!=================================================

!=================================================
SUBROUTINE w2pi(var_w,var_pi)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_pi
INTEGER :: i, k
CALL set_area_pi
CALL set_area_expand(expand)
!IF (ANY(var_w(imin:imax,kmin:kmax+1) == undef)) STOP "w2pi WRONG!!!"
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_pi(i,k) = (var_w(i,k) + var_w(i,k+1))/2.
	END DO
END DO
!OMP END PARALLEL DO
END SUBROUTINE w2pi
!=================================================

!=================================================
SUBROUTINE w2u(var_w,var_u)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: var_u
INTEGER :: i, k
CALL set_area_u
CALL set_area_expand(expand)
!IF (ANY(var_w(imin:imax+1,kmin:kmax+1) == undef)) STOP "w2u WRONG!!!"
!OMP PARALLEL DO
DO k = kmin, kmax
	DO i = imin, imax
		var_u(i,k) = (var_w(i,k) + var_w(i,k+1) + var_w(i+1,k) + var_w(i+1,k+1))/4.
	END DO
END DO
!OMP END PARALLEL DO
END SUBROUTINE w2u
!=================================================

!=================================================
SUBROUTINE ppx_u(var_pi,output)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_pi
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: output
INTEGER :: i, k
CALL set_area_u
!IF (ANY(var_pi(imin:imax+1,kmin:kmax) == undef)) STOP "ppx_u WRONG!!!"
DO k = kmin, kmax
	DO i = imin, imax
		output(i,k) = (var_pi(i+1,k) - var_pi(i,k))/dx
	END DO
END DO
END SUBROUTINE ppx_u
!=================================================

!=================================================
SUBROUTINE ppzeta_u(var_vir,output)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_vir
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: output
INTEGER :: i, k
CALL set_area_u
!IF (ANY(var_vir(imin:imax,kmin:kmax+1) == undef)) STOP "ppzeta_u WRONG!!!"
DO k = kmin, kmax
	DO i = imin, imax
		output(i,k) = (var_vir(i,k+1) - var_vir(i,k))/dz
	END DO
END DO
END SUBROUTINE ppzeta_u
!=================================================

!=================================================
SUBROUTINE ppx_w(var_vir,output)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_vir
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: output
INTEGER :: i, k
CALL set_area_w
!IF (ANY(var_vir(imin-1:imax,kmin:kmax) == undef)) STOP "ppx_w WRONG!!!"
DO k = kmin, kmax
	DO i = imin, imax
		output(i,k) = (var_vir(i,k) - var_vir(i-1,k))/dx
	END DO
END DO
END SUBROUTINE ppx_w
!=================================================

!=================================================
SUBROUTINE ppzeta_w(var_pi,output)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_pi
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: output
INTEGER :: i, k
CALL set_area_w
!IF (ANY(var_pi(imin:imax,kmin-1:kmax) == undef)) STOP "ppzeta_w WRONG!!!"
DO k = kmin, kmax
	DO i = imin, imax
		output(i,k) = (var_pi(i,k) - var_pi(i,k-1))/dz
	END DO
END DO
END SUBROUTINE ppzeta_w
!=================================================

!=================================================
SUBROUTINE ppx_pi(var_u,output)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: output
INTEGER :: i, k
CALL set_area_pi
!IF (ANY(var_u(imin-1:imax,kmin:kmax) == undef)) STOP "ppx_pi WRONG!!!"
DO k = kmin, kmax
	DO i = imin, imax
		output(i,k) = (var_u(i,k) - var_u(i-1,k))/dx
	END DO
END DO
END SUBROUTINE ppx_pi
!=================================================

!=================================================
SUBROUTINE ppzeta_pi(var_w,output)
IMPLICIT NONE
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: output
INTEGER :: i, k
CALL set_area_pi
!IF (ANY(var_w(imin:imax,kmin:kmax+1) == undef)) STOP "ppzeta_pi WRONG!!!"
DO k = kmin, kmax
	DO i = imin, imax
		output(i,k) = (var_w(i,k+1) - var_w(i,k))/dz
	END DO
END DO
END SUBROUTINE ppzeta_pi
!=================================================

!=================================================
END MODULE sp_module_gridvar
!=================================================

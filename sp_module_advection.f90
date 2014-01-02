!=================================================
MODULE sp_module_advection
USE sp_module_constant
USE sp_module_model
USE sp_module_gridvar
USE sp_module_debug
IMPLICIT NONE
!=================================================
CONTAINS

!=================================================
SUBROUTINE calc_advection_u(var_u,A_u,uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE(grid) uGrid, wGrid, piGrid, virGrid
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_u
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: A_u

REAL(kd), DIMENSION(ims:ime,kms:kme) :: rhou_pi, rhouvar_pi
REAL(kd), DIMENSION(ims:ime,kms:kme) :: rhou_vir, rhouvar_vir
REAL(kd), DIMENSION(ims:ime,kms:kme) :: rhow_vir, rhowvar_vir

REAL(kd) :: fa, fb, fc, fd, fe, ff

REAL(kd), DIMENSION(ims:ime,kms:kme) :: PrhouPx_u, PrhouvarPx_u
REAL(kd), DIMENSION(ims:ime,kms:kme) :: PrhouPzeta_u, PrhouvarPzeta_u
REAL(kd), DIMENSION(ims:ime,kms:kme) :: PrhowPzeta_u, PrhowvarPzeta_u

REAL(kd) :: temp_a, temp_b, temp_c
INTEGER :: i, k
!=================================================
CALL set_area_pi
CALL set_area_expand(expand)
SELECT CASE (AdvectionScheme)
CASE(2)
	!$OMP PARALLEL DO PRIVATE(fa)
	DO k = kmin, kmax
		DO i = imin, imax
			rhou_pi(i,k) = piGrid%rho_0(i,k)*piGrid%u(i,k)
			fa = var_u(i-1,k) + var_u(i,k)
			rhouvar_pi(i,k) = rhou_pi(i,k)*fa/2.
		END DO
	END DO
	!$OMP END PARALLEL DO
CASE(5)
	!$OMP PARALLEL DO PRIVATE(fa,fb,fc,fd,fe,ff)
	DO k = kmin, kmax
		DO i = imin, imax
			rhou_pi(i,k) = piGrid%rho_0(i,k)*piGrid%u(i,k)
			fa = var_u(i-1,k) + var_u(i,k)
			fb = var_u(i-2,k) + var_u(i+1,k)
			fc = var_u(i-3,k) + var_u(i+2,k)
			rhouvar_pi(i,k) = rhou_pi(i,k)/60.*(37*fa - 8*fb + fc)
			fd = - var_u(i-1,k) + var_u(i,k)
			fe = - var_u(i-2,k) + var_u(i+1,k)
			ff = - var_u(i-3,k) + var_u(i+2,k)
			rhouvar_pi(i,k) = rhouvar_pi(i,k) - ABS(piGrid%u(i,k))/60.*(10*fd - 5*fe + ff)
		END DO
	END DO
	!$OMP END PARALLEL DO
CASE DEFAULT
	STOP "Wrong AdvectionScheme!!!"
END SELECT

CALL set_area_vir
CALL set_area_expand(expand)
SELECT CASE (AdvectionScheme)
CASE(2)
	!$OMP PARALLEL DO PRIVATE(fa)
	DO k = kmin, kmax
		DO i = imin, imax
			rhou_vir(i,k) = virGrid%rho_0(i,k)*virGrid%u(i,k)
			rhow_vir(i,k) = virGrid%rho_0(i,k)*virGrid%w(i,k)
			fa = var_u(i,k-1) + var_u(i,k)
			rhouvar_vir(i,k) = rhou_vir(i,k)*fa/2.
			rhowvar_vir(i,k) = rhow_vir(i,k)*fa/2.
		END DO
	END DO
	!$OMP END PARALLEL DO
CASE(5)
	!$OMP PARALLEL DO PRIVATE(fa,fb,fc,fd,fe,ff)
	DO k = kmin, kmax
		DO i = imin, imax
			rhou_vir(i,k) = virGrid%rho_0(i,k)*virGrid%u(i,k)
			rhow_vir(i,k) = virGrid%rho_0(i,k)*virGrid%w(i,k)
			fa = var_u(i,k-1) + var_u(i,k)
			fb = var_u(i,k-2) + var_u(i,k+1)
			fc = var_u(i,k-3) + var_u(i,k+2)
			rhouvar_vir(i,k) = rhou_vir(i,k)/60.*(37*fa - 8*fb + fc)
			rhowvar_vir(i,k) = rhow_vir(i,k)/60.*(37*fa - 8*fb + fc)
			fd = - var_u(i,k-1) + var_u(i,k)
			fe = - var_u(i,k-2) + var_u(i,k+1)
			ff = - var_u(i,k-3) + var_u(i,k+2)
			rhouvar_vir(i,k) = rhouvar_vir(i,k) - ABS(virGrid%u(i,k))/60.*(10*fd - 5*fe + ff)
			rhowvar_vir(i,k) = rhowvar_vir(i,k) - ABS(virGrid%w(i,k))/60.*(10*fd - 5*fe + ff)
			IF (k <= kmin+2 .OR. k >= kmax-2) THEN
				rhouvar_vir(i,k) = rhou_vir(i,k)*fa/2.
				rhowvar_vir(i,k) = rhow_vir(i,k)*fa/2.
			END IF
		END DO
	END DO
	!$OMP END PARALLEL DO
CASE DEFAULT
	STOP "Wrong AdvectionScheme!!!"
END SELECT

CALL ppx_u(rhou_pi,PrhouPx_u)
CALL ppx_u(rhouvar_pi,PrhouvarPx_u)

CALL ppzeta_u(rhou_vir,PrhouPzeta_u)
CALL ppzeta_u(rhouvar_vir,PrhouvarPzeta_u)

CALL ppzeta_u(rhow_vir,PrhowPzeta_u)
CALL ppzeta_u(rhowvar_vir,PrhowvarPzeta_u)

CALL set_area_u
!$OMP PARALLEL DO PRIVATE(temp_a,temp_b,temp_c)
DO k = kmin, kmax
	DO i = imin, imax
		temp_a = - 1./uGrid%rho_0(i,k)*(PrhouvarPx_u(i,k) - var_u(i,k)*PrhouPx_u(i,k))
		temp_b = - uGrid%G(i,k)/uGrid%rho_0(i,k)*(PrhouvarPzeta_u(i,k) - var_u(i,k)*PrhouPzeta_u(i,k))
		temp_c = - uGrid%H(i)/uGrid%rho_0(i,k)*(PrhowvarPzeta_u(i,k) - var_u(i,k)*PrhowPzeta_u(i,k))
		A_u(i,k) = temp_a + temp_b + temp_c
	END DO
END DO
!$OMP END PARALLEL DO
!=================================================
END SUBROUTINE calc_advection_u
!=================================================

!=================================================
SUBROUTINE calc_advection_w(var_w,A_w,uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE(grid) uGrid, wGrid, piGrid, virGrid
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_w
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: A_w

REAL(kd), DIMENSION(ims:ime,kms:kme) :: rhou_vir, rhouvar_vir
REAL(kd), DIMENSION(ims:ime,kms:kme) :: rhou_pi, rhouvar_pi
REAL(kd), DIMENSION(ims:ime,kms:kme) :: rhow_pi, rhowvar_pi

REAL(kd) :: fa, fb, fc, fd, fe, ff

REAL(kd), DIMENSION(ims:ime,kms:kme) :: PrhouPx_w, PrhouvarPx_w
REAL(kd), DIMENSION(ims:ime,kms:kme) :: PrhouPzeta_w, PrhouvarPzeta_w
REAL(kd), DIMENSION(ims:ime,kms:kme) :: PrhowPzeta_w, PrhowvarPzeta_w

REAL(kd) :: temp_a, temp_b, temp_c
INTEGER :: i, k
!=================================================
CALL set_area_vir
CALL set_area_expand(expand)
SELECT CASE (AdvectionScheme)
CASE(2)
	!$OMP PARALLEL DO PRIVATE(fa)
	DO k = kmin, kmax
		DO i = imin, imax
			rhou_vir(i,k) = virGrid%rho_0(i,k)*virGrid%u(i,k)
			fa = var_w(i,k) + var_w(i+1,k)
			rhouvar_vir(i,k) = rhou_vir(i,k)*fa/2.
		END DO
	END DO
	!$OMP END PARALLEL DO
CASE(5)
	!$OMP PARALLEL DO PRIVATE(fa,fb,fc,fd,fe,ff)
	DO k = kmin, kmax
		DO i = imin, imax
			rhou_vir(i,k) = virGrid%rho_0(i,k)*virGrid%u(i,k)
			fa = var_w(i,k) + var_w(i+1,k)
			fb = var_w(i-1,k) + var_w(i+2,k)
			fc = var_w(i-2,k) + var_w(i+3,k)
			rhouvar_vir(i,k) = rhou_vir(i,k)/60.*(37*fa - 8*fb + fc)
			fd = - var_w(i,k) + var_w(i+1,k)
			fe = - var_w(i-1,k) + var_w(i+2,k)
			ff = - var_w(i-2,k) + var_w(i+3,k)
			rhouvar_vir(i,k) = rhouvar_vir(i,k) - ABS(virGrid%u(i,k))/60.*(10*fd - 5*fe + ff)
		END DO
	END DO
	!$OMP END PARALLEL DO
CASE DEFAULT
	STOP "Wrong AdvectionScheme!!!"
END SELECT

CALL set_area_pi
CALL set_area_expand(expand)
SELECT CASE (AdvectionScheme)
CASE(2)
	!$OMP PARALLEL DO PRIVATE(fa)
	DO k = kmin, kmax
		DO i = imin, imax
			rhou_pi(i,k) = piGrid%rho_0(i,k)*piGrid%u(i,k)
			rhow_pi(i,k) = piGrid%rho_0(i,k)*piGrid%w(i,k)
			fa = var_w(i,k) + var_w(i,k+1)
			rhouvar_pi(i,k) = rhou_pi(i,k)*fa/2.
			rhowvar_pi(i,k) = rhow_pi(i,k)*fa/2.
		END DO
	END DO
	!$OMP END PARALLEL DO
CASE(5)
	!$OMP PARALLEL DO PRIVATE(fa,fb,fc,fd,fe,ff)
	DO k = kmin, kmax
		DO i = imin, imax
			rhou_pi(i,k) = piGrid%rho_0(i,k)*piGrid%u(i,k)
			rhow_pi(i,k) = piGrid%rho_0(i,k)*piGrid%w(i,k)
			fa = var_w(i,k) + var_w(i,k+1)
			fb = var_w(i,k-1) + var_w(i,k+2)
			fc = var_w(i,k-2) + var_w(i,k+3)
			rhouvar_pi(i,k) = rhou_pi(i,k)/60.*(37*fa - 8*fb + fc)
			rhowvar_pi(i,k) = rhow_pi(i,k)/60.*(37*fa - 8*fb + fc)
			fd = - var_w(i,k) + var_w(i,k+1)
			fe = - var_w(i,k-1) + var_w(i,k+2)
			ff = - var_w(i,k-2) + var_w(i,k+3)
			rhouvar_pi(i,k) = rhouvar_pi(i,k) - ABS(piGrid%u(i,k))/60.*(10*fd - 5*fe + ff)
			rhowvar_pi(i,k) = rhowvar_pi(i,k) - ABS(piGrid%w(i,k))/60.*(10*fd - 5*fe + ff)
			IF (k <= kmin+1 .OR. k >= kmax-1) THEN
				rhouvar_pi(i,k) = rhou_pi(i,k)*fa/2.
				rhowvar_pi(i,k) = rhow_pi(i,k)*fa/2.
			END IF
		END DO
	END DO
	!$OMP END PARALLEL DO
CASE DEFAULT
	STOP "Wrong AdvectionScheme!!!"
END SELECT

CALL ppx_w(rhou_vir,PrhouPx_w)
CALL ppx_w(rhouvar_vir,PrhouvarPx_w)

CALL ppzeta_w(rhou_pi,PrhouPzeta_w)
CALL ppzeta_w(rhouvar_pi,PrhouvarPzeta_w)

CALL ppzeta_w(rhow_pi,PrhowPzeta_w)
CALL ppzeta_w(rhowvar_pi,PrhowvarPzeta_w)

CALL set_area_w
!$OMP PARALLEL DO PRIVATE(temp_a,temp_b,temp_c)
DO k = kmin, kmax
	DO i = imin, imax
		temp_a = - 1./wGrid%rho_0(i,k)*(PrhouvarPx_w(i,k) - var_w(i,k)*PrhouPx_w(i,k))
		temp_b = - wGrid%G(i,k)/wGrid%rho_0(i,k)*(PrhouvarPzeta_w(i,k) - var_w(i,k)*PrhouPzeta_w(i,k))
		temp_c = - wGrid%H(i)/wGrid%rho_0(i,k)*(PrhowvarPzeta_w(i,k) - var_w(i,k)*PrhowPzeta_w(i,k))
		A_w(i,k) = temp_a + temp_b + temp_c
	END DO
END DO
!$OMP END PARALLEL DO
!=================================================
END SUBROUTINE calc_advection_w
!=================================================

!=================================================
SUBROUTINE calc_advection_pi(var_pi,A_pi,uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE(grid) uGrid, wGrid, piGrid, virGrid
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var_pi
REAL(kd), DIMENSION(ims:ime,kms:kme), INTENT(OUT) :: A_pi

REAL(kd), DIMENSION(ims:ime,kms:kme) :: rhou_u, rhouvar_u
REAL(kd), DIMENSION(ims:ime,kms:kme) :: rhou_w, rhouvar_w
REAL(kd), DIMENSION(ims:ime,kms:kme) :: rhow_w, rhowvar_w

REAL(kd) :: fa, fb, fc, fd, fe, ff

REAL(kd), DIMENSION(ims:ime,kms:kme) :: PrhouPx_pi, PrhouvarPx_pi
REAL(kd), DIMENSION(ims:ime,kms:kme) :: PrhouPzeta_pi, PrhouvarPzeta_pi
REAL(kd), DIMENSION(ims:ime,kms:kme) :: PrhowPzeta_pi, PrhowvarPzeta_pi

REAL(kd) :: temp_a, temp_b, temp_c
INTEGER :: i, k
!=================================================
CALL set_area_u
CALL set_area_expand(expand)
SELECT CASE (AdvectionScheme)
CASE(2)
	!$OMP PARALLEL DO PRIVATE(fa)
	DO k = kmin, kmax
		DO i = imin, imax
			rhou_u(i,k) = uGrid%rho_0(i,k)*uGrid%u(i,k)
			fa = var_pi(i,k) + var_pi(i+1,k)
			rhouvar_u(i,k) = rhou_u(i,k)*fa/2.
		END DO
	END DO
	!$OMP END PARALLEL DO
CASE(5)
	!$OMP PARALLEL DO PRIVATE(fa,fb,fc,fd,fe,ff)
	DO k = kmin, kmax
		DO i = imin, imax
			rhou_u(i,k) = uGrid%rho_0(i,k)*uGrid%u(i,k)
			fa = var_pi(i,k) + var_pi(i+1,k)
			fb = var_pi(i-1,k) + var_pi(i+2,k)
			fc = var_pi(i-2,k) + var_pi(i+3,k)
			rhouvar_u(i,k) = rhou_u(i,k)/60.*(37*fa - 8*fb + fc)
			fd = - var_pi(i,k) + var_pi(i+1,k)
			fe = - var_pi(i-1,k) + var_pi(i+2,k)
			ff = - var_pi(i-2,k) + var_pi(i+3,k)
			rhouvar_u(i,k) = rhouvar_u(i,k) - ABS(uGrid%u(i,k))/60.*(10*fd - 5*fe + ff)
		END DO
	END DO
	!$OMP END PARALLEL DO
CASE DEFAULT
	STOP "Wrong AdvectionScheme!!!"
END SELECT

CALL set_area_w
CALL set_area_expand(expand)
SELECT CASE (AdvectionScheme)
CASE(2)
	!$OMP PARALLEL DO PRIVATE(fa)
	DO k = kmin, kmax
		DO i = imin, imax
			rhou_w(i,k) = wGrid%rho_0(i,k)*wGrid%u(i,k)
			rhow_w(i,k) = wGrid%rho_0(i,k)*wGrid%w(i,k)
			fa = var_pi(i,k-1) + var_pi(i,k)
			rhouvar_w(i,k) = rhou_w(i,k)*fa/2.
			rhowvar_w(i,k) = rhow_w(i,k)*fa/2.
		END DO
	END DO
	!$OMP END PARALLEL DO
CASE(5)
	!$OMP PARALLEL DO PRIVATE(fa,fb,fc,fd,fe,ff)
	DO k = kmin, kmax
		DO i = imin, imax
			rhou_w(i,k) = wGrid%rho_0(i,k)*wGrid%u(i,k)
			rhow_w(i,k) = wGrid%rho_0(i,k)*wGrid%w(i,k)
			fa = var_pi(i,k-1) + var_pi(i,k)
			fb = var_pi(i,k-2) + var_pi(i,k+1)
			fc = var_pi(i,k-3) + var_pi(i,k+2)
			rhouvar_w(i,k) = rhou_w(i,k)/60.*(37*fa - 8*fb + fc)
			rhowvar_w(i,k) = rhow_w(i,k)/60.*(37*fa - 8*fb + fc)
			fd = - var_pi(i,k-1) + var_pi(i,k)
			fe = - var_pi(i,k-2) + var_pi(i,k+1)
			ff = - var_pi(i,k-3) + var_pi(i,k+2)
			rhouvar_w(i,k) = rhouvar_w(i,k) - ABS(wGrid%u(i,k))/60.*(10*fd - 5*fe + ff)
			rhowvar_w(i,k) = rhowvar_w(i,k) - ABS(wGrid%w(i,k))/60.*(10*fd - 5*fe + ff)
			IF (k <= kmin+2 .OR. k >= kmax-2) THEN
				rhouvar_w(i,k) = rhou_w(i,k)*fa/2.
				rhowvar_w(i,k) = rhow_w(i,k)*fa/2.
			END IF
		END DO
	END DO
	!$OMP END PARALLEL DO
CASE DEFAULT
	STOP "Wrong AdvectionScheme!!!"
END SELECT

CALL ppx_pi(rhou_u,PrhouPx_pi)
CALL ppx_pi(rhouvar_u,PrhouvarPx_pi)

CALL ppzeta_pi(rhou_w,PrhouPzeta_pi)
CALL ppzeta_pi(rhouvar_w,PrhouvarPzeta_pi)

CALL ppzeta_pi(rhow_w,PrhowPzeta_pi)
CALL ppzeta_pi(rhowvar_w,PrhowvarPzeta_pi)

CALL set_area_pi
!$OMP PARALLEL DO PRIVATE(temp_a,temp_b,temp_c)
DO k = kmin, kmax
	DO i = imin, imax
		temp_a = - 1./piGrid%rho_0(i,k)*(PrhouvarPx_pi(i,k) - var_pi(i,k)*PrhouPx_pi(i,k))
		temp_b = - piGrid%G(i,k)/piGrid%rho_0(i,k)*(PrhouvarPzeta_pi(i,k) - var_pi(i,k)*PrhouPzeta_pi(i,k))
		temp_c = - piGrid%H(i)/piGrid%rho_0(i,k)*(PrhowvarPzeta_pi(i,k) - var_pi(i,k)*PrhowPzeta_pi(i,k))
		A_pi(i,k) = temp_a + temp_b + temp_c
	END DO
END DO
!$OMP END PARALLEL DO
!=================================================
END SUBROUTINE calc_advection_pi
!=================================================

!=================================================
END MODULE sp_module_advection
!=================================================

!=================================================
! Super-Parametertization System (SPS)
!-------------------------------------------------
! Version: 0.2
! Author: Feng Zhu
! Email: zhuf.atmos@gmail.com
! Date: 2014-06-12 18:18:45
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!=================================================
MODULE sp_module_interpolate
USE sp_module_constant
USE sp_module_model
USE sp_module_gridvar
USE sp_module_debug
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================
!/////////////////////////////////////////////////////////////////////
!=================================================
! Basic interpolate.
!=================================================
SUBROUTINE basic_interpolate(Main,uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
!=================================================
TYPE(mainvar), INTENT(IN) :: Main
TYPE(grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!-------------------------------------------------
INTEGER :: i, k
!=================================================
CALL u2w(Main%u,wGrid%u)
CALL u2pi(Main%u,piGrid%u)
CALL u2vir(Main%u,virGrid%u)

CALL w2pi(Main%w,piGrid%w)
CALL w2vir(Main%w,virGrid%w)

CALL pi2vir(Main%pi_1,virGrid%pi_1)

CALL w2u(wGrid%theta_M,uGrid%theta_M)
CALL w2pi(wGrid%theta_M,piGrid%theta_M)

!CALL w2pi(Main%theta,piGrid%theta)

CALL w2pi(Main%qv,piGrid%qv)
CALL w2pi(wGrid%Mqv,piGrid%Mqv)
CALL w2pi(wGrid%Dqv,piGrid%Dqv)

CALL w2pi(wGrid%theta,piGrid%theta)
CALL w2pi(wGrid%Mtheta,piGrid%Mtheta)
CALL w2pi(wGrid%Dtheta,piGrid%Dtheta)

!CALL w2u(Main%qv,uGrid%qv)
!CALL w2u(Main%qc,uGrid%qc)
!CALL w2u(Main%qr,uGrid%qr)

!CALL w2vir(Main%qv,virGrid%qv)
!CALL w2vir(Main%qc,virGrid%qc)
!CALL w2vir(Main%qr,virGrid%qr)
!=================================================
END SUBROUTINE basic_interpolate
!=================================================

!=================================================
! Calculate virtual theta
!=================================================
SUBROUTINE calc_virTheta(uGrid,wGrid,piGrid,virGrid)
IMPLICIT NONE
TYPE (grid), INTENT(INOUT) :: uGrid, wGrid, piGrid, virGrid
!=================================================
! theta_v, theta_M_0, theta_M_1
!-------------------------------------------------
wGrid%theta_v = wGrid%theta*(1. + 0.61*wGrid%qv)
CALL w2pi(wGrid%theta_v,piGrid%theta_v)
wGrid%theta_M = wGrid%theta_v*(1. - wGrid%qc - wGrid%qi)
CALL w2pi(wGrid%theta_M,piGrid%theta_M)
wGrid%theta_1 = wGrid%theta - wGrid%theta_0 ! for output
wGrid%theta_M_1 = wGrid%theta_M - wGrid%theta_M_0

piGrid%pi = piGrid%pi_1 + piGrid%pi_0
CALL pi2w(piGrid%pi,wGrid%pi)

!CALL debug_ascii_output(uGrid%u,"u")
!CALL debug_ascii_output(wGrid%w,"w")
!CALL debug_ascii_output(piGrid%pi_1,"pi_1")
!CALL debug_ascii_output(piGrid%pi_0,"pi_0_pi")
!CALL debug_ascii_output(uGrid%pi_0,"pi_0_u")
!CALL debug_ascii_output(wGrid%pi_0,"pi_0_w")
!CALL debug_ascii_output(virGrid%pi_0,"pi_0_vir")
!CALL debug_ascii_output(piGrid%pi,"pi")
!CALL debug_ascii_output(wGrid%theta,"theta")
!CALL debug_ascii_output(wGrid%theta_0,"theta_0_w")
!CALL debug_ascii_output(uGrid%theta_0,"theta_0_u")
!CALL debug_ascii_output(piGrid%theta_0,"theta_0_pi")
!CALL debug_ascii_output(virGrid%theta_0,"theta_0_vir")
!CALL debug_ascii_output(wGrid%theta_1,"theta_1")
!CALL debug_ascii_output(wGrid%qv,"qv")
!CALL debug_ascii_output(wGrid%qc,"qc")
!CALL debug_ascii_output(wGrid%qr,"qr")
!CALL debug_ascii_output(wGrid%qi,"qi")
!CALL debug_ascii_output(wGrid%qs,"qs")
!CALL debug_ascii_output(wGrid%qg,"qg")
!CALL debug_ascii_output(uGrid%rho_0,"rho_0_u")
!CALL debug_ascii_output(wGrid%rho_0,"rho_0_w")
!CALL debug_ascii_output(piGrid%rho_0,"rho_0_pi")
!CALL debug_ascii_output(virGrid%rho_0,"rho_0_vir")

!CALL debug_SFSG
!=================================================
END SUBROUTINE calc_virTheta
!=================================================

!/////////////////////////////////////////////////////////////////////
!=================================================
END MODULE sp_module_interpolate
!=================================================

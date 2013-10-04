!=================================================
! The debug module of SPS-dynamic-integrate.
!-------------------------------------------------
! Version: 0.11
! Author: Zhu F.
! Email: lyricorpse@gmail.com
! Date: 2013-05-16 19:08:57 
! Copyright: This software is provided under a CC BY-NC-SA 3.0 License(http://creativecommons.org/licenses/by-nc-sa/3.0/deed.zh)
!=================================================
MODULE sp_module_debug
USE sp_module_constant
USE sp_module_model
IMPLICIT NONE
!=================================================
CONTAINS
!=================================================

!=================================================
! Output ASCII file.
!=================================================
SUBROUTINE debug_ascii_output(var)
IMPLICIT NONE
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var
INTEGER :: i, k
!=================================================
OPEN(1, FILE = 'debug/debug_field.txt')
DO k = kms, kme
	DO i = ims, ime - 1
		WRITE(1, "(F11.5\)") var(i,k)
	END DO
	WRITE(1, "(F11.5)") var(ime,k)
END DO
CLOSE(1)
!=================================================
END SUBROUTINE debug_ascii_output
!=================================================

!=================================================
! Initiate value.
!=================================================
SUBROUTINE debug_undef_all( var00,var01,var02,var03,var04,var05,var06,var07,var08,var09,  &
                            var10,var11,var12,var13,var14,var15,var16,var17,var18,var19,  &
                            var20,var21,var22,var23,var24,var25,var26,var27,var28,var29,  &
                            var30,var31,var32,var33,var34,var35,var36,var37,var38,var39,  &
                            var40,var41,var42,var43,var44,var45,var46,var47,var48,var49   )
IMPLICIT NONE
!=================================================
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(INOUT), OPTIONAL :: &
      var00,var01,var02,var03,var04,var05,var06,var07,var08,var09,  &
      var10,var11,var12,var13,var14,var15,var16,var17,var18,var19,  &
      var20,var21,var22,var23,var24,var25,var26,var27,var28,var29,  &
      var30,var31,var32,var33,var34,var35,var36,var37,var38,var39,  &
      var40,var41,var42,var43,var44,var45,var46,var47,var48,var49   
!=================================================
IF (PRESENT(var00)) var00 = undef
IF (PRESENT(var01)) var01 = undef
IF (PRESENT(var02)) var02 = undef
IF (PRESENT(var03)) var03 = undef
IF (PRESENT(var04)) var04 = undef
IF (PRESENT(var05)) var05 = undef
IF (PRESENT(var06)) var06 = undef
IF (PRESENT(var07)) var07 = undef
IF (PRESENT(var08)) var08 = undef
IF (PRESENT(var09)) var09 = undef
IF (PRESENT(var10)) var10 = undef
IF (PRESENT(var11)) var11 = undef
IF (PRESENT(var12)) var12 = undef
IF (PRESENT(var13)) var13 = undef
IF (PRESENT(var14)) var14 = undef
IF (PRESENT(var15)) var15 = undef
IF (PRESENT(var16)) var16 = undef
IF (PRESENT(var17)) var17 = undef
IF (PRESENT(var18)) var18 = undef
IF (PRESENT(var19)) var19 = undef
IF (PRESENT(var20)) var20 = undef
IF (PRESENT(var21)) var21 = undef
IF (PRESENT(var22)) var22 = undef
IF (PRESENT(var23)) var23 = undef
IF (PRESENT(var24)) var24 = undef
IF (PRESENT(var25)) var25 = undef
IF (PRESENT(var26)) var26 = undef
IF (PRESENT(var27)) var27 = undef
IF (PRESENT(var28)) var28 = undef
IF (PRESENT(var29)) var29 = undef
IF (PRESENT(var30)) var30 = undef
IF (PRESENT(var31)) var31 = undef
IF (PRESENT(var32)) var32 = undef
IF (PRESENT(var33)) var33 = undef
IF (PRESENT(var34)) var34 = undef
IF (PRESENT(var35)) var35 = undef
IF (PRESENT(var36)) var36 = undef
IF (PRESENT(var37)) var37 = undef
IF (PRESENT(var38)) var38 = undef
IF (PRESENT(var39)) var39 = undef
IF (PRESENT(var40)) var40 = undef
IF (PRESENT(var41)) var41 = undef
IF (PRESENT(var42)) var42 = undef
IF (PRESENT(var43)) var43 = undef
IF (PRESENT(var44)) var44 = undef
IF (PRESENT(var45)) var45 = undef
IF (PRESENT(var46)) var46 = undef
IF (PRESENT(var47)) var47 = undef
IF (PRESENT(var48)) var48 = undef
IF (PRESENT(var49)) var49 = undef
!=================================================
END SUBROUTINE debug_undef_all
!=================================================

!=================================================
! Test boundary.
!=================================================
SUBROUTINE debug_test_boundary(var)
IMPLICIT NONE
REAL(preci), DIMENSION(ims:ime,kms:kme), INTENT(IN) :: var
INTEGER :: i, k
INTEGER :: i_left, i_right
INTEGER :: k_bot, k_top
!=================================================
DO i = ims, ime - 1
	IF (var(i,kme/2) == undef .AND. var(i+1,kme/2) /= undef) THEN
		i_left = i + 1
	END IF
	IF (var(i,kme/2) /= undef .AND. var(i+1,kme/2) == undef) THEN
		i_right = i
	END IF
END DO

DO k = kms, kme - 1
	IF (var(ite/2,k) == undef .AND. var(ite/2,k+1) /= undef) THEN
		k_bot = k+1
	END IF
	IF (var(ite/2,k) /= undef .AND. var(ite/2,k+1) == undef) THEN
		k_top = k
	END IF
END DO
WRITE(*,*) "The boundary of var:"
WRITE(*,*) "i_left: ", i_left
WRITE(*,*) "i_right: ", i_right
WRITE(*,*) "k_bot ", k_bot
WRITE(*,*) "k_top ", k_top
!=================================================
END SUBROUTINE debug_test_boundary
!=================================================

!=================================================
END MODULE sp_module_debug
!=================================================

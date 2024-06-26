SUBROUTINE SRTM_CMBGB23

!     BAND 23:  8050-12850 cm-1 (low - H2O; high - nothing)
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM , JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOESRTM  , ONLY : NGN
USE YOESRTWN , ONLY : NGC, NGS, RWGT
!USE YOESRTWN , ONLY : NGC, NGS, NGN, RWGT
USE YOESRTA23, ONLY : KA, SELFREF, FORREF, SFLUXREF, RAYL, &
                    & KAC, SELFREFC, FORREFC, SFLUXREFC, RAYLC

IMPLICIT NONE

! Local variables
INTEGER(KIND=JPIM) :: JT, JP, IGC, IPR, IPRSM
REAL(KIND=JPRB)    :: ZSUMK, ZSUMF1, ZSUMF2

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_CMBGB23',0,ZHOOK_HANDLE)

DO JT = 1,5
  DO JP = 1,13
    IPRSM = 0
    DO IGC = 1,NGC(8)
      ZSUMK = 0.
      DO IPR = 1, NGN(NGS(7)+IGC)
        IPRSM = IPRSM + 1
        ZSUMK = ZSUMK + KA(JT,JP,IPRSM)*RWGT(IPRSM+112)
      ENDDO
      KAC(JT,JP,IGC) = ZSUMK
    ENDDO
  ENDDO
ENDDO

DO JT = 1,10
  IPRSM = 0
  DO IGC = 1,NGC(8)
    ZSUMK = 0.
    DO IPR = 1, NGN(NGS(7)+IGC)
      IPRSM = IPRSM + 1
      ZSUMK = ZSUMK + SELFREF(JT,IPRSM)*RWGT(IPRSM+112)
    ENDDO
    SELFREFC(JT,IGC) = ZSUMK
  ENDDO
ENDDO

DO JT = 1,3
  IPRSM = 0
  DO IGC = 1,NGC(8)
    ZSUMK = 0.
    DO IPR = 1, NGN(NGS(7)+IGC)
      IPRSM = IPRSM + 1
      ZSUMK = ZSUMK + FORREF(JT,IPRSM)*RWGT(IPRSM+112)
    ENDDO
    FORREFC(JT,IGC) = ZSUMK
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(8)
  ZSUMF1 = 0.
  ZSUMF2 = 0.
  DO IPR = 1, NGN(NGS(7)+IGC)
    IPRSM = IPRSM + 1
    ZSUMF1 = ZSUMF1 + SFLUXREF(IPRSM)
    ZSUMF2 = ZSUMF2 + RAYL(IPRSM)*RWGT(IPRSM+112)
  ENDDO
  SFLUXREFC(IGC) = ZSUMF1
  RAYLC(IGC) = ZSUMF2
ENDDO

!     -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_CMBGB23',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_CMBGB23


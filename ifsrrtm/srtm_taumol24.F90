SUBROUTINE SRTM_TAUMOL24 &
 & ( KIDIA   , KFDIA    , KLEV,&
 & P_FAC00   , P_FAC01  , P_FAC10   , P_FAC11,&
 & K_JP      , K_JT     , K_JT1     , P_ONEMINUS,&
 & P_COLH2O  , P_COLMOL , P_COLO2   , P_COLO3,&
 & K_LAYTROP , P_SELFFAC, P_SELFFRAC, K_INDSELF  , P_FORFAC, P_FORFRAC, K_INDFOR,&
 & P_SFLUXZEN, P_TAUG   , P_TAUR    , PRMU0   &
 & )  

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

!     BAND 24:  12850-16000 cm-1 (low - H2O,O2; high - O2)

! Modifications
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     JJMorcrette 2003-02-24 adapted to ECMWF environment
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!     JJMorcrette 20110610 Flexible configuration for number of g-points

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARSRTM  , ONLY : JPG
USE YOESRTM  , ONLY : NG24
USE YOESRTA24, ONLY : ABSA, ABSB, FORREFC, SELFREFC, SFLUXREFC, &
 & ABSO3AC, ABSO3BC, RAYLAC, RAYLBC, LAYREFFR, STRRAT  
USE YOESRTWN , ONLY : NSPA, NSPB

IMPLICIT NONE

!-- Output
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC00(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC01(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC10(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC11(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JP(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT1(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_ONEMINUS(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLH2O(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLMOL(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLO2(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLO3(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFRAC(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_INDSELF(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FORFAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FORFRAC(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_INDFOR(KIDIA:KFDIA,KLEV) 

REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SFLUXZEN(KIDIA:KFDIA,JPG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUG(KIDIA:KFDIA,KLEV,JPG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAUR(KIDIA:KFDIA,KLEV,JPG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA)
!- from INTFAC      
!- from INTIND
!- from PRECISE             
!- from PROFDATA             
!- from SELF             
INTEGER(KIND=JPIM) :: IG, IND0, IND1, INDS, INDF, JS, I_LAY, I_LAYSOLFR(KIDIA:KFDIA), I_NLAYERS, IPLON

INTEGER(KIND=JPIM) :: I_LAY_NEXT

REAL(KIND=JPRB) :: Z_FAC000, Z_FAC001, Z_FAC010, Z_FAC011, Z_FAC100, Z_FAC101,&
 & Z_FAC110, Z_FAC111, Z_FS, Z_SPECCOMB, Z_SPECMULT, Z_SPECPARM, &
 & Z_TAURAY  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SRTM_TAUMOL24',0,ZHOOK_HANDLE)

I_NLAYERS = KLEV

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

I_LAYSOLFR(KIDIA:KFDIA) = K_LAYTROP(KIDIA:KFDIA)

DO I_LAY = 1, I_NLAYERS
  I_LAY_NEXT = MIN(I_NLAYERS, I_LAY+1)
  DO IPLON = KIDIA, KFDIA
    IF (PRMU0(IPLON) > 0.0_JPRB) THEN
      IF (I_LAY <= K_LAYTROP(IPLON)) THEN
        IF (K_JP(IPLON,I_LAY) < LAYREFFR .AND. K_JP(IPLON,I_LAY_NEXT) >= LAYREFFR) &
         & I_LAYSOLFR(IPLON) = MIN(I_LAY+1,K_LAYTROP(IPLON))  
        Z_SPECCOMB = P_COLH2O(IPLON,I_LAY) + STRRAT*P_COLO2(IPLON,I_LAY)
        Z_SPECPARM = P_COLH2O(IPLON,I_LAY)/Z_SPECCOMB 
        IF (Z_SPECPARM >= P_ONEMINUS(IPLON)) Z_SPECPARM = P_ONEMINUS(IPLON)
        Z_SPECMULT = 8.*(Z_SPECPARM)
        JS = 1 + INT(Z_SPECMULT)
        Z_FS = MOD(Z_SPECMULT, 1.0_JPRB )
        ! Z_FAC000 = (1. - Z_FS) * P_FAC00(I_LAY)
        ! Z_FAC010 = (1. - Z_FS) * P_FAC10(I_LAY)
        ! Z_FAC100 = Z_FS * P_FAC00(I_LAY)
        ! Z_FAC110 = Z_FS * P_FAC10(I_LAY)
        ! Z_FAC001 = (1. - Z_FS) * P_FAC01(I_LAY)
        ! Z_FAC011 = (1. - Z_FS) * P_FAC11(I_LAY)
        ! Z_FAC101 = Z_FS * P_FAC01(I_LAY)
        ! Z_FAC111 = Z_FS * P_FAC11(I_LAY)
        IND0 = ((K_JP(IPLON,I_LAY)-1)*5+(K_JT(IPLON,I_LAY)-1))*NSPA(24) + JS
        IND1 = (K_JP(IPLON,I_LAY)*5+(K_JT1(IPLON,I_LAY)-1))*NSPA(24) + JS
        INDS = K_INDSELF(IPLON,I_LAY)
        INDF = K_INDFOR(IPLON,I_LAY)

        !  DO IG = 1, NG(24)
!CDIR UNROLL=NG24
        DO IG = 1 , NG24
          Z_TAURAY = P_COLMOL(IPLON,I_LAY) * (RAYLAC(IG,JS) + &
           & Z_FS * (RAYLAC(IG,JS+1) - RAYLAC(IG,JS)))  
          P_TAUG(IPLON,I_LAY,IG) = Z_SPECCOMB * &
           !    & (Z_FAC000 * ABSA(IND0,IG) + &
           !    & Z_FAC100 * ABSA(IND0+1,IG) + &
           !    & Z_FAC010 * ABSA(IND0+9,IG) + &
           !    & Z_FAC110 * ABSA(IND0+10,IG) + &
           !    & Z_FAC001 * ABSA(IND1,IG) + &
           !    & Z_FAC101 * ABSA(IND1+1,IG) + &
           !    & Z_FAC011 * ABSA(IND1+9,IG) + &
           !    & Z_FAC111 * ABSA(IND1+10,IG)) + &
           & (&
           & (1. - Z_FS) * ( ABSA(IND0,IG) * P_FAC00(IPLON,I_LAY) + &
           &                 ABSA(IND0+9,IG) * P_FAC10(IPLON,I_LAY) + &
           &                 ABSA(IND1,IG) * P_FAC01(IPLON,I_LAY) + &
           &                 ABSA(IND1+9,IG) * P_FAC11(IPLON,I_LAY) ) + &
           & Z_FS        * ( ABSA(IND0+1,IG) * P_FAC00(IPLON,I_LAY) + &
           &                 ABSA(IND0+10,IG) * P_FAC10(IPLON,I_LAY) + &
           &                 ABSA(IND1+1,IG) * P_FAC01(IPLON,I_LAY) + &
           &                 ABSA(IND1+10,IG) * P_FAC11(IPLON,I_LAY) ) &
           & ) + &
           & P_COLO3(IPLON,I_LAY) * ABSO3AC(IG) + &
           & P_COLH2O(IPLON,I_LAY) *  &
           & (P_SELFFAC(IPLON,I_LAY) * (SELFREFC(INDS,IG) + &
           & P_SELFFRAC(IPLON,I_LAY) * &
           & (SELFREFC(INDS+1,IG) - SELFREFC(INDS,IG))) + &
           & P_FORFAC(IPLON,I_LAY) * (FORREFC(INDF,IG) +  &
           & P_FORFRAC(IPLON,I_LAY) * &
           & (FORREFC(INDF+1,IG) - FORREFC(INDF,IG))))  
          !     &           + TAURAY
          !    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
          IF (I_LAY == I_LAYSOLFR(IPLON)) P_SFLUXZEN(IPLON,IG) = SFLUXREFC(IG,JS) &
           & + Z_FS * (SFLUXREFC(IG,JS+1) - SFLUXREFC(IG,JS))  
          P_TAUR(IPLON,I_LAY,IG) =  Z_TAURAY
        ENDDO
      ENDIF
    ENDIF
  ENDDO
ENDDO

DO I_LAY = 1, I_NLAYERS
  DO IPLON = KIDIA, KFDIA
    IF (PRMU0(IPLON) > 0.0_JPRB) THEN
      IF (I_LAY >= K_LAYTROP(IPLON)+1) THEN
        IND0 = ((K_JP(IPLON,I_LAY)-13)*5+(K_JT(IPLON,I_LAY)-1))*NSPB(24) + 1
        IND1 = ((K_JP(IPLON,I_LAY)-12)*5+(K_JT1(IPLON,I_LAY)-1))*NSPB(24) + 1

        !  DO IG = 1, NG(24)
!CDIR UNROLL=NG24
        DO IG = 1 , NG24
          Z_TAURAY = P_COLMOL(IPLON,I_LAY) * RAYLBC(IG)
          P_TAUG(IPLON,I_LAY,IG) = P_COLO2(IPLON,I_LAY) * &
           & (P_FAC00(IPLON,I_LAY) * ABSB(IND0,IG) + &
           & P_FAC10(IPLON,I_LAY) * ABSB(IND0+1,IG) + &
           & P_FAC01(IPLON,I_LAY) * ABSB(IND1,IG) + &
           & P_FAC11(IPLON,I_LAY) * ABSB(IND1+1,IG)) + &
           & P_COLO3(IPLON,I_LAY) * ABSO3BC(IG)  
          !     &          + TAURAY
          !    SSA(LAY,IG) = TAURAY/TAUG(LAY,IG)
          P_TAUR(IPLON,I_LAY,IG) = Z_TAURAY
        ENDDO
      ENDIF
    ENDIF
  ENDDO
ENDDO

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_TAUMOL24',1,ZHOOK_HANDLE)

END SUBROUTINE SRTM_TAUMOL24

! radintg_zrgp_mod.fypp - Wrap the block-allocated radiation fields from RADINTG in a FIELD API stack
!
! (C) Copyright 2022- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.


MODULE RADINTG_ZRGP_MOD

USE PARKIND1     , ONLY : JPRB, JPIM
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK

#ifdef HAVE_FIELD_API
USE FIELD_MODULE
USE FIELD_BASIC_MODULE
#endif

IMPLICIT NONE

PRIVATE
PUBLIC :: RADINTG_ZRGP_TYPE

INTEGER, PARAMETER :: NUNDEFLD = -99999999

TYPE RADINTG_ZRGP_TYPE
  ! Field counts and offset indices for ZRGP
  INTEGER(KIND=JPIM) :: IFLDSIN, IFLDSOUT, IFLDSTOT
  INTEGER(KIND=JPIM) :: IINBEG, IINEND, IOUTBEG, IOUTEND
  INTEGER(KIND=JPIM) :: igi
  INTEGER(KIND=JPIM) :: imu0
  INTEGER(KIND=JPIM) :: iamu0
  INTEGER(KIND=JPIM) :: iemiss
  INTEGER(KIND=JPIM) :: its
  INTEGER(KIND=JPIM) :: islm
  INTEGER(KIND=JPIM) :: iccnl
  INTEGER(KIND=JPIM) :: iccno
  INTEGER(KIND=JPIM) :: ibas
  INTEGER(KIND=JPIM) :: itop
  INTEGER(KIND=JPIM) :: igelam
  INTEGER(KIND=JPIM) :: igemu
  INTEGER(KIND=JPIM) :: iclon
  INTEGER(KIND=JPIM) :: islon
  INTEGER(KIND=JPIM) :: iald
  INTEGER(KIND=JPIM) :: ialp
  INTEGER(KIND=JPIM) :: iti
  INTEGER(KIND=JPIM) :: ipr
  INTEGER(KIND=JPIM) :: iqs
  INTEGER(KIND=JPIM) :: iwv
  INTEGER(KIND=JPIM) :: iclc
  INTEGER(KIND=JPIM) :: ilwa
  INTEGER(KIND=JPIM) :: iiwa
  INTEGER(KIND=JPIM) :: iswa
  INTEGER(KIND=JPIM) :: irwa
  INTEGER(KIND=JPIM) :: irra
  INTEGER(KIND=JPIM) :: idp
  INTEGER(KIND=JPIM) :: ifsd
  INTEGER(KIND=JPIM) :: iecpo3
  INTEGER(KIND=JPIM) :: ihpr
  INTEGER(KIND=JPIM) :: iaprs
  INTEGER(KIND=JPIM) :: ihti
  INTEGER(KIND=JPIM) :: ipert
  INTEGER(KIND=JPIM) :: iprogaero
  INTEGER(KIND=JPIM) :: ire_liq
  INTEGER(KIND=JPIM) :: ire_ice
  INTEGER(KIND=JPIM) :: ioverlap
  INTEGER(KIND=JPIM) :: iaero
  INTEGER(KIND=JPIM) :: ifrsod
  INTEGER(KIND=JPIM) :: ifrted
  INTEGER(KIND=JPIM) :: ifrsodc
  INTEGER(KIND=JPIM) :: ifrtedc
  INTEGER(KIND=JPIM) :: iemit
  INTEGER(KIND=JPIM) :: isudu
  INTEGER(KIND=JPIM) :: iuvdf
  INTEGER(KIND=JPIM) :: iparf
  INTEGER(KIND=JPIM) :: iparcf
  INTEGER(KIND=JPIM) :: itincf
  INTEGER(KIND=JPIM) :: ifdir
  INTEGER(KIND=JPIM) :: ifdif
  INTEGER(KIND=JPIM) :: icdir
  INTEGER(KIND=JPIM) :: ilwderivative
  INTEGER(KIND=JPIM) :: iswdirectband
  INTEGER(KIND=JPIM) :: iswdiffuseband
  INTEGER(KIND=JPIM) :: ifrso
  INTEGER(KIND=JPIM) :: iswfc
  INTEGER(KIND=JPIM) :: ifrth
  INTEGER(KIND=JPIM) :: ilwfc
  INTEGER(KIND=JPIM) :: iaer
  INTEGER(KIND=JPIM) :: ioz
  INTEGER(KIND=JPIM) :: iico2
  INTEGER(KIND=JPIM) :: iich4
  INTEGER(KIND=JPIM) :: iin2o
  INTEGER(KIND=JPIM) :: ino2
  INTEGER(KIND=JPIM) :: ic11
  INTEGER(KIND=JPIM) :: ic12
  INTEGER(KIND=JPIM) :: ic22
  INTEGER(KIND=JPIM) :: icl4
  INTEGER(KIND=JPIM) :: igix

#ifdef HAVE_FIELD_API
  ! Field stack wrapper for ZRGP
  CLASS(FIELD_3RB), POINTER :: FIELD_WRAPPER
  TYPE(FIELD_BASIC_PTR), ALLOCATABLE :: MEMBERS(:)
#endif

CONTAINS
  PROCEDURE :: SETUP => RADINTG_ZRGP_SETUP
#ifdef HAVE_FIELD_API
  PROCEDURE :: SETUP_FIELD => RADINTG_ZRGP_SETUP_FIELD
#endif
  PROCEDURE :: COPY_INPUTS => IFS_COPY_INPUTS_TO_BLOCKED
  PROCEDURE :: COPY_FLUXES => IFS_COPY_FLUXES_FROM_BLOCKED

END TYPE RADINTG_ZRGP_TYPE

CONTAINS

INTEGER(KIND=JPIM) FUNCTION INDRAD(KNEXT,KFLDS,LDUSE)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KNEXT
INTEGER(KIND=JPIM),INTENT(IN) :: KFLDS
LOGICAL,INTENT(IN) :: LDUSE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('RADINTG:INDRAD',0,ZHOOK_HANDLE)

IF( LDUSE )THEN
  INDRAD=KNEXT
  KNEXT=KNEXT+KFLDS
ELSE
  INDRAD=NUNDEFLD
ENDIF

IF (LHOOK) CALL DR_HOOK('RADINTG:INDRAD',1,ZHOOK_HANDLE)

END FUNCTION INDRAD


SUBROUTINE RADINTG_ZRGP_SETUP( &
  & SELF, NLEV, NLWEMISS, &
  & NLWOUT, NSW, NFSD, NRFTOTAL_RADGRID, &
  & NPROGAER, NRADAER, &
  & LDEBUG, LSPPRAD, LRAYFM, &
  & LAPPROXLWUPDATE, LAPPROXSWUPDATE, &
  & LEPO3RA, LDIAGFORCING)

  USE YOMLUN_ECRAD       , ONLY : NULOUT

  IMPLICIT NONE

  CLASS(RADINTG_ZRGP_TYPE), INTENT(INOUT):: SELF
  INTEGER, INTENT(IN)                    :: NLEV
  INTEGER, INTENT(IN)                    :: NLWEMISS, NLWOUT, NSW
  INTEGER, INTENT(IN)                    :: NFSD, NRFTOTAL_RADGRID
  INTEGER, INTENT(IN)                    :: NPROGAER, NRADAER
  LOGICAL, INTENT(IN)                    :: LDEBUG, LSPPRAD, LRAYFM
  LOGICAL, INTENT(IN)                    :: LAPPROXLWUPDATE, LAPPROXSWUPDATE
  LOGICAL, INTENT(IN)                    :: LEPO3RA, LDIAGFORCING

  INTEGER(KIND=JPIM), ALLOCATABLE     :: MEMBER_MAP(:)
  INTEGER(KIND=JPIM)                  :: INEXT
  REAL(KIND=JPHOOK)                   :: ZHOOK_HANDLE

#ifdef BITIDENTITY_TESTING
  LOGICAL, PARAMETER :: LTEST_BITID = .TRUE.
#else
  LOGICAL, PARAMETER :: LTEST_BITID = .FALSE.
#endif

#include "abor1.intfb.h"

  IF (LHOOK) CALL DR_HOOK('RADINTG_ZRGP_SETUP',0,ZHOOK_HANDLE)

  INEXT = 1
  SELF%IINBEG=1
  SELF%igi = INDRAD( INEXT, 1, ldebug)
  SELF%imu0 = INDRAD( INEXT, 1, .true.)
  SELF%iamu0 = INDRAD( INEXT, 1, .true.)
  SELF%iemiss = INDRAD( INEXT, nlwemiss, .true.)
  SELF%its = INDRAD( INEXT, 1, .true.)
  SELF%islm = INDRAD( INEXT, 1, .true.)
  SELF%iccnl = INDRAD( INEXT, 1, .true.)
  SELF%iccno = INDRAD( INEXT, 1, .true.)
  SELF%ibas = INDRAD( INEXT, 1, .true.)
  SELF%itop = INDRAD( INEXT, 1, .true.)
  SELF%igelam = INDRAD( INEXT, 1, .true.)
  SELF%igemu = INDRAD( INEXT, 1, .true.)
  SELF%iclon = INDRAD( INEXT, 1, .true.)
  SELF%islon = INDRAD( INEXT, 1, .true.)
  SELF%iald = INDRAD( INEXT, nsw, .true.)
  SELF%ialp = INDRAD( INEXT, nsw, .true.)
  SELF%iti = INDRAD( INEXT, nlev, .true.)
  SELF%ipr = INDRAD( INEXT, nlev, .true.)
  SELF%iqs = INDRAD( INEXT, nlev, .true.)
  SELF%iwv = INDRAD( INEXT, nlev, .true.)
  SELF%iclc = INDRAD( INEXT, nlev, .true.)
  SELF%ilwa = INDRAD( INEXT, nlev, .true.)
  SELF%iiwa = INDRAD( INEXT, nlev, .true.)
  SELF%iswa = INDRAD( INEXT, nlev, .true.)
  SELF%irwa = INDRAD( INEXT, nlev, .true.)
  SELF%irra = INDRAD( INEXT, nlev, .true.)
  SELF%idp = INDRAD( INEXT, nlev, .true.)
  SELF%ifsd = INDRAD( INEXT, nlev, nfsd>0)
  SELF%iecpo3 = INDRAD( INEXT, nlev, .not.lrayfm.and.lepo3ra)
  SELF%ihpr = INDRAD( INEXT, nlev+1, .true.)
  SELF%iaprs = INDRAD( INEXT, nlev+1, .true.)
  SELF%ihti = INDRAD( INEXT, nlev+1, .true.)
  SELF%ipert = INDRAD( INEXT, nrftotal_radgrid, lspprad)
  SELF%iprogaero = INDRAD( INEXT, nprogaer*nlev, .not.lrayfm.and.nradaer>0)
  SELF%ire_liq = INDRAD( INEXT, nlev, ltest_bitid)
  SELF%ire_ice = INDRAD( INEXT, nlev, ltest_bitid)
  SELF%ioverlap = INDRAD( INEXT, nlev, ltest_bitid)
  SELF%IINEND = INEXT-1
  SELF%IOUTBEG = INEXT
  SELF%iaero = INDRAD( INEXT, nradaer*nlev, .true.)
  SELF%ifrsod = INDRAD( INEXT, 1, .true.)
  SELF%ifrted = INDRAD( INEXT, nlwout, .true.)
  SELF%ifrsodc = INDRAD( INEXT, 1, .true.)
  SELF%ifrtedc = INDRAD( INEXT, 1, .true.)
  SELF%iemit = INDRAD( INEXT, 1, .true.)
  SELF%isudu = INDRAD( INEXT, 1, .true.)
  SELF%iuvdf = INDRAD( INEXT, 1, .true.)
  SELF%iparf = INDRAD( INEXT, 1, .true.)
  SELF%iparcf = INDRAD( INEXT, 1, .true.)
  SELF%itincf = INDRAD( INEXT, 1, .true.)
  SELF%ifdir = INDRAD( INEXT, 1, .true.)
  SELF%ifdif = INDRAD( INEXT, 1, .true.)
  SELF%icdir = INDRAD( INEXT, 1, .true.)
  SELF%ilwderivative = INDRAD( INEXT, nlev+1, lapproxlwupdate)
  SELF%iswdirectband = INDRAD( INEXT, nsw, lapproxswupdate)
  SELF%iswdiffuseband = INDRAD( INEXT, nsw, lapproxswupdate)
  SELF%ifrso = INDRAD( INEXT, nlev+1, .true.)
  SELF%iswfc = INDRAD( INEXT, nlev+1, .true.)
  SELF%ifrth = INDRAD( INEXT, nlev+1, .true.)
  SELF%ilwfc = INDRAD( INEXT, nlev+1, .true.)
  SELF%iaer = INDRAD( INEXT, 6*nlev, ldiagforcing)
  SELF%ioz = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%iico2 = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%iich4 = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%iin2o = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%ino2 = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%ic11 = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%ic12 = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%ic22 = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%icl4 = INDRAD( INEXT, nlev, ldiagforcing)
  SELF%igix = INDRAD( INEXT, 1, ldebug)
  SELF%IOUTEND = INEXT-1
  SELF%iaer = INDRAD( INEXT, 6*nlev, .not.ldiagforcing)
  SELF%ioz = INDRAD( INEXT, nlev, .not.(ldiagforcing.or.lrayfm))
  SELF%iico2 = INDRAD( INEXT, nlev, .not.ldiagforcing)
  SELF%iich4 = INDRAD( INEXT, nlev, .not.ldiagforcing)
  SELF%iin2o = INDRAD( INEXT, nlev, .not.ldiagforcing)
  SELF%ino2 = INDRAD( INEXT, nlev, .not.ldiagforcing)
  SELF%ic11 = INDRAD( INEXT, nlev, .not.ldiagforcing)
  SELF%ic12 = INDRAD( INEXT, nlev, .not.ldiagforcing)
  SELF%ic22 = INDRAD( INEXT, nlev, .not.ldiagforcing)
  SELF%icl4 = INDRAD( INEXT, nlev, .not.ldiagforcing)

  IF (LDEBUG) THEN
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IGI',SELF%igi
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IMU0',SELF%imu0
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IAMU0',SELF%iamu0
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IEMISS',SELF%iemiss
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ITS',SELF%its
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ISLM',SELF%islm
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ICCNL',SELF%iccnl
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ICCNO',SELF%iccno
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IBAS',SELF%ibas
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ITOP',SELF%itop
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IGELAM',SELF%igelam
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IGEMU',SELF%igemu
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ICLON',SELF%iclon
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ISLON',SELF%islon
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IALD',SELF%iald
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IALP',SELF%ialp
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ITI',SELF%iti
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IPR',SELF%ipr
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IQS',SELF%iqs
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IWV',SELF%iwv
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ICLC',SELF%iclc
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ILWA',SELF%ilwa
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IIWA',SELF%iiwa
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ISWA',SELF%iswa
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IRWA',SELF%irwa
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IRRA',SELF%irra
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IDP',SELF%idp
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IFSD',SELF%ifsd
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IECPO3',SELF%iecpo3
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IHPR',SELF%ihpr
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IAPRS',SELF%iaprs
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IHTI',SELF%ihti
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IPERT',SELF%ipert
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IPROGAERO',SELF%iprogaero
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IRE_LIQ',SELF%ire_liq
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IRE_ICE',SELF%ire_ice
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IOVERLAP',SELF%ioverlap
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IAERO',SELF%iaero
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IFRSOD',SELF%ifrsod
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IFRTED',SELF%ifrted
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IFRSODC',SELF%ifrsodc
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IFRTEDC',SELF%ifrtedc
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IEMIT',SELF%iemit
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ISUDU',SELF%isudu
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IUVDF',SELF%iuvdf
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IPARF',SELF%iparf
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IPARCF',SELF%iparcf
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ITINCF',SELF%itincf
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IFDIR',SELF%ifdir
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IFDIF',SELF%ifdif
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ICDIR',SELF%icdir
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ILWDERIVATIVE',SELF%ilwderivative
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ISWDIRECTBAND',SELF%iswdirectband
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ISWDIFFUSEBAND',SELF%iswdiffuseband
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IFRSO',SELF%ifrso
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ISWFC',SELF%iswfc
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IFRTH',SELF%ifrth
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ILWFC',SELF%ilwfc
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IAER',SELF%iaer
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IOZ',SELF%ioz
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IICO2',SELF%iico2
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IICH4',SELF%iich4
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IIN2O',SELF%iin2o
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'INO2',SELF%ino2
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IC11',SELF%ic11
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IC12',SELF%ic12
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IC22',SELF%ic22
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'ICL4',SELF%icl4
    WRITE(NULOUT,'("RADINTG: ",A7,"=",I12)') 'IGIX',SELF%igix
  ENDIF

  SELF%IFLDSIN = SELF%IINEND - SELF%IINBEG + 1
  SELF%IFLDSOUT = SELF%IOUTEND - SELF%IOUTBEG + 1
  SELF%IFLDSTOT = INEXT - 1

  WRITE(NULOUT,'("RADINTG: IFLDSIN   =",I12)')SELF%IFLDSIN
  WRITE(NULOUT,'("RADINTG: IFLDSOUT  =",I12)')SELF%IFLDSOUT
  WRITE(NULOUT,'("RADINTG: IFLDSTOT  =",I12)')SELF%IFLDSTOT

  IF (LHOOK) CALL DR_HOOK('RADINTG_ZRGP_SETUP',1,ZHOOK_HANDLE)

END SUBROUTINE RADINTG_ZRGP_SETUP

#ifdef HAVE_FIELD_API

SUBROUTINE FIELD_INDRAD(MEMBER_MAP, KIDX, KNEXT, KFLDS, LDUSE)
  INTEGER(KIND=JPIM), INTENT(INOUT) :: MEMBER_MAP(:)
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDX
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KNEXT
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLDS
  LOGICAL, INTENT(IN) :: LDUSE
  INTEGER(KIND=JPIM) :: ISTART, IEND
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('RADINTG:FIELD_INDRAD',0,ZHOOK_HANDLE)

  ISTART = KNEXT
  IF( LDUSE .AND. KFLDS > 0 ) THEN
      IEND = ISTART + KFLDS - 1
      KNEXT = IEND + 1
  ELSE
      IEND = ISTART
  ENDIF
  MEMBER_MAP(2*KIDX-1) = ISTART
  MEMBER_MAP(2*KIDX) = IEND

  IF (LHOOK) CALL DR_HOOK('RADINTG:FIELD_INDRAD',1,ZHOOK_HANDLE)

END SUBROUTINE FIELD_INDRAD

SUBROUTINE RADINTG_ZRGP_SETUP_FIELD( &
  & SELF, ZRGP, NLEV, NLWEMISS, &
  & NLWOUT, NSW, NFSD, NRFTOTAL_RADGRID, &
  & NPROGAER, NRADAER, &
  & LDEBUG, LSPPRAD, LRAYFM, &
  & LAPPROXLWUPDATE, LAPPROXSWUPDATE, &
  & LEPO3RA, LDIAGFORCING)

  USE FIELD_FACTORY_MODULE
  USE YOMLUN_ECRAD, ONLY: NULOUT

  IMPLICIT NONE

  CLASS(RADINTG_ZRGP_TYPE), INTENT(INOUT):: SELF
  REAL(KIND=JPRB), INTENT(INOUT), TARGET :: ZRGP(:,:,:)
  INTEGER, INTENT(IN)                    :: NLEV
  INTEGER, INTENT(IN)                    :: NLWEMISS, NLWOUT, NSW
  INTEGER, INTENT(IN)                    :: NFSD, NRFTOTAL_RADGRID
  INTEGER, INTENT(IN)                    :: NPROGAER, NRADAER
  LOGICAL, INTENT(IN)                    :: LDEBUG, LSPPRAD, LRAYFM
  LOGICAL, INTENT(IN)                    :: LAPPROXLWUPDATE, LAPPROXSWUPDATE
  LOGICAL, INTENT(IN)                    :: LEPO3RA, LDIAGFORCING

  INTEGER(KIND=JPIM), ALLOCATABLE     :: MEMBER_MAP(:)
  INTEGER(KIND=JPIM), ALLOCATABLE     :: MEMBER_RANKS(:)
  INTEGER(KIND=JPIM)                  :: INEXT
  REAL(KIND=JPHOOK)                   :: ZHOOK_HANDLE

#ifdef BITIDENTITY_TESTING
  LOGICAL, PARAMETER :: LTEST_BITID = .TRUE.
#else
  LOGICAL, PARAMETER :: LTEST_BITID = .FALSE.
#endif

#include "abor1.intfb.h"

  IF (LHOOK) CALL DR_HOOK('RADINTG_ZRGP_SETUP_FIELD',0,ZHOOK_HANDLE)

  ALLOCATE(MEMBER_MAP(158))
  ALLOCATE(MEMBER_RANKS(79))

  INEXT = 1
  ! igi
  CALL FIELD_INDRAD( MEMBER_MAP, 1, INEXT, 1, ldebug)
  MEMBER_RANKS(1) = 2
  ! imu0
  CALL FIELD_INDRAD( MEMBER_MAP, 2, INEXT, 1, .true.)
  MEMBER_RANKS(2) = 2
  ! iamu0
  CALL FIELD_INDRAD( MEMBER_MAP, 3, INEXT, 1, .true.)
  MEMBER_RANKS(3) = 2
  ! iemiss
  CALL FIELD_INDRAD( MEMBER_MAP, 4, INEXT, nlwemiss, .true.)
  MEMBER_RANKS(4) = 3
  ! its
  CALL FIELD_INDRAD( MEMBER_MAP, 5, INEXT, 1, .true.)
  MEMBER_RANKS(5) = 2
  ! islm
  CALL FIELD_INDRAD( MEMBER_MAP, 6, INEXT, 1, .true.)
  MEMBER_RANKS(6) = 2
  ! iccnl
  CALL FIELD_INDRAD( MEMBER_MAP, 7, INEXT, 1, .true.)
  MEMBER_RANKS(7) = 2
  ! iccno
  CALL FIELD_INDRAD( MEMBER_MAP, 8, INEXT, 1, .true.)
  MEMBER_RANKS(8) = 2
  ! ibas
  CALL FIELD_INDRAD( MEMBER_MAP, 9, INEXT, 1, .true.)
  MEMBER_RANKS(9) = 2
  ! itop
  CALL FIELD_INDRAD( MEMBER_MAP, 10, INEXT, 1, .true.)
  MEMBER_RANKS(10) = 2
  ! igelam
  CALL FIELD_INDRAD( MEMBER_MAP, 11, INEXT, 1, .true.)
  MEMBER_RANKS(11) = 2
  ! igemu
  CALL FIELD_INDRAD( MEMBER_MAP, 12, INEXT, 1, .true.)
  MEMBER_RANKS(12) = 2
  ! iclon
  CALL FIELD_INDRAD( MEMBER_MAP, 13, INEXT, 1, .true.)
  MEMBER_RANKS(13) = 2
  ! islon
  CALL FIELD_INDRAD( MEMBER_MAP, 14, INEXT, 1, .true.)
  MEMBER_RANKS(14) = 2
  ! iald
  CALL FIELD_INDRAD( MEMBER_MAP, 15, INEXT, nsw, .true.)
  MEMBER_RANKS(15) = 3
  ! ialp
  CALL FIELD_INDRAD( MEMBER_MAP, 16, INEXT, nsw, .true.)
  MEMBER_RANKS(16) = 3
  ! iti
  CALL FIELD_INDRAD( MEMBER_MAP, 17, INEXT, nlev, .true.)
  MEMBER_RANKS(17) = 3
  ! ipr
  CALL FIELD_INDRAD( MEMBER_MAP, 18, INEXT, nlev, .true.)
  MEMBER_RANKS(18) = 3
  ! iqs
  CALL FIELD_INDRAD( MEMBER_MAP, 19, INEXT, nlev, .true.)
  MEMBER_RANKS(19) = 3
  ! iwv
  CALL FIELD_INDRAD( MEMBER_MAP, 20, INEXT, nlev, .true.)
  MEMBER_RANKS(20) = 3
  ! iclc
  CALL FIELD_INDRAD( MEMBER_MAP, 21, INEXT, nlev, .true.)
  MEMBER_RANKS(21) = 3
  ! ilwa
  CALL FIELD_INDRAD( MEMBER_MAP, 22, INEXT, nlev, .true.)
  MEMBER_RANKS(22) = 3
  ! iiwa
  CALL FIELD_INDRAD( MEMBER_MAP, 23, INEXT, nlev, .true.)
  MEMBER_RANKS(23) = 3
  ! iswa
  CALL FIELD_INDRAD( MEMBER_MAP, 24, INEXT, nlev, .true.)
  MEMBER_RANKS(24) = 3
  ! irwa
  CALL FIELD_INDRAD( MEMBER_MAP, 25, INEXT, nlev, .true.)
  MEMBER_RANKS(25) = 3
  ! irra
  CALL FIELD_INDRAD( MEMBER_MAP, 26, INEXT, nlev, .true.)
  MEMBER_RANKS(26) = 3
  ! idp
  CALL FIELD_INDRAD( MEMBER_MAP, 27, INEXT, nlev, .true.)
  MEMBER_RANKS(27) = 3
  ! ifsd
  CALL FIELD_INDRAD( MEMBER_MAP, 28, INEXT, nlev, nfsd>0)
  MEMBER_RANKS(28) = 3
  ! iecpo3
  CALL FIELD_INDRAD( MEMBER_MAP, 29, INEXT, nlev, .not.lrayfm.and.lepo3ra)
  MEMBER_RANKS(29) = 3
  ! ihpr
  CALL FIELD_INDRAD( MEMBER_MAP, 30, INEXT, nlev+1, .true.)
  MEMBER_RANKS(30) = 3
  ! iaprs
  CALL FIELD_INDRAD( MEMBER_MAP, 31, INEXT, nlev+1, .true.)
  MEMBER_RANKS(31) = 3
  ! ihti
  CALL FIELD_INDRAD( MEMBER_MAP, 32, INEXT, nlev+1, .true.)
  MEMBER_RANKS(32) = 3
  ! ipert
  CALL FIELD_INDRAD( MEMBER_MAP, 33, INEXT, nrftotal_radgrid, lspprad)
  MEMBER_RANKS(33) = 3
  ! iprogaero
  CALL FIELD_INDRAD( MEMBER_MAP, 34, INEXT, nprogaer*nlev, .not.lrayfm.and.nradaer>0)
  MEMBER_RANKS(34) = 3
  ! ire_liq
  CALL FIELD_INDRAD( MEMBER_MAP, 35, INEXT, nlev, ltest_bitid)
  MEMBER_RANKS(35) = 3
  ! ire_ice
  CALL FIELD_INDRAD( MEMBER_MAP, 36, INEXT, nlev, ltest_bitid)
  MEMBER_RANKS(36) = 3
  ! ioverlap
  CALL FIELD_INDRAD( MEMBER_MAP, 37, INEXT, nlev, ltest_bitid)
  MEMBER_RANKS(37) = 3
  ! iaero
  CALL FIELD_INDRAD( MEMBER_MAP, 38, INEXT, nradaer*nlev, .true.)
  MEMBER_RANKS(38) = 3
  ! ifrsod
  CALL FIELD_INDRAD( MEMBER_MAP, 39, INEXT, 1, .true.)
  MEMBER_RANKS(39) = 2
  ! ifrted
  CALL FIELD_INDRAD( MEMBER_MAP, 40, INEXT, nlwout, .true.)
  MEMBER_RANKS(40) = 3
  ! ifrsodc
  CALL FIELD_INDRAD( MEMBER_MAP, 41, INEXT, 1, .true.)
  MEMBER_RANKS(41) = 2
  ! ifrtedc
  CALL FIELD_INDRAD( MEMBER_MAP, 42, INEXT, 1, .true.)
  MEMBER_RANKS(42) = 2
  ! iemit
  CALL FIELD_INDRAD( MEMBER_MAP, 43, INEXT, 1, .true.)
  MEMBER_RANKS(43) = 2
  ! isudu
  CALL FIELD_INDRAD( MEMBER_MAP, 44, INEXT, 1, .true.)
  MEMBER_RANKS(44) = 2
  ! iuvdf
  CALL FIELD_INDRAD( MEMBER_MAP, 45, INEXT, 1, .true.)
  MEMBER_RANKS(45) = 2
  ! iparf
  CALL FIELD_INDRAD( MEMBER_MAP, 46, INEXT, 1, .true.)
  MEMBER_RANKS(46) = 2
  ! iparcf
  CALL FIELD_INDRAD( MEMBER_MAP, 47, INEXT, 1, .true.)
  MEMBER_RANKS(47) = 2
  ! itincf
  CALL FIELD_INDRAD( MEMBER_MAP, 48, INEXT, 1, .true.)
  MEMBER_RANKS(48) = 2
  ! ifdir
  CALL FIELD_INDRAD( MEMBER_MAP, 49, INEXT, 1, .true.)
  MEMBER_RANKS(49) = 2
  ! ifdif
  CALL FIELD_INDRAD( MEMBER_MAP, 50, INEXT, 1, .true.)
  MEMBER_RANKS(50) = 2
  ! icdir
  CALL FIELD_INDRAD( MEMBER_MAP, 51, INEXT, 1, .true.)
  MEMBER_RANKS(51) = 2
  ! ilwderivative
  CALL FIELD_INDRAD( MEMBER_MAP, 52, INEXT, nlev+1, lapproxlwupdate)
  MEMBER_RANKS(52) = 3
  ! iswdirectband
  CALL FIELD_INDRAD( MEMBER_MAP, 53, INEXT, nsw, lapproxswupdate)
  MEMBER_RANKS(53) = 3
  ! iswdiffuseband
  CALL FIELD_INDRAD( MEMBER_MAP, 54, INEXT, nsw, lapproxswupdate)
  MEMBER_RANKS(54) = 3
  ! ifrso
  CALL FIELD_INDRAD( MEMBER_MAP, 55, INEXT, nlev+1, .true.)
  MEMBER_RANKS(55) = 3
  ! iswfc
  CALL FIELD_INDRAD( MEMBER_MAP, 56, INEXT, nlev+1, .true.)
  MEMBER_RANKS(56) = 3
  ! ifrth
  CALL FIELD_INDRAD( MEMBER_MAP, 57, INEXT, nlev+1, .true.)
  MEMBER_RANKS(57) = 3
  ! ilwfc
  CALL FIELD_INDRAD( MEMBER_MAP, 58, INEXT, nlev+1, .true.)
  MEMBER_RANKS(58) = 3
  ! iaer
  CALL FIELD_INDRAD( MEMBER_MAP, 59, INEXT, 6*nlev, ldiagforcing)
  MEMBER_RANKS(59) = 3
  ! ioz
  CALL FIELD_INDRAD( MEMBER_MAP, 60, INEXT, nlev, ldiagforcing)
  MEMBER_RANKS(60) = 3
  ! iico2
  CALL FIELD_INDRAD( MEMBER_MAP, 61, INEXT, nlev, ldiagforcing)
  MEMBER_RANKS(61) = 3
  ! iich4
  CALL FIELD_INDRAD( MEMBER_MAP, 62, INEXT, nlev, ldiagforcing)
  MEMBER_RANKS(62) = 3
  ! iin2o
  CALL FIELD_INDRAD( MEMBER_MAP, 63, INEXT, nlev, ldiagforcing)
  MEMBER_RANKS(63) = 3
  ! ino2
  CALL FIELD_INDRAD( MEMBER_MAP, 64, INEXT, nlev, ldiagforcing)
  MEMBER_RANKS(64) = 3
  ! ic11
  CALL FIELD_INDRAD( MEMBER_MAP, 65, INEXT, nlev, ldiagforcing)
  MEMBER_RANKS(65) = 3
  ! ic12
  CALL FIELD_INDRAD( MEMBER_MAP, 66, INEXT, nlev, ldiagforcing)
  MEMBER_RANKS(66) = 3
  ! ic22
  CALL FIELD_INDRAD( MEMBER_MAP, 67, INEXT, nlev, ldiagforcing)
  MEMBER_RANKS(67) = 3
  ! icl4
  CALL FIELD_INDRAD( MEMBER_MAP, 68, INEXT, nlev, ldiagforcing)
  MEMBER_RANKS(68) = 3
  ! igix
  CALL FIELD_INDRAD( MEMBER_MAP, 69, INEXT, 1, ldebug)
  MEMBER_RANKS(69) = 2
  ! iaer
  CALL FIELD_INDRAD( MEMBER_MAP, 70, INEXT, 6*nlev, .not.ldiagforcing)
  MEMBER_RANKS(70) = 3
  ! ioz
  CALL FIELD_INDRAD( MEMBER_MAP, 71, INEXT, nlev, .not.(ldiagforcing.or.lrayfm))
  MEMBER_RANKS(71) = 3
  ! iico2
  CALL FIELD_INDRAD( MEMBER_MAP, 72, INEXT, nlev, .not.ldiagforcing)
  MEMBER_RANKS(72) = 3
  ! iich4
  CALL FIELD_INDRAD( MEMBER_MAP, 73, INEXT, nlev, .not.ldiagforcing)
  MEMBER_RANKS(73) = 3
  ! iin2o
  CALL FIELD_INDRAD( MEMBER_MAP, 74, INEXT, nlev, .not.ldiagforcing)
  MEMBER_RANKS(74) = 3
  ! ino2
  CALL FIELD_INDRAD( MEMBER_MAP, 75, INEXT, nlev, .not.ldiagforcing)
  MEMBER_RANKS(75) = 3
  ! ic11
  CALL FIELD_INDRAD( MEMBER_MAP, 76, INEXT, nlev, .not.ldiagforcing)
  MEMBER_RANKS(76) = 3
  ! ic12
  CALL FIELD_INDRAD( MEMBER_MAP, 77, INEXT, nlev, .not.ldiagforcing)
  MEMBER_RANKS(77) = 3
  ! ic22
  CALL FIELD_INDRAD( MEMBER_MAP, 78, INEXT, nlev, .not.ldiagforcing)
  MEMBER_RANKS(78) = 3
  ! icl4
  CALL FIELD_INDRAD( MEMBER_MAP, 79, INEXT, nlev, .not.ldiagforcing)
  MEMBER_RANKS(79) = 3

  CALL FIELD_NEW(SELF%FIELD_WRAPPER, SELF%MEMBERS, DATA=ZRGP, MEMBER_MAP=MEMBER_MAP, MEMBER_RANKS=MEMBER_RANKS)

  IF (LDEBUG) THEN
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IGI',MEMBER_MAP(1),MEMBER_MAP(2)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IMU0',MEMBER_MAP(3),MEMBER_MAP(4)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IAMU0',MEMBER_MAP(5),MEMBER_MAP(6)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IEMISS',MEMBER_MAP(7),MEMBER_MAP(8)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ITS',MEMBER_MAP(9),MEMBER_MAP(10)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ISLM',MEMBER_MAP(11),MEMBER_MAP(12)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ICCNL',MEMBER_MAP(13),MEMBER_MAP(14)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ICCNO',MEMBER_MAP(15),MEMBER_MAP(16)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IBAS',MEMBER_MAP(17),MEMBER_MAP(18)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ITOP',MEMBER_MAP(19),MEMBER_MAP(20)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IGELAM',MEMBER_MAP(21),MEMBER_MAP(22)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IGEMU',MEMBER_MAP(23),MEMBER_MAP(24)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ICLON',MEMBER_MAP(25),MEMBER_MAP(26)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ISLON',MEMBER_MAP(27),MEMBER_MAP(28)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IALD',MEMBER_MAP(29),MEMBER_MAP(30)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IALP',MEMBER_MAP(31),MEMBER_MAP(32)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ITI',MEMBER_MAP(33),MEMBER_MAP(34)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IPR',MEMBER_MAP(35),MEMBER_MAP(36)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IQS',MEMBER_MAP(37),MEMBER_MAP(38)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IWV',MEMBER_MAP(39),MEMBER_MAP(40)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ICLC',MEMBER_MAP(41),MEMBER_MAP(42)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ILWA',MEMBER_MAP(43),MEMBER_MAP(44)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IIWA',MEMBER_MAP(45),MEMBER_MAP(46)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ISWA',MEMBER_MAP(47),MEMBER_MAP(48)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IRWA',MEMBER_MAP(49),MEMBER_MAP(50)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IRRA',MEMBER_MAP(51),MEMBER_MAP(52)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IDP',MEMBER_MAP(53),MEMBER_MAP(54)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IFSD',MEMBER_MAP(55),MEMBER_MAP(56)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IECPO3',MEMBER_MAP(57),MEMBER_MAP(58)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IHPR',MEMBER_MAP(59),MEMBER_MAP(60)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IAPRS',MEMBER_MAP(61),MEMBER_MAP(62)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IHTI',MEMBER_MAP(63),MEMBER_MAP(64)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IPERT',MEMBER_MAP(65),MEMBER_MAP(66)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IPROGAERO',MEMBER_MAP(67),MEMBER_MAP(68)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IRE_LIQ',MEMBER_MAP(69),MEMBER_MAP(70)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IRE_ICE',MEMBER_MAP(71),MEMBER_MAP(72)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IOVERLAP',MEMBER_MAP(73),MEMBER_MAP(74)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IAERO',MEMBER_MAP(75),MEMBER_MAP(76)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IFRSOD',MEMBER_MAP(77),MEMBER_MAP(78)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IFRTED',MEMBER_MAP(79),MEMBER_MAP(80)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IFRSODC',MEMBER_MAP(81),MEMBER_MAP(82)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IFRTEDC',MEMBER_MAP(83),MEMBER_MAP(84)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IEMIT',MEMBER_MAP(85),MEMBER_MAP(86)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ISUDU',MEMBER_MAP(87),MEMBER_MAP(88)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IUVDF',MEMBER_MAP(89),MEMBER_MAP(90)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IPARF',MEMBER_MAP(91),MEMBER_MAP(92)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IPARCF',MEMBER_MAP(93),MEMBER_MAP(94)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ITINCF',MEMBER_MAP(95),MEMBER_MAP(96)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IFDIR',MEMBER_MAP(97),MEMBER_MAP(98)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IFDIF',MEMBER_MAP(99),MEMBER_MAP(100)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ICDIR',MEMBER_MAP(101),MEMBER_MAP(102)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ILWDERIVATIVE',MEMBER_MAP(103),MEMBER_MAP(104)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ISWDIRECTBAND',MEMBER_MAP(105),MEMBER_MAP(106)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ISWDIFFUSEBAND',MEMBER_MAP(107),MEMBER_MAP(108)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IFRSO',MEMBER_MAP(109),MEMBER_MAP(110)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ISWFC',MEMBER_MAP(111),MEMBER_MAP(112)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IFRTH',MEMBER_MAP(113),MEMBER_MAP(114)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ILWFC',MEMBER_MAP(115),MEMBER_MAP(116)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IAER',MEMBER_MAP(117),MEMBER_MAP(118)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IOZ',MEMBER_MAP(119),MEMBER_MAP(120)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IICO2',MEMBER_MAP(121),MEMBER_MAP(122)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IICH4',MEMBER_MAP(123),MEMBER_MAP(124)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IIN2O',MEMBER_MAP(125),MEMBER_MAP(126)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'INO2',MEMBER_MAP(127),MEMBER_MAP(128)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IC11',MEMBER_MAP(129),MEMBER_MAP(130)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IC12',MEMBER_MAP(131),MEMBER_MAP(132)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IC22',MEMBER_MAP(133),MEMBER_MAP(134)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ICL4',MEMBER_MAP(135),MEMBER_MAP(136)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IGIX',MEMBER_MAP(137),MEMBER_MAP(138)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IAER',MEMBER_MAP(139),MEMBER_MAP(140)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IOZ',MEMBER_MAP(141),MEMBER_MAP(142)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IICO2',MEMBER_MAP(143),MEMBER_MAP(144)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IICH4',MEMBER_MAP(145),MEMBER_MAP(146)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IIN2O',MEMBER_MAP(147),MEMBER_MAP(148)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'INO2',MEMBER_MAP(149),MEMBER_MAP(150)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IC11',MEMBER_MAP(151),MEMBER_MAP(152)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IC12',MEMBER_MAP(153),MEMBER_MAP(154)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'IC22',MEMBER_MAP(155),MEMBER_MAP(156)
    WRITE(NULOUT,'("RADINTG_ZRGP_SETUP_FIELD: ",A7,"=",I8,":",I8)') 'ICL4',MEMBER_MAP(157),MEMBER_MAP(158)
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RADINTG_ZRGP_SETUP_FIELD_API',1,ZHOOK_HANDLE)

END SUBROUTINE RADINTG_ZRGP_SETUP_FIELD

#endif

subroutine ifs_copy_inputs_to_blocked ( &
  & zrgp_fields, yradiation, ncol, nlev, nproma, &
  & single_level, thermodynamics, gas, cloud, aerosol, &
  & sin_latitude, longitude_rad, land_frac, pressure_fl, temperature_fl, &
  & zrgp, thermodynamics_out, iseed)

  use radiation_single_level,   only : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_gas,            only : gas_type, IMassMixingRatio, &
        &   IH2O, ICO2, IO3, IN2O, ICH4, ICFC11, ICFC12, IHCFC22, ICCL4
  use radiation_cloud,          only : cloud_type
  use radiation_aerosol,        only : aerosol_type
  use radiation_setup,          only : tradiation

  implicit none

  class(RADINTG_ZRGP_TYPE), intent(in)    :: zrgp_fields

  ! Configuration for the radiation scheme, IFS style
  type(tradiation), intent(in)          :: yradiation

  integer, intent(in) :: ncol, nlev, nproma  ! Number of columns and levels

  ! Derived types for the inputs to the radiation scheme
  type(single_level_type), intent(in)   :: single_level
  type(thermodynamics_type), intent(in) :: thermodynamics
  type(gas_type), intent(in)            :: gas
  type(cloud_type), intent(in)          :: cloud
  type(aerosol_type), intent(in)        :: aerosol

  ! Additional input data, required for effective radii calculation
  real(jprb), dimension(:), intent(in)   :: sin_latitude, longitude_rad, land_frac
  real(jprb), dimension(:,:), intent(in) :: pressure_fl, temperature_fl

  ! monolithic IFS data structure to pass to radiation scheme
  real(kind=jprb), intent(out), allocatable :: zrgp(:,:,:)

  ! Empty thermodynamics type to store pressure_hl for output at the end
  type(thermodynamics_type), intent(inout), optional  :: thermodynamics_out

  ! Seed for random number generator
  integer, intent(out), allocatable, optional :: iseed(:,:)

  ! number of column blocks, block size
  integer :: ngpblks

  integer :: jrl, ibeg, iend, il, ib, ifld, jemiss, jalb, jlev, joff, jaer

  ! Extract some config values
  ngpblks=(ncol-1)/nproma+1              ! number of column blocks

  ! Allocate blocked data structure
  allocate(zrgp(nproma,zrgp_fields%ifldstot,ngpblks))
  if(present(thermodynamics_out)) allocate(thermodynamics_out%pressure_hl(ncol,nlev+1))
  if(present(iseed)) allocate(iseed(nproma,ngpblks))

  ! First touch
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)&
  !$OMP&PRIVATE(IB,IFLD)
  do ib=1,ngpblks
    do ifld=1,zrgp_fields%ifldstot
      zrgp(:,ifld,ib) = 0._jprb
    enddo
    if(present(iseed)) iseed(:,ib) = 0
  enddo
  !$OMP END PARALLEL DO

  associate(yderad=>yradiation%yrerad, rad_config=>yradiation%rad_config)

    ! REPLACED ich4 with iich4 due to clash
    ! REPLACED in2o with iin2o due to clash
    ! REPLACED ico2 with iico2 due to clash

    !  -------------------------------------------------------
    !
    !  INPUT LOOP
    !
    !  -------------------------------------------------------

    !$OMP PARALLEL DO SCHEDULE(RUNTIME)&
    !$OMP&PRIVATE(JRL,IBEG,IEND,IL,IB,JAER,JOFF,JLEV,JALB)
    do jrl=1,ncol,nproma

      ibeg=jrl
      iend=min(ibeg+nproma-1,ncol)
      il=iend-ibeg+1
      ib=(jrl-1)/nproma+1

      !* RADINTG:  3.      PREPARE INPUT ARRAYS

      ! zrgp(1:il,imu0,ib)  = ???
      zrgp(1:il,zrgp_fields%iamu0,ib)  =  single_level%cos_sza(ibeg:iend)   ! cosine of solar zenith ang (mu0)

      do jemiss=1,yderad%nlwemiss
        zrgp(1:il,zrgp_fields%iemiss+jemiss-1,ib)  =  single_level%lw_emissivity(ibeg:iend,jemiss)
      enddo

      zrgp(1:il,zrgp_fields%its,ib)      = single_level%skin_temperature(ibeg:iend)  ! skin temperature
      zrgp(1:il,zrgp_fields%islm,ib)     = land_frac(ibeg:iend) ! land-sea mask
      zrgp(1:il,zrgp_fields%iccnl,ib)    = yderad%rccnlnd ! CCN over land
      zrgp(1:il,zrgp_fields%iccno,ib)    = yderad%rccnsea ! CCN over sea
      ! zrgp(1:il,ibas,ib)     = ???
      ! zrgp(1:il,itop,ib)     = ???
      zrgp(1:il,zrgp_fields%igelam,ib)   = longitude_rad(ibeg:iend) ! longitude
      zrgp(1:il,zrgp_fields%igemu,ib)    = sin_latitude(ibeg:iend) ! sine of latitude
      ! zrgp(1:il,iclon,ib)    = ???
      ! zrgp(1:il,islon,ib)    = ???

      do jalb=1,yderad%nsw
        zrgp(1:il,zrgp_fields%iald+jalb-1,ib)  =  single_level%sw_albedo(ibeg:iend,jalb)
      enddo

      if (allocated(single_level%sw_albedo_direct)) then
        do jalb=1,yderad%nsw
          zrgp(1:il,zrgp_fields%ialp+jalb-1,ib)  =  single_level%sw_albedo_direct(ibeg:iend,jalb)
        end do
      else
        do jalb=1,yderad%nsw
          zrgp(1:il,zrgp_fields%ialp+jalb-1,ib)  =  single_level%sw_albedo(ibeg:iend,jalb)
        end do
      end if

      do jlev=1,nlev
        zrgp(1:il,zrgp_fields%iti+jlev-1,ib)   = temperature_fl(ibeg:iend,jlev) ! full level temperature
        zrgp(1:il,zrgp_fields%ipr+jlev-1,ib)   = pressure_fl(ibeg:iend,jlev) ! full level pressure
        ! zrgp(1:il,iqs+jlev-1,ib)   = ???
      enddo

      do jlev=1,nlev
        zrgp(1:il,zrgp_fields%iwv+jlev-1,ib)   = gas%mixing_ratio(ibeg:iend,jlev,IH2O) ! this is already in MassMixingRatio units
        if (rad_config%do_clouds) then
          zrgp(1:il,zrgp_fields%iclc+jlev-1,ib)  = cloud%fraction(ibeg:iend,jlev)
          zrgp(1:il,zrgp_fields%ilwa+jlev-1,ib)  = cloud%q_liq(ibeg:iend,jlev)
          zrgp(1:il,zrgp_fields%iiwa+jlev-1,ib)  = cloud%q_ice(ibeg:iend,jlev)
        else
          zrgp(1:il,zrgp_fields%iclc+jlev-1,ib)  = 0._jprb
          zrgp(1:il,zrgp_fields%ilwa+jlev-1,ib)  = 0._jprb
          zrgp(1:il,zrgp_fields%iiwa+jlev-1,ib)  = 0._jprb
        endif
        zrgp(1:il,zrgp_fields%iswa+jlev-1,ib)  = 0._jprb  ! snow
        zrgp(1:il,zrgp_fields%irwa+jlev-1,ib)  = 0._jprb  ! rain

        ! zrgp(1:il,irra+jlev-1,ib)  = ???
        ! zrgp(1:il,idp+jlev-1,ib)   = ???
        ! zrgp(1:il,ifsd+jlev-1,ib)   = ???
        ! zrgp(1:il,iecpo3+jlev-1,ib) = ???
      enddo

      zrgp(1:il,zrgp_fields%iaer:zrgp_fields%iaer+nlev,ib)  =  0._jprb ! old aerosol, not used
      if (yderad%naermacc == 1) then
        joff=zrgp_fields%iaero
        do jaer=1,rad_config%n_aerosol_types
          do jlev=1,nlev
            zrgp(1:il,joff,ib) = aerosol%mixing_ratio(ibeg:iend,jlev,jaer)
            joff=joff+1
          enddo
        enddo
      endif

      do jlev=1,nlev+1
        ! zrgp(1:il,ihpr+jlev-1,ib)  = ???
        zrgp(1:il,zrgp_fields%iaprs+jlev-1,ib) = thermodynamics%pressure_hl(ibeg:iend,jlev)
        zrgp(1:il,zrgp_fields%ihti+jlev-1,ib)  = thermodynamics%temperature_hl(ibeg:iend,jlev)
      enddo

      ! -- by default, globally averaged concentrations (mmr)
      call gas%get(ICO2, IMassMixingRatio, zrgp(1:il,zrgp_fields%iico2:zrgp_fields%iico2+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICH4, IMassMixingRatio, zrgp(1:il,zrgp_fields%iich4:zrgp_fields%iich4+nlev-1,ib), istartcol=ibeg)
      call gas%get(IN2O, IMassMixingRatio, zrgp(1:il,zrgp_fields%iin2o:zrgp_fields%iin2o+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICFC11, IMassMixingRatio, zrgp(1:il,zrgp_fields%ic11:zrgp_fields%ic11+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICFC12, IMassMixingRatio, zrgp(1:il,zrgp_fields%ic12:zrgp_fields%ic12+nlev-1,ib), istartcol=ibeg)
      call gas%get(IHCFC22,IMassMixingRatio, zrgp(1:il,zrgp_fields%ic22:zrgp_fields%ic22+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICCL4,  IMassMixingRatio, zrgp(1:il,zrgp_fields%icl4:zrgp_fields%icl4+nlev-1,ib), istartcol=ibeg)
      call gas%get(IO3, IMassMixingRatio, zrgp(1:il,zrgp_fields%ioz:zrgp_fields%ioz+nlev-1,ib), istartcol=ibeg)
      ! convert ozone kg/kg to Pa*kg/kg
      ! do jlev=1,nlev
      !   zrgp(1:il,zrgp_fields%ioz+jlev-1,ib)  = zrgp(1:il,zrgp_fields%ioz+jlev-1,ib) &
      !         &                       * (thermodynamics%pressure_hl(ibeg:iend,jlev+1) &
      !         &                         - thermodynamics%pressure_hl(ibeg:iend,jlev))
      ! enddo

      ! local workaround variables for standalone input files
#ifdef BITIDENTITY_TESTING
      ! To validate results against standalone ecrad, we overwrite effective
      ! radii, cloud overlap and seed with input values
      if (rad_config%do_clouds) then
        do jlev=1,nlev
          ! missing full-level temperature and pressure as well as land-sea-mask
          zrgp(1:il,zrgp_fields%ire_liq+jlev-1,ib) = cloud%re_liq(ibeg:iend,jlev)
          zrgp(1:il,zrgp_fields%ire_ice+jlev-1,ib) = cloud%re_ice(ibeg:iend,jlev)
        enddo
        do jlev=1,nlev-1
          ! for the love of it, I can't figure this one out. Probably to do with
          ! my crude approach of setting PGEMU?
          zrgp(1:il,zrgp_fields%ioverlap+jlev-1,ib) = cloud%overlap_param(ibeg:iend,jlev)
        enddo
        if(present(iseed)) iseed(1:il,ib) = single_level%iseed(ibeg:iend)
      else
        do jlev=1,nlev
          ! missing full-level temperature and pressure as well as land-sea-mask
          zrgp(1:il,zrgp_fields%ire_liq+jlev-1,ib) = 0._jprb
          zrgp(1:il,zrgp_fields%ire_ice+jlev-1,ib) = 0._jprb
        enddo
        do jlev=1,nlev-1
          zrgp(1:il,zrgp_fields%ioverlap+jlev-1,ib) = 0._jprb
        enddo
        if(present(iseed)) iseed(1:il,ib) = 0
      endif ! do_clouds
#endif
    enddo
    !$OMP END PARALLEL DO

    ! Store pressure for output
    if(present(thermodynamics_out)) thermodynamics_out%pressure_hl(:,:) = thermodynamics%pressure_hl(:,:)

  end associate

end subroutine ifs_copy_inputs_to_blocked

subroutine ifs_copy_fluxes_from_blocked(&
    & zrgp_fields, yradiation, ncol, nlev, nproma, &
    & zrgp, flux, flux_sw_direct_normal, flux_uv, flux_par, flux_par_clear,&
    & emissivity_out, flux_diffuse_band, flux_direct_band)
  use radiation_setup,          only : tradiation
  use radiation_flux,           only : flux_type

  class(radintg_zrgp_type), intent(in)     :: zrgp_fields

  ! Configuration for the radiation scheme, IFS style
  type(tradiation), intent(in)          :: yradiation

  integer, intent(in) :: ncol, nlev, nproma  ! Number of columns and levels

  ! monolithic IFS data structure passed to radiation scheme
  real(kind=jprb), intent(inout), allocatable :: zrgp(:,:,:)

  ! Derived type containing outputs from the radiation scheme
  type(flux_type), intent(inout)              :: flux

  ! Additional output fluxes as arrays
  real(jprb), dimension(:), intent(inout)     :: flux_sw_direct_normal, flux_uv, flux_par,&
                                                 & flux_par_clear, emissivity_out
  real(jprb), dimension(:,:), intent(inout) :: flux_diffuse_band, flux_direct_band

  ! number of column blocks, block size
  integer :: ngpblks

  integer :: jrl, ibeg, iend, il, ib, jlev, jg

  ! Extract some config values
  ngpblks=(ncol-1)/nproma+1              ! number of column blocks

    !  -------------------------------------------------------
    !
    !  OUTPUT LOOP
    !
    !  -------------------------------------------------------

    !$OMP PARALLEL DO SCHEDULE(RUNTIME)&
    !$OMP&PRIVATE(JRL,IBEG,IEND,IL,IB,JLEV,JG)
    do jrl=1,ncol,nproma
      ibeg=jrl
      iend=min(ibeg+nproma-1,ncol)
      il=iend-ibeg+1
      ib=(jrl-1)/nproma+1

      do jlev=1,nlev+1
        flux%sw_up(ibeg:iend,jlev) = zrgp(1:il,zrgp_fields%ifrso+jlev-1,ib)
        flux%lw_up(ibeg:iend,jlev) = zrgp(1:il,zrgp_fields%ifrth+jlev-1,ib)
        flux%sw_up_clear(ibeg:iend,jlev) = zrgp(1:il,zrgp_fields%iswfc+jlev-1,ib)
        flux%lw_up_clear(ibeg:iend,jlev) = zrgp(1:il,zrgp_fields%ilwfc+jlev-1,ib)
        if (yradiation%yrerad%lapproxlwupdate) then
          flux%lw_derivatives(ibeg:iend,jlev) = zrgp(1:il,zrgp_fields%ilwderivative+jlev-1,ib)
        else
          flux%lw_derivatives(ibeg:iend,jlev) = 0.0_jprb
        endif
      end do
      flux%sw_dn(ibeg:iend,nlev+1) = zrgp(1:il,zrgp_fields%ifrsod,ib)
      flux%lw_dn(ibeg:iend,nlev+1) = zrgp(1:il,zrgp_fields%ifrted,ib)
      flux%sw_dn_clear(ibeg:iend,nlev+1) = zrgp(1:il,zrgp_fields%ifrsodc,ib)
      flux%lw_dn_clear(ibeg:iend,nlev+1) = zrgp(1:il,zrgp_fields%ifrtedc,ib)
      flux%sw_dn_direct(ibeg:iend,nlev+1) = zrgp(1:il,zrgp_fields%ifdir,ib)
      flux%sw_dn_direct_clear(ibeg:iend,nlev+1) = zrgp(1:il,zrgp_fields%icdir,ib)
      flux_sw_direct_normal(ibeg:iend) = zrgp(1:il,zrgp_fields%isudu,ib)
      flux_uv(ibeg:iend) = zrgp(1:il,zrgp_fields%iuvdf,ib)
      flux_par(ibeg:iend) = zrgp(1:il,zrgp_fields%iparf,ib)
      flux_par_clear(ibeg:iend) = zrgp(1:il,zrgp_fields%iparcf,ib)
      flux%sw_dn(ibeg:iend,1) = zrgp(1:il,zrgp_fields%itincf,ib)
      emissivity_out(ibeg:iend) = zrgp(1:il,zrgp_fields%iemit,ib)
      if (yradiation%yrerad%lapproxswupdate) then
        do jg=1,yradiation%yrerad%nsw
          flux_diffuse_band(ibeg:iend,jg) = zrgp(1:il,zrgp_fields%iswdiffuseband+jg-1,ib)
          flux_direct_band(ibeg:iend,jg) = zrgp(1:il,zrgp_fields%iswdirectband+jg-1,ib)
        end do
      else
        flux_diffuse_band(ibeg:iend,:) = 0.0_jprb
        flux_direct_band(ibeg:iend,:) = 0.0_jprb
      endif
    end do

    deallocate(zrgp)

end subroutine ifs_copy_fluxes_from_blocked

END MODULE RADINTG_ZRGP_MOD

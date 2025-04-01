! ifs_blocking.F90 - Reshuffle ecRad data into an NPROMA-blocked data structure
!
! (C) Copyright 2022- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Balthasar Reuter
! Email:   balthasar.reuter@ecmwf.int
!

module ifs_blocking

  use parkind1,                 only : jprb, jpim ! Working precision, integer type

#ifdef HAVE_FIELD_API
  use field_module
  use field_stack_module
#endif

  implicit none

  public

! Variable offsets in ZRGP
  type ifs_config_type
  integer :: ifldstot
  integer :: igi
  integer :: iamu0
  integer :: iemiss
  integer :: its
  integer :: islm
  integer :: iccnl
  integer :: ibas
  integer :: itop
  integer :: igelam
  integer :: igemu
  integer :: iclon
  integer :: islon
  integer :: iald
  integer :: ialp
  integer :: iti
  integer :: ipr
  integer :: iqs
  integer :: iwv
  integer :: iclc
  integer :: ilwa
  integer :: iiwa
  integer :: iswa
  integer :: irwa
  integer :: irra
  integer :: idp
  integer :: ioz
  integer :: ihpr
  integer :: iaprs
  integer :: ihti
  integer :: iaero
  integer :: ifrsod
  integer :: ifrted
  integer :: ifrsodc
  integer :: ifrtedc
  integer :: iemit
  integer :: isudu
  integer :: iuvdf
  integer :: iparf
  integer :: iparcf
  integer :: itincf
  integer :: ifdir
  integer :: ifdif
  integer :: icdir
  integer :: ilwderivative
  integer :: iswdirectband
  integer :: iswdiffuseband
  integer :: ifrso
  integer :: iswfc
  integer :: ifrth
  integer :: ilwfc
  integer :: iaer
  integer :: iico2
  integer :: iich4
  integer :: iin2o
  integer :: ino2
  integer :: ic11
  integer :: ic12
  integer :: ic22
  integer :: icl4
  integer :: igix
  integer :: iccno
  integer :: ire_liq
  integer :: ire_ice
  integer :: ioverlap
end type ifs_config_type

#ifdef HAVE_FIELD_API
! Field stack wrapper for ZRGP
  TYPE RADINTG_ZRGP_TYPE
  CLASS(FIELD_3RB), POINTER :: FIELD_WRAPPER
  TYPE(FIELD_3RB_STACK_MEMBER_PTR), ALLOCATABLE :: MEMBERS(:)
  END TYPE RADINTG_ZRGP_TYPE
#endif

contains

integer(kind=jpim) function indrad(knext,kflds,lduse)

  integer(kind=jpim), intent(inout) :: knext
  integer(kind=jpim), intent(in) :: kflds
  logical, intent(in) :: lduse

  if( lduse ) then
    indrad=knext
    knext=knext+kflds
  else
    indrad=-99999999
  endif

end function indrad

subroutine ifs_setup_indices (driver_config, ifs_config, yradiation, nlev)

  use radiation_io,             only : nulout
  use ecrad_driver_config,      only : driver_config_type
  use radiation_setup,          only : tradiation

  ! Configuration specific to this driver
  type(driver_config_type), intent(in)     :: driver_config
  type(ifs_config_type), intent(inout)     :: ifs_config

  ! Configuration for the radiation scheme, IFS style
  type(tradiation), intent(inout)          :: yradiation

  integer, intent(inout) :: nlev

  integer :: inext

  !  INITIALISE INDICES FOR VARIABLE

  ! INDRAD is a CONTAIN'd function (now a module function)

  inext = 1
  ifs_config%igi = indrad( &
      & inext, 1, &
      & .false.)
  ifs_config%iamu0 = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%iemiss = indrad( &
      & inext, yradiation%yrerad%nlwemiss, &
      & .true.)
  ifs_config%its = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%islm = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%iccnl = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%ibas = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%itop = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%igelam = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%igemu = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%iclon = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%islon = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%iald = indrad( &
      & inext, yradiation%yrerad%nsw, &
      & .true.)
  ifs_config%ialp = indrad( &
      & inext, yradiation%yrerad%nsw, &
      & .true.)
  ifs_config%iti = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%ipr = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%iqs = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%iwv = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%iclc = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%ilwa = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%iiwa = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%iswa = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%irwa = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%irra = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%idp = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%ioz = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%ihpr = indrad( &
      & inext, nlev+1, &
      & .true.)
  ifs_config%iaprs = indrad( &
      & inext, nlev+1, &
      & .true.)
  ifs_config%ihti = indrad( &
      & inext, nlev+1, &
      & .true.)
  ifs_config%iaero = indrad( &
      & inext, yradiation%rad_config%n_aerosol_types*nlev, &
      & .true.)
  ifs_config%ifrsod = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%ifrted = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%ifrsodc = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%ifrtedc = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%iemit = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%isudu = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%iuvdf = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%iparf = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%iparcf = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%itincf = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%ifdir = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%ifdif = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%icdir = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%ilwderivative = indrad( &
      & inext, nlev+1, &
      & yradiation%yrerad%lapproxlwupdate)
  ifs_config%iswdirectband = indrad( &
      & inext, yradiation%yrerad%nsw, &
      & yradiation%yrerad%lapproxswupdate)
  ifs_config%iswdiffuseband = indrad( &
      & inext, yradiation%yrerad%nsw, &
      & yradiation%yrerad%lapproxswupdate)
  ifs_config%ifrso = indrad( &
      & inext, nlev+1, &
      & .true.)
  ifs_config%iswfc = indrad( &
      & inext, nlev+1, &
      & .true.)
  ifs_config%ifrth = indrad( &
      & inext, nlev+1, &
      & .true.)
  ifs_config%ilwfc = indrad( &
      & inext, nlev+1, &
      & .true.)
  ifs_config%iaer = indrad( &
      & inext, 6*nlev, &
      & .true.)
  ifs_config%iico2 = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%iich4 = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%iin2o = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%ino2 = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%ic11 = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%ic12 = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%ic22 = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%icl4 = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%igix = indrad( &
      & inext, 1, &
      & .false.)
  ifs_config%iccno = indrad( &
      & inext, 1, &
      & .true.)
  ifs_config%ire_liq = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%ire_ice = indrad( &
      & inext, nlev, &
      & .true.)
  ifs_config%ioverlap = indrad( &
      & inext, nlev-1, &
      & .true.)
  ifs_config%ifldstot = inext  - 1

  if( driver_config%iverbose > 4 )then
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%igi', &
        & ifs_config%igi, &
        & 1, &
        & '.false.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iamu0', &
        & ifs_config%iamu0, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iemiss', &
        & ifs_config%iemiss, &
        & yradiation%yrerad%nlwemiss, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%its', &
        & ifs_config%its, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%islm', &
        & ifs_config%islm, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iccnl', &
        & ifs_config%iccnl, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ibas', &
        & ifs_config%ibas, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%itop', &
        & ifs_config%itop, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%igelam', &
        & ifs_config%igelam, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%igemu', &
        & ifs_config%igemu, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iclon', &
        & ifs_config%iclon, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%islon', &
        & ifs_config%islon, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iald', &
        & ifs_config%iald, &
        & yradiation%yrerad%nsw, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ialp', &
        & ifs_config%ialp, &
        & yradiation%yrerad%nsw, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iti', &
        & ifs_config%iti, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ipr', &
        & ifs_config%ipr, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iqs', &
        & ifs_config%iqs, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iwv', &
        & ifs_config%iwv, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iclc', &
        & ifs_config%iclc, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ilwa', &
        & ifs_config%ilwa, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iiwa', &
        & ifs_config%iiwa, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iswa', &
        & ifs_config%iswa, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%irwa', &
        & ifs_config%irwa, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%irra', &
        & ifs_config%irra, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%idp', &
        & ifs_config%idp, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ioz', &
        & ifs_config%ioz, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ihpr', &
        & ifs_config%ihpr, &
        & nlev+1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iaprs', &
        & ifs_config%iaprs, &
        & nlev+1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ihti', &
        & ifs_config%ihti, &
        & nlev+1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iaero', &
        & ifs_config%iaero, &
        & yradiation%rad_config%n_aerosol_types*nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ifrsod', &
        & ifs_config%ifrsod, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ifrted', &
        & ifs_config%ifrted, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ifrsodc', &
        & ifs_config%ifrsodc, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ifrtedc', &
        & ifs_config%ifrtedc, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iemit', &
        & ifs_config%iemit, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%isudu', &
        & ifs_config%isudu, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iuvdf', &
        & ifs_config%iuvdf, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iparf', &
        & ifs_config%iparf, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iparcf', &
        & ifs_config%iparcf, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%itincf', &
        & ifs_config%itincf, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ifdir', &
        & ifs_config%ifdir, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ifdif', &
        & ifs_config%ifdif, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%icdir', &
        & ifs_config%icdir, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ilwderivative', &
        & ifs_config%ilwderivative, &
        & nlev+1, &
        & 'yradiation%yrerad%lapproxlwupdate'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iswdirectband', &
        & ifs_config%iswdirectband, &
        & yradiation%yrerad%nsw, &
        & 'yradiation%yrerad%lapproxswupdate'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iswdiffuseband', &
        & ifs_config%iswdiffuseband, &
        & yradiation%yrerad%nsw, &
        & 'yradiation%yrerad%lapproxswupdate'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ifrso', &
        & ifs_config%ifrso, &
        & nlev+1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iswfc', &
        & ifs_config%iswfc, &
        & nlev+1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ifrth', &
        & ifs_config%ifrth, &
        & nlev+1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ilwfc', &
        & ifs_config%ilwfc, &
        & nlev+1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iaer', &
        & ifs_config%iaer, &
        & 6*nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iico2', &
        & ifs_config%iico2, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iich4', &
        & ifs_config%iich4, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iin2o', &
        & ifs_config%iin2o, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ino2', &
        & ifs_config%ino2, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ic11', &
        & ifs_config%ic11, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ic12', &
        & ifs_config%ic12, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ic22', &
        & ifs_config%ic22, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%icl4', &
        & ifs_config%icl4, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%igix', &
        & ifs_config%igix, &
        & 1, &
        & '.false.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%iccno', &
        & ifs_config%iccno, &
        & 1, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ire_liq', &
        & ifs_config%ire_liq, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ire_ice', &
        & ifs_config%ire_ice, &
        & nlev, &
        & '.true.'
    write(nulout,'(a," = ",i0,", dim = ",i0,", condition = ",a)') &
        & 'ifs_config%ioverlap', &
        & ifs_config%ioverlap, &
        & nlev-1, &
        & '.true.'
    write(nulout,'("ifldstot=",i0)')ifs_config%ifldstot
  endif

end subroutine ifs_setup_indices

#ifdef HAVE_FIELD_API

SUBROUTINE FIELD_INDRAD(MEMBER_MAP, KIDX, KNEXT, KFLDS, LDUSE)
  INTEGER(KIND=JPIM), INTENT(INOUT) :: MEMBER_MAP(:)
  INTEGER(KIND=JPIM), INTENT(IN) :: KIDX
  INTEGER(KIND=JPIM), INTENT(INOUT) :: KNEXT
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLDS
  LOGICAL, INTENT(IN) :: LDUSE
  INTEGER(KIND=JPIM) :: ISTART, IEND

  ISTART = KNEXT
  IF( LDUSE ) THEN
      IEND = ISTART + KFLDS - 1
      KNEXT = IEND + 1
  ELSE
      IEND = ISTART
  ENDIF
  MEMBER_MAP(2*KIDX-1) = ISTART
  MEMBER_MAP(2*KIDX) = IEND
END SUBROUTINE FIELD_INDRAD

SUBROUTINE IFS_SETUP_FAPI(ZRGP_FIELDS, ZRGP, DRIVER_CONFIG, YRADIATION, NLEV)

  USE FIELD_FACTORY_MODULE
  USE ECRAD_DRIVER_CONFIG,      ONLY : DRIVER_CONFIG_TYPE
  USE RADIATION_SETUP,          ONLY : TRADIATION
  USE RADIATION_IO,             ONLY : NULERR

  IMPLICIT NONE

  TYPE(RADINTG_ZRGP_TYPE), INTENT(INOUT) :: ZRGP_FIELDS
  REAL(KIND=JPRB), INTENT(INOUT), TARGET :: ZRGP(:,:,:)
  TYPE(DRIVER_CONFIG_TYPE), INTENT(IN)   :: DRIVER_CONFIG
  TYPE(TRADIATION), INTENT(IN)           :: YRADIATION
  INTEGER, INTENT(IN)                     :: NLEV

  INTEGER(KIND=JPIM), ALLOCATABLE     :: MEMBER_MAP(:)
  INTEGER(KIND=JPIM)                  :: INEXT

#include "abor1.intfb.h"

  ALLOCATE(MEMBER_MAP(128))

  INEXT = 1
  ! igi
  CALL FIELD_INDRAD( MEMBER_MAP, 1, INEXT, 1, .false.)
  ! iamu0
  CALL FIELD_INDRAD( MEMBER_MAP, 2, INEXT, 1, .true.)
  ! iemiss
  CALL FIELD_INDRAD( MEMBER_MAP, 3, INEXT, yradiation%yrerad%nlwemiss, .true.)
  ! its
  CALL FIELD_INDRAD( MEMBER_MAP, 4, INEXT, 1, .true.)
  ! islm
  CALL FIELD_INDRAD( MEMBER_MAP, 5, INEXT, 1, .true.)
  ! iccnl
  CALL FIELD_INDRAD( MEMBER_MAP, 6, INEXT, 1, .true.)
  ! ibas
  CALL FIELD_INDRAD( MEMBER_MAP, 7, INEXT, 1, .true.)
  ! itop
  CALL FIELD_INDRAD( MEMBER_MAP, 8, INEXT, 1, .true.)
  ! igelam
  CALL FIELD_INDRAD( MEMBER_MAP, 9, INEXT, 1, .true.)
  ! igemu
  CALL FIELD_INDRAD( MEMBER_MAP, 10, INEXT, 1, .true.)
  ! iclon
  CALL FIELD_INDRAD( MEMBER_MAP, 11, INEXT, 1, .true.)
  ! islon
  CALL FIELD_INDRAD( MEMBER_MAP, 12, INEXT, 1, .true.)
  ! iald
  CALL FIELD_INDRAD( MEMBER_MAP, 13, INEXT, yradiation%yrerad%nsw, .true.)
  ! ialp
  CALL FIELD_INDRAD( MEMBER_MAP, 14, INEXT, yradiation%yrerad%nsw, .true.)
  ! iti
  CALL FIELD_INDRAD( MEMBER_MAP, 15, INEXT, nlev, .true.)
  ! ipr
  CALL FIELD_INDRAD( MEMBER_MAP, 16, INEXT, nlev, .true.)
  ! iqs
  CALL FIELD_INDRAD( MEMBER_MAP, 17, INEXT, nlev, .true.)
  ! iwv
  CALL FIELD_INDRAD( MEMBER_MAP, 18, INEXT, nlev, .true.)
  ! iclc
  CALL FIELD_INDRAD( MEMBER_MAP, 19, INEXT, nlev, .true.)
  ! ilwa
  CALL FIELD_INDRAD( MEMBER_MAP, 20, INEXT, nlev, .true.)
  ! iiwa
  CALL FIELD_INDRAD( MEMBER_MAP, 21, INEXT, nlev, .true.)
  ! iswa
  CALL FIELD_INDRAD( MEMBER_MAP, 22, INEXT, nlev, .true.)
  ! irwa
  CALL FIELD_INDRAD( MEMBER_MAP, 23, INEXT, nlev, .true.)
  ! irra
  CALL FIELD_INDRAD( MEMBER_MAP, 24, INEXT, nlev, .true.)
  ! idp
  CALL FIELD_INDRAD( MEMBER_MAP, 25, INEXT, nlev, .true.)
  ! ioz
  CALL FIELD_INDRAD( MEMBER_MAP, 26, INEXT, nlev, .true.)
  ! ihpr
  CALL FIELD_INDRAD( MEMBER_MAP, 27, INEXT, nlev+1, .true.)
  ! iaprs
  CALL FIELD_INDRAD( MEMBER_MAP, 28, INEXT, nlev+1, .true.)
  ! ihti
  CALL FIELD_INDRAD( MEMBER_MAP, 29, INEXT, nlev+1, .true.)
  ! iaero
  CALL FIELD_INDRAD( MEMBER_MAP, 30, INEXT, yradiation%rad_config%n_aerosol_types*nlev, .true.)
  ! ifrsod
  CALL FIELD_INDRAD( MEMBER_MAP, 31, INEXT, 1, .true.)
  ! ifrted
  CALL FIELD_INDRAD( MEMBER_MAP, 32, INEXT, 1, .true.)
  ! ifrsodc
  CALL FIELD_INDRAD( MEMBER_MAP, 33, INEXT, 1, .true.)
  ! ifrtedc
  CALL FIELD_INDRAD( MEMBER_MAP, 34, INEXT, 1, .true.)
  ! iemit
  CALL FIELD_INDRAD( MEMBER_MAP, 35, INEXT, 1, .true.)
  ! isudu
  CALL FIELD_INDRAD( MEMBER_MAP, 36, INEXT, 1, .true.)
  ! iuvdf
  CALL FIELD_INDRAD( MEMBER_MAP, 37, INEXT, 1, .true.)
  ! iparf
  CALL FIELD_INDRAD( MEMBER_MAP, 38, INEXT, 1, .true.)
  ! iparcf
  CALL FIELD_INDRAD( MEMBER_MAP, 39, INEXT, 1, .true.)
  ! itincf
  CALL FIELD_INDRAD( MEMBER_MAP, 40, INEXT, 1, .true.)
  ! ifdir
  CALL FIELD_INDRAD( MEMBER_MAP, 41, INEXT, 1, .true.)
  ! ifdif
  CALL FIELD_INDRAD( MEMBER_MAP, 42, INEXT, 1, .true.)
  ! icdir
  CALL FIELD_INDRAD( MEMBER_MAP, 43, INEXT, 1, .true.)
  ! ilwderivative
  CALL FIELD_INDRAD( MEMBER_MAP, 44, INEXT, nlev+1, yradiation%yrerad%lapproxlwupdate)
  ! iswdirectband
  CALL FIELD_INDRAD( MEMBER_MAP, 45, INEXT, yradiation%yrerad%nsw, yradiation%yrerad%lapproxswupdate)
  ! iswdiffuseband
  CALL FIELD_INDRAD( MEMBER_MAP, 46, INEXT, yradiation%yrerad%nsw, yradiation%yrerad%lapproxswupdate)
  ! ifrso
  CALL FIELD_INDRAD( MEMBER_MAP, 47, INEXT, nlev+1, .true.)
  ! iswfc
  CALL FIELD_INDRAD( MEMBER_MAP, 48, INEXT, nlev+1, .true.)
  ! ifrth
  CALL FIELD_INDRAD( MEMBER_MAP, 49, INEXT, nlev+1, .true.)
  ! ilwfc
  CALL FIELD_INDRAD( MEMBER_MAP, 50, INEXT, nlev+1, .true.)
  ! iaer
  CALL FIELD_INDRAD( MEMBER_MAP, 51, INEXT, 6*nlev, .true.)
  ! iico2
  CALL FIELD_INDRAD( MEMBER_MAP, 52, INEXT, nlev, .true.)
  ! iich4
  CALL FIELD_INDRAD( MEMBER_MAP, 53, INEXT, nlev, .true.)
  ! iin2o
  CALL FIELD_INDRAD( MEMBER_MAP, 54, INEXT, nlev, .true.)
  ! ino2
  CALL FIELD_INDRAD( MEMBER_MAP, 55, INEXT, nlev, .true.)
  ! ic11
  CALL FIELD_INDRAD( MEMBER_MAP, 56, INEXT, nlev, .true.)
  ! ic12
  CALL FIELD_INDRAD( MEMBER_MAP, 57, INEXT, nlev, .true.)
  ! ic22
  CALL FIELD_INDRAD( MEMBER_MAP, 58, INEXT, nlev, .true.)
  ! icl4
  CALL FIELD_INDRAD( MEMBER_MAP, 59, INEXT, nlev, .true.)
  ! igix
  CALL FIELD_INDRAD( MEMBER_MAP, 60, INEXT, 1, .false.)
  ! iccno
  CALL FIELD_INDRAD( MEMBER_MAP, 61, INEXT, 1, .true.)
  ! ire_liq
  CALL FIELD_INDRAD( MEMBER_MAP, 62, INEXT, nlev, .true.)
  ! ire_ice
  CALL FIELD_INDRAD( MEMBER_MAP, 63, INEXT, nlev, .true.)
  ! ioverlap
  CALL FIELD_INDRAD( MEMBER_MAP, 64, INEXT, nlev-1, .true.)

  CALL FIELD_NEW(ZRGP_FIELDS%FIELD_WRAPPER, ZRGP_FIELDS%MEMBERS, DATA=ZRGP, MEMBER_MAP=MEMBER_MAP)

END SUBROUTINE IFS_SETUP_FAPI

#endif

subroutine ifs_copy_inputs_to_blocked ( &
  & driver_config, ifs_config, yradiation, ncol, nlev, &
  & single_level, thermodynamics, gas, cloud, aerosol, &
  & sin_latitude, longitude_rad, land_frac, pressure_fl, temperature_fl, &
  & zrgp, thermodynamics_out, iseed)

  use radiation_single_level,   only : single_level_type
  use radiation_thermodynamics, only : thermodynamics_type
  use radiation_gas,            only : gas_type, IMassMixingRatio, &
        &   IH2O, ICO2, IO3, IN2O, ICH4, ICFC11, ICFC12, IHCFC22, ICCL4
  use radiation_cloud,          only : cloud_type
  use radiation_aerosol,        only : aerosol_type
  use ecrad_driver_config,      only : driver_config_type
  use radiation_setup,          only : tradiation

  implicit none

  ! Configuration specific to this driver
  type(driver_config_type), intent(in)     :: driver_config

  type(ifs_config_type), intent(in)     :: ifs_config

  ! Configuration for the radiation scheme, IFS style
  type(tradiation), intent(in)          :: yradiation

  integer, intent(in) :: ncol, nlev         ! Number of columns and levels

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
  integer :: ngpblks, nproma

  integer :: jrl, ibeg, iend, il, ib, ifld, jemiss, jalb, jlev, joff, jaer

  ! Extract some config values
  nproma=driver_config%nblocksize        ! nproma size
  ngpblks=(ncol-1)/nproma+1              ! number of column blocks

  ! Allocate blocked data structure
  allocate(zrgp(nproma,ifs_config%ifldstot,ngpblks))
  if(present(thermodynamics_out)) allocate(thermodynamics_out%pressure_hl(ncol,nlev+1))
  if(present(iseed)) allocate(iseed(nproma,ngpblks))

  ! First touch
  !$OMP PARALLEL DO SCHEDULE(RUNTIME)&
  !$OMP&PRIVATE(IB,IFLD)
  do ib=1,ngpblks
    do ifld=1,ifs_config%ifldstot
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
      zrgp(1:il,ifs_config%iamu0,ib)  =  single_level%cos_sza(ibeg:iend)   ! cosine of solar zenith ang (mu0)

      do jemiss=1,yderad%nlwemiss
        zrgp(1:il,ifs_config%iemiss+jemiss-1,ib)  =  single_level%lw_emissivity(ibeg:iend,jemiss)
      enddo

      zrgp(1:il,ifs_config%its,ib)      = single_level%skin_temperature(ibeg:iend)  ! skin temperature
      zrgp(1:il,ifs_config%islm,ib)     = land_frac(ibeg:iend) ! land-sea mask
      zrgp(1:il,ifs_config%iccnl,ib)    = yderad%rccnlnd ! CCN over land
      zrgp(1:il,ifs_config%iccno,ib)    = yderad%rccnsea ! CCN over sea
      ! zrgp(1:il,ibas,ib)     = ???
      ! zrgp(1:il,itop,ib)     = ???
      zrgp(1:il,ifs_config%igelam,ib)   = longitude_rad(ibeg:iend) ! longitude
      zrgp(1:il,ifs_config%igemu,ib)    = sin_latitude(ibeg:iend) ! sine of latitude
      ! zrgp(1:il,iclon,ib)    = ???
      ! zrgp(1:il,islon,ib)    = ???

      do jalb=1,yderad%nsw
        zrgp(1:il,ifs_config%iald+jalb-1,ib)  =  single_level%sw_albedo(ibeg:iend,jalb)
      enddo

      if (allocated(single_level%sw_albedo_direct)) then
        do jalb=1,yderad%nsw
          zrgp(1:il,ifs_config%ialp+jalb-1,ib)  =  single_level%sw_albedo_direct(ibeg:iend,jalb)
        end do
      else
        do jalb=1,yderad%nsw
          zrgp(1:il,ifs_config%ialp+jalb-1,ib)  =  single_level%sw_albedo(ibeg:iend,jalb)
        end do
      end if

      do jlev=1,nlev
        zrgp(1:il,ifs_config%iti+jlev-1,ib)   = temperature_fl(ibeg:iend,jlev) ! full level temperature
        zrgp(1:il,ifs_config%ipr+jlev-1,ib)   = pressure_fl(ibeg:iend,jlev) ! full level pressure
        ! zrgp(1:il,iqs+jlev-1,ib)   = ???
      enddo

      do jlev=1,nlev
        zrgp(1:il,ifs_config%iwv+jlev-1,ib)   = gas%mixing_ratio(ibeg:iend,jlev,IH2O) ! this is already in MassMixingRatio units
        if (rad_config%do_clouds) then
          zrgp(1:il,ifs_config%iclc+jlev-1,ib)  = cloud%fraction(ibeg:iend,jlev)
          zrgp(1:il,ifs_config%ilwa+jlev-1,ib)  = cloud%q_liq(ibeg:iend,jlev)
          zrgp(1:il,ifs_config%iiwa+jlev-1,ib)  = cloud%q_ice(ibeg:iend,jlev)
        else
          zrgp(1:il,ifs_config%iclc+jlev-1,ib)  = 0._jprb
          zrgp(1:il,ifs_config%ilwa+jlev-1,ib)  = 0._jprb
          zrgp(1:il,ifs_config%iiwa+jlev-1,ib)  = 0._jprb
        endif
        zrgp(1:il,ifs_config%iswa+jlev-1,ib)  = 0._jprb  ! snow
        zrgp(1:il,ifs_config%irwa+jlev-1,ib)  = 0._jprb  ! rain

        ! zrgp(1:il,irra+jlev-1,ib)  = ???
        ! zrgp(1:il,idp+jlev-1,ib)   = ???
        ! zrgp(1:il,ifsd+jlev-1,ib)   = ???
        ! zrgp(1:il,iecpo3+jlev-1,ib) = ???
      enddo

      zrgp(1:il,ifs_config%iaer:ifs_config%iaer+nlev,ib)  =  0._jprb ! old aerosol, not used
      if (yderad%naermacc == 1) then
        joff=ifs_config%iaero
        do jaer=1,rad_config%n_aerosol_types
          do jlev=1,nlev
            zrgp(1:il,joff,ib) = aerosol%mixing_ratio(ibeg:iend,jlev,jaer)
            joff=joff+1
          enddo
        enddo
      endif

      do jlev=1,nlev+1
        ! zrgp(1:il,ihpr+jlev-1,ib)  = ???
        zrgp(1:il,ifs_config%iaprs+jlev-1,ib) = thermodynamics%pressure_hl(ibeg:iend,jlev)
        zrgp(1:il,ifs_config%ihti+jlev-1,ib)  = thermodynamics%temperature_hl(ibeg:iend,jlev)
      enddo

      ! -- by default, globally averaged concentrations (mmr)
      call gas%get(ICO2, IMassMixingRatio, zrgp(1:il,ifs_config%iico2:ifs_config%iico2+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICH4, IMassMixingRatio, zrgp(1:il,ifs_config%iich4:ifs_config%iich4+nlev-1,ib), istartcol=ibeg)
      call gas%get(IN2O, IMassMixingRatio, zrgp(1:il,ifs_config%iin2o:ifs_config%iin2o+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICFC11, IMassMixingRatio, zrgp(1:il,ifs_config%ic11:ifs_config%ic11+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICFC12, IMassMixingRatio, zrgp(1:il,ifs_config%ic12:ifs_config%ic12+nlev-1,ib), istartcol=ibeg)
      call gas%get(IHCFC22,IMassMixingRatio, zrgp(1:il,ifs_config%ic22:ifs_config%ic22+nlev-1,ib), istartcol=ibeg)
      call gas%get(ICCL4,  IMassMixingRatio, zrgp(1:il,ifs_config%icl4:ifs_config%icl4+nlev-1,ib), istartcol=ibeg)
      call gas%get(IO3, IMassMixingRatio, zrgp(1:il,ifs_config%ioz:ifs_config%ioz+nlev-1,ib), istartcol=ibeg)
      ! convert ozone kg/kg to Pa*kg/kg
      ! do jlev=1,nlev
      !   zrgp(1:il,ifs_config%ioz+jlev-1,ib)  = zrgp(1:il,ifs_config%ioz+jlev-1,ib) &
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
          zrgp(1:il,ifs_config%ire_liq+jlev-1,ib) = cloud%re_liq(ibeg:iend,jlev)
          zrgp(1:il,ifs_config%ire_ice+jlev-1,ib) = cloud%re_ice(ibeg:iend,jlev)
        enddo
        do jlev=1,nlev-1
          ! for the love of it, I can't figure this one out. Probably to do with
          ! my crude approach of setting PGEMU?
          zrgp(1:il,ifs_config%ioverlap+jlev-1,ib) = cloud%overlap_param(ibeg:iend,jlev)
        enddo
        if(present(iseed)) iseed(1:il,ib) = single_level%iseed(ibeg:iend)
      else
        do jlev=1,nlev
          ! missing full-level temperature and pressure as well as land-sea-mask
          zrgp(1:il,ifs_config%ire_liq+jlev-1,ib) = 0._jprb
          zrgp(1:il,ifs_config%ire_ice+jlev-1,ib) = 0._jprb
        enddo
        do jlev=1,nlev-1
          zrgp(1:il,ifs_config%ioverlap+jlev-1,ib) = 0._jprb
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
    & driver_config, ifs_config, yradiation, ncol, nlev,&
    & zrgp, flux, flux_sw_direct_normal, flux_uv, flux_par, flux_par_clear,&
    & emissivity_out, flux_diffuse_band, flux_direct_band)
  use ecrad_driver_config,      only : driver_config_type
  use radiation_setup,          only : tradiation
  use radiation_flux,           only : flux_type

  ! Configuration specific to this driver
  type(driver_config_type), intent(in)     :: driver_config

  type(ifs_config_type), intent(in)     :: ifs_config

  ! Configuration for the radiation scheme, IFS style
  type(tradiation), intent(in)          :: yradiation

  integer, intent(in) :: ncol, nlev         ! Number of columns and levels

  ! monolithic IFS data structure passed to radiation scheme
  real(kind=jprb), intent(inout), allocatable :: zrgp(:,:,:)

  ! Derived type containing outputs from the radiation scheme
  type(flux_type), intent(inout)              :: flux

  ! Additional output fluxes as arrays
  real(jprb), dimension(:), intent(inout)     :: flux_sw_direct_normal, flux_uv, flux_par,&
                                                 & flux_par_clear, emissivity_out
  real(jprb), dimension(:,:), intent(inout) :: flux_diffuse_band, flux_direct_band

  ! number of column blocks, block size
  integer :: ngpblks, nproma

  integer :: jrl, ibeg, iend, il, ib, jlev, jg

  ! Extract some config values
  nproma=driver_config%nblocksize        ! nproma size
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
        flux%sw_up(ibeg:iend,jlev) = zrgp(1:il,ifs_config%ifrso+jlev-1,ib)
        flux%lw_up(ibeg:iend,jlev) = zrgp(1:il,ifs_config%ifrth+jlev-1,ib)
        flux%sw_up_clear(ibeg:iend,jlev) = zrgp(1:il,ifs_config%iswfc+jlev-1,ib)
        flux%lw_up_clear(ibeg:iend,jlev) = zrgp(1:il,ifs_config%ilwfc+jlev-1,ib)
        if (yradiation%yrerad%lapproxlwupdate) then
          flux%lw_derivatives(ibeg:iend,jlev) = zrgp(1:il,ifs_config%ilwderivative+jlev-1,ib)
        else
          flux%lw_derivatives(ibeg:iend,jlev) = 0.0_jprb
        endif
      end do
      flux%sw_dn(ibeg:iend,nlev+1) = zrgp(1:il,ifs_config%ifrsod,ib)
      flux%lw_dn(ibeg:iend,nlev+1) = zrgp(1:il,ifs_config%ifrted,ib)
      flux%sw_dn_clear(ibeg:iend,nlev+1) = zrgp(1:il,ifs_config%ifrsodc,ib)
      flux%lw_dn_clear(ibeg:iend,nlev+1) = zrgp(1:il,ifs_config%ifrtedc,ib)
      flux%sw_dn_direct(ibeg:iend,nlev+1) = zrgp(1:il,ifs_config%ifdir,ib)
      flux%sw_dn_direct_clear(ibeg:iend,nlev+1) = zrgp(1:il,ifs_config%icdir,ib)
      flux_sw_direct_normal(ibeg:iend) = zrgp(1:il,ifs_config%isudu,ib)
      flux_uv(ibeg:iend) = zrgp(1:il,ifs_config%iuvdf,ib)
      flux_par(ibeg:iend) = zrgp(1:il,ifs_config%iparf,ib)
      flux_par_clear(ibeg:iend) = zrgp(1:il,ifs_config%iparcf,ib)
      flux%sw_dn(ibeg:iend,1) = zrgp(1:il,ifs_config%itincf,ib)
      emissivity_out(ibeg:iend) = zrgp(1:il,ifs_config%iemit,ib)
      if (yradiation%yrerad%lapproxswupdate) then
        do jg=1,yradiation%yrerad%nsw
          flux_diffuse_band(ibeg:iend,jg) = zrgp(1:il,ifs_config%iswdiffuseband+jg-1,ib)
          flux_direct_band(ibeg:iend,jg) = zrgp(1:il,ifs_config%iswdirectband+jg-1,ib)
        end do
      else
        flux_diffuse_band(ibeg:iend,:) = 0.0_jprb
        flux_direct_band(ibeg:iend,:) = 0.0_jprb
      endif
    end do

    deallocate(zrgp)

end subroutine ifs_copy_fluxes_from_blocked

end module ifs_blocking

! radiation_field_type_mod.F90 - FIELD API backed data structures
! that own the memory used by ecrads derived types
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
! Author:  Johan Ericsson
! Email:   johan.ericsson@ecmwf.int
!


module radiation_field_type_module

  use parkind1,                 only : jprb, jpim ! Working precision, integer type

  ! FIELD API imports
  use field_module, only: field_2rb, field_3rb, field_4rb, field_2im
  use field_factory_module

  ! radiation imports
  use radiation_single_level,     only: single_level_type
  use radiation_thermodynamics,   only: thermodynamics_type
  use radiation_gas,              only: gas_type
  use radiation_gas_constants,    only: NMaxGases
  use radiation_cloud,            only: cloud_type
  use radiation_aerosol,          only: aerosol_type
  use radiation_flux,             only: flux_type

  implicit none

  public

  type single_level_field_type

    integer :: nblocks, ncol, nalbedobands, nemissbands, n_bands_sw
    ! single level type access pointers
    real(jprb), pointer, dimension(:) :: &
         &   cos_sza=>null(), &          ! (ncol) Cosine of solar zenith angle
         &   skin_temperature=>null()    ! (ncol) Skin temperature (K)
    real(jprb), pointer, dimension(:,:) :: &
         &   sw_albedo=>null(), &        ! (ncol,nalbedobands)
         &   sw_albedo_direct=>null()    ! (ncol,nalbedobands)
    real(jprb), pointer, dimension(:,:) :: &
         &   lw_emissivity=>null()       ! (ncol,nemissbands) If
    real(jprb), pointer, dimension(:,:) :: &
         &   lw_emission=>null()         ! (ncol,nemissbands)
    real(jprb), pointer, dimension(:) :: &
         &   spectral_solar_scaling=>null() ! (n_bands_sw)
    integer, pointer, dimension(:) :: iseed=>null() ! (ncol)

    ! single level type fields
    class(field_2rb), pointer :: &
         &   f_cos_sza=>null(), &          ! (ncol,nblocks) Cosine of solar zenith angle
         &   f_skin_temperature=>null()    ! (ncol,nblocks) Skin temperature (K)
    class(field_3rb), pointer :: &
         &   f_sw_albedo=>null(), &        ! (ncol,nalbedobands,nblocks)
         &   f_sw_albedo_direct=>null()    ! (ncol,nalbedobands,nblocks)
    class(field_3rb), pointer :: &
         &   f_lw_emissivity=>null()       ! (ncol,nemissbands,nblocks) If
    class(field_3rb), pointer :: &
         &   f_lw_emission=>null()         ! (ncol,nemissbands,nblocks)
    class(field_2rb), pointer :: &
         &   f_spectral_solar_scaling=>null() ! (n_bands_sw,nblocks)
    class(field_2im), pointer :: f_iseed=>null() ! (ncol,nblocks)

    logical                   :: is_simple_surface = .true.

  contains

    procedure :: init => single_level_field_init
    procedure :: final => single_level_field_final
    procedure :: update_view => single_level_field_update_view
    procedure :: update_single_level => single_level_field_update_single_level

  end type single_level_field_type

  type thermodynamics_field_type

     real(jprb), pointer, dimension(:,:) :: &
          &  pressure_hl=>null(), &   ! (ncol,nlev+1) pressure (Pa)
          &  temperature_hl=>null()   ! (ncol,nlev+1) temperature (K)
     real(jprb), pointer, dimension(:,:) :: &
          &  h2o_sat_liq=>null() ! (ncol,nlev) specific humidity at liquid
                         ! saturation (kg/kg)
     class(field_3rb), pointer :: &
          &  f_pressure_hl=>null(), &   ! (ncol,nlev+1,nblocks) pressure (Pa)
          &  f_temperature_hl=>null()   ! (ncol,nlev+1,nblocks) temperature (K)
     class(field_3rb), pointer :: &
          &  f_h2o_sat_liq=>null() ! (ncol,nlev,nblocks) specific humidity at liquid
                         ! saturation (kg/kg)

   contains

    procedure :: init => thermodynamics_field_init
    procedure :: final => thermodynamics_field_final
    procedure :: update_view => thermodynamics_field_update_view
    procedure :: update_thermodynamics => thermodynamics_field_update_thermodynamics

  end type thermodynamics_field_type

  type gas_field_type
    integer :: ncol           = 0 ! Number of columns in mixing_ratio
    integer :: nlev           = 0 ! Number of levels  in mixing_ratio

    real(jprb), pointer, dimension(:,:,:) :: mixing_ratio=>null()

    class(field_4rb), pointer :: f_mixing_ratio=>null() ! (ncol, nlev, NMaxGases, nblks)

  contains

    procedure :: init => gas_field_init
    procedure :: final => gas_field_final
    procedure :: update_view => gas_field_update_view
    procedure :: update_gas => gas_field_update_gas

  end type gas_field_type

  type cloud_field_type
    integer                                   :: ntype = 0
    logical                                   :: ntype_present = .false.
    real(jprb), pointer, dimension(:,:,:) :: & ! (ncol,nlev,ntype)
         &  mixing_ratio=>null(), &
         &  effective_radius=>null()

    real(jprb), pointer, dimension(:,:) :: & ! (ncol,nlev)
         &  q_liq=>null(),  q_ice=>null(),  &
         &  re_liq=>null(), re_ice=>null()
    real(jprb), pointer, dimension(:,:) :: fraction=>null()
    real(jprb), pointer, dimension(:,:) :: fractional_std=>null() ! (ncol,nlev)
    real(jprb), pointer, dimension(:,:) :: & ! (ncol,nlev).
         &  inv_cloud_effective_size=>null() ! (ncol,nlev)
    real(jprb), pointer, dimension(:,:) :: &
         &  inv_inhom_effective_size=>null()
    real(jprb), pointer, dimension(:,:) :: overlap_param=>null() ! (ncol,nlev-1)

    class(field_4rb), pointer :: & ! (ncol,nlev,ntype,nblocks)
         &  f_mixing_ratio=>null(), &
         &  f_effective_radius=>null()

    class(field_3rb), pointer :: f_fraction=>null()
    class(field_3rb), pointer :: f_fractional_std=>null() ! (ncol,nlev,nblocks)
    class(field_3rb), pointer :: & ! (ncol,nlev,nblocks).
         &  f_inv_cloud_effective_size=>null() ! (ncol,nlev,nblocks)
    class(field_3rb), pointer :: &
         &  f_inv_inhom_effective_size=>null()
    class(field_3rb), pointer :: f_overlap_param=>null() ! (ncol,nlev-1,nblocks)

  contains

    procedure :: init => cloud_field_init
    procedure :: final => cloud_field_final
    procedure :: update_view => cloud_field_update_view
    procedure :: update_cloud => cloud_field_update_cloud

  end type cloud_field_type

  type aerosol_field_type
    class(field_4rb), pointer :: & ! (ncol,istartlev:iendlev,config%n_aerosol_types, nblocks)
          &  f_mixing_ratio=>null()
    real(jprb), pointer, dimension(:,:,:) :: & ! (ncol,istartlev:iendlev,config%n_aerosol_types)
          &  mixing_ratio=>null()

     integer :: istartlev, iendlev
     logical :: is_direct = .false.

  contains
    procedure :: init           => aerosol_field_init
    procedure :: final          => aerosol_field_final
    procedure :: update_view    => aerosol_field_update_view
    procedure :: update_aerosol => aerosol_field_update_aerosol
  end type aerosol_field_type


  type flux_field_type
    real(jprb), pointer, dimension(:,:) :: &  ! (ncol,nlev+1)
         &  lw_up=>null(), lw_dn=>null(), &
         &  sw_up=>null(), sw_dn=>null(), &
         &  sw_dn_direct=>null(), &
         &  lw_up_clear=>null(), lw_dn_clear=>null(), &
         &  sw_up_clear=>null(), sw_dn_clear=>null(), &
         &  sw_dn_direct_clear=>null()

    class(field_3rb), pointer :: &  ! (ncol,nlev+1,nblocks)
         &  f_lw_up=>null(), f_lw_dn=>null(), &
         &  f_sw_up=>null(), f_sw_dn=>null(), &
         &  f_sw_dn_direct=>null(), &
         &  f_lw_up_clear=>null(), f_lw_dn_clear=>null(), &
         &  f_sw_up_clear=>null(), f_sw_dn_clear=>null(), &
         &  f_sw_dn_direct_clear=>null()

    
    real(jprb), pointer, dimension(:,:,:) :: & ! (nband,ncol,nlev+1)
         &  lw_up_band=>null(), lw_dn_band=>null(), &
         &  sw_up_band=>null(), sw_dn_band=>null(), &
         &  sw_dn_direct_band=>null(), &
         &  lw_up_clear_band=>null(), lw_dn_clear_band=>null(), &
         &  sw_up_clear_band=>null(), sw_dn_clear_band=>null(), &
         &  sw_dn_direct_clear_band=>null()

    class(field_4rb), pointer :: & ! (nband,ncol,nlev+1,nblocks)
         &  f_lw_up_band=>null(), f_lw_dn_band=>null(), &
         &  f_sw_up_band=>null(), f_sw_dn_band=>null(), &
         &  f_sw_dn_direct_band=>null(), &
         &  f_lw_up_clear_band=>null(), f_lw_dn_clear_band=>null(), &
         &  f_sw_up_clear_band=>null(), f_sw_dn_clear_band=>null(), &
         &  f_sw_dn_direct_clear_band=>null()


    real(jprb), pointer, dimension(:,:) :: &  ! (ng,ncol)
         &  lw_dn_surf_g=>null(), lw_dn_surf_clear_g=>null(), &
         &  sw_dn_diffuse_surf_g=>null(), sw_dn_direct_surf_g=>null(), &
         &  sw_dn_diffuse_surf_clear_g=>null(), sw_dn_direct_surf_clear_g=>null()

    class(field_3rb), pointer :: &  ! (ng,ncol,nblocks)
         &  f_lw_dn_surf_g=>null(), f_lw_dn_surf_clear_g=>null(), &
         &  f_sw_dn_diffuse_surf_g=>null(), f_sw_dn_direct_surf_g=>null(), &
         &  f_sw_dn_diffuse_surf_clear_g=>null(), f_sw_dn_direct_surf_clear_g=>null()

    
    real(jprb), pointer, dimension(:,:) :: &  !(ng,ncol)
         &  lw_up_toa_g=>null(), lw_up_toa_clear_g=>null(), &
         &  sw_dn_toa_g=>null(), sw_up_toa_g=>null(), sw_up_toa_clear_g=>null()
    
    class(field_3rb), pointer :: &  !(ng,ncol,nblocks)
         &  f_lw_up_toa_g=>null(), f_lw_up_toa_clear_g=>null(), &
         &  f_sw_dn_toa_g=>null(), f_sw_up_toa_g=>null(), f_sw_up_toa_clear_g=>null()


    real(jprb), pointer, dimension(:,:) :: &  ! (nband,ncol)
         &  sw_dn_surf_band=>null(), sw_dn_direct_surf_band=>null(), &
         &  sw_dn_surf_clear_band=>null(), sw_dn_direct_surf_clear_band=>null()
    
    class(field_3rb), pointer :: &  ! (nband,ncol,nblocks)
         &  f_sw_dn_surf_band=>null(), f_sw_dn_direct_surf_band=>null(), &
         &  f_sw_dn_surf_clear_band=>null(), f_sw_dn_direct_surf_clear_band=>null()


    real(jprb), pointer, dimension(:,:) :: &  ! (nband,ncol)
         &  lw_up_toa_band=>null(), lw_up_toa_clear_band=>null(), &
         &  sw_dn_toa_band=>null(), sw_up_toa_band=>null(), sw_up_toa_clear_band=>null()

    class(field_3rb), pointer :: &
         &  f_lw_up_toa_band=>null(), f_lw_up_toa_clear_band=>null(), &
         &  f_sw_dn_toa_band=>null(), f_sw_up_toa_band=>null(), f_sw_up_toa_clear_band=>null()


    real(jprb), pointer, dimension(:,:) :: &
         &  lw_dn_surf_canopy=>null(), &
         &  sw_dn_diffuse_surf_canopy=>null(), sw_dn_direct_surf_canopy=>null()

    class(field_3rb), pointer :: &
         &  f_lw_dn_surf_canopy=>null(), &
         &  f_sw_dn_diffuse_surf_canopy=>null(), f_sw_dn_direct_surf_canopy=>null()


    real(jprb), pointer, dimension(:) :: &
         &  cloud_cover_lw=>null(), cloud_cover_sw=>null()

    class(field_2rb), pointer :: &
         &  f_cloud_cover_lw=>null(), f_cloud_cover_sw=>null()


    real(jprb), pointer, dimension(:,:) :: &  ! (ncol,nlev+1)
          &  lw_derivatives=>null()
    class(field_3rb), pointer :: &  ! (ncol,nlev+1,nblocks)
          &  f_lw_derivatives=>null()
   
  contains
    procedure :: init           => flux_field_init
    procedure :: final          => flux_field_final
    procedure :: update_view    => flux_field_update_view
    procedure :: update_flux => flux_field_update_flux
 
  end type flux_field_type


contains


!-----------------------------------------------------------------------
! single_level_field_type procedures

  !---------------------------------------------------------------------
  ! Initialise single_field_type
  subroutine single_level_field_init(this, nblocks, ncol, nalbedobands, nemisbands, &
       &                           use_sw_albedo_direct, is_simple_surface)

    use yomhook, only : lhook, dr_hook, jphook

    class(single_level_field_type), intent(inout) :: this
    integer,                  intent(in)    :: nblocks, ncol, nalbedobands, nemisbands
    logical,        optional, intent(in)    :: use_sw_albedo_direct
    logical,        optional, intent(in)    :: is_simple_surface

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type_module:single_level_field_init',0,hook_handle)

    if (present(is_simple_surface)) this%is_simple_surface = is_simple_surface

    call field_new(this%f_cos_sza, ubounds=[ncol, nblocks], persistent=.true.)

    if (this%is_simple_surface) then
      call field_new(this%f_skin_temperature, ubounds=[ncol, nblocks], persistent=.true.)
    else
      call field_new(this%f_lw_emission, ubounds=[ncol, nemisbands, nblocks], persistent=.true.)
    end if
    call field_new(this%f_lw_emissivity, ubounds=[ncol, nemisbands, nblocks], persistent=.true.)

    call field_new(this%f_sw_albedo, ubounds=[ncol, nalbedobands, nblocks], persistent=.true.)

    if (present(use_sw_albedo_direct)) then
      if (use_sw_albedo_direct) then
        call field_new(this%f_sw_albedo_direct, ubounds=[ncol, nalbedobands, nblocks], persistent=.true.)
      end if
    end if

    call field_new(this%f_iseed, ubounds=[ncol, nblocks], persistent=.true.)

    if (lhook) call dr_hook('radiation_field_type_module:single_level_field_init',1,hook_handle)
  end subroutine single_level_field_init


  !---------------------------------------------------------------------
  ! Finalise single_level_field_type
  subroutine single_level_field_final(this)

    use yomhook, only : lhook, dr_hook, jphook

    class(single_level_field_type), intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:single_level_field_final',0,hook_handle)

    if (associated(this%f_cos_sza)) then
      call field_delete(this%f_cos_sza)
    end if
    this%f_cos_sza => null()
    this%cos_sza => null()

    if (associated(this%f_skin_temperature)) then
      call field_delete(this%f_skin_temperature)
    end if
    this%f_skin_temperature => null()
    this%skin_temperature => null()

    if (associated(this%f_sw_albedo)) then
      call field_delete(this%f_sw_albedo)
    end if
    this%f_sw_albedo => null()
    this%sw_albedo => null()

    if (associated(this%f_sw_albedo_direct)) then
      call field_delete(this%f_sw_albedo_direct)
    end if
    this%f_sw_albedo_direct => null()
    this%sw_albedo_direct => null()

    if (associated(this%f_lw_emissivity)) then
      call field_delete(this%f_lw_emissivity)
    end if
    this%f_lw_emissivity => null()
    this%lw_emissivity => null()

    if (associated(this%f_lw_emission)) then
      call field_delete(this%f_lw_emission)
    end if
    this%f_lw_emission => null()
    this%lw_emission => null()

    if (associated(this%f_spectral_solar_scaling)) then
      call field_delete(this%f_spectral_solar_scaling)
    end if
    this%f_spectral_solar_scaling => null()
    this%spectral_solar_scaling => null()

    if (associated(this%f_iseed)) then
      call field_delete(this%f_iseed)
    end if
    this%f_iseed => null()
    this%iseed => null()

    if (lhook) call dr_hook('radiation_field_type:single_level_field_final',1,hook_handle)

  end subroutine single_level_field_final

  !---------------------------------------------------------------------
  ! Update view pointers with block level views
  subroutine single_level_field_update_view(this, block_index)

    use yomhook, only : lhook, dr_hook, jphook

    class(single_level_field_type), intent(inout) :: this
    integer,                        intent(in)    :: block_index
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:single_level_field_update_view',0,hook_handle)

    if (associated(this%f_cos_sza)) then
      this%cos_sza => this%f_cos_sza%get_view(block_index)
    end if

    if (associated(this%f_skin_temperature)) then
      this%skin_temperature => this%f_skin_temperature%get_view(block_index)
    end if

    if (associated(this%f_sw_albedo)) then
      this%sw_albedo => this%f_sw_albedo%get_view(block_index)
    end if

    if (associated(this%f_sw_albedo_direct)) then
      this%sw_albedo_direct => this%f_sw_albedo_direct%get_view(block_index)
    end if

    if (associated(this%f_lw_emissivity)) then
      this%lw_emissivity => this%f_lw_emissivity%get_view(block_index)
    end if

    if (associated(this%f_lw_emission)) then
      this%lw_emission => this%f_lw_emission%get_view(block_index)
    end if

    if (associated(this%f_spectral_solar_scaling)) then
      this%spectral_solar_scaling => this%f_spectral_solar_scaling%get_view(block_index)
    end if

    if (associated(this%f_iseed)) then
      this%iseed => this%f_iseed%get_view(block_index)
    end if

    if (lhook) call dr_hook('radiation_field_type:single_level_field_update_view',1,hook_handle)

  end subroutine single_level_field_update_view

  !---------------------------------------------------------------------
  ! Update single_level_type with the view pointers of this object
  subroutine single_level_field_update_single_level(this, single_level)

    use yomhook, only : lhook, dr_hook, jphook

    class(single_level_field_type), intent(inout) :: this
    class(single_level_type),       intent(inout) :: single_level

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:single_level_field_update_single_level',0,hook_handle)

    if (associated(this%cos_sza)) then
      single_level%cos_sza => this%cos_sza
    end if

    if (associated(this%skin_temperature)) then
      single_level%skin_temperature => this%skin_temperature
    end if

    if (associated(this%sw_albedo)) then
      single_level%sw_albedo => this%sw_albedo
    end if

    if (associated(this%sw_albedo_direct)) then
      single_level%sw_albedo_direct => this%sw_albedo_direct
    end if

    if (associated(this%lw_emissivity)) then
      single_level%lw_emissivity => this%lw_emissivity
    end if

    if (associated(this%lw_emission)) then
      single_level%lw_emission => this%lw_emission
    end if

    if (associated(this%spectral_solar_scaling)) then
      single_level%spectral_solar_scaling => this%spectral_solar_scaling
    end if

    single_level%is_simple_surface = this%is_simple_surface

    if (associated(this%f_iseed)) then
      single_level%iseed => this%iseed
    end if

    if (lhook) call dr_hook('radiation_field_type:single_level_field_update_single_level',1,hook_handle)

  end subroutine single_level_field_update_single_level


!-----------------------------------------------------------------------
! thermodynamics_field_type procedures

  !---------------------------------------------------------------------
  ! Allocate variables with specified dimensions
  subroutine thermodynamics_field_init(this, nblocks, ncol, nlev, use_h2o_sat)

    use yomhook,  only : lhook, dr_hook, jphook

    class(thermodynamics_field_type), intent(inout) :: this
    integer, intent(in)           :: nblocks  ! Number of blocks
    integer, intent(in)           :: ncol  ! Number of columns
    integer, intent(in)           :: nlev  ! Number of levels
    logical, intent(in), optional :: use_h2o_sat ! Allocate h2o_sat_liq?

    logical :: use_h2o_sat_local

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:thermodynamics_field_init',0,hook_handle)

    call field_new(this%f_pressure_hl, ubounds=[ncol,nlev+1, nblocks], persistent=.true.)
    call field_new(this%f_temperature_hl, ubounds=[ncol,nlev+1, nblocks], persistent=.true.)

    use_h2o_sat_local = .false.
    if (present(use_h2o_sat)) then
      use_h2o_sat_local = use_h2o_sat
    end if

    if (use_h2o_sat_local) then
      call field_new(this%f_h2o_sat_liq, ubounds=[ncol,nlev, nblocks], persistent=.true.)
    end if

    if (lhook) call dr_hook('radiation_field_type:thermodynamics_field_init',1,hook_handle)

  end subroutine thermodynamics_field_init

  !---------------------------------------------------------------------
  ! Deallocate variables
  subroutine thermodynamics_field_final(this)

    use yomhook,  only : lhook, dr_hook, jphook

    class(thermodynamics_field_type), intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:thermodynamics_field_final',0,hook_handle)

    if (associated(this%f_pressure_hl)) then
      call field_delete(this%f_pressure_hl)
    end if
    this%f_pressure_hl=>null()
    this%pressure_hl=>null()
    if (associated(this%f_temperature_hl)) then
      call field_delete(this%f_temperature_hl)
    end if
    this%f_temperature_hl=>null()
    this%temperature_hl=>null()
    if (associated(this%f_h2o_sat_liq)) then
      call field_delete(this%f_h2o_sat_liq)
    end if
    this%f_h2o_sat_liq=>null()
    this%h2o_sat_liq=>null()

    if (lhook) call dr_hook('radiation_field_type:thermodynamics_field_final',1,hook_handle)

  end subroutine thermodynamics_field_final

  !---------------------------------------------------------------------
  ! Update view pointers with block level views
  subroutine thermodynamics_field_update_view(this, block_index)

    use yomhook,  only : lhook, dr_hook, jphook

    class(thermodynamics_field_type), intent(inout) :: this
    integer,                          intent(in)    :: block_index

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:thermodynamics_field_update_view',0,hook_handle)

    if (associated(this%f_pressure_hl)) then
      this%pressure_hl => this%f_pressure_hl%get_view(block_index)
    end if
    if (associated(this%f_temperature_hl)) then
      this%temperature_hl => this%f_temperature_hl%get_view(block_index)
    end if
    if (associated(this%f_h2o_sat_liq)) then
      this%h2o_sat_liq => this%f_h2o_sat_liq%get_view(block_index)
    end if

    if (lhook) call dr_hook('radiation_field_type:thermodynamics_field_update_view',1,hook_handle)

  end subroutine thermodynamics_field_update_view

  !---------------------------------------------------------------------
  ! Update thermodynamics type with the view pointers of this object
  subroutine thermodynamics_field_update_thermodynamics(this, thermodynamics)

    use yomhook,  only : lhook, dr_hook, jphook

    class(thermodynamics_field_type), intent(inout) :: this
    class(thermodynamics_type),       intent(inout) :: thermodynamics

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:thermodynamics_field_update_thermodynamics',0,hook_handle)

    if (associated(this%f_pressure_hl)) then
      thermodynamics%pressure_hl => this%pressure_hl
    end if
    if (associated(this%f_temperature_hl)) then
      thermodynamics%temperature_hl => this%temperature_hl
    end if
    if (associated(this%f_h2o_sat_liq)) then
      thermodynamics%h2o_sat_liq => this%h2o_sat_liq
    end if

    if (lhook) call dr_hook('radiation_field_type:thermodynamics_field_update_thermodynamics',1,hook_handle)

  end subroutine thermodynamics_field_update_thermodynamics


!-----------------------------------------------------------------------
! Gas field type procedures

  !---------------------------------------------------------------------
  ! Initialise gas field type
  subroutine gas_field_init(this, nblocks, ncol, nlev)

    use yomhook, only : lhook, dr_hook, jphook

    class(gas_field_type), intent(inout) :: this
    integer,         intent(in)    :: nblocks, ncol, nlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:gas_field_init',0,hook_handle)

    call field_new(this%f_mixing_ratio, ubounds=[ncol, nlev, NMaxGases, nblocks], persistent=.true., init_value=0.0_jprb)

    this%ncol = ncol
    this%nlev = nlev

    if (lhook) call dr_hook('radiation_field_type:gas_field_init',1,hook_handle)

  end subroutine gas_field_init

  !---------------------------------------------------------------------
  ! Finalise gas field type
  subroutine gas_field_final(this)

    use yomhook, only : lhook, dr_hook, jphook

    class(gas_field_type),  intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:gas_field_final',0,hook_handle)

    call field_delete(this%f_mixing_ratio)
    this%f_mixing_ratio => null()

    this%ncol = 0
    this%nlev = 0

    if (lhook) call dr_hook('radiation_field_type:gas_field_final',1,hook_handle)

  end subroutine gas_field_final

  !---------------------------------------------------------------------
  ! Update view pointer
  subroutine gas_field_update_view(this, block_index)

    use yomhook, only : lhook, dr_hook, jphook

    class(gas_field_type),  intent(inout) :: this
    integer,                intent(in)    :: block_index

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:gas_field_update_view',0,hook_handle)

    if (associated(this%f_mixing_ratio)) then
      this%mixing_ratio => this%f_mixing_ratio%get_view(block_index)
    end if

    if (lhook) call dr_hook('radiation_field_type:gas_field_update_view',1,hook_handle)

  end subroutine gas_field_update_view

  !---------------------------------------------------------------------
  ! Update gas field type
  subroutine gas_field_update_gas(this, gas)

    use yomhook, only : lhook, dr_hook, jphook

    class(gas_field_type),  intent(inout) :: this
    class(gas_type),        intent(inout) :: gas

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:gas_field_update_gas',0,hook_handle)
    if (associated(this%mixing_ratio)) gas%mixing_ratio => this%mixing_ratio
    gas%ncol = this%ncol
    gas%nlev = this%nlev

    if (lhook) call dr_hook('radiation_field_type:gas_field_update_gas',1,hook_handle)

  end subroutine gas_field_update_gas


!-----------------------------------------------------------------------
! cloud_field_type procedures

  !---------------------------------------------------------------------
  ! Initialise cloud_field_type
  subroutine cloud_field_init(this, nblocks, ncol, nlev, ntype, use_inhom_effective_size, frac_std)

    use yomhook,     only : lhook, dr_hook, jphook

    class(cloud_field_type), intent(inout), target :: this
    integer, intent(in)              :: nblocks   ! Total number of blocks
    integer, intent(in)              :: ncol   ! Number of columns
    integer, intent(in)              :: nlev   ! Number of levels
    integer, intent(in), optional    :: ntype
    logical, intent(in), optional    :: use_inhom_effective_size
    real(jprb), intent(in), optional :: frac_std ! Fractional std

    real(jprb)   :: frac_std_local = 1.0_jprb
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:cloud_field_init',0,hook_handle)

    if (present(ntype)) then
      this%ntype = ntype
      this%ntype_present = .true.
    else
      this%ntype = 2
      this%ntype_present = .false.
    end if

    if (present(frac_std)) frac_std_local = frac_std

    call field_new(this%f_mixing_ratio, ubounds=[ncol,nlev,this%ntype, nblocks], persistent=.true.)
    call field_new(this%f_effective_radius, ubounds=[ncol,nlev,this%ntype, nblocks], persistent=.true.)

    call field_new(this%f_fraction, ubounds=[ncol,nlev, nblocks], persistent=.true.)
    call field_new(this%f_overlap_param, ubounds=[ncol,nlev-1, nblocks], persistent=.true.)
    call field_new(this%f_fractional_std, ubounds=[ncol,nlev, nblocks], persistent=.true., init_value=frac_std_local)
    call field_new(this%f_inv_cloud_effective_size, ubounds=[ncol,nlev, nblocks], persistent=.true.)
    call field_new(this%f_inv_inhom_effective_size, ubounds=[ncol,nlev, nblocks], persistent=.true.)

    if (lhook) call dr_hook('radiation_radiation_field_type:cloud_field_init',1,hook_handle)

  end subroutine cloud_field_init

  !---------------------------------------------------------------------
  ! cloud_field_type finalisation
  subroutine cloud_field_final(this)

    use yomhook,     only : lhook, dr_hook, jphook

    class(cloud_field_type), intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:cloud_field_final',0,hook_handle)

    this%q_liq => null()
    this%q_ice => null()
    this%re_liq => null()
    this%re_ice => null()

    if (associated(this%f_mixing_ratio)) then
      call field_delete(this%f_mixing_ratio)
    end if
    this%f_mixing_ratio=>null()
    this%mixing_ratio=>null()
    if (associated(this%f_effective_radius)) then
      call field_delete(this%f_effective_radius)
    end if
    this%f_effective_radius=>null()
    this%effective_radius=>null()
    if (associated(this%f_fraction)) then
      call field_delete(this%f_fraction)
    end if
    this%f_fraction=>null()
    this%fraction=>null()
    if (associated(this%f_overlap_param)) then
      call field_delete(this%f_overlap_param)
    end if
    this%f_overlap_param=>null()
    this%overlap_param=>null()
    if (associated(this%f_fractional_std)) then
      call field_delete(this%f_fractional_std)
    end if
    this%f_fractional_std=>null()
    this%fractional_std=>null()
    if (associated(this%f_inv_cloud_effective_size)) then
      call field_delete(this%f_inv_cloud_effective_size)
    end if
    this%f_inv_cloud_effective_size=>null()
    this%inv_cloud_effective_size=>null()
    if (associated(this%f_inv_inhom_effective_size)) then
      call field_delete(this%f_inv_inhom_effective_size)
    end if
    this%f_inv_inhom_effective_size=>null()
    this%inv_inhom_effective_size=>null()

    this%ntype_present = .false.
    this%ntype = 0

    this%q_liq=>null()
    this%q_ice=>null()
    this%re_liq=>null()
    this%re_ice=>null()

    if (lhook) call dr_hook('radiation_field_type:cloud_field_final',1,hook_handle)

  end subroutine cloud_field_final

  !---------------------------------------------------------------------
  ! Update view pointers of cloud_field_type
  subroutine cloud_field_update_view(this, block_index)

    use yomhook,     only : lhook, dr_hook, jphook

    class(cloud_field_type), intent(inout)  :: this
    integer, intent(in)                     :: block_index

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:cloud_field_update_view',0,hook_handle)


    if (associated(this%f_mixing_ratio)) then
      this%mixing_ratio => this%f_mixing_ratio%get_view(block_index)
    end if
    if (associated(this%f_mixing_ratio)) then
      this%mixing_ratio => this%f_mixing_ratio%get_view(block_index)
    end if
    if (associated(this%f_effective_radius)) then
      this%effective_radius => this%f_effective_radius%get_view(block_index)
    end if

    if (associated(this%f_fraction)) then
      this%fraction => this%f_fraction%get_view(block_index)
    end if
    if (associated(this%f_overlap_param)) then
      this%overlap_param => this%f_overlap_param%get_view(block_index)
    end if
    if (associated(this%f_fractional_std)) then
      this%fractional_std => this%f_fractional_std%get_view(block_index)
    end if
    if (associated(this%f_inv_cloud_effective_size)) then
      this%inv_cloud_effective_size => this%f_inv_cloud_effective_size%get_view(block_index)
    end if

    if (associated(this%f_inv_inhom_effective_size)) then
      this%inv_inhom_effective_size => this%f_inv_inhom_effective_size%get_view(block_index)
    end if

    if (.not. this%ntype_present) then
      ! Older interface in which only liquid and ice are supported
      if (associated(this%mixing_ratio)) then
        this%q_liq  => this%mixing_ratio(:,:,1)
        this%q_ice  => this%mixing_ratio(:,:,2)
      end if
      if (associated(this%effective_radius)) then
        this%re_liq => this%effective_radius(:,:,1)
        this%re_ice => this%effective_radius(:,:,2)
      end if
    end if

    if (lhook) call dr_hook('radiation_radiation_field_type:cloud_field_update_view',1,hook_handle)

  end subroutine cloud_field_update_view

  !---------------------------------------------------------------------
  ! Update cloud pointers of cloud_field_type
  subroutine cloud_field_update_cloud(this, ylcloud)

    use yomhook,     only : lhook, dr_hook, jphook

    class(cloud_field_type), intent(inout)  :: this
    class(cloud_type), intent(inout)        :: ylcloud

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:cloud_field_update_cloud',0,hook_handle)


    if (associated(this%mixing_ratio)) then
      ylcloud%mixing_ratio => this%mixing_ratio
    end if
    if (associated(this%effective_radius)) then
      ylcloud%effective_radius => this%effective_radius
    end if

    if (associated(this%fraction)) then
      ylcloud%fraction => this%fraction
    end if
    if (associated(this%overlap_param)) then
      ylcloud%overlap_param => this%overlap_param
    end if
    if (associated(this%fractional_std)) then
      ylcloud%fractional_std => this%fractional_std
    end if
    if (associated(this%inv_cloud_effective_size)) then
      ylcloud%inv_cloud_effective_size => this%inv_cloud_effective_size
    end if

    if (associated(this%inv_inhom_effective_size)) then
      ylcloud%inv_inhom_effective_size => this%inv_inhom_effective_size
    end if

    if (.not. this%ntype_present) then
      if (associated(this%mixing_ratio)) then
        ylcloud%q_liq  => this%q_liq
        ylcloud%q_ice  => this%q_ice
      end if
      if (associated(this%effective_radius)) then
        ylcloud%re_liq => this%re_liq
        ylcloud%re_ice => this%re_ice
      end if
    end if

    if (lhook) call dr_hook('radiation_radiation_field_type:cloud_field_update_cloud',1,hook_handle)

  end subroutine cloud_field_update_cloud


!-----------------------------------------------------------------------
! aerosol_field_type procedures

  !---------------------------------------------------------------------
  ! Initialise aerosol_field_type
  subroutine aerosol_field_init(this, nblocks, ncol, istartlev, iendlev, ntype)

    use yomhook,     only : lhook, dr_hook, jphook

    class(aerosol_field_type), intent(inout)  :: this
    integer, intent(in)                       :: nblocks  ! Total number of blocks
    integer, intent(in)                       :: ncol  ! Number of columns
    integer, intent(in)                       :: istartlev, iendlev ! Level range
    integer, intent(in)                       :: ntype ! Number of aerosol types

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:aerosol_field_init',0,hook_handle)

    this%is_direct = .false.
    this%istartlev = istartlev
    this%iendlev   = iendlev

    call field_new(this%f_mixing_ratio, lbounds=[1,istartlev,1,1], &
                   &                    ubounds=[ncol,iendlev,ntype,nblocks], persistent=.true.)

    if (lhook) call dr_hook('radiation_radiation_field_type:aerosol_field_init',1,hook_handle)

  end subroutine aerosol_field_init

  !---------------------------------------------------------------------
  ! aerosol_field_type finalisation
  subroutine aerosol_field_final(this)

    use yomhook,     only : lhook, dr_hook, jphook

    class(aerosol_field_type), intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:aerosol_field_final',0,hook_handle)
    if (associated(this%f_mixing_ratio)) then
      call field_delete(this%f_mixing_ratio)
    end if
    this%f_mixing_ratio=>null()
    this%mixing_ratio=>null()

    if (lhook) call dr_hook('radiation_field_type:aerosol_field_final',1,hook_handle)

  end subroutine aerosol_field_final

  !---------------------------------------------------------------------
  ! Update view pointers of aerosol_field_type
  subroutine aerosol_field_update_view(this, block_index)

    use yomhook,     only : lhook, dr_hook, jphook

    class(aerosol_field_type), intent(inout)  :: this
    integer, intent(in)                     :: block_index

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:aerosol_field_update_view',0,hook_handle)

    if (associated(this%f_mixing_ratio)) then
      this%mixing_ratio => this%f_mixing_ratio%get_view(block_index)
    end if

    if (lhook) call dr_hook('radiation_radiation_field_type:aerosol_field_update_view',1,hook_handle)

  end subroutine aerosol_field_update_view

  !---------------------------------------------------------------------
  ! Update aerosol pointers of aerosol_field_type
  subroutine aerosol_field_update_aerosol(this, ylaerosol)

    use yomhook,     only : lhook, dr_hook, jphook

    class(aerosol_field_type), intent(inout)  :: this
    class(aerosol_type), intent(inout)        :: ylaerosol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:aerosol_field_update_aerosol',0,hook_handle)

    if (associated(this%mixing_ratio)) then
      ylaerosol%mixing_ratio => this%mixing_ratio
    end if

    ylaerosol%is_direct=this%is_direct
    ylaerosol%istartlev=this%istartlev
    ylaerosol%iendlev=this%iendlev

    if (lhook) call dr_hook('radiation_radiation_field_type:aerosol_field_update_aerosol',1,hook_handle)

  end subroutine aerosol_field_update_aerosol


!-----------------------------------------------------------------------
! flux_field_type procedures

  !---------------------------------------------------------------------
  ! Initialise flux_field_type
  subroutine flux_field_init(this, nblocks, config, istartcol, iendcol, nlev)

    use yomhook,          only : lhook, dr_hook, jphook
    use radiation_io,     only : nulerr, radiation_abort
    use radiation_config, only: config_type

    class(flux_field_type), intent(inout) :: this
    type(config_type), intent(in)         :: config
    integer, intent(in)                   :: nblocks  ! Total number of blocks
    integer, intent(in)                   :: istartcol, iendcol, nlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:flux_field_init',0,hook_handle)

    ! Allocate longwave arrays
    if (config%do_lw) then
      call field_new(this%f_lw_up,lbounds=[istartcol,1,1], ubounds=[iendcol,nlev+1,nblocks], persistent=.true.)
      call field_new(this%f_lw_dn,lbounds=[istartcol,1,1], ubounds=[iendcol,nlev+1,nblocks], persistent=.true.)
      if (config%do_clear) then
        call field_new(this%f_lw_up_clear,lbounds=[istartcol,1,1], &
            &                             ubounds=[iendcol,nlev+1,nblocks], persistent=.true.)
        call field_new(this%f_lw_dn_clear,lbounds=[istartcol,1,1], &
            &                             ubounds=[iendcol,nlev+1,nblocks], persistent=.true.)
      end if

      if (config%do_save_spectral_flux) then
        if (config%n_spec_lw == 0) then
          write(nulerr,'(a)') '*** Error: number of LW spectral points to save not yet defined ' &
               & // 'so cannot allocate spectral flux arrays'
          call radiation_abort()
        end if

        call field_new(this%f_lw_up_band,lbounds=[1,istartcol,1,1], &
            &                            ubounds=[config%n_spec_lw,iendcol,nlev+1,nblocks], persistent=.true.)
        call field_new(this%f_lw_dn_band,lbounds=[1,istartcol,1,1], &
            &                            ubounds=[config%n_spec_lw,iendcol,nlev+1,nblocks], persistent=.true.)
        if (config%do_clear) then
          call field_new(this%f_lw_up_clear_band,lbounds=[1,istartcol,1,1], &
              & ubounds=[config%n_spec_lw,iendcol,nlev+1,nblocks], persistent=.true.)
          call field_new(this%f_lw_dn_clear_band,lbounds=[1,istartcol,1,1], &
              & ubounds=[config%n_spec_lw,iendcol,nlev+1,nblocks], persistent=.true.)
        end if
      end if

      if (config%do_lw_derivatives) then
        call field_new(this%f_lw_derivatives,lbounds=[istartcol,1,1], &
            &                                ubounds=[iendcol,nlev+1,nblocks], persistent=.true.)
      end if

      if (config%do_toa_spectral_flux) then
        if (config%n_bands_lw == 0) then
          write(nulerr,'(a)') '*** Error: number of LW bands not yet defined ' &
               & // 'so cannot allocate TOA spectral flux arrays'
          call radiation_abort()
        end if
        call field_new(this%f_lw_up_toa_band,lbounds=[1,istartcol,1], &
            & ubounds=[config%n_bands_lw, iendcol,nblocks], persistent=.true.)
        if (config%do_clear) then
          call field_new(this%f_lw_up_toa_clear_band,lbounds=[1,istartcol,1], &
            & ubounds=[config%n_bands_lw, iendcol,nblocks], persistent=.true.)
        end if
      end if

      ! Allocate g-point downwelling fluxes at surface, and TOA fluxes
      call field_new(this%f_lw_dn_surf_g,lbounds=[1,istartcol,1], &
          & ubounds=[config%n_g_lw,iendcol,nblocks], persistent=.true.)
      call field_new(this%f_lw_up_toa_g ,lbounds=[1,istartcol,1], &
          & ubounds=[config%n_g_lw,iendcol,nblocks], persistent=.true.)
      if (config%do_clear) then
        call field_new(this%f_lw_dn_surf_clear_g,lbounds=[1,istartcol,1], &
            & ubounds=[config%n_g_lw,iendcol,nblocks], persistent=.true.)
        call field_new(this%f_lw_up_toa_clear_g ,lbounds=[1,istartcol,1], &
            & ubounds=[config%n_g_lw,iendcol,nblocks], persistent=.true.)
      end if

      if (config%do_canopy_fluxes_lw) then
        ! Downward fluxes at top of canopy at the spectral resolution
        ! used in the canopy radiative transfer scheme
        call field_new(this%f_lw_dn_surf_canopy,lbounds=[1,istartcol,1], &
            & ubounds=[config%n_canopy_bands_lw,iendcol,nblocks], persistent=.true.)
      end if
    end if

    ! Allocate shortwave arrays
    if (config%do_sw) then
      call field_new(this%f_sw_up,lbounds=[istartcol,1,1], ubounds=[iendcol,nlev+1,nblocks], persistent=.true.)
      call field_new(this%f_sw_dn,lbounds=[istartcol,1,1], ubounds=[iendcol,nlev+1,nblocks], persistent=.true.)
      if (config%do_sw_direct) then
        call field_new(this%f_sw_dn_direct,lbounds=[istartcol,1,1], &
            & ubounds=[iendcol,nlev+1,nblocks], persistent=.true.)
      end if
      if (config%do_clear) then
        call field_new(this%f_sw_up_clear,lbounds=[istartcol,1,1], &
            & ubounds=[iendcol,nlev+1,nblocks], persistent=.true.)
        call field_new(this%f_sw_dn_clear,lbounds=[istartcol,1,1], &
            & ubounds=[iendcol,nlev+1,nblocks], persistent=.true.)
        if (config%do_sw_direct) then
          call field_new(this%f_sw_dn_direct_clear,lbounds=[istartcol,1,1], &
              & ubounds=[iendcol,nlev+1,nblocks], persistent=.true.)
        end if
      end if

      if (config%do_save_spectral_flux) then
        if (config%n_spec_sw == 0) then
          write(nulerr,'(a)') '*** Error: number of SW spectral points to save not yet defined ' &
               & // 'so cannot allocate spectral flux arrays'
          call radiation_abort()
        end if

        call field_new(this%f_sw_up_band,lbounds=[1,istartcol,1,1], &
            & ubounds=[config%n_spec_sw,iendcol,nlev+1,nblocks], persistent=.true.)
        call field_new(this%f_sw_dn_band,lbounds=[1,istartcol,1,1], &
            & ubounds=[config%n_spec_sw,iendcol,nlev+1,nblocks], persistent=.true.)

        if (config%do_sw_direct) then
          call field_new(this%f_sw_dn_direct_band,lbounds=[1,istartcol,1,1], &
                & ubounds=[config%n_spec_sw,iendcol,nlev+1,nblocks], persistent=.true.)
        end if
        if (config%do_clear) then
          call field_new(this%f_sw_up_clear_band,lbounds=[1,istartcol,1,1], &
                & ubounds=[config%n_spec_sw,iendcol,nlev+1,nblocks], persistent=.true.)
          call field_new(this%f_sw_dn_clear_band,lbounds=[1,istartcol,1,1], &
                & ubounds=[config%n_spec_sw,iendcol,nlev+1,nblocks], persistent=.true.)
          if (config%do_sw_direct) then
            call field_new(this%f_sw_dn_direct_clear_band,lbounds=[1,istartcol,1,1], &
                & ubounds=[config%n_spec_sw,iendcol, nlev+1,nblocks], persistent=.true.)
          end if
        end if
      end if

      if (config%do_surface_sw_spectral_flux) then
        if (config%n_bands_sw == 0) then
          write(nulerr,'(a)') '*** Error: number of SW bands not yet defined ' &
               & // 'so cannot allocate TOA spectral flux arrays'
          call radiation_abort()
        end if
        call field_new(this%f_sw_dn_surf_band,lbounds=[1,istartcol,1], &
                & ubounds=[config%n_bands_sw,iendcol,nblocks], persistent=.true.)
        call field_new(this%f_sw_dn_direct_surf_band,lbounds=[1,istartcol,1], &
                & ubounds=[config%n_bands_sw,iendcol,nblocks], persistent=.true.)
        if (config%do_clear) then
          call field_new(this%f_sw_dn_surf_clear_band,lbounds=[1,istartcol,1], &
                  & ubounds=[config%n_bands_sw,iendcol,nblocks], persistent=.true.)
          call field_new(this%f_sw_dn_direct_surf_clear_band,lbounds=[1,istartcol,1], &
                  & ubounds=[config%n_bands_sw,iendcol,nblocks], persistent=.true.)
        end if
      end if

      if (config%do_toa_spectral_flux) then
        if (config%n_bands_sw == 0) then
          write(nulerr,'(a)') '*** Error: number of SW bands not yet defined ' &
               & // 'so cannot allocate surface spectral flux arrays'
          call radiation_abort()
        end if
        call field_new(this%f_sw_dn_toa_band,lbounds=[1,istartcol,1], &
            & ubounds=[config%n_bands_sw, iendcol,nblocks], persistent=.true.)
        call field_new(this%f_sw_up_toa_band,lbounds=[1,istartcol,1], &
            & ubounds=[config%n_bands_sw, iendcol,nblocks], persistent=.true.)
        if (config%do_clear) then
          call field_new(this%f_sw_up_toa_clear_band,lbounds=[1,istartcol,1], &
              & ubounds=[config%n_bands_sw, iendcol,nblocks], persistent=.true.)
        end if
      end if

      ! Allocate g-point downwelling fluxes at surface, and TOA fluxes
      call field_new(this%f_sw_dn_diffuse_surf_g,lbounds=[1,istartcol,1], &
          & ubounds=[config%n_g_sw,iendcol,nblocks], persistent=.true.)
      call field_new(this%f_sw_dn_direct_surf_g ,lbounds=[1,istartcol,1], &
          & ubounds=[config%n_g_sw,iendcol,nblocks], persistent=.true.)
      call field_new(this%f_sw_dn_toa_g         ,lbounds=[1,istartcol,1], &
          & ubounds=[config%n_g_sw,iendcol,nblocks], persistent=.true.)
      call field_new(this%f_sw_up_toa_g         ,lbounds=[1,istartcol,1], &
          & ubounds=[config%n_g_sw,iendcol,nblocks], persistent=.true.)
      if (config%do_clear) then
        call field_new(this%f_sw_dn_diffuse_surf_clear_g,lbounds=[1,istartcol,1], &
            & ubounds=[config%n_g_sw,iendcol,nblocks], persistent=.true.)
        call field_new(this%f_sw_dn_direct_surf_clear_g ,lbounds=[1,istartcol,1], &
            & ubounds=[config%n_g_sw,iendcol,nblocks], persistent=.true.)
        call field_new(this%f_sw_up_toa_clear_g         ,lbounds=[1,istartcol,1], &
            & ubounds=[config%n_g_sw,iendcol,nblocks], persistent=.true.)
      end if

      if (config%do_canopy_fluxes_sw) then
        ! Downward fluxes at top of canopy at the spectral resolution
        ! used in the canopy radiative transfer scheme
        call field_new(this%f_sw_dn_diffuse_surf_canopy,lbounds=[1,istartcol,1], &
            & ubounds=[config%n_canopy_bands_sw,iendcol,nblocks], persistent=.true.)
        call field_new(this%f_sw_dn_direct_surf_canopy ,lbounds=[1,istartcol,1], &
            & ubounds=[config%n_canopy_bands_sw,iendcol,nblocks], persistent=.true.)
      end if
    end if

    ! Allocate cloud cover arrays
    call field_new(this%f_cloud_cover_lw,lbounds=[istartcol,1], &
          & ubounds=[iendcol,nblocks], persistent=.true., init_value=-1.0_jprb)
    call field_new(this%f_cloud_cover_sw,lbounds=[istartcol,1], &
          & ubounds=[iendcol,nblocks], persistent=.true., init_value=-1.0_jprb)

    if (lhook) call dr_hook('radiation_flux:allocate',1,hook_handle)

    if (lhook) call dr_hook('radiation_radiation_field_type:flux_field_init',1,hook_handle)

  end subroutine flux_field_init



  !---------------------------------------------------------------------
  ! flux_field_type finalisation
  subroutine flux_field_final(this)

    use yomhook,     only : lhook, dr_hook, jphook

    class(flux_field_type), intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:flux_field_final',0,hook_handle)

    if (associated(this%f_lw_up)) then
      call field_delete(this%f_lw_up)
    end if
    this%f_lw_up=>null()
    this%lw_up=>null()

    if (associated(this%f_lw_dn)) then
      call field_delete(this%f_lw_dn)
    end if
    this%f_lw_dn=>null()
    this%lw_dn=>null()

    if (associated(this%f_sw_up)) then
      call field_delete(this%f_sw_up)
    end if
    this%f_sw_up=>null()
    this%sw_up=>null()

    if (associated(this%f_sw_dn)) then
      call field_delete(this%f_sw_dn)
    end if
    this%f_sw_dn=>null()
    this%sw_dn=>null()

    if (associated(this%f_sw_dn_direct)) then
      call field_delete(this%f_sw_dn_direct)
    end if
    this%f_sw_dn_direct=>null()
    this%sw_dn_direct=>null()

    if (associated(this%f_lw_up_clear)) then
      call field_delete(this%f_lw_up_clear)
    end if
    this%f_lw_up_clear=>null()
    this%lw_up_clear=>null()

    if (associated(this%f_lw_dn_clear)) then
      call field_delete(this%f_lw_dn_clear)
    end if
    this%f_lw_dn_clear=>null()
    this%lw_dn_clear=>null()

    if (associated(this%f_sw_up_clear)) then
      call field_delete(this%f_sw_up_clear)
    end if
    this%f_sw_up_clear=>null()
    this%sw_up_clear=>null()

    if (associated(this%f_sw_dn_clear)) then
      call field_delete(this%f_sw_dn_clear)
    end if
    this%f_sw_dn_clear=>null()
    this%sw_dn_clear=>null()

    if (associated(this%f_sw_dn_direct_clear)) then
      call field_delete(this%f_sw_dn_direct_clear)
    end if
    this%f_sw_dn_direct_clear=>null()
    this%sw_dn_direct_clear=>null()

    if (associated(this%f_lw_up_band)) then
      call field_delete(this%f_lw_up_band)
    end if
    this%f_lw_up_band=>null()
    this%lw_up_band=>null()

    if (associated(this%f_lw_dn_band)) then
      call field_delete(this%f_lw_dn_band)
    end if
    this%f_lw_dn_band=>null()
    this%lw_dn_band=>null()

    if (associated(this%f_sw_up_band)) then
      call field_delete(this%f_sw_up_band)
    end if
    this%f_sw_up_band=>null()
    this%sw_up_band=>null()

    if (associated(this%f_sw_dn_band)) then
      call field_delete(this%f_sw_dn_band)
    end if
    this%f_sw_dn_band=>null()
    this%sw_dn_band=>null()

    if (associated(this%f_sw_dn_direct_band)) then
      call field_delete(this%f_sw_dn_direct_band)
    end if
    this%f_sw_dn_direct_band=>null()
    this%sw_dn_direct_band=>null()

    if (associated(this%f_lw_up_clear_band)) then
      call field_delete(this%f_lw_up_clear_band)
    end if
    this%f_lw_up_clear_band=>null()
    this%lw_up_clear_band=>null()

    if (associated(this%f_lw_dn_clear_band)) then
      call field_delete(this%f_lw_dn_clear_band)
    end if
    this%f_lw_dn_clear_band=>null()
    this%lw_dn_clear_band=>null()

    if (associated(this%f_sw_up_clear_band)) then
      call field_delete(this%f_sw_up_clear_band)
    end if
    this%f_sw_up_clear_band=>null()
    this%sw_up_clear_band=>null()

    if (associated(this%f_sw_dn_clear_band)) then
      call field_delete(this%f_sw_dn_clear_band)
    end if
    this%f_sw_dn_clear_band=>null()
    this%sw_dn_clear_band=>null()

    if (associated(this%f_sw_dn_direct_clear_band)) then
      call field_delete(this%f_sw_dn_direct_clear_band)
    end if
    this%f_sw_dn_direct_clear_band=>null()
    this%sw_dn_direct_clear_band=>null()

    if (associated(this%f_lw_dn_surf_g)) then
      call field_delete(this%f_lw_dn_surf_g)
    end if
    this%f_lw_dn_surf_g=>null()
    this%lw_dn_surf_g=>null()

    if (associated(this%f_lw_dn_surf_clear_g)) then
      call field_delete(this%f_lw_dn_surf_clear_g)
    end if
    this%f_lw_dn_surf_clear_g=>null()
    this%lw_dn_surf_clear_g=>null()

    if (associated(this%f_sw_dn_diffuse_surf_g)) then
      call field_delete(this%f_sw_dn_diffuse_surf_g)
    end if
    this%f_sw_dn_diffuse_surf_g=>null()
    this%sw_dn_diffuse_surf_g=>null()

    if (associated(this%f_sw_dn_direct_surf_g)) then
      call field_delete(this%f_sw_dn_direct_surf_g)
    end if
    this%f_sw_dn_direct_surf_g=>null()
    this%sw_dn_direct_surf_g=>null()

    if (associated(this%f_sw_dn_diffuse_surf_clear_g)) then
      call field_delete(this%f_sw_dn_diffuse_surf_clear_g)
    end if
    this%f_sw_dn_diffuse_surf_clear_g=>null()
    this%sw_dn_diffuse_surf_clear_g=>null()

    if (associated(this%f_sw_dn_direct_surf_clear_g)) then
      call field_delete(this%f_sw_dn_direct_surf_clear_g)
    end if
    this%f_sw_dn_direct_surf_clear_g=>null()
    this%sw_dn_direct_surf_clear_g=>null()

    if (associated(this%f_lw_up_toa_g)) then
      call field_delete(this%f_lw_up_toa_g)
    end if
    this%f_lw_up_toa_g=>null()
    this%lw_up_toa_g=>null()

    if (associated(this%f_lw_up_toa_clear_g)) then
      call field_delete(this%f_lw_up_toa_clear_g)
    end if
    this%f_lw_up_toa_clear_g=>null()
    this%lw_up_toa_clear_g=>null()

    if (associated(this%f_sw_dn_toa_g)) then
      call field_delete(this%f_sw_dn_toa_g)
    end if
    this%f_sw_dn_toa_g=>null()
    this%sw_dn_toa_g=>null()

    if (associated(this%f_sw_up_toa_g)) then
      call field_delete(this%f_sw_up_toa_g)
    end if
    this%f_sw_up_toa_g=>null()
    this%sw_up_toa_g=>null()

    if (associated(this%f_sw_up_toa_clear_g)) then
      call field_delete(this%f_sw_up_toa_clear_g)
    end if
    this%f_sw_up_toa_clear_g=>null()
    this%sw_up_toa_clear_g=>null()

    if (associated(this%f_sw_dn_surf_band)) then
      call field_delete(this%f_sw_dn_surf_band)
    end if
    this%f_sw_dn_surf_band=>null()
    this%sw_dn_surf_band=>null()

    if (associated(this%f_sw_dn_direct_surf_band)) then
      call field_delete(this%f_sw_dn_direct_surf_band)
    end if
    this%f_sw_dn_direct_surf_band=>null()
    this%sw_dn_direct_surf_band=>null()

    if (associated(this%f_sw_dn_surf_clear_band)) then
      call field_delete(this%f_sw_dn_surf_clear_band)
    end if
    this%f_sw_dn_surf_clear_band=>null()
    this%sw_dn_surf_clear_band=>null()

    if (associated(this%f_sw_dn_direct_surf_clear_band)) then
      call field_delete(this%f_sw_dn_direct_surf_clear_band)
    end if
    this%f_sw_dn_direct_surf_clear_band=>null()
    this%sw_dn_direct_surf_clear_band=>null()

    if (associated(this%f_lw_up_toa_band)) then
      call field_delete(this%f_lw_up_toa_band)
    end if
    this%f_lw_up_toa_band=>null()
    this%lw_up_toa_band=>null()

    if (associated(this%f_lw_up_toa_clear_band)) then
      call field_delete(this%f_lw_up_toa_clear_band)
    end if
    this%f_lw_up_toa_clear_band=>null()
    this%lw_up_toa_clear_band=>null()

    if (associated(this%f_sw_dn_toa_band)) then
      call field_delete(this%f_sw_dn_toa_band)
    end if
    this%f_sw_dn_toa_band=>null()
    this%sw_dn_toa_band=>null()

    if (associated(this%f_sw_up_toa_band)) then
      call field_delete(this%f_sw_up_toa_band)
    end if
    this%f_sw_up_toa_band=>null()
    this%sw_up_toa_band=>null()

    if (associated(this%f_sw_up_toa_clear_band)) then
      call field_delete(this%f_sw_up_toa_clear_band)
    end if
    this%f_sw_up_toa_clear_band=>null()
    this%sw_up_toa_clear_band=>null()

    if (associated(this%f_lw_dn_surf_canopy)) then
      call field_delete(this%f_lw_dn_surf_canopy)
    end if
    this%f_lw_dn_surf_canopy=>null()
    this%lw_dn_surf_canopy=>null()

    if (associated(this%f_sw_dn_diffuse_surf_canopy)) then
      call field_delete(this%f_sw_dn_diffuse_surf_canopy)
    end if
    this%f_sw_dn_diffuse_surf_canopy=>null()
    this%sw_dn_diffuse_surf_canopy=>null()

    if (associated(this%f_sw_dn_direct_surf_canopy)) then
      call field_delete(this%f_sw_dn_direct_surf_canopy)
    end if
    this%f_sw_dn_direct_surf_canopy=>null()
    this%sw_dn_direct_surf_canopy=>null()

    if (associated(this%f_cloud_cover_lw)) then
      call field_delete(this%f_cloud_cover_lw)
    end if
    this%f_cloud_cover_lw=>null()
    this%cloud_cover_lw=>null()

    if (associated(this%f_cloud_cover_sw)) then
      call field_delete(this%f_cloud_cover_sw)
    end if
    this%f_cloud_cover_sw=>null()
    this%cloud_cover_sw=>null()

    if (associated(this%f_lw_derivatives)) then
      call field_delete(this%f_lw_derivatives)
    end if
    this%f_lw_derivatives=>null()
    this%lw_derivatives=>null()

    if (lhook) call dr_hook('radiation_field_type:flux_field_final',1,hook_handle)

  end subroutine flux_field_final

  !---------------------------------------------------------------------
  ! Update view pointers of flux_field_type
  subroutine flux_field_update_view(this, block_index)

    use yomhook,     only : lhook, dr_hook, jphook

    class(flux_field_type), intent(inout)  :: this
    integer, intent(in)                     :: block_index

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:flux_field_update_view',0,hook_handle)

    if (associated(this%f_lw_up)) then
      this%lw_up => this%f_lw_up%get_view(block_index)
    end if

    if (associated(this%f_lw_dn)) then
      this%lw_dn => this%f_lw_dn%get_view(block_index)
    end if

    if (associated(this%f_sw_up)) then
      this%sw_up => this%f_sw_up%get_view(block_index)
    end if

    if (associated(this%f_sw_dn)) then
      this%sw_dn => this%f_sw_dn%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_direct)) then
      this%sw_dn_direct => this%f_sw_dn_direct%get_view(block_index)
    end if

    if (associated(this%f_lw_up_clear)) then
      this%lw_up_clear => this%f_lw_up_clear%get_view(block_index)
    end if

    if (associated(this%f_lw_dn_clear)) then
      this%lw_dn_clear => this%f_lw_dn_clear%get_view(block_index)
    end if

    if (associated(this%f_sw_up_clear)) then
      this%sw_up_clear => this%f_sw_up_clear%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_clear)) then
      this%sw_dn_clear => this%f_sw_dn_clear%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_direct_clear)) then
      this%sw_dn_direct_clear => this%f_sw_dn_direct_clear%get_view(block_index)
    end if

    if (associated(this%f_lw_up_band)) then
      this%lw_up_band => this%f_lw_up_band%get_view(block_index)
    end if

    if (associated(this%f_lw_dn_band)) then
      this%lw_dn_band => this%f_lw_dn_band%get_view(block_index)
    end if

    if (associated(this%f_sw_up_band)) then
      this%sw_up_band => this%f_sw_up_band%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_band)) then
      this%sw_dn_band => this%f_sw_dn_band%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_direct_band)) then
      this%sw_dn_direct_band => this%f_sw_dn_direct_band%get_view(block_index)
    end if

    if (associated(this%f_lw_up_clear_band)) then
      this%lw_up_clear_band => this%f_lw_up_clear_band%get_view(block_index)
    end if

    if (associated(this%f_lw_dn_clear_band)) then
      this%lw_dn_clear_band => this%f_lw_dn_clear_band%get_view(block_index)
    end if

    if (associated(this%f_sw_up_clear_band)) then
      this%sw_up_clear_band => this%f_sw_up_clear_band%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_clear_band)) then
      this%sw_dn_clear_band => this%f_sw_dn_clear_band%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_direct_clear_band)) then
      this%sw_dn_direct_clear_band => this%f_sw_dn_direct_clear_band%get_view(block_index)
    end if

    if (associated(this%f_lw_dn_surf_g)) then
      this%lw_dn_surf_g => this%f_lw_dn_surf_g%get_view(block_index)
    end if

    if (associated(this%f_lw_dn_surf_clear_g)) then
      this%lw_dn_surf_clear_g => this%f_lw_dn_surf_clear_g%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_diffuse_surf_g)) then
      this%sw_dn_diffuse_surf_g => this%f_sw_dn_diffuse_surf_g%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_direct_surf_g)) then
      this%sw_dn_direct_surf_g => this%f_sw_dn_direct_surf_g%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_diffuse_surf_clear_g)) then
      this%sw_dn_diffuse_surf_clear_g => this%f_sw_dn_diffuse_surf_clear_g%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_direct_surf_clear_g)) then
      this%sw_dn_direct_surf_clear_g => this%f_sw_dn_direct_surf_clear_g%get_view(block_index)
    end if

    if (associated(this%f_lw_up_toa_g)) then
      this%lw_up_toa_g => this%f_lw_up_toa_g%get_view(block_index)
    end if

    if (associated(this%f_lw_up_toa_clear_g)) then
      this%lw_up_toa_clear_g => this%f_lw_up_toa_clear_g%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_toa_g)) then
      this%sw_dn_toa_g => this%f_sw_dn_toa_g%get_view(block_index)
    end if

    if (associated(this%f_sw_up_toa_g)) then
      this%sw_up_toa_g => this%f_sw_up_toa_g%get_view(block_index)
    end if

    if (associated(this%f_sw_up_toa_clear_g)) then
      this%sw_up_toa_clear_g => this%f_sw_up_toa_clear_g%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_surf_band)) then
      this%sw_dn_surf_band => this%f_sw_dn_surf_band%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_direct_surf_band)) then
      this%sw_dn_direct_surf_band => this%f_sw_dn_direct_surf_band%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_surf_clear_band)) then
      this%sw_dn_surf_clear_band => this%f_sw_dn_surf_clear_band%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_direct_surf_clear_band)) then
      this%sw_dn_direct_surf_clear_band => this%f_sw_dn_direct_surf_clear_band%get_view(block_index)
    end if

    if (associated(this%f_lw_up_toa_band)) then
      this%lw_up_toa_band => this%f_lw_up_toa_band%get_view(block_index)
    end if

    if (associated(this%f_lw_up_toa_clear_band)) then
      this%lw_up_toa_clear_band => this%f_lw_up_toa_clear_band%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_toa_band)) then
      this%sw_dn_toa_band => this%f_sw_dn_toa_band%get_view(block_index)
    end if

    if (associated(this%f_sw_up_toa_band)) then
      this%sw_up_toa_band => this%f_sw_up_toa_band%get_view(block_index)
    end if

    if (associated(this%f_sw_up_toa_clear_band)) then
      this%sw_up_toa_clear_band => this%f_sw_up_toa_clear_band%get_view(block_index)
    end if

    if (associated(this%f_lw_dn_surf_canopy)) then
      this%lw_dn_surf_canopy => this%f_lw_dn_surf_canopy%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_diffuse_surf_canopy)) then
      this%sw_dn_diffuse_surf_canopy => this%f_sw_dn_diffuse_surf_canopy%get_view(block_index)
    end if

    if (associated(this%f_sw_dn_direct_surf_canopy)) then
      this%sw_dn_direct_surf_canopy => this%f_sw_dn_direct_surf_canopy%get_view(block_index)
    end if

    if (associated(this%f_cloud_cover_lw)) then
      this%cloud_cover_lw => this%f_cloud_cover_lw%get_view(block_index)
    end if

    if (associated(this%f_cloud_cover_sw)) then
      this%cloud_cover_sw => this%f_cloud_cover_sw%get_view(block_index)
    end if

    if (associated(this%f_lw_derivatives)) then
      this%lw_derivatives => this%f_lw_derivatives%get_view(block_index)
    end if

    if (lhook) call dr_hook('radiation_radiation_field_type:flux_field_update_view',1,hook_handle)

  end subroutine flux_field_update_view

  !---------------------------------------------------------------------
  ! Update flux pointers of flux_field_type
  subroutine flux_field_update_flux(this, flux)

    use yomhook,     only : lhook, dr_hook, jphook

    class(flux_field_type), intent(inout) :: this
    class(flux_type), intent(inout)       :: flux

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:flux_field_update_flux',0,hook_handle)

    if (associated(this%lw_up)) then
      flux%lw_up => this%lw_up
    end if

    if (associated(this%lw_dn)) then
      flux%lw_dn => this%lw_dn
    end if

    if (associated(this%sw_up)) then
      flux%sw_up => this%sw_up
    end if

    if (associated(this%sw_dn)) then
      flux%sw_dn => this%sw_dn
    end if

    if (associated(this%sw_dn_direct)) then
      flux%sw_dn_direct => this%sw_dn_direct
    end if

    if (associated(this%lw_up_clear)) then
      flux%lw_up_clear => this%lw_up_clear
    end if

    if (associated(this%lw_dn_clear)) then
      flux%lw_dn_clear => this%lw_dn_clear
    end if

    if (associated(this%sw_up_clear)) then
      flux%sw_up_clear => this%sw_up_clear
    end if

    if (associated(this%sw_dn_clear)) then
      flux%sw_dn_clear => this%sw_dn_clear
    end if

    if (associated(this%sw_dn_direct_clear)) then
      flux%sw_dn_direct_clear => this%sw_dn_direct_clear
    end if

    if (associated(this%lw_up_band)) then
      flux%lw_up_band => this%lw_up_band
    end if

    if (associated(this%lw_dn_band)) then
      flux%lw_dn_band => this%lw_dn_band
    end if

    if (associated(this%sw_up_band)) then
      flux%sw_up_band => this%sw_up_band
    end if

    if (associated(this%sw_dn_band)) then
      flux%sw_dn_band => this%sw_dn_band
    end if

    if (associated(this%sw_dn_direct_band)) then
      flux%sw_dn_direct_band => this%sw_dn_direct_band
    end if

    if (associated(this%lw_up_clear_band)) then
      flux%lw_up_clear_band => this%lw_up_clear_band
    end if

    if (associated(this%lw_dn_clear_band)) then
      flux%lw_dn_clear_band => this%lw_dn_clear_band
    end if

    if (associated(this%sw_up_clear_band)) then
      flux%sw_up_clear_band => this%sw_up_clear_band
    end if

    if (associated(this%sw_dn_clear_band)) then
      flux%sw_dn_clear_band => this%sw_dn_clear_band
    end if

    if (associated(this%sw_dn_direct_clear_band)) then
      flux%sw_dn_direct_clear_band => this%sw_dn_direct_clear_band
    end if

    if (associated(this%lw_dn_surf_g)) then
      flux%lw_dn_surf_g => this%lw_dn_surf_g
    end if

    if (associated(this%lw_dn_surf_clear_g)) then
      flux%lw_dn_surf_clear_g => this%lw_dn_surf_clear_g
    end if

    if (associated(this%sw_dn_diffuse_surf_g)) then
      flux%sw_dn_diffuse_surf_g => this%sw_dn_diffuse_surf_g
    end if

    if (associated(this%sw_dn_direct_surf_g)) then
      flux%sw_dn_direct_surf_g => this%sw_dn_direct_surf_g
    end if

    if (associated(this%sw_dn_diffuse_surf_clear_g)) then
      flux%sw_dn_diffuse_surf_clear_g => this%sw_dn_diffuse_surf_clear_g
    end if

    if (associated(this%sw_dn_direct_surf_clear_g)) then
      flux%sw_dn_direct_surf_clear_g => this%sw_dn_direct_surf_clear_g
    end if

    if (associated(this%lw_up_toa_g)) then
      flux%lw_up_toa_g => this%lw_up_toa_g
    end if

    if (associated(this%lw_up_toa_clear_g)) then
      flux%lw_up_toa_clear_g => this%lw_up_toa_clear_g
    end if

    if (associated(this%sw_dn_toa_g)) then
      flux%sw_dn_toa_g => this%sw_dn_toa_g
    end if

    if (associated(this%sw_up_toa_g)) then
      flux%sw_up_toa_g => this%sw_up_toa_g
    end if

    if (associated(this%sw_up_toa_clear_g)) then
      flux%sw_up_toa_clear_g => this%sw_up_toa_clear_g
    end if

    if (associated(this%sw_dn_surf_band)) then
      flux%sw_dn_surf_band => this%sw_dn_surf_band
    end if

    if (associated(this%sw_dn_direct_surf_band)) then
      flux%sw_dn_direct_surf_band => this%sw_dn_direct_surf_band
    end if

    if (associated(this%sw_dn_surf_clear_band)) then
      flux%sw_dn_surf_clear_band => this%sw_dn_surf_clear_band
    end if

    if (associated(this%sw_dn_direct_surf_clear_band)) then
      flux%sw_dn_direct_surf_clear_band => this%sw_dn_direct_surf_clear_band
    end if

    if (associated(this%lw_up_toa_band)) then
      flux%lw_up_toa_band => this%lw_up_toa_band
    end if

    if (associated(this%lw_up_toa_clear_band)) then
      flux%lw_up_toa_clear_band => this%lw_up_toa_clear_band
    end if

    if (associated(this%sw_dn_toa_band)) then
      flux%sw_dn_toa_band => this%sw_dn_toa_band
    end if

    if (associated(this%sw_up_toa_band)) then
      flux%sw_up_toa_band => this%sw_up_toa_band
    end if

    if (associated(this%sw_up_toa_clear_band)) then
      flux%sw_up_toa_clear_band => this%sw_up_toa_clear_band
    end if

    if (associated(this%lw_dn_surf_canopy)) then
      flux%lw_dn_surf_canopy => this%lw_dn_surf_canopy
    end if

    if (associated(this%sw_dn_diffuse_surf_canopy)) then
      flux%sw_dn_diffuse_surf_canopy => this%sw_dn_diffuse_surf_canopy
    end if

    if (associated(this%sw_dn_direct_surf_canopy)) then
      flux%sw_dn_direct_surf_canopy => this%sw_dn_direct_surf_canopy
    end if

    if (associated(this%cloud_cover_lw)) then
      flux%cloud_cover_lw => this%cloud_cover_lw
    end if

    if (associated(this%cloud_cover_sw)) then
      flux%cloud_cover_sw => this%cloud_cover_sw
    end if

    if (associated(this%lw_derivatives)) then
      flux%lw_derivatives => this%lw_derivatives
    end if

    if (lhook) call dr_hook('radiation_radiation_field_type:flux_field_update_flux',1,hook_handle)

  end subroutine flux_field_update_flux

end module radiation_field_type_module

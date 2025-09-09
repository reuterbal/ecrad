# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

include(ecrad_fypp_macros)

set( ifs_SOURCES
        ice_effective_radius.F90
        liquid_effective_radius.F90
        radiation_scheme.F90
        radiation_setup.F90
        yoerdu.F90
        yomrip.F90
        yoephy.F90
        yoecld.F90
        yoe_spectral_planck.F90
        cloud_overlap_decorr_len.F90
        yoerad.F90
        ifs_blocking.F90
)

# Export the list of sources for Loki processing
list( TRANSFORM ifs_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/" )
set( ifs_SOURCES ${ifs_SOURCES} PARENT_SCOPE )

# Generate radintg with fypp (if enabled)
ecrad_ifs_process_fypp( radintg_zrgp_mod )

ecbuild_add_library(
    TARGET ifs.${PREC}
    TYPE OBJECT
    SOURCES
        ${ifs_SOURCES}
        radintg_zrgp_mod.F90
    PUBLIC_INCLUDES
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/module_ifs>"
    PRIVATE_LIBS
        ecrad.${PREC}
        ecrad_base.${PREC}
)
set_target_properties( ifs.${PREC}
    PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module_ifs"
)


if( HAVE_FIELD_API )

    # fypp radintg and add it to the source list
    ecrad_ifs_process_fypp( radiation_scheme_layer_mod )

    ecbuild_add_library(
        TARGET ifs_field_api.${PREC}
        TYPE OBJECT
        SOURCES
            ${ifs_SOURCES}
            radintg_zrgp_mod.F90
            radiation_scheme_layer_mod.F90
        PUBLIC_DEFINITIONS
            HAVE_FIELD_API
        PUBLIC_INCLUDES
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/module_ifs_field_api>"
        PRIVATE_LIBS
            ecrad.${PREC}
            ecrad_base.${PREC}
        PUBLIC_LIBS
            field_api_${PREC}
            $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
    )

    set_target_properties( ifs_field_api.${PREC}
        PROPERTIES
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module_ifs_field_api"
    )

endif()

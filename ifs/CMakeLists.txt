# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

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

# If FYPP has been enabled, regenerate the F90 source files first
if( HAVE_FYPP AND NOT TARGET ecrad_fypp_files )
    set( fypp_files radiation_scheme_layer ifs_blocking )
    foreach( fypp_file ${fypp_files} )
        add_custom_command(
            OUTPUT
                ${CMAKE_CURRENT_SOURCE_DIR}/${fypp_file}.F90
            COMMAND
                ${FYPP} -m os -m field_config -M ${CMAKE_CURRENT_SOURCE_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}/${fypp_file}.fypp
                ${CMAKE_CURRENT_SOURCE_DIR}/${fypp_file}.F90
            DEPENDS
                ${CMAKE_CURRENT_SOURCE_DIR}/${fypp_file}.fypp
                ${CMAKE_CURRENT_SOURCE_DIR}/radintg_config.yaml
        )
    endforeach()
    list( TRANSFORM fypp_files APPEND .F90 )
    add_custom_target( ecrad_fypp_files ALL DEPENDS ${fypp_files} )
endif()

ecbuild_add_library(
    TARGET ifs.${PREC}
    TYPE OBJECT
    SOURCES ${ifs_SOURCES}
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

    ecbuild_add_library(
        TARGET ifs_field_api.${PREC}
        TYPE OBJECT
        SOURCES
            ${ifs_SOURCES}
            radiation_scheme_layer.F90
        DEFINITIONS
            $<$<NOT:${HAVE_OMP}>:NO_OPENMP>
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

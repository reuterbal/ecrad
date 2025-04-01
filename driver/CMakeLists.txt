# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

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
    TARGET driver_lib.${PREC}
    TYPE OBJECT
    SOURCES
        ecrad_driver_read_input.F90
        ecrad_driver_config.F90
        print_matrix_mod.F90
    PUBLIC_DEFINITIONS
        $<$<NOT:${HAVE_OMP}>:NO_OPENMP>
    PUBLIC_LIBS
        ecrad.${PREC}
        ecrad_base.${PREC}
        $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
)

ecbuild_add_executable(
    TARGET ecrad_${PREC}
    SOURCES
        ecrad_driver.F90
    LIBS
        driver_lib.${PREC}
)

ecbuild_add_executable(
    TARGET ecrad_ifs_${PREC}
    SOURCES
        ecrad_ifs_driver.F90
    LIBS
        ifs.${PREC}
        driver_lib.${PREC}
        $<${HAVE_NVTX}:${NVTX_TARGET}>
)

ecbuild_add_executable(
    TARGET ecrad_ifs_blocked_${PREC}
    SOURCES
        ecrad_ifs_driver_blocked.F90
        ifs_blocking.F90
    DEFINITIONS
        $<$<NOT:${HAVE_OMP}>:NO_OPENMP>
    LIBS
        driver_lib.${PREC}
        ifs.${PREC}
        ecrad.${PREC}
        $<${HAVE_NVTX}:${NVTX_TARGET}>
)
set_target_properties( ecrad_ifs_blocked_${PREC}
    PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module-blocked"
)

# For ecrad_ifs_driver_blocked we have to disable bounds checking
# (which is enabled in Debug builds) because sequence association
# and values initialized to -9999 for unused arguments is used in
# the call to RADIATION_SCHEME.
if( CMAKE_BUILD_TYPE MATCHES "Debug" )
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        set_source_files_properties(
            "ecrad_ifs_driver_blocked.F90"
            PROPERTIES COMPILE_FLAGS "-fcheck=no-bounds"
        )
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set_source_files_properties(
            "ecrad_ifs_driver_blocked.F90"
            PROPERTIES COMPILE_FLAGS "-check nobounds"
        )
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
        set_source_files_properties(
            "ecrad_ifs_driver_blocked.F90"
            PROPERTIES COMPILE_FLAGS "-Mnobounds"
        )
    endif()
endif()

if( HAVE_FIELD_API )
    ecbuild_add_executable(
        TARGET ecrad_ifs_field_api_${PREC}
        SOURCES
            ecrad_ifs_driver_field_api.F90
            ifs_blocking.F90
            radiation_scheme_layer.F90
        DEFINITIONS
            $<$<NOT:${HAVE_OMP}>:NO_OPENMP>
            HAVE_FIELD_API
        LIBS
            driver_lib.${PREC}
            ifs.${PREC}
            ecrad.${PREC}
            field_api_${PREC}
    )
    set_target_properties( ecrad_ifs_field_api_${PREC}
        PROPERTIES
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module-field_api"
    )
endif()

# Create a symlink for each driver without the precision suffix.
# The default is to link to the double precision version if it is built,
# and single precision otherwise.
if ( "${PREC}" STREQUAL "dp" OR NOT HAVE_DOUBLE_PRECISION )
    foreach( tgt ecrad ecrad_ifs ecrad_ifs_blocked ecrad_ifs_field_api )
        if( TARGET ${tgt}_${PREC} )
            add_custom_command(
                TARGET ${tgt}_${PREC}
                POST_BUILD
                COMMAND
                    "${CMAKE_COMMAND}" -E create_symlink ${tgt}_${PREC} ${tgt}
                BYPRODUCTS ${tgt}
                WORKING_DIRECTORY
                    "${CMAKE_BINARY_DIR}/bin"
                COMMENT "Creating symbolic link from ${tgt} to ${tgt}_${PREC}"
            )
            install(
                FILES ${CMAKE_BINARY_DIR}/bin/${tgt}
                DESTINATION bin
            )
        endif()
    endforeach()
endif()

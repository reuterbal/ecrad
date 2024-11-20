# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

macro( ecrad_ifs_process_fypp fypp_file )
    # If FYPP has been enabled, regenerate the F90 source files otherwise simply link
    # or copy commited versions
    if( HAVE_FYPP )
        add_custom_command(
            OUTPUT
                ${CMAKE_CURRENT_BINARY_DIR}/${fypp_file}.F90
            COMMAND
                ${FYPP} -m os -m field_config -M ${CMAKE_CURRENT_SOURCE_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}/${fypp_file}.fypp
                ${CMAKE_CURRENT_BINARY_DIR}/${fypp_file}.F90
            DEPENDS
                ${CMAKE_CURRENT_SOURCE_DIR}/${fypp_file}.fypp
                ${CMAKE_CURRENT_SOURCE_DIR}/radiation_fields_config.yaml
        )
    else()
        file( CREATE_LINK
            ${CMAKE_CURRENT_SOURCE_DIR}/${fypp_file}.F90
            ${CMAKE_CURRENT_BINARY_DIR}/${fypp_file}.F90
            SYMBOLIC COPY_ON_ERROR
        )
    endif()
endmacro()


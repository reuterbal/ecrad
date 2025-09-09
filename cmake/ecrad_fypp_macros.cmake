# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

macro( ecrad_ifs_process_fypp fypp_file )

    get_filename_component(input_dir ${fypp_file} DIRECTORY)
    if(input_dir STREQUAL "")
        # Check for file in current source dir if only a file name is given
        set(fypp_file_dir ${CMAKE_CURRENT_SOURCE_DIR})
        set(file_name ${fypp_file})
        set(file_path ${CMAKE_CURRENT_SOURCE_DIR}/${fypp_file})
    else()
        # Assume full path is given
        get_filename_component(file_path ${fypp_file} ABSOLUTE
                              BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
        get_filename_component(fypp_file_dir ${file_path} DIRECTORY)
        get_filename_component(file_name ${file_path} NAME)
    endif()
    if (NOT EXISTS ${file_path}.fypp)
        message (FATAL_ERROR "The fypp file: ${file_name}.fypp does not exist in the current source directory: ${fypp_file_dir}")
    endif()
    # If FYPP has been enabled, regenerate the F90 source files otherwise simply link
    # or copy commited versions
    if( HAVE_FYPP )
        add_custom_command(
            OUTPUT
                ${CMAKE_CURRENT_BINARY_DIR}/${file_name}.F90
            COMMAND
                ${FYPP} -m os -m field_config -M ${fypp_file_dir}
                ${file_path}.fypp
                ${file_path}.F90
            DEPENDS
                ${file_path}.fypp
                ${fypp_file_dir}/radiation_fields_config.yaml
        )
    else()
        file( CREATE_LINK
            ${file_path}.F90
            ${CMAKE_CURRENT_BINARY_DIR}/${file_name}.F90
            SYMBOLIC COPY_ON_ERROR
        )
    endif()
endmacro()


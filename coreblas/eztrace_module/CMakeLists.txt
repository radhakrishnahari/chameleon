###
#
# @file CMakeLists.txt
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
###
#
#  @project CHAMELEON
#  CHAMELEON is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
# @version 1.3.0
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2024-02-18
#
###

if (NOT EZTRACE_FOUND)
    find_package(EZTRACE)
endif()

if (EZTRACE_FOUND)

    set(EZTRACE_eztrace_create_plugin_DIR "EZTRACE_eztrace_create_plugin_DIR-NOTFOUND")
    find_path(EZTRACE_eztrace_create_plugin_DIR
      NAMES eztrace_create_plugin
      HINTS ${EZTRACE_DIR_FOUND}/bin)
    mark_as_advanced(EZTRACE_eztrace_create_plugin_DIR)

    if (EZTRACE_eztrace_create_plugin_DIR)

        set(EZTRACE_CREATE_PLUGIN "${EZTRACE_eztrace_create_plugin_DIR}/eztrace_create_plugin")

        set( COREBLAS_EZTRACE_MODULES core tcore )

        foreach( _module ${COREBLAS_EZTRACE_MODULES} )
            add_custom_command(
                OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_module}
                COMMAND ${EZTRACE_CREATE_PLUGIN} -o ${CMAKE_CURRENT_BINARY_DIR}/${_module} > /dev/null
                ARGS ${CMAKE_CURRENT_SOURCE_DIR}/coreblas_${_module}_eztrace_module
                DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/coreblas_${_module}_eztrace_module
                )
            add_custom_target(
                eztrace-module-chameleon_${_module}-dir ALL
                DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${_module}
                )
            add_custom_command(
                OUTPUT libeztrace-convert-chameleon_${_module}.so
                COMMAND make
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_module}
                DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${_module}
                )
            add_custom_target(
                eztrace-module-chameleon_${_module}-libs ALL
                DEPENDS libeztrace-convert-chameleon_${_module}.so
                )
            # installation
            # ------------
            install(
                FILES
                ${CMAKE_CURRENT_BINARY_DIR}/${_module}/libeztrace-autostart-chameleon_${_module}.so
                ${CMAKE_CURRENT_BINARY_DIR}/${_module}/libeztrace-chameleon_${_module}.so
                ${CMAKE_CURRENT_BINARY_DIR}/${_module}/libeztrace-convert-chameleon_${_module}.so
                DESTINATION ${EZTRACE_LIBRARY_DIRS}
                )
        endforeach()
    endif (EZTRACE_eztrace_create_plugin_DIR)

endif (EZTRACE_FOUND)

###
### END CMakeLists.txt
###

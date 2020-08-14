###
#
# @file GenPkgConfig.cmake
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
# @version 1.0.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2020-03-03
#
###

###
#
# CONVERT_LIBSTYLE_TO_PKGCONFIG: convert a libraries list to follow the pkg-config style
#                                used in CLEAN_LIB_LIST
#
###
MACRO(CONVERT_LIBSTYLE_TO_PKGCONFIG _liblist)
    set(${_liblist}_CPY "${${_liblist}}")
    set(${_liblist} "")
    foreach(_dep ${${_liblist}_CPY})
        if (${_dep} MATCHES "^/")
            get_filename_component(dep_libname ${_dep} NAME)
            get_filename_component(dep_libdir  ${_dep} DIRECTORY)
            STRING(REPLACE "lib"    "" dep_libname "${dep_libname}")
            STRING(REPLACE ".so"    "" dep_libname "${dep_libname}")
            STRING(REPLACE ".a"     "" dep_libname "${dep_libname}")
            STRING(REPLACE ".dylib" "" dep_libname "${dep_libname}")
            STRING(REPLACE ".dll"   "" dep_libname "${dep_libname}")
            list(APPEND ${_liblist} -L${dep_libdir} -l${dep_libname})
        elseif(NOT ${_dep} MATCHES "^-")
            list(APPEND ${_liblist} "-l${_dep}")
        else()
            list(APPEND ${_liblist} ${_dep})
        endif()
    endforeach()
ENDMACRO(CONVERT_LIBSTYLE_TO_PKGCONFIG)

###
#
# CLEAN_LIB_LIST: clean libraries lists to follow the pkg-config style
#                 used in GENERATE_PKGCONFIG_FILE
#
###
MACRO(CLEAN_LIB_LIST _package)
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_LIBS)
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_LIBS_PRIVATE)
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_REQUIRED)
    list(REMOVE_DUPLICATES ${_package}_PKGCONFIG_REQUIRED_PRIVATE)
    CONVERT_LIBSTYLE_TO_PKGCONFIG(${_package}_PKGCONFIG_LIBS)
    CONVERT_LIBSTYLE_TO_PKGCONFIG(${_package}_PKGCONFIG_LIBS_PRIVATE)
    STRING(REPLACE ";" " " ${_package}_PKGCONFIG_LIBS "${${_package}_PKGCONFIG_LIBS}")
    STRING(REPLACE ";" " " ${_package}_PKGCONFIG_LIBS_PRIVATE "${${_package}_PKGCONFIG_LIBS_PRIVATE}")
    STRING(REPLACE ";" " " ${_package}_PKGCONFIG_REQUIRED "${${_package}_PKGCONFIG_REQUIRED}")
    STRING(REPLACE ";" " " ${_package}_PKGCONFIG_REQUIRED_PRIVATE "${${_package}_PKGCONFIG_REQUIRED_PRIVATE}")
ENDMACRO(CLEAN_LIB_LIST)

###
#
# GENERATE_PKGCONFIG_FILE: generate files chameleon.pc, coreblas.pc and cudablas.pc
#
###
MACRO(GENERATE_PKGCONFIG_FILE)

    # The definitions that should be given to users (change the API)
    set(CHAMELEON_PKGCONFIG_DEFINITIONS "")
    set(COREBLAS_PKGCONFIG_DEFINITIONS "")
    set(CUDABLAS_PKGCONFIG_DEFINITIONS "")

    # The link flags specific to this package and any required libraries
    # that don't support PkgConfig
    set(CHAMELEON_PKGCONFIG_LIBS "-lchameleon")
    set(COREBLAS_PKGCONFIG_LIBS  "-lcoreblas")
    set(CUDABLAS_PKGCONFIG_LIBS  "-lcudablas")

    # The link flags for private libraries required by this package but not
    # exposed to applications
    set(CHAMELEON_PKGCONFIG_LIBS_PRIVATE "")
    set(COREBLAS_PKGCONFIG_LIBS_PRIVATE  "")
    set(CUDABLAS_PKGCONFIG_LIBS_PRIVATE  "")

    # A list of packages required by this package
    set(CHAMELEON_PKGCONFIG_REQUIRED "hqr")
    set(COREBLAS_PKGCONFIG_REQUIRED  "")
    set(CUDABLAS_PKGCONFIG_REQUIRED  "")

    # A list of private packages required by this package but not exposed to
    # applications
    set(CHAMELEON_PKGCONFIG_REQUIRED_PRIVATE "")
    set(COREBLAS_PKGCONFIG_REQUIRED_PRIVATE  "")
    set(CUDABLAS_PKGCONFIG_REQUIRED_PRIVATE  "")

    if(CHAMELEON_SCHED_OPENMP)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon_openmp)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS_PRIVATE "${OpenMP_C_LIBRARIES}")
    endif()
    if(CHAMELEON_SCHED_PARSEC)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon_parsec)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS_PRIVATE "${PARSEC_LIBRARIES_DEP}")
    endif()
    if(CHAMELEON_SCHED_QUARK)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon_quark)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS_PRIVATE "${QUARK_LIBRARIES_DEP}")
    endif()
    if(CHAMELEON_SCHED_STARPU)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon_starpu)
        if ( CHAMELEON_USE_MPI )
            list(APPEND CHAMELEON_PKGCONFIG_REQUIRED_PRIVATE starpumpi-${STARPU_VERSION})
        else()
            list(APPEND CHAMELEON_PKGCONFIG_REQUIRED_PRIVATE starpu-${STARPU_VERSION})
        endif()
    endif()

    if(NOT CHAMELEON_SIMULATION)

        list(APPEND COREBLAS_PKGCONFIG_LIBS_PRIVATE
        ${LAPACKE_LIBRARIES_DEP}
        ${CBLAS_LIBRARIES_DEP}
        )
        list(APPEND CHAMELEON_PKGCONFIG_LIBS_PRIVATE
        ${EXTRA_LIBRARIES}
        )
        list(APPEND CHAMELEON_PKGCONFIG_REQUIRED "coreblas")

        if(CHAMELEON_USE_CUDA)
            list(APPEND CUDABLAS_PKGCONFIG_LIBS_PRIVATE ${EXTRA_LIBRARIES_CUDA})
            list(APPEND CHAMELEON_PKGCONFIG_REQUIRED "cudablas")
        endif()

    else(NOT CHAMELEON_SIMULATION)

        if(CHAMELEON_USE_CUDA)
            list(APPEND CHAMELEON_PKGCONFIG_LIBS -lcudablas)
        endif()
        list(APPEND CHAMELEON_PKGCONFIG_LIBS
        -lcoreblas
        ${EXTRA_LIBRARIES}
        )

    endif(NOT CHAMELEON_SIMULATION)

    # Define required package
    # -----------------------
    CLEAN_LIB_LIST(CHAMELEON)
    CLEAN_LIB_LIST(COREBLAS)
    if(CHAMELEON_USE_CUDA)
        CLEAN_LIB_LIST(CUDABLAS)
    endif()

    # Create .pc file
    # ---------------
    SET(_output_chameleon_file "${CMAKE_BINARY_DIR}/chameleon.pc")
    SET(_output_coreblas_file "${CMAKE_BINARY_DIR}/coreblas.pc")
    if(CHAMELEON_USE_CUDA)
        SET(_output_cudablas_file "${CMAKE_BINARY_DIR}/cudablas.pc")
    endif()

    # TODO: add url of CHAMELEON releases in .pc file
    CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/lib/pkgconfig/chameleon.pc.in" "${_output_chameleon_file}" @ONLY)
    CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/lib/pkgconfig/coreblas.pc.in"  "${_output_coreblas_file}" @ONLY)
    if(CHAMELEON_USE_CUDA)
        CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/lib/pkgconfig/cudablas.pc.in"  "${_output_cudablas_file}" @ONLY)
    endif()

    # installation
    # ------------
    INSTALL(FILES ${_output_chameleon_file} DESTINATION lib/pkgconfig)
    INSTALL(FILES ${_output_coreblas_file}  DESTINATION lib/pkgconfig)
    INSTALL(FILES ${_output_cudablas_file}  DESTINATION lib/pkgconfig)

ENDMACRO(GENERATE_PKGCONFIG_FILE)

##
## @end file GenPkgConfig.cmake
##

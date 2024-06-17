###
#
# @file GenPkgConfig.cmake
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
# @version 1.2.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2022-02-22
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
# GENERATE_PKGCONFIG_FILE: generate files chameleon_lapack.pc, chameleon.pc, coreblas.pc and gpucublas.pc or gpuhipblas.pc
#
###
MACRO(GENERATE_PKGCONFIG_FILE)

    # The definitions that should be given to users (change the API)
    set(CHAMELEON_LAPACK_PKGCONFIG_DEFINITIONS "")
    set(CHAMELEON_PKGCONFIG_DEFINITIONS "")
    set(COREBLAS_PKGCONFIG_DEFINITIONS "")
    set(GPUCUBLAS_PKGCONFIG_DEFINITIONS "")
    set(GPUHIPBLAS_PKGCONFIG_DEFINITIONS "")

    # The link flags specific to this package and any required libraries
    # that don't support PkgConfig
    set(CHAMELEON_LAPACK_PKGCONFIG_LIBS "-lchameleon_lapack")
    set(CHAMELEON_PKGCONFIG_LIBS "-lchameleon")
    set(COREBLAS_PKGCONFIG_LIBS  "-lcoreblas")
    set(GPUCUBLAS_PKGCONFIG_LIBS  "-lgpucublas")
    set(GPUHIPBLAS_PKGCONFIG_LIBS  "-lgpuhipblas")

    # The link flags for private libraries required by this package but not
    # exposed to applications
    set(CHAMELEON_LAPACK_PKGCONFIG_LIBS_PRIVATE "")
    set(CHAMELEON_PKGCONFIG_LIBS_PRIVATE "")
    set(COREBLAS_PKGCONFIG_LIBS_PRIVATE  "")
    set(GPUCUBLAS_PKGCONFIG_LIBS_PRIVATE  "")
    set(GPUHIPBLAS_PKGCONFIG_LIBS_PRIVATE  "")

    # A list of packages required by this package
    set(CHAMELEON_LAPACK_PKGCONFIG_REQUIRED "chameleon")
    set(CHAMELEON_PKGCONFIG_REQUIRED "hqr")
    set(COREBLAS_PKGCONFIG_REQUIRED  "")
    set(GPUCUBLAS_PKGCONFIG_REQUIRED  "")
    set(GPUHIPBLAS_PKGCONFIG_REQUIRED  "")

    # A list of private packages required by this package but not exposed to
    # applications
    set(CHAMELEON_LAPACK_PKGCONFIG_REQUIRED_PRIVATE "")
    set(CHAMELEON_PKGCONFIG_REQUIRED_PRIVATE "")
    set(COREBLAS_PKGCONFIG_REQUIRED_PRIVATE  "")
    set(GPUCUBLAS_PKGCONFIG_REQUIRED_PRIVATE  "")
    set(GPUHIPBLAS_PKGCONFIG_REQUIRED_PRIVATE  "")

    if(CHAMELEON_SCHED_OPENMP)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon_openmp)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS_PRIVATE "${OpenMP_C_LIBRARIES}")
    elseif(CHAMELEON_SCHED_PARSEC)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon_parsec)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS_PRIVATE "${PARSEC_LIBRARIES}")
        elseif(CHAMELEON_SCHED_QUARK)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon_quark)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS_PRIVATE "${QUARK_LIBRARIES_DEP}")
    elseif(CHAMELEON_SCHED_STARPU)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon_starpu)
        if ( CHAMELEON_SIMULATION )
            # need PUBLIC here because simulation mode with starpu requires to replace the
            # main function by a starpu main one
            if ( CHAMELEON_USE_MPI )
                list(APPEND CHAMELEON_PKGCONFIG_REQUIRED starpumpi-${CHAMELEON_STARPU_VERSION})
            else()
                list(APPEND CHAMELEON_PKGCONFIG_REQUIRED starpu-${CHAMELEON_STARPU_VERSION})
            endif()
        else()
            if ( CHAMELEON_USE_MPI )
                list(APPEND CHAMELEON_PKGCONFIG_REQUIRED_PRIVATE starpumpi-${CHAMELEON_STARPU_VERSION})
            else()
                list(APPEND CHAMELEON_PKGCONFIG_REQUIRED_PRIVATE starpu-${CHAMELEON_STARPU_VERSION})
            endif()
        endif()
    endif()

    if(NOT CHAMELEON_SIMULATION)

        list(APPEND COREBLAS_PKGCONFIG_LIBS_PRIVATE
        ${LAPACKE_LIBRARIES}
        ${CBLAS_LIBRARIES}
        )
        list(APPEND CHAMELEON_PKGCONFIG_REQUIRED "coreblas")

        if(CHAMELEON_USE_CUDA)
            list(APPEND GPUCUBLAS_PKGCONFIG_LIBS_PRIVATE ${CUDA_CUBLAS_LIBRARIES})
            list(APPEND GPUCUBLAS_PKGCONFIG_REQUIRED "cuda")
            list(APPEND CHAMELEON_PKGCONFIG_REQUIRED "gpucublas")
        endif()

        if(CHAMELEON_USE_HIP)
            list(APPEND GPUHIPBLAS_PKGCONFIG_LIBS_PRIVATE ${HIPBLAS_LIBRARIES})
            list(APPEND GPUHIPBLAS_PKGCONFIG_LIBS_PRIVATE ${HIP_LIBRARIES})
            list(APPEND CHAMELEON_PKGCONFIG_REQUIRED "gpuhipblas")
        endif()

    endif(NOT CHAMELEON_SIMULATION)

    list(APPEND CHAMELEON_PKGCONFIG_LIBS_PRIVATE
    ${M_LIBRARIES}
    )

    if(CHAMELEON_USE_MPI)
        list(APPEND CHAMELEON_PKGCONFIG_REQUIRED "mpi")
    endif()

    # Define required package
    # -----------------------
    CLEAN_LIB_LIST(CHAMELEON)
    CLEAN_LIB_LIST(COREBLAS)
    if(CHAMELEON_USE_CUDA)
        CLEAN_LIB_LIST(GPUCUBLAS)
    endif()
    if(CHAMELEON_USE_HIP)
        CLEAN_LIB_LIST(GPUHIPBLAS)
    endif()

    # Create .pc file
    # ---------------
    SET(_output_chameleon_lapack_file "${CMAKE_BINARY_DIR}/chameleon_lapack.pc")
    SET(_output_chameleon_file "${CMAKE_BINARY_DIR}/chameleon.pc")
    SET(_output_coreblas_file "${CMAKE_BINARY_DIR}/coreblas.pc")
    if(CHAMELEON_USE_CUDA)
        SET(_output_gpucublas_file "${CMAKE_BINARY_DIR}/gpucublas.pc")
    endif()
    if(CHAMELEON_USE_HIP)
        SET(_output_gpuhipblas_file "${CMAKE_BINARY_DIR}/gpuhipblas.pc")
    endif()

    # TODO: add url of CHAMELEON releases in .pc file
    CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/lib/pkgconfig/chameleon_lapack.pc.in" "${_output_chameleon_lapack_file}" @ONLY)
    CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/lib/pkgconfig/chameleon.pc.in" "${_output_chameleon_file}" @ONLY)
    CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/lib/pkgconfig/coreblas.pc.in" "${_output_coreblas_file}" @ONLY)
    if(CHAMELEON_USE_CUDA)
        CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/lib/pkgconfig/gpucublas.pc.in" "${_output_gpucublas_file}" @ONLY)
    endif()
    if(CHAMELEON_USE_HIP)
        CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/lib/pkgconfig/gpuhipblas.pc.in" "${_output_gpuhipblas_file}" @ONLY)
    endif()

    # installation
    # ------------
    INSTALL(FILES ${_output_chameleon_lapack_file} DESTINATION lib/pkgconfig)
    INSTALL(FILES ${_output_chameleon_file} DESTINATION lib/pkgconfig)
    INSTALL(FILES ${_output_coreblas_file} DESTINATION lib/pkgconfig)
    INSTALL(FILES ${_output_gpucublas_file} DESTINATION lib/pkgconfig)
    INSTALL(FILES ${_output_gpuhipblas_file} DESTINATION lib/pkgconfig)

ENDMACRO(GENERATE_PKGCONFIG_FILE)

##
## @end file GenPkgConfig.cmake
##

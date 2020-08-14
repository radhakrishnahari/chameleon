###
#
# @file quark/CMakeLists.txt
#
# @copyright 2009-2015 The University of Tennessee and The University of
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
cmake_minimum_required(VERSION 2.8)

find_package( QUARK REQUIRED COMPONENTS HWLOC )

# Add definition and include_dir if found
if(QUARK_FOUND)
  message("-- ${Blue}Add definition CHAMELEON_SCHED_QUARK"
    " - Activate QUARK in Chameleon${ColourReset}")
  set(CHAMELEON_SCHED "${CHAMELEON_SCHED};QUARK" CACHE STRING "Runtime systems available" FORCE)
  if (QUARK_INCLUDE_DIRS_DEP)
    include_directories(${QUARK_INCLUDE_DIRS_DEP})
  endif()
  if(QUARK_CFLAGS_OTHER_DEP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${QUARK_CFLAGS_OTHER_DEP}")
    string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  endif()
  if(QUARK_LDFLAGS_OTHER_DEP)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${QUARK_LDFLAGS_OTHER_DEP}")
    string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
  endif()
  if(QUARK_LIBRARY_DIRS_DEP)
    list(APPEND CMAKE_INSTALL_RPATH "${QUARK_LIBRARY_DIRS_DEP}")
  endif()
  if (QUARK_LIBRARIES_DEP)
    set(CHAMELEON_DEP  "${QUARK_LIBRARIES_DEP};${CHAMELEON_DEP}" PARENT_SCOPE)
  endif()
else(QUARK_FOUND)
  if(CHAMELEON_VERBOSE_FIND_PACKAGE)
    if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
      if (NOT HWLOC_hwloc.h_DIRS)
        Print_Find_Header_Status(hwloc hwloc.h)
      endif ()
      if (NOT HWLOC_hwloc_LIBRARY)
        Print_Find_Library_Status(hwloc libhwloc)
      endif ()
    endif()
    if (NOT QUARK_quark.h_DIRS)
      Print_Find_Header_Status(quark quark.h)
    endif ()
    if (NOT QUARK_quark_LIBRARY)
      Print_Find_Library_Status(quark libquark)
    endif ()
  else()
    message(WARNING "QUARK library has not been found and CHAMELEON_VERBOSE_FIND_PACKAGE is set to OFF."
      " Try to activate CHAMELEON_VERBOSE_FIND_PACKAGE option (-DCHAMELEON_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
  endif()
  if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
    message(FATAL_ERROR "hwloc library is required but has not been found")
  endif()
  message(FATAL_ERROR "QUARK library is required but has not been found")
  return()
endif(QUARK_FOUND)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )
include_directories( ${CMAKE_CURRENT_BINARY_DIR}/include )

# # Generate headers for all possible precisions
# # --------------------------------------------
# set(RUNTIME_HDRS_GENERATED "")
# set(ZHDR
# )
#
# precisions_rules_py(RUNTIME_HDRS_GENERATED "${ZHDR}"
#   PRECISIONS "s;d;c;z;ds;zc"
#   TARGETDIR  "include")

# Define the list of headers
# --------------------------
set(RUNTIME_HDRS
  include/chameleon_quark.h
  include/core_blas_dag.h
  )

# # Add generated headers
# # ---------------------
# foreach( hdr_file ${RUNTIME_HDRS_GENERATED} )
#   list(APPEND RUNTIME_HDRS ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
# endforeach()

# Force generation of headers
# ---------------------------
add_custom_target(
  runtime_quark_include
  ALL SOURCES ${RUNTIME_HDRS})

# Installation
# ------------
install(
  FILES ${RUNTIME_HDRS}
  DESTINATION include/runtime/quark )

# Generate the Chameleon common for all possible precisions
# ---------------------------------------------------------
set(RUNTIME_COMMON_GENERATED "")
set(ZSRC
  control/runtime_zprofiling.c
  control/runtime_zlocality.c
  )

precisions_rules_py(
  RUNTIME_COMMON_GENERATED "${ZSRC}"
  PRECISIONS "${CHAMELEON_PRECISION}"
  TARGETDIR "control")

set(RUNTIME_COMMON
  control/runtime_async.c
  control/runtime_context.c
  control/runtime_control.c
  control/runtime_descriptor.c
  control/runtime_options.c
  control/runtime_profiling.c
  ${RUNTIME_COMMON_GENERATED}
  )

# Generate the Chameleon sources for all possible precisions
# ----------------------------------------------------------
set(RUNTIME_SRCS_GENERATED "")
set(ZSRC
  ${CODELETS_ZSRC}
  )

precisions_rules_py(RUNTIME_SRCS_GENERATED "${ZSRC}"
  PRECISIONS "${CHAMELEON_PRECISION}"
  TARGETDIR "codelets")

set(RUNTIME_SRCS
  ${RUNTIME_COMMON}
  ${RUNTIME_SRCS_GENERATED}
  ${CODELETS_SRC}
  )

# Force generation of sources
# ---------------------------
add_custom_target(quark_sources ALL SOURCES ${RUNTIME_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};runtime_quark_include;quark_sources" CACHE INTERNAL "List of targets of sources")

# Add library
# -----------
add_library(chameleon_quark ${RUNTIME_SRCS})
set_property(TARGET chameleon_quark PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET chameleon_quark PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

target_link_libraries(chameleon_quark coreblas ${QUARK_LIBRARIES_DEP})
if(CHAMELEON_USE_CUDA)
  target_link_libraries(chameleon_quark cudablas)
endif(CHAMELEON_USE_CUDA)

add_dependencies(chameleon_quark
  chameleon_include
  control_include
  runtime_quark_include
  quark_sources
)
if (NOT CHAMELEON_SIMULATION)
    add_dependencies(chameleon_quark coreblas_include)
endif()

# installation
# ------------
install(TARGETS chameleon_quark
  DESTINATION lib)

###
### END CMakeLists.txt
###

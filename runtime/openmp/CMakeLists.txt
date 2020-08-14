###
#
# @file openmp/CMakeLists.txt
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

find_package(OpenMP REQUIRED)

if ( OPENMP_FOUND )
  message("-- ${Blue}Add definition CHAMELEON_SCHED_OPENMP"
    " - Activate OpenMP in Chameleon${ColourReset}")
  include_directories(${OpenMP_C_INCLUDE_DIRS})
  if (OpenMP_C_LIBRARIES)
    set(CHAMELEON_DEP  "${OpenMP_C_LIBRARIES};${CHAMELEON_DEP}" PARENT_SCOPE)
  endif()
else ( OPENMP_FOUND )
  message(FATAL_ERROR "Something went wrong when finding OpenMP")
  return()
endif ( OPENMP_FOUND )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )
include_directories( ${CMAKE_CURRENT_BINARY_DIR}/include )

# Define the list of headers
# --------------------------
set(RUNTIME_HDRS
  include/chameleon_openmp.h
  )

# Force generation of headers
# ---------------------------
add_custom_target(
  runtime_openmp_include
  ALL SOURCES ${RUNTIME_HDRS})

# Installation
# ------------
install(
  FILES ${RUNTIME_HDRS}
  DESTINATION include/runtime/openmp )

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
add_custom_target(openmp_sources ALL SOURCES ${RUNTIME_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};runtime_openmp_include;openmp_sources" CACHE INTERNAL "List of targets of sources")

# Add library
# -----------
add_library(chameleon_openmp ${RUNTIME_SRCS})
set_property(TARGET chameleon_openmp PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET chameleon_openmp PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

target_link_libraries(chameleon_openmp ${OpenMP_C_LIBRARIES})
if (NOT CHAMELEON_SIMULATION)
  target_link_libraries(chameleon_openmp coreblas)
endif()

add_dependencies(chameleon_openmp
  chameleon_include
  control_include
  runtime_openmp_include
  openmp_sources
  )

add_dependencies(chameleon_openmp coreblas_include)

# installation
# ------------
install(TARGETS chameleon_openmp
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

###
### END CMakeLists.txt
###

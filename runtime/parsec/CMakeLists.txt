###
#
# @file parsec/CMakeLists.txt
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

# create list of components in order to make a single call to find_package(starpu...)
set(PARSEC_COMPONENT_LIST "HWLOC")
if(CHAMELEON_USE_MPI)
  list(APPEND PARSEC_COMPONENT_LIST "MPI")
endif()
if(CHAMELEON_USE_CUDA)
  list(APPEND PARSEC_COMPONENT_LIST "CUDA")
endif()

find_package( PARSEC REQUIRED
  COMPONENTS ${PARSEC_COMPONENT_LIST} )

if(PARSEC_FOUND)
  message("-- ${Blue}Add definition CHAMELEON_SCHED_PARSEC"
    " - Activate PaRSEC in Chameleon${ColourReset}")
  set(CHAMELEON_SCHED "${CHAMELEON_SCHED};PARSEC" CACHE STRING "Runtime systems available" FORCE)
  if (PARSEC_INCLUDE_DIRS_DEP)
    include_directories(${PARSEC_INCLUDE_DIRS_DEP})
  endif()
  if(PARSEC_CFLAGS_OTHER_DEP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PARSEC_CFLAGS_OTHER_DEP}")
    string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  endif()
  if(PARSEC_LDFLAGS_OTHER_DEP)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${PARSEC_LDFLAGS_OTHER_DEP}")
    string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
  endif()
  if(PARSEC_LIBRARY_DIRS_DEP)
    list(APPEND CMAKE_INSTALL_RPATH "${PARSEC_LIBRARY_DIRS_DEP}")
  endif()
  if (PARSEC_LIBRARIES_DEP)
    set(CHAMELEON_DEP  "${PARSEC_LIBRARIES_DEP};${CHAMELEON_DEP}" PARENT_SCOPE)
  endif()
  string(REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
  # to test if symbols exist, e.g.
  # set(CMAKE_REQUIRED_INCLUDES "${PARSEC_INCLUDE_DIRS_DEP}")
  # set(CMAKE_REQUIRED_FLAGS "${PARSEC_CFLAGS_OTHER_DEP}")
  # set(CMAKE_REQUIRED_LDFLAGS "${PARSEC_LDFLAGS_OTHER_DEP}")
  # set(CMAKE_REQUIRED_LIBRARIES "${PARSEC_LIBRARIES_DEP}")
  # check_function_exists(parsec_init HAVE_PARSEC_INIT)
else ( PARSEC_FOUND )
  if(CHAMELEON_VERBOSE_FIND_PACKAGE)
    if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
      if (NOT HWLOC_hwloc.h_DIRS)
        Print_Find_Header_Status(hwloc hwloc.h)
      endif ()
      if (NOT HWLOC_hwloc_LIBRARY)
        Print_Find_Library_Status(hwloc libhwloc)
      endif ()
    endif()
    if( (NOT PARSEC_FOUND) OR (NOT PARSEC_LIBRARIES_DEP) )
      foreach(parsec_hdr ${PARSEC_hdrs_to_find})
        if (NOT PARSEC_${parsec_hdr}_INCLUDE_DIRS)
          Print_Find_Header_Status(parsec ${parsec_hdr})
        endif ()
      endforeach()
      if(PARSEC_VERSION_STRING)
        foreach(parsec_lib ${PARSEC_libs_to_find})
          if (NOT PARSEC_${parsec_lib}_LIBRARY)
            Print_Find_Library_Status(parsec ${parsec_lib})
          endif ()
        endforeach()
      endif()
    endif()
  else(CHAMELEON_VERBOSE_FIND_PACKAGE)
    message(WARNING "PaRSEC library has not been found and CHAMELEON_VERBOSE_FIND_PACKAGE is set to OFF."
      " Try to activate CHAMELEON_VERBOSE_FIND_PACKAGE option (-DCHAMELEON_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
  endif(CHAMELEON_VERBOSE_FIND_PACKAGE)
  if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
    message(FATAL_ERROR "hwloc library is required but has not been found")
  endif()
  if( (NOT PARSEC_FOUND) OR (NOT PARSEC_LIBRARIES_DEP) )
    message(FATAL_ERROR "PaRSEC library is required but has not been found")
  endif()
  return()
endif ( PARSEC_FOUND )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )
include_directories( ${CMAKE_CURRENT_BINARY_DIR}/include )

# Generate headers for all possible precisions
# --------------------------------------------
set(RUNTIME_HDRS_GENERATED "")
#set(ZHDR
#)
#
#precisions_rules_py(RUNTIME_HDRS_GENERATED "${ZHDR}"
#                    PRECISIONS "s;d;c;z;ds;zc"
#                    TARGETDIR  "include")

# Define the list of headers
# --------------------------
set(RUNTIME_HDRS
  include/chameleon_parsec.h
  )

# Add generated headers
# ---------------------
foreach( hdr_file ${RUNTIME_HDRS_GENERATED} )
  list(APPEND RUNTIME_HDRS ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()

# Force generation of headers
# ---------------------------
add_custom_target(
  runtime_parsec_include
  ALL SOURCES ${RUNTIME_HDRS})

# Installation
# ------------
install(
  FILES ${RUNTIME_HDRS}
  DESTINATION include/runtime/parsec )

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
add_custom_target(parsec_sources ALL SOURCES ${RUNTIME_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};runtime_parsec_include;parsec_sources" CACHE INTERNAL "List of targets of sources")

# Add library
# -----------
add_library(chameleon_parsec ${RUNTIME_SRCS})
set_property(TARGET chameleon_parsec PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET chameleon_parsec PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

target_link_libraries(chameleon_parsec
  ${PARSEC_LIBRARIES_DEP})
if (NOT CHAMELEON_SIMULATION)
  target_link_libraries(chameleon_parsec coreblas)
  if(CHAMELEON_USE_CUDA)
    target_link_libraries(chameleon_parsec cudablas)
  endif(CHAMELEON_USE_CUDA)
endif(NOT CHAMELEON_SIMULATION)

add_dependencies(chameleon_parsec
  chameleon_include
  control_include
  runtime_parsec_include
  parsec_sources
)

if (NOT CHAMELEON_SIMULATION)
  add_dependencies(chameleon_parsec coreblas_include)
  if (CHAMELEON_USE_CUDA)
    add_dependencies(chameleon_parsec cudablas_include)
  endif()
endif()

# installation
# ------------
install(TARGETS chameleon_parsec
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

###
### END CMakeLists.txt
###

###
#
#  @copyright 2013-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.3.0
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2024-06-13
#
###
cmake_minimum_required(VERSION 3.5)

project(LINK_CHAMELEON C CXX Fortran)

# look for CHAMELEON on the system
# Hint: use CMAKE_PREFIX_PATH (env. var. or cmake var.) to the installation
# directory of CHAMELEON if not installed in a standard path
find_package(CHAMELEON REQUIRED)

if (TARGET CHAMELEON::chameleon)
  message(STATUS "CHAMELEON_BIN_DIR ${CHAMELEON_BIN_DIR}")
  message(STATUS "CHAMELEON_INC_DIR ${CHAMELEON_INC_DIR}")
  message(STATUS "CHAMELEON_LIB_DIR ${CHAMELEON_LIB_DIR}")
else()
  message(FATAL_ERROR "Target CHAMELEON::chameleon is not found, check your CHAMELEONConfig.cmake.")
endif()

set(CMAKE_BUILD_RPATH ${CHAMELEON_LIB_DIR})

add_executable(link_chameleon_c ${CMAKE_CURRENT_SOURCE_DIR}/../../example/link_chameleon/link_chameleon.c)
target_link_libraries(link_chameleon_c PRIVATE CHAMELEON::chameleon)

# Launch executables with ctest
enable_testing()
include(CTest)
add_test(link_chameleon_c link_chameleon_c)

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
#  @version 1.3.0
#  @author Florent Pruvost
#  @date 2025-07-23
#
###
cmake_minimum_required(VERSION 3.12)
project(LINK_CHAMELEON Fortran C CXX)

option(CHAMELEON_PREC_S "Enable single precision example" ON)
option(CHAMELEON_PREC_D "Enable double precision example" ON)
option(CHAMELEON_PREC_C "Enable complex precision example" ON)
option(FORTRAN "Enable Fortran example" OFF)
option(CPP "Enable C++ example" OFF)
option(MPI "Set to ON if Chameleon is configured with MPI" OFF)

find_package(CHAMELEON REQUIRED)

if (TARGET CHAMELEON::chameleon)
  message(STATUS "CHAMELEON_BIN_DIR ${CHAMELEON_BIN_DIR}")
  message(STATUS "CHAMELEON_INC_DIR ${CHAMELEON_INC_DIR}")
  message(STATUS "CHAMELEON_LIB_DIR ${CHAMELEON_LIB_DIR}")
else()
  message(FATAL_ERROR "Target CHAMELEON::chameleon is not found, check your CHAMELEONConfig.cmake.")
endif()

set(CMAKE_BUILD_RPATH ${CHAMELEON_LIB_DIR})
enable_testing()
include(CTest)

if (CHAMELEON_PREC_S)
  add_executable(link_chameleon_c link_chameleon.c)
  target_link_libraries(link_chameleon_c PRIVATE CHAMELEON::chameleon)
  add_test(link_chameleon_c link_chameleon_c)
endif()
if (CHAMELEON_PREC_D AND FORTRAN)
  add_executable(link_chameleon_f link_chameleon.f90)
  target_link_libraries(link_chameleon_f PRIVATE CHAMELEON::chameleon)
  add_test(link_chameleon_f link_chameleon_f)
endif()
if (CHAMELEON_PREC_C AND CPP)
  add_executable(link_chameleon_cpp link_chameleon.cpp)
  if (MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    target_link_libraries(link_chameleon_cpp PRIVATE CHAMELEON::chameleon MPI::MPI_CXX)
  else()
    target_link_libraries(link_chameleon_cpp PRIVATE CHAMELEON::chameleon)
  endif()
  target_compile_definitions(link_chameleon_cpp PRIVATE CHAMELEON_COMPLEX_CPP)
  add_test(link_chameleon_cpp link_chameleon_cpp)
endif()
###
### END CMakeLists.txt
###

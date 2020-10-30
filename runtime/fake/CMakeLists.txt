###
#
# @file fake/CMakeLists.txt
#
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
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2020-10-28
#
###
cmake_minimum_required(VERSION 2.8)

set(RUNTIME_SRCS "runtime_fake.c")

# Add library
# -----------
add_library(chameleon_fake ${RUNTIME_SRCS})
set_property(TARGET chameleon_fake PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET chameleon_fake PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

if (NOT CHAMELEON_SIMULATION)
  target_link_libraries(chameleon_fake coreblas)
endif()

add_dependencies(chameleon_fake
  chameleon_include
  control_include
  coreblas_include
  )

# installation
# ------------
install(TARGETS chameleon_fake
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

###
### END CMakeLists.txt
###

###
#
# @file CMakeLists.txt
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
# @version 1.3.0
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2023-07-06
#
###

# Generate header files
# ---------------------
set(GPUCUBLAS_HDRS_GENERATED "")
set(ZHDR
    gpucublas/gpucublas_z.h
    gpucublas/gpucublas_zc.h
)
precisions_rules_py(
  GPUCUBLAS_HDRS_GENERATED "${ZHDR}"
  TARGETDIR gpucublas
  PRECISIONS "s;d;c;z;zc;ds" )

# Define the list of headers
# --------------------------
configure_file("gpucublas.h.in"
               "gpucublas.h"
               @ONLY)

set(GPUCUBLAS_HDRS
  ${CMAKE_CURRENT_BINARY_DIR}/gpucublas.h
)

# Add generated headers
# ---------------------
foreach( hdr_file ${GPUCUBLAS_HDRS_GENERATED} )
  list(APPEND GPUCUBLAS_HDRS ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()

# Force generation of headers
# ---------------------------
add_custom_target(gpucublas_include ALL SOURCES ${GPUCUBLAS_HDRS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};gpucublas_include" CACHE INTERNAL "List of targets of sources")

# Installation
# ------------
install( FILES ${CMAKE_CURRENT_BINARY_DIR}/gpucublas.h
         DESTINATION include )

install( FILES ${GPUCUBLAS_HDRS}
         DESTINATION include/gpucublas )

###
### END CMakeLists.txt
###

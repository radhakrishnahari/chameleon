###
#
# @file CMakeLists.txt
#
# @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
#  @author Abel Calluaud
#  @date 2024-04-27
#
###

set(CHAMELEON_STARPU_VERSION "1.3" CACHE STRING "necessary STARPU API version")

find_package(STARPU ${CHAMELEON_STARPU_VERSION} REQUIRED)
morse_export_imported_target(MORSE STARPU starpu chameleon)


# Set list of starpu_rapack sources
set(
    CHAMELEON_STARPU_RAPACK_SRCS
    starpu/control/runtime_lrmem.c
    starpu/control/runtime_rpk_ctx.c
)

set(
    CHAMELEON_STARPU_RAPACK_ZSRCS
    starpu/control/runtime_rpk_zctx.c
)

# Set list of starpu_rapack headers
set(
	CHAMELEON_STARPU_RAPACK_HDRS
	include/runtime_rpk.h
	include/runtime_lrmem.h
	include/cham_tile_interface.h
)

set(
    CHAMELEON_STARPU_RAPACK_ZHDRS
    include/runtime_rpk_z.h
)

if(CHAMELEON_USE_HMATOSS)
    set(
        CHAMELEON_STARPU_RAPACK_SRCS
        ${CHAMELEON_STARPU_RAPACK_SRCS}
        starpu/control/runtime_hmat_interface.c
    )
    set(
        CHAMELEON_STARPU_RAPACK_ZSRCS
        ${CHAMELEON_STARPU_RAPACK_ZSRCS}
        starpu/control/runtime_hmat_zinterface.c
    )
    set(
        CHAMELEON_STARPU_RAPACK_HDRS
        ${CHAMELEON_STARPU_RAPACK_HDRS}
        include/runtime_hmat.h
    )
    set(
        CHAMELEON_STARPU_RAPACK_ZHDRS
        ${CHAMELEON_STARPU_RAPACK_ZHDRS}
        include/runtime_hmat_z.h
    )
endif()

message(STATUS "rapack srcs ${CHAMELEON_STARPU_RAPACK_SRCS} ${CHAMELEON_STARPU_RAPACK_ZSRCS}")

# Generate precision variants sources
precisions_rules_py(
    CHAMELEON_STARPU_RAPACK_ZSRCS_GENERATED
    "${CHAMELEON_STARPU_RAPACK_ZSRCS}"
    PRECISIONS "${CHAMELEON_PRECISION}"
)

# Generate precision variants headers
precisions_rules_py(
    CHAMELEON_STARPU_RAPACK_ZHDRS_GENERATED
    "${CHAMELEON_STARPU_RAPACK_ZHDRS}"
    TARGETDIR "include/runtime_rpk"
    PRECISIONS "s;d;c;z;zc;ds"
)

# Add a custom target for starpu_rapack headers
add_custom_target(
    chameleon_starpu_rapack_include
	ALL SOURCES ${CHAMELEON_STARPU_RAPACK_HDRS} ${CHAMELEON_STARPU_RAPACK_ZHDRS_GENERATED}
)

# Add a custom target for starpu_rapack sources
add_custom_target(
    chameleon_starpu_rapack_sources
    ALL SOURCES ${CHAMELEON_STARPU_RAPACK_SRCS} ${CHAMELEON_STARPU_RAPACK_ZSRCS_GENERATED}
)

# Add the starpu_rapack library
add_library(
    chameleon_starpu_rapack
    ${CHAMELEON_STARPU_RAPACK_SRCS}
    ${CHAMELEON_STARPU_RAPACK_ZSRCS_GENERATED}
)

add_dependencies(
    chameleon_starpu_rapack
    chameleon_starpu_rapack_include
    chameleon_starpu_rapack_sources
)

set_target_properties(
    chameleon_starpu_rapack
    PROPERTIES VERSION ${CHAMELEON_VERSION}
)

set_target_properties(
    chameleon_starpu_rapack
    PROPERTIES SOVERSION ${CHAMELEON_VERSION_MAJOR}
)

set_property(
    TARGET chameleon_starpu_rapack
    PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib"
)

target_include_directories(
    chameleon_starpu_rapack PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CHAMELEON_SOURCE_DIR}>/include
    $<BUILD_INTERFACE:${CHAMELEON_BINARY_DIR}>/include
)

add_dependencies(
    chameleon_starpu_rapack
    chameleon_include
    control_include
    runtime_starpu_include
)

target_link_libraries(
    chameleon_starpu_rapack
    PUBLIC
    MORSE::STARPU
    rapack
)

if (CHAMELEON_USE_HMATOSS)
    target_link_libraries(
        chameleon_starpu_rapack
        PUBLIC
        HMAT::hmat
    )
endif()

# Library installation

install(
    EXPORT chameleon_starpu_rapackTargets
    NAMESPACE CHAMELEON::
    DESTINATION lib/cmake/chameleon
)

install(
    TARGETS chameleon_starpu_rapack
    EXPORT chameleon_starpu_rapackTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

install(
    FILES include/runtime_rpk.h
    DESTINATION include
)

install(
    FILES ${CHAMELEON_STARPU_RAPACK_ZHDRS_GENERATED}
    DESTINATION include/runtime_rpk
)

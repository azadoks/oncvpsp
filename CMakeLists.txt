cmake_minimum_required(VERSION 3.22...3.31)

project(
    oncvpsp
    VERSION 4.0.1
    DESCRIPTION "Optimized norm-conserving Vanderbilt pseudopotentials"
    HOMEPAGE_URL "https://github.com/oncvpsp/oncvpsp"
    LANGUAGES Fortran C)

if (PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds are not supported; configure in a separate directory, e.g.\n"
        "  cmake -S . -B build")
endif ()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

##########################################################
# Options
##########################################################
option(ONCVPSP_WITH_Libxc "Enable Libxc support" ON)
option(ONCVPSP_WITH_TOML "Enable TOML input support" ON)
option(ONCVPSP_WITH_FR "Build the fully-relativistic executables (oncvpspr*.x)" ON)
option(ONCVPSP_BLAS_DOT "Pseudo-atom overlaps by the linked BLAS (faster; results depend on it at round-off)" OFF)
option(ONCVPSP_WITH_C_API "Build the shared library liboncvpsp with the C interface (include/oncvpsp.h)" ON)

# The shared library links the static one, and the dependencies built here, which must then be
# position-independent
if (ONCVPSP_WITH_C_API)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif ()

# Where each dependency comes from: AUTO uses the system's if it is found and works and otherwise
# downloads and builds it (with CPM), SYSTEM requires the system's and FETCH always builds it
set(ONCVPSP_DEPENDENCY_SOURCES AUTO SYSTEM FETCH)
set(ONCVPSP_LAPACK_SOURCE AUTO CACHE STRING "Where LAPACK and BLAS come from: AUTO, SYSTEM or FETCH")
set(ONCVPSP_Libxc_SOURCE AUTO CACHE STRING "Where Libxc comes from: AUTO, SYSTEM or FETCH")
set(ONCVPSP_TOML_SOURCE AUTO CACHE STRING "Where toml-f comes from: AUTO, SYSTEM or FETCH")
foreach (dep IN ITEMS LAPACK Libxc TOML)
    set_property(CACHE ONCVPSP_${dep}_SOURCE PROPERTY STRINGS ${ONCVPSP_DEPENDENCY_SOURCES})
    if (NOT ONCVPSP_${dep}_SOURCE IN_LIST ONCVPSP_DEPENDENCY_SOURCES)
        message(FATAL_ERROR "ONCVPSP_${dep}_SOURCE must be one of ${ONCVPSP_DEPENDENCY_SOURCES}, "
                            "not '${ONCVPSP_${dep}_SOURCE}'")
    endif ()
endforeach ()

# RelWithDebInfo (-O2 -g with gfortran) is the build the test references are checked against
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
    message(STATUS "No build type given, using ${CMAKE_BUILD_TYPE}")
endif ()

include(GNUInstallDirs)
# The installed executables keep the paths to the shared libraries they link outside the system
# directories (such as a Libxc under Libxc_ROOT)
if (NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
endif ()
if (PROJECT_IS_TOP_LEVEL)
    include(CTest)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
endif ()

###########################################################
# Dependencies
###########################################################
# The dependencies' option() commands take the values given to them in CPMAddPackage OPTIONS (as
# normal variables) even where their cmake_minimum_required predates that behaviour
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
include(CPM)
include(OncvpspCompilerOptions)
include(OncvpspLapack)
if (ONCVPSP_WITH_Libxc)
    include(OncvpspLibxc)
endif ()

# toml-f, for the TOML input
if (ONCVPSP_WITH_TOML)
    if (NOT ONCVPSP_TOML_SOURCE STREQUAL "FETCH")
        find_package(toml-f 0.4 CONFIG QUIET)
        if (toml-f_FOUND)
            set(ONCVPSP_TOML_DESCRIPTION "toml-f ${toml-f_VERSION} (${toml-f_DIR})")
        elseif (ONCVPSP_TOML_SOURCE STREQUAL "SYSTEM")
            message(FATAL_ERROR "toml-f not found; set toml-f_ROOT or CMAKE_PREFIX_PATH to its "
                                "installation, or ONCVPSP_TOML_SOURCE to AUTO or FETCH")
        endif ()
    endif ()
    if (NOT TARGET toml-f::toml-f)
        CPMAddPackage(
            NAME toml-f
            GITHUB_REPOSITORY toml-f/toml-f
            VERSION 0.4.3
            EXCLUDE_FROM_ALL YES
            SYSTEM YES
            OPTIONS "BUILD_TESTING OFF" "TOML_F_BUILD_TESTING OFF" "CMAKE_WARN_DEPRECATED OFF"
        )
        add_library(toml-f::toml-f ALIAS toml-f-lib)
        set(ONCVPSP_TOML_DESCRIPTION "toml-f 0.4.3 (built)")
    endif ()
endif ()

# xmlf90, whose wxml writes the PSML output. It is always built, with a patch that adds xml_OpenUnit:
# it writes a document to a unit that is already open (the standard output, where the PSML file goes
# like the psp8 and upf files) and leaves it open at xml_Close, which frees the document's element
# stack. The patch's hash is in the patch command, so that a changed patch is applied again
set(xmlf90_patch "${PROJECT_SOURCE_DIR}/cmake/patches/xmlf90-1.6.3.patch")
file(SHA256 "${xmlf90_patch}" xmlf90_patch_sha256)
CPMAddPackage(
    NAME xmlf90
    GIT_REPOSITORY "https://gitlab.com/siesta-project/libraries/xmlf90.git"
    GIT_TAG 1.6.3
    PATCH_COMMAND "${CMAKE_COMMAND}" "-DPATCH=${xmlf90_patch}" "-DPATCH_SHA256=${xmlf90_patch_sha256}"
                  -P "${PROJECT_SOURCE_DIR}/cmake/ApplyPatch.cmake"
    EXCLUDE_FROM_ALL YES
    SYSTEM YES
)

###########################################################
# ONCVPSP
###########################################################
add_subdirectory(src)
if (PROJECT_IS_TOP_LEVEL AND BUILD_TESTING)
    add_subdirectory(tests)
endif ()

# The executables, and the shared library with its header; the static library is internal
install(TARGETS ${ONCVPSP_EXECUTABLES} RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
if (ONCVPSP_WITH_C_API)
    install(TARGETS oncvpsp_c LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
    install(FILES include/oncvpsp.h DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
endif ()
# The Python package's (pyproject.toml, scikit-build-core): the shared library beside its modules
if (SKBUILD AND ONCVPSP_WITH_C_API)
    install(TARGETS oncvpsp_c LIBRARY DESTINATION oncvpsp NAMELINK_SKIP COMPONENT python)
endif ()

message(STATUS "")
message(STATUS "ONCVPSP ${PROJECT_VERSION}")
message(STATUS "  Build type:  ${CMAKE_BUILD_TYPE}")
message(STATUS "  Fortran:     ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "  LAPACK:      ${ONCVPSP_LAPACK_DESCRIPTION}")
if (ONCVPSP_WITH_Libxc)
    message(STATUS "  Libxc:       ${ONCVPSP_Libxc_DESCRIPTION}")
else ()
    message(STATUS "  Libxc:       disabled")
endif ()
if (ONCVPSP_WITH_TOML)
    message(STATUS "  TOML input:  ${ONCVPSP_TOML_DESCRIPTION}")
else ()
    message(STATUS "  TOML input:  disabled")
endif ()
if (NOT ONCVPSP_WITH_FR)
    message(STATUS "  Fully-relativistic executables: disabled")
endif ()
if (ONCVPSP_BLAS_DOT)
    message(STATUS "  Pseudo-atom overlaps: BLAS (dgemm)")
endif ()
if (ONCVPSP_TEST_PYTHON_DESCRIPTION)
    message(STATUS "  Test Python: ${ONCVPSP_TEST_PYTHON_DESCRIPTION}")
endif ()
message(STATUS "  Install to:  ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")

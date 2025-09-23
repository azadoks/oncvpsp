cmake_minimum_required(VERSION 3.15)

# Include Guard
set(CONFIGURED_ONCE TRUE CACHE INTERNAL
    "A flag showing, that CMake has configured at least once.")

# Prevent in-source builds
if(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL ERROR "Source folder cannot be safely used as build folder!")
endif()

##########################################################
# oncvpsp
##########################################################
project(oncvpsp
    VERSION 4.0.1
    DESCRIPTION "Optimized norm-conserving Vanderbilt pseudopotentials"
    LANGUAGES Fortran C)

##########################################################
# Build types
##########################################################
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building in debug mode")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0 -Wall -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow")
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Building in release mode")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -ffast-math -funroll-loops -DNDEBUG")
endif (CMAKE_BUILD_TYPE STREQUAL "Release")

##########################################################
# Options
##########################################################
option(ONCVPSP_ENABLE_EXTERN_FALLBACKS
    "Enable fallbacks to git submodules if external dependencies are not found"
    ON)

##########################################################
# Define the paths for static libraries and executables
##########################################################
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    CACHE
    PATH "Single output directory for building all libraries."
)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    CACHE
    PATH "Single output directory for building all libraries."
)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    CACHE
    PATH "Single output directory for building all executables."
)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
    CACHE
    PATH "Single output directory for building all Fortran modules."
)

###########################################################
# Build helpers
###########################################################
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

###########################################################
# Git submodules
###########################################################
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

###########################################################
# Common variables
###########################################################
set(COMMON_LIBRARIES "")
set(COMMON_INCLUDE_DIRECTORIES "")

###########################################################
# LAPACK
###########################################################
find_package(LAPACK HINTS ${LAPACK_ROOT} ${BLAS_ROOT})
if (LAPACK_FOUND)
    message(STATUS "LAPACK found.")
    message(STATUS "  LAPACK libraries: ${LAPACK_LIBRARIES}")
    message(STATUS "  LAPACK linker flags: ${LAPACK_LINKER_FLAGS}")
elseif (ONCVPSP_ENABLE_EXTERN_FALLBACKS)
    message(STATUS "Installing LAPACK via submodule")
    add_subdirectory(extern/lapack)
else ()
    message(FATAL ERROR "Failed to find LAPACK.")
endif()
list(APPEND COMMON_LIBRARIES ${LAPACK_LIBRARIES})

###########################################################
# HDF5
###########################################################
if (HDF5_USE_STATIC_LIBRARIES)
    find_package(ZLIB REQUIRED)
    find_package(LIBAEC REQUIRED)
endif()
find_package(HDF5 COMPONENTS Fortran HL REQUIRED)
if (NOT HDF5_FOUND)
    message(FATAL ERROR "Failed to find HDF5.")
else()
    message(STATUS "HDF5 found.")
    message(STATUS "  HDF5 version: ${HDF5_VERSION}")
    message(STATUS "  HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
    message(STATUS "  HDF5 libraries: ${HDF5_LIBRARIES}")
    message(STATUS "  HDF5 high-level libraries: ${HDF5_HL_LIBRARIES}")
endif()
list(APPEND COMMON_LIBRARIES ${HDF5_HL_LIBRARIES})
list(APPEND COMMON_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})

###########################################################
# LIBXC
###########################################################
find_package(Libxc COMPONENTS Fortran)
if (NOT Libxc_FOUND)
    message(STATUS "Libxc search failed in CMake Module mode, trying Config mode")
    find_package(Libxc COMPONENTS Fortran)
endif()
if (Libxc_FOUND)
    message(STATUS "Libxc version ${Libxc_VERSION} found.")
    if (${Libxc_VERSION} VERSION_GREATER_EQUAL 5 AND ${Libxc_VERSION} VERSION_LESS 7)
    else()
        message(STATUS "Only Libxc major versions 5 and 6 are supported. Setting Libxc_FOUND to FALSE.")
        set(Libxc_FOUND FALSE)
    endif()
    set(LIBXC_LIBRARIES Libxc::xcf03)
endif()
if (NOT Libxc_FOUND)
    if (ONCVPSP_ENABLE_EXTERN_FALLBACKS)
        message(STATUS "Installing Libxc via submodule")
        set(ENABLE_FORTRAN ON CACHE BOOL "Enable Fortran bindings" FORCE)
        add_subdirectory(extern/libxc)
        set(LIBXC_LIBRARIES Libxc::xcf03)
    else()
        message(FATAL_ERROR "Failed to find Libxc package (5.0.0 <= x <= 7.0.0) with Fortran enabled.")
    endif()
else ()
endif()
list(APPEND COMMON_LIBRARIES ${LIBXC_LIBRARIES})

###########################################################
# CMake Configuration
###########################################################
# Save compile commands to `compile_commands.json` file
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
link_directories(${CMAKE_LIBRARY_PATH} ${COMMON_INCLUDE_DIRECTORIES})
include_directories(${CMAKE_INCLUDE_PATH} ${COMMON_INCLUDE_DIRECTORIES})

###########################################################
# Subdirectories
###########################################################
add_subdirectory(src)
add_subdirectory(app)

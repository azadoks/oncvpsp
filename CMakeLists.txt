cmake_minimum_required(VERSION 3.15)

# Include Guard
set(CONFIGURED_ONCE TRUE CACHE INTERNAL
    "A flag showing, that CMake has configured at least once.")

# Prevent in-source builds
if(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL ERROR "Source folder cannot be safely used as build folder!")
endif()

##########################################################
# ONCVPSP
##########################################################
project(oncvpsp
    VERSION 4.0.1
    DESCRIPTION "Optimized norm-conserving Vanderbilt pseudopotentials"
    LANGUAGES Fortran)

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
# BLAS
###########################################################
if (NOT BLA_VENDOR)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64.*")
        message(STATUS "Trying to find BLAS from Intel MKL")
        SET(BLA_VENDOR Intel10_64lp_seq)
        find_package(BLAS)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "i686.*|i386.*|x86.*")
        message(STATUS "Trying to find BLAS from Intel MKL - 32bit")
        SET(BLA_VENDOR Intel10_32)
        find_package(BLAS)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64.*|AARCH64.*|arm64.*|ARM64.*)")
        message(STATUS "Trying to find BLAS from ARM Performance Library")
        SET(BLA_VENDOR Arm)
        find_package(BLAS)
    endif()
    if (NOT BLAS_FOUND)
        message(STATUS "Trying to find alternative BLAS libraries")
        SET(BLA_VENDOR All)
        find_package(BLAS)
    endif()
else()
    find_package(BLAS)
endif()
if (BLAS_FOUND)
    include_directories(${BLAS_INCLUDE_DIRS})
else()
    message(FATAL ERROR "Failed to find BLAS.")
endif()

###########################################################
# LAPACK
###########################################################
if (NOT BLA_VENDOR)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64.*")
        message(STATUS "Trying to find LAPACK from Intel MKL")
        SET(BLA_VENDOR Intel10_64lp_seq)
        find_package(LAPACK)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "i686.*|i386.*|x86.*")
        message(STATUS "Trying to find LAPACK from Intel MKL - 32bit")
        SET(BLA_VENDOR Intel10_32)
        find_package(LAPACK)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64.*|AARCH64.*|arm64.*|ARM64.*)")
        message(STATUS "Trying to find LAPACK from ARM Performance Library")
        SET(BLA_VENDOR Arm)
        find_package(LAPACK)
    endif()
    if (NOT LAPACK_FOUND)
        message(STATUS "Trying to find alternative LAPACK libraries")
        SET(BLA_VENDOR All)
        find_package(LAPACK)
    endif()
else()
    find_package(LAPACK)
endif()
if (LAPACK_FOUND)
    include_directories(${LAPACK_INCLUDE_DIRS})
else()
    message(FATAL ERROR "Failed to find LAPACK.")
endif()

###########################################################
# LIBXC
###########################################################
# target_compile_definitions(oncvpsp_libxc INTERFACE "__LIBXC")
find_package(Libxc COMPONENTS Fortran)
if (NOT Libxc_FOUND)
    message(STATUS "Libxc search failed in CMake Module mode, trying Config mode")
    find_package(Libxc COMPONENTS Fortran REQUIRED)
endif()
if (${Libxc_FOUND})
    message(STATUS "Libxc version ${Libxc_VERSION} found.")
    if (${Libxc_VERSION} VERSION_GREATER_EQUAL 5 AND ${Libxc_VERSION} VERSION_LESS 7)
    else()
        message(FATAL_ERROR "Only Libxc major versions 5 and 6 are supported.")
    endif()
else()
    message(FATAL_ERROR "Failed to find Libxc package (5.0.0 <= x <= 7.0.0) with Fortran enabled.")
endif()

###########################################################
# CMake Configuration
###########################################################
# Save compile commands to `compile_commands.json` file
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
link_directories(${CMAKE_LIBRARY_PATH})
include_directories(${CMAKE_INCLUDE_PATH})
###########################################################
# Subdirectories
###########################################################
add_subdirectory(src)
add_subdirectory(app)

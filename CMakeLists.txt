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
    LANGUAGES Fortran)

##########################################################
# Build types
##########################################################
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)
if (CMAKE_BUILD_TYPE STREQUAL "debug")
    message(STATUS "Building in debug mode")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0 -Wall -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow")
endif (CMAKE_BUILD_TYPE STREQUAL "debug")
if (CMAKE_BUILD_TYPE STREQUAL "release")
    message(STATUS "Building in release mode")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -ffast-math -funroll-loops -DNDEBUG")
endif (CMAKE_BUILD_TYPE STREQUAL "release")

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
# LAPACK
###########################################################
find_package(LAPACK)
if (LAPACK_FOUND)
elseif (ONCVPSP_ENABLE_EXTERN_FALLBACKS)
    message(STATUS "Installing LAPACK via submodule")
    add_subdirectory(extern/lapack)
    set(LAPACK_LIBRARIES lapack)
else ()
    message(FATAL ERROR "Failed to find LAPACK.")
endif()

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

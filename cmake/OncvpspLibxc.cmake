# Libxc, with its Fortran interface (xcf03, or xcf90 before Libxc 5), as the target oncvpsp::libxc.
#
# With ONCVPSP_Libxc_SOURCE AUTO or SYSTEM, the system's is taken from, in order,
#   1. Libxc's CMake package (LibxcConfig.cmake, installed by CMake builds of Libxc 5 and later)
#   2. pkg-config (libxc and libxcf03 or libxcf90), which searches PKG_CONFIG_PATH
#   3. the libraries and xc_version.h under the usual prefixes
# all of which search Libxc_ROOT first if it is set. It must be version 3 or later, but not 5.0.0, and
# a program using its Fortran module must compile and link, which it does not if the module was
# compiled by another compiler. With AUTO, Libxc 7.0.0 is built if none is found or works; with FETCH
# it is always built.

include(CheckFortranSourceCompiles)

set(_libxc_fetch_version 7.0.0)

# Whether Libxc version VERSION is one ONCVPSP can use; the reason it cannot goes to REASON
function(_oncvpsp_libxc_version_ok version reason)
    if (version VERSION_LESS 3)
        set(${reason} "version ${version} is older than 3" PARENT_SCOPE)
    elseif (version VERSION_EQUAL 5.0.0)
        set(${reason} "version 5.0.0 has serious bugs" PARENT_SCOPE)
    else ()
        set(${reason} "" PARENT_SCOPE)
    endif ()
endfunction()

# Whether a program using the Fortran interface of Libxc VERSION and xc_version.h compiles and links
# with the libraries or targets in ARGN; the result goes to the cache variable VAR
function(_oncvpsp_libxc_links var version)
    if (version VERSION_LESS 5)
        set(f "f90")
    else ()
        set(f "f03")
    endif ()
    unset(${var} CACHE)
    set(CMAKE_REQUIRED_LIBRARIES ${ARGN})
    set(CMAKE_REQUIRED_QUIET ON)
    check_fortran_source_compiles("
program libxc_check
  use xc_${f}_lib_m
  implicit none
#include \"xc_version.h\"
  integer :: vmajor, vminor, vmicro
  call xc_${f}_version(vmajor, vminor, vmicro)
end program libxc_check
" ${var} SRC_EXT F90)
endfunction()

# The version in the xc_version.h under the directories in ARGN, to VAR
function(_oncvpsp_libxc_header_version var)
    find_file(_xc_version_h xc_version.h PATHS ${ARGN} NO_DEFAULT_PATH NO_CACHE)
    set(version "")
    if (_xc_version_h)
        file(STRINGS "${_xc_version_h}" defines REGEX "^#define XC_(MAJOR|MINOR|MICRO)_VERSION")
        foreach (part IN ITEMS MAJOR MINOR MICRO)
            string(REGEX MATCH "XC_${part}_VERSION +([0-9]+)" _ "${defines}")
            list(APPEND version "${CMAKE_MATCH_1}")
        endforeach ()
        list(JOIN version "." version)
    endif ()
    set(${var} "${version}" PARENT_SCOPE)
endfunction()

set(_libxc_libraries "")
if (NOT ONCVPSP_Libxc_SOURCE STREQUAL "FETCH")
    # each candidate is "where|version|library,library,..."
    set(_libxc_candidates "")

    # 1. Libxc's CMake package
    find_package(Libxc CONFIG QUIET)
    if (Libxc_FOUND)
        if (TARGET Libxc::xcf03)
            list(APPEND _libxc_candidates "CMake package ${Libxc_DIR}|${Libxc_VERSION}|Libxc::xcf03,Libxc::xc")
        elseif (TARGET Libxc::xcf90)
            list(APPEND _libxc_candidates "CMake package ${Libxc_DIR}|${Libxc_VERSION}|Libxc::xcf90,Libxc::xc")
        else ()
            message(STATUS "Libxc: ${Libxc_DIR} has no Fortran interface")
        endif ()
    endif ()

    # 2. pkg-config
    find_package(PkgConfig QUIET)
    if (PkgConfig_FOUND)
        set(_pkg_config_path "$ENV{PKG_CONFIG_PATH}")
        if (DEFINED Libxc_ROOT)
            set(ENV{PKG_CONFIG_PATH} "${Libxc_ROOT}/lib/pkgconfig:${Libxc_ROOT}/lib64/pkgconfig:${_pkg_config_path}")
        endif ()
        pkg_check_modules(ONCVPSP_PC_Libxc QUIET IMPORTED_TARGET libxc)
        if (ONCVPSP_PC_Libxc_FOUND)
            if (ONCVPSP_PC_Libxc_VERSION VERSION_LESS 5)
                pkg_check_modules(ONCVPSP_PC_LibxcFortran QUIET IMPORTED_TARGET libxcf90)
            else ()
                pkg_check_modules(ONCVPSP_PC_LibxcFortran QUIET IMPORTED_TARGET libxcf03)
            endif ()
            if (ONCVPSP_PC_LibxcFortran_FOUND)
                # pkg-config leaves out system include directories such as /usr/include, which
                # Fortran compilers do not search for modules
                set_property(TARGET PkgConfig::ONCVPSP_PC_LibxcFortran APPEND PROPERTY
                             INTERFACE_INCLUDE_DIRECTORIES "${ONCVPSP_PC_LibxcFortran_INCLUDEDIR}")
                list(APPEND _libxc_candidates
                    "pkg-config ${ONCVPSP_PC_Libxc_PREFIX}|${ONCVPSP_PC_Libxc_VERSION}|PkgConfig::ONCVPSP_PC_LibxcFortran,PkgConfig::ONCVPSP_PC_Libxc")
            endif ()
        endif ()
        set(ENV{PKG_CONFIG_PATH} "${_pkg_config_path}")
    endif ()

    # 3. the libraries and header
    find_path(ONCVPSP_Libxc_INCLUDE_DIR xc_version.h HINTS ${Libxc_ROOT} ENV Libxc_ROOT PATH_SUFFIXES include)
    foreach (_name IN ITEMS xc xcf03 xcf90)
        string(TOUPPER "${_name}" _var)
        find_library(ONCVPSP_Libxc_${_var}_LIBRARY ${_name} HINTS ${Libxc_ROOT} ENV Libxc_ROOT PATH_SUFFIXES lib lib64)
    endforeach ()
    mark_as_advanced(ONCVPSP_Libxc_INCLUDE_DIR ONCVPSP_Libxc_XC_LIBRARY ONCVPSP_Libxc_XCF03_LIBRARY
                     ONCVPSP_Libxc_XCF90_LIBRARY)
    if (ONCVPSP_Libxc_INCLUDE_DIR AND ONCVPSP_Libxc_XC_LIBRARY)
        _oncvpsp_libxc_header_version(_version "${ONCVPSP_Libxc_INCLUDE_DIR}")
        if (_version VERSION_LESS 5)
            set(_fortran_library "${ONCVPSP_Libxc_XCF90_LIBRARY}")
        else ()
            set(_fortran_library "${ONCVPSP_Libxc_XCF03_LIBRARY}")
        endif ()
        if (_version AND _fortran_library)
            if (NOT TARGET oncvpsp_libxc_found)
                add_library(oncvpsp_libxc_found INTERFACE IMPORTED)
                target_include_directories(oncvpsp_libxc_found INTERFACE "${ONCVPSP_Libxc_INCLUDE_DIR}")
                target_link_libraries(oncvpsp_libxc_found INTERFACE "${_fortran_library}" "${ONCVPSP_Libxc_XC_LIBRARY}")
            endif ()
            list(APPEND _libxc_candidates "${ONCVPSP_Libxc_XC_LIBRARY}|${_version}|oncvpsp_libxc_found")
        endif ()
    endif ()

    # the first candidate that ONCVPSP can use
    foreach (_candidate IN LISTS _libxc_candidates)
        string(REPLACE "|" ";" _candidate "${_candidate}")
        list(POP_FRONT _candidate _where _version)
        string(REPLACE "," ";" _candidate "${_candidate}")
        _oncvpsp_libxc_version_ok("${_version}" _reason)
        if (NOT _reason)
            string(MAKE_C_IDENTIFIER "ONCVPSP_Libxc_LINKS_${_where}" _check)
            _oncvpsp_libxc_links(${_check} "${_version}" ${_candidate})
            if (NOT ${_check})
                string(CONCAT _reason "a program using its Fortran module does not compile and link "
                                      "(was the module compiled by another compiler?)")
            endif ()
        endif ()
        if (_reason)
            message(STATUS "Libxc: not using ${_version} from ${_where}: ${_reason}")
        else ()
            set(_libxc_libraries ${_candidate})
            set(ONCVPSP_Libxc_DESCRIPTION "Libxc ${_version} (${_where})")
            break ()
        endif ()
    endforeach ()

    if (NOT _libxc_libraries AND ONCVPSP_Libxc_SOURCE STREQUAL "SYSTEM")
        message(FATAL_ERROR
            "Libxc: no usable Libxc with a Fortran interface found. Set Libxc_ROOT to its installation, "
            "set ONCVPSP_Libxc_SOURCE to AUTO or FETCH to build Libxc ${_libxc_fetch_version}, or "
            "ONCVPSP_WITH_Libxc to OFF")
    endif ()
endif ()

if (NOT _libxc_libraries)
    # Only the energies and their first and second derivatives (vxc, fxc) are used. The release archive
    # is much faster to download than a clone of the repository.
    CPMAddPackage(
        NAME libxc
        VERSION ${_libxc_fetch_version}
        URL "https://gitlab.com/libxc/libxc/-/archive/${_libxc_fetch_version}/libxc-${_libxc_fetch_version}.tar.bz2"
        URL_HASH SHA256=e9ae69f8966d8de6b7585abd9fab588794ada1fab8f689337959a35abbf9527d
        EXCLUDE_FROM_ALL YES
        SYSTEM YES
        OPTIONS "BUILD_TESTING OFF" "BUILD_SHARED_LIBS OFF" "ENABLE_FORTRAN ON" "ENABLE_PYTHON OFF"
                "ENABLE_XHOST OFF" "DISABLE_KXC ON" "DISABLE_LXC ON"
    )
    # Libxc's targets give neither xc_version.h nor the Fortran modules to the targets that use them
    set_target_properties(xcf03 PROPERTIES Fortran_MODULE_DIRECTORY "${libxc_BINARY_DIR}/modules")
    add_library(oncvpsp_libxc_built INTERFACE)
    target_include_directories(oncvpsp_libxc_built INTERFACE "${libxc_BINARY_DIR}" "${libxc_BINARY_DIR}/modules")
    target_link_libraries(oncvpsp_libxc_built INTERFACE xcf03 xc)
    set(_libxc_libraries oncvpsp_libxc_built)
    set(ONCVPSP_Libxc_DESCRIPTION "Libxc ${_libxc_fetch_version} (built)")
endif ()

add_library(oncvpsp_libxc INTERFACE)
add_library(oncvpsp::libxc ALIAS oncvpsp_libxc)
target_link_libraries(oncvpsp_libxc INTERFACE ${_libxc_libraries})
unset(_libxc_libraries)
unset(_libxc_candidates)

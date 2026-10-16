# LAPACK and BLAS, as the target oncvpsp::lapack.
#
# With ONCVPSP_LAPACK_SOURCE AUTO or SYSTEM, the system's are taken from, in order,
#   1. ONCVPSP_LAPACK_LIBRARIES, a list of libraries or linker flags, if it is set
#   2. CMake's FindLAPACK, which honours BLA_VENDOR (e.g. OpenBLAS, FlexiBLAS, Intel10_64lp,
#      Apple), BLA_STATIC, and LAPACK_ROOT and BLAS_ROOT
#   3. pkg-config (the lapack and blas modules), which searches PKG_CONFIG_PATH, unless BLA_VENDOR
#      is set
# and must link a program calling the routines ONCVPSP uses. With AUTO, Reference-LAPACK is built if
# none does; with FETCH it is always built.

include(CheckFortranSourceCompiles)

# Whether a program calling the LAPACK routines ONCVPSP uses links with the libraries or targets
# in ARGN; the result goes to the cache variable VAR
function(_oncvpsp_lapack_links var)
    # checked at every configuration, since what is found depends on BLA_VENDOR and the like
    unset(${var} CACHE)
    set(CMAKE_REQUIRED_LIBRARIES ${ARGN})
    set(CMAKE_REQUIRED_QUIET ON)
    check_fortran_source_compiles([=[
program lapack_check
  implicit none
  integer :: info, ipiv(2)
  double precision :: a(2,2), b(2), w(2), s(2), u(2,2), vt(2,2), work(20)
  a = reshape([2d0, 1d0, 1d0, 2d0], [2,2]); b = 1d0
  call dsyev('V', 'U', 2, a, 2, w, work, 20, info)
  call dgesv(2, 1, a, 2, ipiv, b, 2, info)
  call dgesvd('A', 'A', 2, 2, a, 2, s, u, 2, vt, 2, work, 20, info)
end program lapack_check
]=] ${var} SRC_EXT F90)
endfunction()

set(_lapack_libraries "")
if (NOT ONCVPSP_LAPACK_SOURCE STREQUAL "FETCH")
    set(ONCVPSP_LAPACK_LIBRARIES "" CACHE STRING
        "LAPACK and BLAS libraries or linker flags, used instead of searching for them")
    if (ONCVPSP_LAPACK_LIBRARIES)
        _oncvpsp_lapack_links(ONCVPSP_LAPACK_LINKS_USER ${ONCVPSP_LAPACK_LIBRARIES})
        if (NOT ONCVPSP_LAPACK_LINKS_USER)
            message(FATAL_ERROR "LAPACK: a program calling dsyev, dgesv and dgesvd does not link with "
                                "ONCVPSP_LAPACK_LIBRARIES=${ONCVPSP_LAPACK_LIBRARIES}")
        endif ()
        set(_lapack_libraries ${ONCVPSP_LAPACK_LIBRARIES})
        set(ONCVPSP_LAPACK_DESCRIPTION "${ONCVPSP_LAPACK_LIBRARIES}")
    else ()
        find_package(LAPACK QUIET)
        if (LAPACK_FOUND)
            _oncvpsp_lapack_links(ONCVPSP_LAPACK_LINKS_FIND LAPACK::LAPACK)
            if (ONCVPSP_LAPACK_LINKS_FIND)
                set(_lapack_libraries LAPACK::LAPACK)
                set(ONCVPSP_LAPACK_DESCRIPTION "${LAPACK_LIBRARIES}")
            else ()
                message(STATUS "LAPACK: ${LAPACK_LIBRARIES} found, but a program calling it does not link")
            endif ()
        endif ()
        # pkg-config cannot be told which vendor BLA_VENDOR asks for
        if (NOT _lapack_libraries AND (NOT BLA_VENDOR OR BLA_VENDOR STREQUAL "All"))
            find_package(PkgConfig QUIET)
            if (PkgConfig_FOUND)
                pkg_check_modules(ONCVPSP_PC_LAPACK QUIET IMPORTED_TARGET lapack blas)
                if (ONCVPSP_PC_LAPACK_FOUND)
                    _oncvpsp_lapack_links(ONCVPSP_LAPACK_LINKS_PKGCONFIG PkgConfig::ONCVPSP_PC_LAPACK)
                    if (ONCVPSP_LAPACK_LINKS_PKGCONFIG)
                        set(_lapack_libraries PkgConfig::ONCVPSP_PC_LAPACK)
                        set(ONCVPSP_LAPACK_DESCRIPTION "${ONCVPSP_PC_LAPACK_LINK_LIBRARIES} (pkg-config)")
                    endif ()
                endif ()
            endif ()
        endif ()
    endif ()
    if (NOT _lapack_libraries AND ONCVPSP_LAPACK_SOURCE STREQUAL "SYSTEM")
        message(FATAL_ERROR
            "LAPACK: no working LAPACK found. Set BLA_VENDOR, LAPACK_ROOT or ONCVPSP_LAPACK_LIBRARIES, "
            "or ONCVPSP_LAPACK_SOURCE to AUTO or FETCH to build Reference-LAPACK")
    endif ()
endif ()

if (NOT _lapack_libraries)
    CPMAddPackage(
        NAME lapack
        GITHUB_REPOSITORY Reference-LAPACK/lapack
        GIT_TAG v3.12.1
        EXCLUDE_FROM_ALL YES
        SYSTEM YES
        OPTIONS "BUILD_TESTING OFF" "BUILD_SHARED_LIBS OFF" "CBLAS OFF" "LAPACKE OFF"
                "BUILD_INDEX64 OFF" "BUILD_INDEX64_EXT_API OFF" "BUILD_SINGLE OFF" "BUILD_COMPLEX OFF" "BUILD_COMPLEX16 OFF"
    )
    set(_lapack_libraries lapack blas)
    set(ONCVPSP_LAPACK_DESCRIPTION "Reference-LAPACK 3.12.1 (built)")
endif ()

add_library(oncvpsp_lapack INTERFACE)
add_library(oncvpsp::lapack ALIAS oncvpsp_lapack)
target_link_libraries(oncvpsp_lapack INTERFACE ${_lapack_libraries})
unset(_lapack_libraries)

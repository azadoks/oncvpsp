# Compiler options for the ONCVPSP targets (not the dependencies), on top of CMake's per-build-type
# defaults (CMAKE_Fortran_FLAGS_<CONFIG>), which stay the user's to change:
#   Debug    run-time checks, and traps on invalid operations, division by zero and overflow; with
#            gfortran, standard Fortran 2008 (-std=f2008) and warnings of conversions that may
#            change a value (-Wconversion), except in the imported files, which keep their form
#   Release  fast math and unrolled loops
# and, with gfortran in RelWithDebInfo (-O2), a higher inlining limit, so the Adams formulas of
# src/adams.inc are inlined into the radial integrators (at -O2 they are calls, about 9 % of the run
# time; the output is unchanged). It is the -O3 default, so Release needs no change
# With gfortran in every build type, a call of a procedure without an explicit interface is an
# error (-Werror=implicit-interface), so every call is checked against its definition: every
# procedure is in a module, and the external libraries are called through interfaces (lapack_m).
# It is the check of Fortran 2018's implicit none (type, external), without that syntax.
# In every build type the procedures are compiled recursive (-frecursive, -recursive, -Mrecursive):
# local arrays are on the stack, never static (gfortran makes those over 64 kB static otherwise),
# and gfortran's run-time check of recursion (-fcheck=all) is off, so that runs can share a process
# in threads (tests/unit/test_run_threads.f90).
# The targets link the interface target oncvpsp::compile_options to get them.

add_library(oncvpsp_compile_options INTERFACE)
add_library(oncvpsp::compile_options ALIAS oncvpsp_compile_options)

set(debug "$<CONFIG:Debug>")
set(release "$<CONFIG:Release>")
set(relwithdebinfo "$<CONFIG:RelWithDebInfo>")
target_compile_options(oncvpsp_compile_options INTERFACE
    # gfortran
    "$<$<Fortran_COMPILER_ID:GNU>:-Werror=implicit-interface;-frecursive>"
    "$<$<AND:${debug},$<Fortran_COMPILER_ID:GNU>>:-O0;-Wall;-fcheck=all;-fbacktrace;-ffpe-trap=invalid,zero,overflow;-std=f2008;-Wconversion>"
    "$<$<AND:${release},$<Fortran_COMPILER_ID:GNU>>:-ffast-math;-funroll-loops>"
    "$<$<AND:${relwithdebinfo},$<Fortran_COMPILER_ID:GNU>>:--param=max-inline-insns-auto=30>"
    # ifx and ifort
    "$<$<Fortran_COMPILER_ID:IntelLLVM,Intel>:-recursive>"
    "$<$<AND:${debug},$<Fortran_COMPILER_ID:IntelLLVM,Intel>>:-O0;-warn;all;-check;all;-traceback;-fpe0>"
    "$<$<AND:${release},$<Fortran_COMPILER_ID:IntelLLVM,Intel>>:-fp-model=fast>"
    # nvfortran
    "$<$<Fortran_COMPILER_ID:NVHPC>:-Mrecursive>"
    "$<$<AND:${debug},$<Fortran_COMPILER_ID:NVHPC>>:-O0;-Mbounds;-traceback;-Ktrap=inv,divz,ovf>"
    "$<$<AND:${release},$<Fortran_COMPILER_ID:NVHPC>>:-fast>"
)
unset(debug)
unset(release)
unset(relwithdebinfo)

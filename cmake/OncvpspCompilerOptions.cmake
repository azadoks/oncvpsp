# Compiler options for the ONCVPSP targets (not the dependencies), on top of CMake's per-build-type
# defaults (CMAKE_Fortran_FLAGS_<CONFIG>), which stay the user's to change:
#   Debug    run-time checks, and traps on invalid operations, division by zero and overflow
#   Release  fast math and unrolled loops
# The targets link the interface target oncvpsp::compile_options to get them.

add_library(oncvpsp_compile_options INTERFACE)
add_library(oncvpsp::compile_options ALIAS oncvpsp_compile_options)

set(debug "$<CONFIG:Debug>")
set(release "$<CONFIG:Release>")
target_compile_options(oncvpsp_compile_options INTERFACE
    # gfortran
    "$<$<AND:${debug},$<Fortran_COMPILER_ID:GNU>>:-O0;-Wall;-fcheck=all;-fbacktrace;-ffpe-trap=invalid,zero,overflow>"
    "$<$<AND:${release},$<Fortran_COMPILER_ID:GNU>>:-ffast-math;-funroll-loops>"
    # ifx and ifort
    "$<$<AND:${debug},$<Fortran_COMPILER_ID:IntelLLVM,Intel>>:-O0;-warn;all;-check;all;-traceback;-fpe0>"
    "$<$<AND:${release},$<Fortran_COMPILER_ID:IntelLLVM,Intel>>:-fp-model=fast>"
    # nvfortran
    "$<$<AND:${debug},$<Fortran_COMPILER_ID:NVHPC>>:-O0;-Mbounds;-traceback;-Ktrap=inv,divz,ovf>"
    "$<$<AND:${release},$<Fortran_COMPILER_ID:NVHPC>>:-fast>"
)
unset(debug)
unset(release)

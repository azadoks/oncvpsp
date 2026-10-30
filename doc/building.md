# Building and Installing

How to build ONCVPSP with CMake: the executables, and the library with its
C interface. The last two sections cover the Python and Julia packages:
installing each, pointing it at the library, and an example that uses most
of its interface. [Python and Julia](wrappers.md) describes the interface
itself.

## Requirements

- CMake 3.22 or later, and a build tool (Make or Ninja).
- A Fortran compiler: gfortran (the references are checked with
  gfortran 11), ifx or ifort, or nvfortran; and a C compiler, for the C
  interface and its test.
- LAPACK and BLAS, Libxc and toml-f, from the system or downloaded and
  built by the configuration (see [Dependencies](#dependencies)). xmlf90,
  for the PSML output, is always downloaded and built.
- git, so that the build records the commit it was built from (otherwise
  the commit is `unknown`), and network access at configuration time for
  whatever is downloaded.

## Quick start

```bash
cmake -S . -B build                    # configure (RelWithDebInfo)
cmake --build build --parallel         # compile
ctest --test-dir build --parallel 8    # test
cmake --install build --prefix ~/.local
```

The executables are then in `build/bin` (`oncvpsp.x`, `oncvpspnr.x`,
`oncvpspr.x`, and the `_ae.x`, `_ps.x`, `oncvpsp_ld.x` and
`oncvpsp_confine.x` tools), and the library in `build/lib`
(`liboncvpsp.so`, with its header `include/oncvpsp.h`).

The same with a preset:

```bash
cmake --preset default
cmake --build --preset default --parallel
ctest --preset default
```

## Build types and compiler flags

The build type (`CMAKE_BUILD_TYPE`) picks the optimization. Each type
starts from CMake's flags for it, `CMAKE_Fortran_FLAGS_<TYPE>`, which stay
yours to change; with gfortran these are `-O2 -g -DNDEBUG` for
RelWithDebInfo, `-O3 -DNDEBUG` for Release and `-g` for Debug. ONCVPSP adds
its own flags on top (`cmake/OncvpspCompilerOptions.cmake`), to its own
targets only, not to the dependencies:

`RelWithDebInfo` (the default)
:   The build the test references are checked against.

    - gfortran: `--param=max-inline-insns-auto=30`, so that the radial
      integrators' Adams formulas are inlined as they are at `-O3` (the
      output is unchanged)
    - ifx, ifort, nvfortran: nothing

`Debug`
:   Run-time checks and floating-point traps, for when something goes
    wrong.

    - gfortran: `-O0 -Wall -fcheck=all -fbacktrace
      -ffpe-trap=invalid,zero,overflow -std=f2008 -Wconversion`
    - ifx, ifort: `-O0 -warn all -check all -traceback -fpe0`
    - nvfortran: `-O0 -Mbounds -traceback -Ktrap=inv,divz,ovf`

`Release`
:   Fast math: faster, and the results move at round-off (about 1e-8 in
    the nonlocal kernels). With `-march=native` as well (the `study`
    preset), about 35 % faster than RelWithDebInfo.

    - gfortran: `-ffast-math -funroll-loops`
    - ifx, ifort: `-fp-model=fast`
    - nvfortran: `-fast`

In every build type:

- gfortran makes a call without an explicit interface an error
  (`-Werror=implicit-interface`), so that every call is checked against its
  definition.
- The procedures are compiled recursive (`-frecursive`, `-recursive`,
  `-Mrecursive`), so local arrays live on the stack and several runs can
  share a process, in threads. Don't remove this flag if you set your own
  flags.

The compiler is chosen as usual, on a fresh build directory:

```bash
cmake -S . -B build-intel \
    -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_C_COMPILER=icx
FC=ifx CC=icx cmake -S . -B build-intel    # the same
```

### Presets

The presets in `CMakePresets.json` each build in their own directory,
`build-<preset>/`, except `default` (in `build/`) and `coverage` (in
`build-cov/`):

| Preset | Build |
|---|---|
| `default` | RelWithDebInfo |
| `debug` | Debug |
| `release` | Release |
| `study` | Release with `-O3 -march=native`; runs only on CPUs like the one that built it |
| `trap` | `-O3 -march=native -funroll-loops` with floating-point traps and no fast math |
| `fetch` | RelWithDebInfo, downloading and building LAPACK, Libxc and toml-f |
| `profile` | RelWithDebInfo with frame pointers, for `perf` |
| `coverage` | `-O0 --coverage`, for `scripts/coverage.py` |

Personal presets (compiler paths, `Libxc_ROOT`) go in
`CMakeUserPresets.json`, which git ignores.

## ONCVPSP's options

`ONCVPSP_WITH_Libxc` (`ON`)
:   The Libxc functionals (`iexc < 0`). Off, only the built-in ones.

`ONCVPSP_WITH_TOML` (`ON`)
:   The TOML input (`--toml-input`), which the C interface and the Python
    and Julia packages need.

`ONCVPSP_WITH_FR` (`ON`)
:   The fully-relativistic executables (`oncvpspr*.x`) and their tests.

`ONCVPSP_WITH_C_API` (`ON`)
:   The shared library `liboncvpsp` with the C interface. It also makes
    everything position-independent (`-fPIC`), which leaves the results
    unchanged.

`ONCVPSP_BLAS_DOT` (`OFF`)
:   Computes the pseudo-atom's overlaps with the linked BLAS (`dgemm`):
    faster with an optimized BLAS, but the results then depend on the BLAS
    at round-off.

`BUILD_TESTING` (`ON`)
:   The tests (`tests/`).

## Dependencies

Each of LAPACK, Libxc and toml-f comes from where
`ONCVPSP_LAPACK_SOURCE`, `ONCVPSP_Libxc_SOURCE` and `ONCVPSP_TOML_SOURCE`
say:

- `AUTO` (the default): the system's copy, if one is found and a program
  using it compiles and links; otherwise it is downloaded and built
  (Reference-LAPACK 3.12.1, Libxc 7.0.0, toml-f 0.4.3).
- `SYSTEM`: the system's copy, or an error.
- `FETCH`: always downloaded and built.

`CPM_SOURCE_CACHE` (e.g. `~/.cache/CPM`) keeps the downloads between build
directories. The configuration ends with a summary of what was found.

### LAPACK and BLAS

The configuration tries these, in order:

1. `ONCVPSP_LAPACK_LIBRARIES`, if given.
2. CMake's FindLAPACK, which follows `BLA_VENDOR`, `BLA_STATIC`,
   `LAPACK_ROOT` and `BLAS_ROOT`. Without `BLA_VENDOR` it takes the first
   BLAS it finds, which on a machine with several can be MKL.
3. pkg-config.

```bash
cmake -S . -B build -DBLA_VENDOR=OpenBLAS
cmake -S . -B build -DBLA_VENDOR=FlexiBLAS
cmake -S . -B build -DBLA_VENDOR=Intel10_64_dyn     # MKL through mkl_rt
cmake -S . -B build -DBLA_VENDOR=Intel10_64lp_seq   # MKL, sequential
cmake -S . -B build -DBLA_VENDOR=Apple              # Accelerate
cmake -S . -B build \
    -DONCVPSP_LAPACK_LIBRARIES="-L/opt/lapack/lib;-llapack;-lblas"
```

A threaded BLAS changes the last digits of the results with its number of
threads. ONCVPSP's matrices are small, so the threads gain little; set
`OPENBLAS_NUM_THREADS=1` (or `MKL_NUM_THREADS=1`) when several runs share
the machine, and for results that repeat bit for bit.

### Libxc and toml-f

**Libxc.** The configuration looks for a CMake package, then pkg-config,
then the library's files, under `Libxc_ROOT` first if it is set. It needs
version 3 or later (not 5.0.0), with its Fortran interface built by the
same compiler: `-DLibxc_ROOT=/opt/libxc`.

**toml-f.** The configuration looks for a CMake package, version 0.4 or
later, under `toml-f_ROOT` or `CMAKE_PREFIX_PATH`.

## Tests

`ctest` runs every test; `-L` selects a kind of test (`ref`, `run`, `unit`,
`split`, `mesh`, `scf`, `toml`, `well`, `output`, `psml`, `psp8`, `pipe`,
`echo`, `errors`, `nist`, `lint`) or a relativistic treatment (`nr`,
`sr`, `fr`), and `-R` selects tests by name.

The output comparisons need a Python with
[pyoncvpsp](https://github.com/azadoks/pyoncvpsp) and upf-tools. The
configuration makes a virtual environment for them in the build tree
(`test-venv/`), unless `ONCVPSP_TEST_PYTHON` names a Python that already
has them. `ONCVPSP_PYONCVPSP` and `ONCVPSP_UPF_TOOLS` choose what that
environment installs. The tests of the Python package also need NumPy in
that Python, and those of the Julia package a `julia` on the `PATH`; each
set is left out, with a message, when what it needs is missing.

## Installing

`cmake --install build --prefix PREFIX` installs:

- the executables, in `PREFIX/bin`;
- the library, in `PREFIX/lib` (`liboncvpsp.so.1`, with the link
  `liboncvpsp.so`);
- the header, in `PREFIX/include/oncvpsp.h`.

The installed executables keep the paths of the shared libraries they link
outside the system directories (a Libxc under `Libxc_ROOT`, say), so they
need no `LD_LIBRARY_PATH`.

## Python

The package `oncvpsp` (`python/oncvpsp`) calls `liboncvpsp` through ctypes
and returns NumPy arrays. It needs Python 3.11 or later and NumPy, and a
library built with `ONCVPSP_WITH_C_API` and `ONCVPSP_WITH_TOML` (both on by
default).

### Installing

**With pip.** This compiles the library and installs it inside the
package, from the source tree:

```bash
pip install .
```

pip builds through scikit-build-core (`pyproject.toml`), without the tests.
CMake options go through `-C` (`--config-settings`):

```bash
# the reference build type
pip install . -C cmake.build-type=RelWithDebInfo
# a BLAS and a Libxc
pip install . -C cmake.define.BLA_VENDOR=OpenBLAS \
    -C cmake.define.Libxc_ROOT=/opt/libxc
# a compiler
FC=ifx CC=icx pip install .
```

!!! note "pip builds Release by default"
    Without `-C cmake.build-type=...`, pip makes a `Release` build, with
    fast math, whose results differ from the RelWithDebInfo executables'
    at round-off.

!!! tip "Results identical to the executables'"
    For results identical to the executables', bit for bit, the package
    needs the executables' build type and the same BLAS. pip configures
    from scratch, and FindLAPACK takes the first BLAS it finds, which can
    be a different one from your CMake build's: give the same `BLA_VENDOR`
    to both. `ldd` on the library shows which BLAS it links.

The wheel (`pip wheel .`) is `py3-none-<platform>`: the library is called
through ctypes, so the wheel doesn't depend on the Python version. It links
the system's LAPACK, Libxc and Fortran runtime, so it is for the machine
that built it, not for distribution as is.

**From a CMake build.** This uses the library of a build directory, with
no installation:

```bash
cmake -S . -B build && cmake --build build --parallel
export ONCVPSP_LIBRARY=$PWD/build/lib/liboncvpsp.so
export PYTHONPATH=$PWD/python
```

The package looks for the library in this order:

1. `ONCVPSP_LIBRARY`;
2. the package's own directory (where pip puts it);
3. the system's library search path (e.g. after `cmake --install` into a
   prefix on `LD_LIBRARY_PATH`).

It checks the library's interface version (`oncv_api_version`) when it is
imported. A library and a package of different interface versions refuse
to load together: rebuild the library from the same commit as the
package.

**A conda Python.** A conda Python may load conda's libquadmath before the
compiler's. Its quad precision rounds differently, and the results then
differ from the executables' at about 1e-9. Preloading the compiler's
library restores them:
{: #conda-libquadmath }

```bash
export LD_PRELOAD=$(gfortran -print-file-name=libquadmath.so.0)
```

### Example

The steps below make one script, run from the top of the source tree (it
reads `tests/data/toml/14_Si.toml`).

The library's build, and a whole run from a TOML file. `irel` is 0, 1 or 2,
or `"nr"`, `"sr"`, `"fr"`. The run is a `Mapping` of its fields to NumPy
arrays in the Fortran shapes, and `field_info` says what each one is:

```python
import numpy as np

import oncvpsp

print(oncvpsp.build_info())  # {'commit': ..., 'source': ..., 'compiled': ...}

run = oncvpsp.run("tests/data/toml/14_Si.toml", irel="sr")
dims = run.dims  # mmax, lmax, lloc, nproj, nkap, ...
print(dims.mmax, dims.lmax, dims.nproj[: dims.lmax + 1])
print(len(run), "fields:", list(run)[:3])
print(oncvpsp.field_info("pot.vkb"))  # its axes, units and what it is
vkb = run["pot.vkb"]  # (r, iproj, l, ikap): (mmax, mxprj, mxl+1, nkap)
```

The physics, without the indices: the mesh, the local and semi-local
potentials, each channel's projectors (only the `nproj` it has) and its
nonlocal operator, the densities in e/bohr³ (the fields are 4π times
them), and a pseudo-wave's convergence with the cutoff:

```python
rr = run.mesh  # bohr
vloc = run.local_potential  # Ha
energies, beta = run.projectors(0)  # l=0: (nproj,), (mmax, nproj)
kk = run.kernel(0)  # sum_i beta_i(r) e_i beta_i(r'), (mmax, mmax)
rho = run.valence_density
print("valence charge:", np.trapezoid(4 * np.pi * rr**2 * rho, rr))
print(run.convergence(0))  # rows (Ecut, err), l=0's first pseudo-wave
```

The tests, `tests.*`, and the all-electron atom, `ae.*`:

```python
# (Ecut, err) per level, projector, l, ikap
cross = run["tests.convergence.cross"]
print("Ecut at 1e-3 Ha, l=0, first projector:", cross[0, 2, 0, 0, 0])
# per test configuration: dE_ae, dE_ps, error
configs = run["tests.configs.de"]
print("excitation errors:", configs[2, : run["tests.configs.ncfg"]])
ee, ae_ld, ps_ld = (run[f"tests.logder.{k}"] for k in ("ee", "ae", "ps"))
print("ghost bound states per l:", run["tests.ghosts.nrow"][:, 0])
print("all-electron eigenvalues:", run["ae.ea"][:, 0])
```

The text output, and the pseudopotential files (psp8, UPF or PSML, by the
extension):

```python
assert "Begin DIAGNOSTICS" in run.text
run.write("Si.psp8")
run.write("Si.upf")
```

The input as a dictionary, with tables merged over it key by key, here for
a scan of `qcut` (the table's other keys stay `base`'s). `check=False`
leaves a failure in the handle instead of raising it, and `text=False`
discards the text output (`text="run.out"` would write it to that file):

```python
base = oncvpsp.read_toml("tests/data/toml/14_Si.toml")
for qcut in (5.0, 6.0, 7.0):
    pp = {"qcut": [qcut, 5.0, 5.0]}
    scan = oncvpsp.run(
        base, irel=1, pseudopotentials=pp, check=False, text=False
    )
    if scan.ifail:  # at 7.0 the s pseudo-wave has a node
        print(qcut, "failed:", scan.error())
    else:
        print(qcut, scan["tests.convergence.opt"][0, 0, 0, 0])
```

The all-electron atom once, then the stages from it: a run from the atom,
bit for bit the whole run's; each channel, with its projectors; and the
node test, the cheapest way to scan for a node:

```python
ae = oncvpsp.all_electron(base, irel=1)
print("ae total energy:", ae["ae.etot"])
same = oncvpsp.run(base, ae=ae, text=False)
assert np.array_equal(same["pot.vkb"], vkb)
for ll in range(dims.lmax + 1):
    ch = oncvpsp.channel(base, ae, ll, mode="projectors", text=False)
    print(ll, ch["channel.pswf"].shape, ch.projectors().energies)
pp = {"qcut": [7.0, 5.0, 5.0]}
print("node at qcut 7:", oncvpsp.has_node(base, ae, 0, pseudopotentials=pp))
```

The model core again on the run's stages, and the confined atom (as
`oncvpsp_confine.x`):

```python
teter = oncvpsp.coupled(run, icmod=3, fcfact=1.5, rcfact=1.2, text=False)
print("pseudo-atom energy, icmod 3:", teter["epstot"])
conf = oncvpsp.confine(run, radii=[3.0, 3.5], walls=())
print("confine:", conf["confine.nmiss"], conf["confine.max_de"])
```

Failures: an `OncvpspError` of the class of its kind (`InputError`,
`ConvergenceError`, `ConstructionError` and its `NodeError`, ...), or with
`check=False` the handle's `ifail` and `failure`. The library's own message
(`Error: missing key z in [oncvpsp]`) goes to stderr. A handle is freed
when it is no longer referenced, at the end of a `with` block, or by
`close()`:

```python
try:
    oncvpsp.run({"oncvpsp": {"atsym": "O"}}, irel=1)
except oncvpsp.InputError as err:
    print(err.kind.name, err.stage.name, err.msg)  # INPUT INPUT 42
bad = oncvpsp.run({"oncvpsp": {"atsym": "O"}}, irel=1, check=False)
print(bad.ifail, bad.failure.stage)

with oncvpsp.run(base, irel=1, text=False) as again:
    print(again.kernel(1).shape)
run.close()
```

## Julia

The package `ONCVPSP` (`julia/ONCVPSP`) calls the same library with
`ccall` and returns Julia arrays, in the Fortran shapes. It needs Julia 1.9
or later and has no dependencies outside the standard library. There is no
registered package or prebuilt library (JLL) yet, so the library comes from
a CMake build.

### Installing

Build the library, then either use the package's own environment:

```bash
cmake -S . -B build && cmake --build build --parallel
export ONCVPSP_LIBRARY=$PWD/build/lib/liboncvpsp.so
julia --project=julia/ONCVPSP
```

or add it, as a development package, to an environment of yours:

```julia
using Pkg
Pkg.develop(path="/path/to/oncvpsp/julia/ONCVPSP")
```

The package opens the library when it is loaded (`using ONCVPSP`):

- If `ONCVPSP_LIBRARY` is set, it opens that file. It must be set before
  loading, in the shell or with `ENV["ONCVPSP_LIBRARY"] = ...` in the
  session (or in `~/.julia/config/startup.jl`).
- Otherwise it looks for `liboncvpsp` on the system's library search path
  (`LD_LIBRARY_PATH`, or an installation prefix the linker searches).

As in Python, it checks the library's interface version, and the
[libquadmath note](#conda-libquadmath) applies to a Julia of its own
toolchain too.

A whole run is `oncvpsp(...)`; `ONCVPSP.run` is the same function, not
exported because `Base` has a `run`. For the same reason the channel's
handle type is `ONCVPSP.Channel`. A run's file is written with Base's
`write("Si.upf", run)`. The package's tests run under ctest
(`ctest --test-dir build -R julia_wrapper`).

### Example

The same steps as the Python example, as one script run from the top of
the source tree. Julia's indices start at 1, so `l` is at index `l+1`, the
first level at 1, and so on: `cross[1, 3, 1, 1, 1]` below is Python's
`cross[0, 2, 0, 0, 0]`, with the same value.

The library's build, and a whole run, an `AbstractDict` of its fields:

```julia
using ONCVPSP

println(build_info())  # Dict("commit" => ..., "source" => ..., ...)

run = oncvpsp("tests/data/toml/14_Si.toml"; irel=:sr)
dims = run.dims                     # mmax, lmax, lloc, nproj, nkap, ...
println(dims.mmax, " ", dims.lmax, " ", dims.nproj[1:dims.lmax+1])
println(length(run), " fields: ", first(keys(run), 3))
display(fieldinfo("pot.vkb"))       # its axes, units and what it is
vkb = run["pot.vkb"]                # (mmax, mxprj, mxl+1, nkap)
```

The physics, without the indices (`l` is the angular momentum here, not
an index):

```julia
rr = mesh(run)                      # bohr
vloc = local_potential(run)         # Ha
energies, beta = projectors(run, 0) # l=0: (nproj,), (mmax, nproj)
kk = kernel(run, 0)                 # sum_i beta_i(r) e_i beta_i(r')
rho = 4π .* rr .^ 2 .* valence_density(run)
println("valence charge: ", sum(diff(rr) .* (rho[1:end-1] .+ rho[2:end]) ./ 2))
println(convergence(run, 0))        # rows (Ecut, err)
```

The tests and the all-electron atom:

```julia
cross = run["tests.convergence.cross"]  # (Ecut, err), level, iproj, l, ikap
println("Ecut at 1e-3 Ha, l=0, first projector: ", cross[1, 3, 1, 1, 1])
de = run["tests.configs.de"]            # dE_ae, dE_ps, error per config.
println("excitation errors: ", de[3, 1:run["tests.configs.ncfg"]])
println("all-electron eigenvalues: ", run["ae.ea"][:, 1])
```

The text output and the pseudopotential files:

```julia
@assert occursin("Begin DIAGNOSTICS", run.text)
write("Si.psp8", run)
write("Si.upf", run)
```

A scan of `qcut`, with tables merged over a dictionary:

```julia
base = ONCVPSP.TOML.parsefile("tests/data/toml/14_Si.toml")
for qcut in (5.0, 6.0, 7.0)
    pp = Dict("qcut" => [qcut, 5.0, 5.0])  # merged over base's table
    scan = oncvpsp(base; irel=1, pseudopotentials=pp, check=false, text=false)
    if scan.ifail != 0  # at 7.0 the s pseudo-wave has a node
        println(qcut, " failed: ", sprint(showerror, ONCVPSP.error_of(scan)))
    else
        println(qcut, " ", scan["tests.convergence.opt"][1, 1, 1, 1])
    end
end
```

The all-electron atom once, then the stages from it:

```julia
ae = all_electron(base; irel=1)
println("ae total energy: ", ae["ae.etot"])
same = oncvpsp(base; ae=ae, text=false)
@assert same["pot.vkb"] == vkb
for ll in 0:dims.lmax
    ch = channel(base, ae, ll; mode=:projectors, text=false)
    println(ll, " ", size(ch["channel.pswf"]), " ", projectors(ch).energies)
end
pp = Dict("qcut" => [7.0, 5.0, 5.0])
println("node at qcut 7: ", has_node(base, ae, 0; pseudopotentials=pp))
```

The model core again, and the confined atom:

```julia
teter = coupled(run, 3; fcfact=1.5, rcfact=1.2, text=false)
println("pseudo-atom energy, icmod 3: ", teter["epstot"])
conf = confine(run; radii=[3.0, 3.5], walls=())
println("confine: ", conf["confine.nmiss"], " ", conf["confine.max_de"])
```

Failures, whose kind and stage are Symbols, and freeing a handle at the
end of a `do` block or by `close` (otherwise the garbage collector does):

```julia
try
    oncvpsp(Dict("oncvpsp" => Dict("atsym" => "O")); irel=1)
catch err
    err isa OncvpspError || rethrow()
    println(err.kind, " ", err.stage, " ", err.msg)  # input input 42
end
bad = oncvpsp(Dict("oncvpsp" => Dict("atsym" => "O")); irel=1, check=false)
println(bad.ifail, " ", bad.failure.stage)

oncvpsp(base; irel=1, text=false) do again
    println(size(kernel(again, 1)))
end
close(run)
```

# oncvpsp

Pseudopotential development

The official repository of the oncvpsp Fortran code to generate
optimized norm-conserving Vanderbilt pseudopotentials.

## Installation and usage

### With CMake

CMake 3.22 or later is needed.

```bash
# Configure the build (example flags in brackets)
cmake -S . -B build [-DLibxc_ROOT=<path>] [-DBLA_VENDOR=OpenBLAS] [-DCMAKE_INSTALL_PREFIX=<path>]
# Compile ONCVPSP
cmake --build build --parallel
# Run the test suite (or one part of it, with -L unit, sr, r, split, toml, psml or pipe)
ctest --test-dir build --parallel 8
# Install the executables
cmake --install build
```

or, with the presets in `CMakePresets.json` (`default`, `debug`, `release`
and `fetch`), `cmake --preset default`, `cmake --build --preset default` and
`ctest --preset default`. Without a build type, the build is
`RelWithDebInfo` (`-O2 -g` with gfortran), the one the test references are
checked against; `Debug` adds run-time checks and floating-point traps, and
`Release` fast math.

The dependencies are LAPACK and BLAS, Libxc (unless `-DONCVPSP_WITH_Libxc=OFF`),
toml-f for the TOML input (unless `-DONCVPSP_WITH_TOML=OFF`) and xmlf90 for
the PSML output. xmlf90 is always downloaded and built, with a small patch
(`cmake/patches/xmlf90-1.6.3.patch`). For the others,
`ONCVPSP_LAPACK_SOURCE`, `ONCVPSP_Libxc_SOURCE` and `ONCVPSP_TOML_SOURCE` say
where they come from:

- `AUTO` (the default): the system's, if one is found and a program using it
  compiles and links, and otherwise Reference-LAPACK 3.12.1, Libxc 7.0.0 or
  toml-f 0.4.3, downloaded and built with [CPM](https://github.com/cpm-cmake/CPM.cmake)
- `SYSTEM`: the system's, or an error
- `FETCH`: always downloaded and built

LAPACK is searched for with CMake's FindLAPACK (which honours `BLA_VENDOR`,
`BLA_STATIC` and `LAPACK_ROOT`) and then pkg-config, or given directly as
`-DONCVPSP_LAPACK_LIBRARIES="-llapack;-lblas"`. Libxc is searched for as a
CMake package, with pkg-config and then by its files, all under `Libxc_ROOT`
first if it is set; it must be version 3 or later (but not 5.0.0) with its
Fortran interface built by the same compiler. The configuration ends with a
summary of what was found.

Downloading needs network access when configuring. Set `CPM_SOURCE_CACHE`
(e.g. `~/.cache/CPM`) to keep the downloads between build directories.

Every build records the git commit of the source tree, whether the source
differs from that commit (with a hash of the differences), and the build
time. The executables print this on one line after the citation notice, and
with `--version`, after their name and relativistic treatment, for example

    build: oncvpspr.x (fully-relativistic), commit 4b47c7355e556fb098724409daaa436d8216790f, source clean, compiled 2026-10-15T10:07:48Z

Outside a git checkout the commit and source state are `unknown`. Both CMake
and make generate this information with `cmake/BuildInfo.cmake`.

### With make

1. Edit `make.inc` (or copy it from older oncvpsp-3*)

2. `make all >& make.log`

3. `cd tests/data; evaluate TEST.report`

The code provided here is intended for a Linux or Unix system with a Fortran 95 compiler.
Your system must have LAPACK and BLAS installed to compile ONCVPSP, and
gnuplot installed and in your `$PATH` to run it as recommended.  If these are
not already on your system, they are available from <www.netlib.org/lapack>
and <www.gnuplot.info>. As of release 3.2.3 and beyond, Libxc is supported
for versions at least to libxc-2.2.1, available
[here](http://www.tddft.org/programs/octopus/wiki/index.php/Libxc).

To take advantage of the improved graphics in release 3.1.1, gnuplot
version 4.6 (or later) is required (see [release_notes.md](doc/release_notes.md)
for a workaround if this is not available).

To compile, edit `makefile.inc` in this directory appropriately for your
libraries, compiler, and optimization flags. In the main directory,

    make all >& make.log

and the executables `oncvpsp.x` and `oncvpsp_nr.x` should be created in `src`.
As of version 2.0.1, the relativistic executable `oncvpsp_r.x` is also made.
I have most recently been using gfortran-4.7.2, but most of the code has
been successfully compiled using an old ifort compiler.

Make runs `./set_path` script. This sets the correct path for your installation
in the various shell scripts in scripts, and duplicates them in the tests directory.

Make then runs `tests/data/Test.sh` testing a full set of features (as of
3.2.1). For each `<prefix>.dat` file, this produces a `<prefix>.out` file.
This is compared with `tests/ref/<prefix>.out` using the `fldiff` utility from
Abinit, and produces `<prefix>.diff`, and summary lines from all diffs
are collected in `TEST.report`.  If a reported summary error appears to be
non-trivial, look at the `<prefix>.diff`.

If `tests/compare.sh` fails, check the first line of `scripts/fldiff.pl`
and make sure it points to your Perl executable. Running `83_Bi` requires
Libxc, and the program will stop with an error message if this was not included.

Once you have successfully built the code it would be a good idea to
follow the steps outlined in [users_guide](doc/users_guide.md).

`oncvpsp.x` is designed to be run by the script `run.sh` found in `tests`.
You will probably want to create a new directory in the oncvpsp main
directory for your own psp input data and results.  The input should be
named `<prefix>.dat`, and the output will be `<prefix>.out`. The command
in your new directory

    ../scripts/run.sh <prefix>

should run the code (5-10 seconds?), create `<prefix>.out` in your directory,
and start gnuplot to lead you through a series of plots ("hit enter to
continue") to evaluate your results.
The command

    ../scripts/replot.sh <prefix>

will do the obvious, and the command

    ../scripts/extract.sh <prefix> <pspdirectory>

will cause an ABINIT-readable pseudopotential file in the psp8 format, or a
Quantum ESPRESSO-readable pseudopotential in the UPF format to be created in the
`~/<pspdirectory>`.  The plotting and pseudopotential data are all saved in the
`<prefix>.out` file.  The last input datum on the first data line determines
which format is produced by ONCVPSP.  A copy of the input data is also
saved at the head of the `<prefix>.out` file.  The remaining text above the
"DATA FOR PLOTTING" line reports a variety of consistency checks, information
about the pseudo wave functions, diagnostic checks, and convergence
information (now also plotted).

To use relativistic psps from release 2.0.2 with abinit, the supplied patch
must be applied to abinit-7.4.3, the current release (01/20/2014).  The
patch should be copied into the abinit-7.4.3 main directory, the command

     patch -p0 <abinit-7.4.3_psp8.patch

executed, and `make` re-run.  The new 7.6.1 patch works with 7.6.1 and 7.6.2.
7.6.3 and beyond no longer need patching.

The non-relativistic executable `src/oncvpsp_nr.x` is run by `tests/run_nr.sh`,
and the fully-relativistic one `src/oncvpsp_r.x` by `tests/run_r.sh`.
Please see the files in `doc` directory for details on the code, input file format, etc.

### Solving the all-electron atom once

When the same atom is regenerated many times with different
pseudopotential parameters (for example in an optimization loop), the
all-electron atom can be solved once and reused. For the
scalar-relativistic code,

    oncvpsp_ae.x --ae-file <prefix>.ae < <prefix>.dat > <prefix>_ae.out
    oncvpsp_ps.x --ae-file <prefix>.ae < <prefix>.dat > <prefix>.out

and likewise `oncvpspnr_ae.x` and `oncvpspnr_ps.x` for the
non-relativistic code and `oncvpspr_ae.x` and `oncvpspr_ps.x` for the
fully-relativistic code. `oncvpsp_ae.x` solves the reference atom and the
all-electron halves of the test configurations and writes them to the
`--ae-file`. `oncvpsp_ps.x` reads that file instead of solving the
all-electron atom, and otherwise gives the same output as `oncvpsp.x`.
Its input may change anything except the reference atom (`atsym`, `z`,
`nc`, `nv`, `iexc` and the reference configuration), and may use any
subset of the test configurations that `oncvpsp_ae.x` was given; it stops
with an error otherwise, or if the file was written for a different
relativistic treatment.

Setting `psfile` to `none` skips writing a pseudopotential file, and
`ncnf` = 0 skips the test configurations other than the reference.

## Documentation

Documentation available at [https://oncvpsp.github.io/oncvpsp/](https://oncvpsp.github.io/oncvpsp/)

## How to cite oncvpsp

If you use oncvpsp in your research, please consider citing the
[following work](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.88.085117):

> Optimized norm-conserving Vanderbilt pseudopotentials
    D. R. Hamann
    Phys. Rev. B 88, 085117 (2013)
    10.1103/PhysRevB.88.085117

```
@article{PhysRevB.88.085117,
  title = {Optimized norm-conserving Vanderbilt pseudopotentials},
  author = {Hamann, D. R.},
  journal = {Phys. Rev. B},
  volume = {88},
  issue = {8},
  pages = {085117},
  numpages = {10},
  year = {2013},
  month = {Aug},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevB.88.085117},
  url = {https://link.aps.org/doi/10.1103/PhysRevB.88.085117}
}
```

## License

See the [License File](./COPYING).

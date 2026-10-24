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
# Run the test suite (or one part of it, with -L unit, sr, r, nr, split, mesh, toml, psml, pipe,
# errors, echo, psp8 or nist)
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

Outside a git checkout the commit and source state are `unknown`. CMake
generates this information with `cmake/BuildInfo.cmake`.

To contribute, install the [pre-commit](https://pre-commit.com) hooks with
`pre-commit install`; they lint every commit, as CI does. See
[Developer Setup](doc/developer_setup.md) for the tools, the editor
integration and the builds, and [Coding Standards](doc/coding_standards.md)
for the rules.

### Running

To run ONCVPSP as recommended, gnuplot (4.6 or later) must be in your
`$PATH`: `scripts/run.sh` leads you through its plots. The scripts find
the executables in `build/bin` of the CMake build of this tree, or in
`$ONCVPSP_BIN_DIR`.

Once you have successfully built the code it would be a good idea to
follow the steps outlined in [users_guide](doc/users_guide.md).

`oncvpsp.x` is designed to be run by the script `run.sh` found in `scripts`.
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

The non-relativistic executable `oncvpspnr.x` is run by `scripts/run_nr.sh`,
and the fully-relativistic one `oncvpspr.x` by `scripts/run_r.sh`.
Please see the files in `doc` directory for details on the code, input file format, etc.

The executables read the text input from stdin, or from a file with
`-i FILE`; `-t FILE` reads the same input as a TOML file instead. The text
input is described in `doc/32_Ge_annotated.dat`, and the TOML input, with
`scripts/input_text_to_toml.py`, which converts one to the other, in the
[users' guide](doc/users_guide.md). Every executable lists its options
with `--help`.

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
relativistic treatment or on a different radial mesh.

### Log derivatives at any radii

The log derivatives (`run_phsft`) can be taken out
of the run and computed afterwards, at any radii and energies:

    oncvpsp_ps.x --ae-file <prefix>.ae --ps-file <prefix>.ps \
        --no-log-derivatives < <prefix>.dat > <prefix>.out
    oncvpsp_ld.x -a <prefix>.ae -p <prefix>.ps -r 2.5,3.0,3.5 > <prefix>_ld.out

`--ps-file` (also on `oncvpsp.x`, `oncvpspnr.x`, `oncvpspr.x` and the
other `_ps.x`) writes the pseudopotential data file, and
`--no-log-derivatives` leaves the log derivatives out of the output.
`oncvpsp_ld.x`, for every relativistic treatment, reads the two files and
prints the `ld_r` and `ld` lines of the DIAGNOSTICS block for each l, j
and radius; one solution per energy serves all the radii. Without `-r`
(radii) and `-e` (`emin,emax,de`) it uses the input's `rxpsh`, `epsh1`,
`epsh2` and `depsh`, and reproduces the run's log derivatives; `-l`
(`l1,l2,...`) computes only those angular momenta.

### The confined atom

The same two files give the pseudopotential and its all-electron atom in a
Wigner-Seitz sphere, a stand-in for the atom in a solid:

    oncvpsp_confine.x -a <prefix>.ae -p <prefix>.ps -r 3.0,3.5 > <prefix>_ws.out

For each l, j and sphere radius `rs` (`-r`, by default 1.25, 1.5 and 2
times the largest core radius) it finds the band edges of the
Wigner-Seitz rule in the energy window (`-e emin,emax,de`, by default the
input's `epsh1` and `epsh2`). The bottom of a band is where r phi'/phi = 0
at `rs`, its centre where it is -l-1, and its top where phi = 0. The
all-electron and pseudopotential energies are printed side by side (`ws`
lines), with n the all-electron state's nodes inside `rs` plus l+1. For
each pseudopotential edge it also prints the cutoff at which the state's
energy in the sphere is converged to 1e-2 ... 1e-6 Ha (`ws_ecut`). This is
the variational energy in the spherical Bessel functions of the same
boundary condition, so channels with no occupied state get a cutoff at the
energies of the valence bands too. `ws_floor` is the error at the largest
cutoff (`-c`, by default 200 Ha), `-t` prints the error at every basis
size (`ws_conv`), and `-l` restricts the angular momenta. The edges come
from the roots of monotonic phases (atan(r phi'/phi) less pi per node
inside `rs`), so they do not depend on the step `de` of the scan that
brackets them. The pseudopotential's phase counts the poles of its
r phi'/phi at `rs` below the energy rather than its nodes, since its wave
function can gain a node inside the projectors. The poles are the
eigenvalues of the sphere with phi(rs) = 0 in the same basis as the cutoffs,
each refined to the root of u(rs). Each all-electron edge's partner is the
pseudopotential edge of the same phase. What does not pair is named: `ws_missing` for an
all-electron edge with no partner, `ws_core` for one that is the band of a
core state, `ws_extra` for a pseudopotential edge with no all-electron
partner (a ghost in the window), and `ws_paired_off` for a pseudopotential
edge of another phase nearer than the partner.

It then compresses the atom. It solves the self-consistent all-electron
atom and pseudo-atom of the reference configuration in a smooth wall,
1/2 k x^2 exp(-0.5/x) with x = r - r0 beyond `r0` (`-R`, by default the
largest core radius), for each strength k (`-w`, by default 0.01, 0.03
and 0.1 Ha/bohr^2). The `wall` lines give the two total energies and
their changes from the free atom, which should agree as closely as the
pseudopotential is transferable. The `wall_eig` lines give the valence
eigenvalues. `-w none` leaves the walls out and `-s` leaves out the
sphere. The last line, `confine_summary max_de ecut4 ecut5 n_missing n_extra
n_paired_off n_wall_failed max_dde`, gives:

- over the edges below 0, the largest |e_ps - e_ae|, and the largest
  cutoffs at 1e-4 and 1e-5 Ha (Infinity if one does not reach it);
- the counts of the `ws_missing`, `ws_extra`, `ws_paired_off` and
  `wall_failed` lines;
- the largest |de_ps - de_ae| of the walls.

A part that is left out gives NaN and -1. The pseudopotential file has carried what the pseudo-atom needs
(format version 2) since this tool was added, so older `.ps` files must
be written again.

Setting `psfile` to `none` skips writing a pseudopotential file, and
`ncnf` = 0 skips the test configurations other than the reference.

### The radial mesh

The calculations use the log radial mesh `rr(i) = rr1 * amesh**(i-1)`
out to `rmax`. By default `amesh` = 1.006, `rmax` = 45 bohr, and
`rr1` = 0.0005/z, capped at 0.00005 except in the fully-relativistic
code. An optional last line of the text input, after the test
configurations, sets them:

    # amesh, rr1, rmax
      1.004  -1  50.0

A value of 0 or less keeps the default, so this line keeps `rr1`. The
TOML input takes them in an optional `[log_mesh]` table, where a missing
key also keeps the default. The output, and the input copied into the
psp8 and UPF files, give the mesh (all three values, as used) only when
the input sets it, so the output for an input without it is unchanged. `rmax` must be larger than `rlmax`, `rr1`
smaller than the largest rc, and `amesh` larger than 1.

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

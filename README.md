# oncvpsp

Pseudopotential development

The official repository of the oncvpsp Fortran code to generate optimized norm-conserving Vanderbilt pseudopotentials.

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

## Installation

### Cloning the repository
This repository contains git submodules for the external dependencies of oncvpsp (BLAS/LAPACK and Libxc).
To clone the submodules with the repository, git requires the `--recursive` flag:

```bash
git clone https://github.com/azadoks/oncvpsp.git --recursive && cd oncvpsp
```

### Building with CMake
You can usually use the following recipe:

> To install in a custom directory, use the `-DCMAKE_INSTALL_PREFIX` flag with `cmake`

> Libxc fallback compilation is currently broken! Compile Libxc (>=v5.0.0 and <7.0.0) and provide the flag `-DLibxc_ROOT="/path/to/libxc"`! 

```bash
mkdir build && cd build
cmake .. -DLibxc_ROOT="/path/to/libxc"
make
make install
```

If you don't want to use the external fallbacks and `cmake` fails to find one of the dependencies, then you can point it to where they are installed by adding one or more flags to the command above:

> You might have to delete the contents of your build directory for the flags to have an effect!

* Libxc: `-DLibxc_ROOT="/path/to/libxc"`
* BLAS: `-DBLAS_ROOT="/path/to/blas"`
* LAPACK: `-DLAPACK_ROOT="/path/to/lapack"`

You can also configure your prefered flavor of BLAS/LAPACK using the `-DBLA_VENDOR` flag.
Supported vendors are listed in the [CMake documentation](https://cmake.org/cmake/help/latest/module/FindBLAS.html#blas-lapack-vendors).

### Running tests
Following a successful build, run the following from the root of the repository:

```bash
./set_path
cd tests/data
./TEST.sh
```


The code provided here is intended for a Linux or Unix system with a Fortran 95 compiler.
Your system must have lapack and blas installed to compile ONCVPSP, and 
gnuplot installed and in your $PATH to run it as recommended.  If these are
not already on your system, they are available from <www.netlib.org/lapack>
and <www.gnuplot.info>. As of release 3.2.3 and beyond, libxc is supported 
for versions at least to libxc-2.2.1, available 
[here](http://www.tddft.org/programs/octopus/wiki/index.php/Libxc).

To take advantage of the improved graphics in release 3.1.1, gnuplot
version 4.6 (or later) is required (see [release_notes.md](doc/release_notes.md)
for a workaround if this is not available).

To compile, edit makefile.inc in this directory appropriately for your 
libraries, compiler, and optimization flags. In the main directory,

    make all >& make.log

and the executables oncvpsp.x and oncvpsp_nr.x should be created in src. 
As of version 2.0.1, the relativistic executable oncvpsp_r.x is also made.
I have most recently been using gfortran-4.7.2, but most of the code has 
been successfully compiled using an old ifort compiler.

Make runs *./set_path* script. This sets the correct path for your installation 
in the various shell scripts in scripts, and duplicates them in the tests directory.

Make then runs *tests/data/Test.sh* testing a full set of features (as of 
3.2.1). For each *<prefix>.dat* file, this produces a *<prefix>.out* file.  
This is compared with *tests/ref/<prefix>.out* using the *fldiff* utility from
Abinit, and produces *<prefix>.diff*, and summary lines from all  diffs
are collected in TEST.report.  If a reported summary error appears to be 
non-trivial, look at the *<prefix>.diff*.

If *tests/compare.sh* fails, check the first line of *scripts/fldiff.pl*
and make sure it points to your perl executable. Running *83_Bi* requires
libxc, and the program will stop with an error message if this was not included.

Once you have successfully built the code it would be a good idea to
follow the steps outlined in [users_guide](doc/users_guide.md).

oncvpsp.x is designed to be run by the script run.sh found in tests.
You will probably want to create a new directory in the oncvpsp main
directory for your own psp input data and results.  The input should be 
named *<prefix>.dat*, and the output will be *<prefix>.out*. The command 
in your new directory

    ../scripts/run.sh <prefix>

should run the code (5-10 seconds?), create *<prefix>.out* in your directory,
and start gnuplot to lead you through a series of plots ("hit enter to
continue") to evaluate your results.
The command

    ../scripts/replot.sh <prefix>

will do the obvious, and the command

    ../scripts/extract.sh <prefix> <pspdirectory>

will cause an ABINIT-readable pseudopotential file in the psp8 format, or a
PWSCF-readable pseudopotential in the UPF format to be created in the 
*~/<pspdirectory>*.  The plotting and pseudopotential data are all saved in the 
*<prefix>.out* file.  The last input datum on the first data line determines
which format is produced by ONCVPSP.  A copy of the input data is also
saved at the head of the *<prefix>.out* file.  The remaining text above the 
"DATA FOR PLOTTING" line reports a variety of consistency checks, information
about the pseudo wave functions, diagnostic checks, and convergence
information (now also plotted).

To use relativistic psps from release 2.0.2 with abinit, the supplied patch
must be applied to abinit-7.4.3, the current release (01/20/2014).  The
patch should be copied into the abinit-7.4.3 main directory, the command

     patch -p0 <abinit-7.4.3_psp8.patch

executed, and "make" re-run.  The new 7.6.1 patch works with 7.6.1 and 7.6.2.
7.6.3 and beyond no longer need patching.

The non-relativistic executable *src/oncvpsp_nr.x* is run by *tests/run_nr.sh*,
and the fully-relativistic one *src/oncvpsp_r.x* by *tests/run_r.sh*.
Please see the files in doc directory for details on the code, input file format, etc.

## Documentation

* [User Guide](./doc/users_guide.md)
* [Ghost States](./doc/ghosts.md)
* [Relativistic](./doc/relativistic.md)
* [Core Correction](./doc/core_correction.md)
* [How to use libxc](./doc/libxc_use.md)
* [How to translate iexc to pwscf names](./doc/pwscf_exc.md)
* [Program Overview](./doc/program_overview.md)
* [Coding Standards](./doc/coding_standards.md)
* [Support](./SUPPORT.md)
* [Release Notes](./doc/release_notes.md)

## License

See the [License File](./COPYING).

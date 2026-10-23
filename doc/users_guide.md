## ONCVPSP USERS' GUIDE

### Introduction

There is no detailed users' manual, and this is primarily an introductory
tutorial for new users.  With the included Phys. Rev. paper, the
program overview, the simple structure of the input data, the examples
provided, and the automated graphics, it should be relatively easy to get
to know and use this code.  As of release 4.0.1, most of the potental
problems disussed in the last few paragraphs are highly improbable.

You should start by reading the paper, and some of the references if
necessary.  The input and output files in tests/data and tests/refs include
all the "OV" pseudopotentials used in the Sec.IV tests in the paper. Next,
I recommend reading the README, and [program_overview](program_overview.md).

The gnuplot graphics are very important, and a good place to go next.
The file 73_Ta_plots.pdf in this directory is a "slide show" which walks you
through a sample of what you will see whenever you run the program.  After
viewing the pdf, go to directory tests/refs and type

    ../../scripts/replot.sh 73_Ta

If this doesn't produce the same sequence of plots, stop here and try to fix
your gnuplot installation and accessibility.  Once this works,  it will pay to
run replot for all the examples in this directory.  This should illustrate what
your own output plots should look like.  If your own data give results that
look odd by these standards, they probably are and should not be used until
adjustments of your input data can produce "nice" plots as exemplified by
these cases.

The input files for oncvpsp, examples of which you'll find in
tests/data, are NOT THE KEYWORD FILES YOU MAY BE USED TO.  The # comment
lines define a strict template in which ALL input parameters must be specified,
in the correct order.  The only defaults are those of optional values at the
end of some lines, and of the optional last line (the log mesh).  (On the
other hand,  there are
not so many input parameters, and you need to understand them all to get good results anyway.)
The same input can be given as a TOML keyword file instead (see "TOML
input" below).

The file `32_Ge_annotated.dat` adds more comments (###) to the usual
`*.dat` file describing each input variable.  Reading this should be your next
assignment.  The code doesn't mind extra comments, and you can in fact
use this as a combined template and abbreviated users' manual for your
own calculations instead of the shorter versions in tests/data.

Next, look at `tests/refs/32_Ge.out`.  The output is pretty-much
self-explanatory.  First, the input data is echoed, adding the all-electron
eigenvalues to the configuration information.  This is formatted so that it
can be copied and re-used as an input file (the eigenvalues will be ignored).
The psp8, UPF and PSML files carry the same echo, without the eigenvalues,
and each of these echoes runs again to the same output. Optional inputs
(`c_light`, the optional values of the lmax line, `rxpsh`, the `icmod=4`
search grid, `slack` and `budget_cutoff` of `icmod` 6 and 7, the log mesh
and `upf_mesh`) are echoed only when the input sets them, so older inputs
echo as they always have.
The code then loops over angular momenta and reports some information and
consistency test results for each pseudo wave function. For unbound states,
the 'half-point radius' is that at which the added barrier has reached half
its asymptotic value. A well state is solved on the log mesh continued
past its end (`rmax`) as far as the state needs to decay to machine
precision, and the output says how far ("mesh continued to r ="). The
well rises only 0.5 Ha above the state, so that is about 36 bohr past
the wall, and far more for a wide well (a small energy); the end of the
mesh is then not a box for it, and nothing depends on `rmax`.  The qq values are the generalized norm-conserving
conditions, Eq.(23).  The constraint-consistency test should be very well
satisfied, or there is an error in the code. The following output is just
information ("Fraction of norm ...").  For the kinetic energy and potential
tests, the ratio should be close to one, and the difference should be small.
For a well state confined by its well, a line compares its residual
kinetic energy at qcut with the well's tail and with the free tail at
the same energy (that of `well_tail = "free"`): inside rc it is the same
wave, so the difference is the well's alone. A WARNING marks a change of
more than 10 %: the cutoffs of that wave then carry the arbitrary well,
which `well_tail = "free"` removes. Another WARNING marks a well state
whose wave vector is not below qcut: it oscillates faster than the
cutoff, and no tail makes its cutoffs meaningful.

`well_tail = "free"` needs every well state at ep >= 0. Below 0 the free
tail grows, and its residual kinetic energy, computed from the transform
of the windowed tail, carries round-off that grows with the tail's peak
squared: 2.5e3 Ha for an I d state at -1.6 Ha, 5e20 Ha for a Cs p state
at -0.49 Ha. Such an input is an error; raise the state's ep (`debl`),
or use `well_tail = "confined"` or `"well"`.

Last comes the cutoff convergence profile for the current wave function.

After the pseudo-wave-function loop is completed, the projectors are
constructed.  The Hermiticity error is discussed in the paragraph following
Eq.(24).  It will be greatly reduced in a non-relativistic run. The data on
the projectors is not of much interest, but if the coefficients are huge or
the overlap is very close to one, it probably pays to change some parameters.
As of release 2.0.1, orthonormal projectors are constructed by default,
so the overlap is no longer reported.

The "Diagnostic tests" results are explained in more detail in
[program_overview](program_overview.md).  Basically, e in and e test should closely match, and
the norm and slope values should be close to one.  For two projectors and a
non-relativistic calculation, the agreement should be within numerical
convergence error for all these quantities.  That error is O(al^4) in the
mesh step: the projectors end at rc with a kink, and the pseudo-atom solver
(vkboutwf) takes their overlaps over their support, as the B matrix is
built, and does not carry the kink into the steps beyond rc, either of which
made it O(al^2) (-40 uHa for the 3d of Fe at amesh 1.006).  Scalar-
relativistic runs keep a mesh-independent difference of a few uHa from the
relativistic correction to the projectors (vrel), which run_vkb also takes
out of a semi-local local potential (lloc<4), as it is vae+vr at rc.  Finally,  the total energies
and eigenvalues of self-consistent pseudo atoms are compared with all-electron
results for the reference and test valence configurations.  You should
always look at these diagnostics, even if the automated graphics all look
good.

Everything below the "DATA FOR PLOTTING"line is for the scripts to read,
not you (plotting data, template for a custom gnuplot input script, and
the psp file).

After you've gotten through the above and have compiled the code
following the README, copy one of the `<prefix>.dat` files in tests/data to a
directory of your own, run it using

    ../scripts/run.sh <prefix>

The sample `<prefix>.dat` files in tests/data should provide a starting
point to produce good psps for new elements. Reading the `<prefix>.out` files,
you should be sure that the various consistency tests look good compared to
the examples. Then, pay attention to the plots.

### TOML input

Every executable but `oncvpsp_ld.x` reads its input from a TOML file with
`-t FILE` (`--toml-input`) instead of the text input on stdin (or of
`-i FILE`). The TOML input is part of the build unless CMake is configured
with `-DONCVPSP_WITH_TOML=OFF`. It holds the same values as the text input,
with the same meanings, limits and defaults, as keys of named tables, so
their order does not matter; a missing table or key that is not optional,
or one of the wrong type, stops the program with a message naming it. So
does a real that is not finite (`nan` or `inf`), here and in the text input,
where it is the read error of its line. The lists of the per-l rows hold one element per l, from 0 to lmax. The output,
and the echo of the input in it and in the pseudopotential files, is that
of the text input. tests/data/toml has the TOML form of every input of
tests/data, and

    scripts/input_text_to_toml.py -i <prefix>.dat -o <prefix>.toml

converts a text input (with `--check`, it also reads the result back and
compares it to the text input). The tables and keys, with the names of the
text input where they differ (see `32_Ge_annotated.dat`):

| table | keys |
|---|---|
| `[oncvpsp]` | `atsym` (string), `z`, `nc`, `nv`, `iexc`; `c_light` optional |
| `[reference_configuration]` | `n`, `l`, `f`, lists of nc + nv |
| `[pseudopotentials]` | `lmax`, and the lists `l`, `rc`, `ep`, `ncon`, `nbas`, `qcut`; optional `ep_as_given` (true or false), `well_onset`, `well_tail` (`"well"`, `"free"` or `"confined"`), `wave_regularization`, `well_norm` (`"space"` or `"rc"`; by default that of `well_tail`) |
| `[local_potential]` | `lloc`, `lpopt`, `rcloc` (the text input's rc(5)), `dvloc0` |
| `[vkb_projectors]` | the lists `l`, `nproj`, `debl`; `well_n` optional, with `well_tail = "confined"` only; `ecut_weight` optional (negative for the default), for channels with no occupied state only |
| `[model_core_charge]` | `icmod`; `fcfact` (required for icmod 1, 2, 3, 6, 7); `rcfact` (required for icmod 3; optional for 4, 6, 7); for icmod 4 the optional grid `fcfact_min`, `fcfact_max`, `fcfact_step`, `rcfact_min`, `rcfact_max`, `rcfact_step`; for icmod 6 and 7 the optional `slack` and `budget_cutoff` (`"max"`, `"weighted"` or `"max_all"`) |
| `[log_derivative_analysis]` | `epsh1`, `epsh2`, `depsh`; `rxpsh` optional |
| `[linear_mesh]` | `a` (drl), `rmax` (rlmax) |
| `[log_mesh]` | optional table: `amesh`, `rr1`, `rmax`, each optional |
| `[pp_output]` | `psfile` (string); `upf_mesh` optional |
| `[[test_configurations]]` | optional, one table per test configuration (at most 4): `n`, `l`, `f`, lists of its valence states only |

For example, tests/data/toml/32_Ge.toml, with two of its four test
configurations:

```toml
[oncvpsp]
atsym = "Ge"
z = 32.0
nc = 5
nv = 3
iexc = 3

[reference_configuration]
n = [1, 2, 2, 3, 3, 3, 4, 4]
l = [0, 0, 1, 0, 1, 2, 0, 1]
f = [2.0, 2.0, 6.0, 2.0, 6.0, 10.0, 2.0, 2.0]

[pseudopotentials]
lmax = 2
l = [0, 1, 2]
rc = [2.6, 2.6, 2.0]
ep = [-0.0, -0.0, 0.0]
ncon = [4, 4, 4]
nbas = [8, 8, 9]
qcut = [5.0, 5.2, 8.4]

[local_potential]
lloc = 4
lpopt = 5
rcloc = 2.0
dvloc0 = 0.0

[vkb_projectors]
l = [0, 1, 2]
nproj = [2, 2, 2]
debl = [1.5, 1.5, 1.5]

[model_core_charge]
icmod = 0

[log_derivative_analysis]
epsh1 = -2.0
epsh2 = 2.0
depsh = 0.02

[linear_mesh]
a = 0.01
rmax = 4.0

[pp_output]
psfile = "both"

[[test_configurations]]
n = [3, 4, 4]
l = [2, 0, 1]
f = [10.0, 1.0, 2.0]

[[test_configurations]]
n = [3, 4, 4]
l = [2, 0, 1]
f = [10.0, 2.0, 1.0]
```

### Potential issues

1) The code will stop if it doesn't find an input data entry it is looking for.
   It will also run some checks to be sure the data is within bounds. The first
   stop will just give you a line number, but the error might be above that.
   The second kind should give more explanation.

2) The code will stop for a number of other error conditions you shouldn't
   encounter if it is built properly and your parameters are reasonable.
   All stops should give some information.  One stop you probably will
   encounter is

        "pspot:  first pseudo wave function has node."

   While in principle the Vanderbilt method can deal with this, it isn't a
   good idea and it's better to have a semi-local potential to look at to make
   sure it looks reasonable.  The fix in this case is usually to increase qcut
   a bit, although changes in nbas, rc or ep for a scattering state may also
   work.

   Another you might encounter is

        "optimize: norm constraint violated (line 260)"

   which will probably only occur if you have an unphysically small rc,
   indicated by a very small value for "Fraction of norm inside rc."  The
   code now gives "First projector wave function outermost peak radius,"
   and rc should be > that value (usually >1.5 * r_peak, see examples).

   Another is

       "wellstate: well potential iteration failed"

   which should be fixable by changing ep (for the first projector) or
   a combination of ep and debl for the second.  rc can also influence
   this.

3) A related problem is an unreasonably large semi-local 1=0 potential at r=0.
   Once again, this bad behavior won't carry over into the projector, and
   might not do any real damage (unless lloc=0,  which is always a bad idea
   anyway).  However, it is best to eliminate it by parameter changes like
   those for the nodes.

4) The next-to-last paragraph of Sec. II of the paper remarks on the
   "strong short-wavelength oscillations" discussed by Troullier and
   Martins as an issue with the Rappe et al. optimization scheme.  While
   I have not observed these for the usual ranges of parameters in the
   potentials I've generated with the present optimization approach,
   it can be made to happen.  It is certainly a good idea to eliminate them
   by adjusting the usual suspects (qcut, nbas, rc, ep).

5) The ncon continuity-constraints may be less important than I initially
   thought, and don't have that much influence on the results (see
   Fig. 3).  I prefer to use ncon = 4 or 5.  When looking at pseudo
   wave functions, a higher the degree of continuity does allows a larger
   rc for comparable agreement with the all-electron wave function.

6) Choosing qcut is a compromise.  While a smaller qcut will typically give
   better convergence up to the corresponding energy, the "foot" of the
   convergence curve  (see Fig. 4) will have an earlier onset and be flatter.
   The cutoff convergence profiles in the output (and in the graphics
   as of release 2.0) can guide you.  Increasing qcut by a small amount,
   say ~0.5, can usually give more uniform convergence at only a small
   sacrifice of low-cutoff convergence.

7) I have had trouble with pseudo atom calculations failing to converge in
   the configuration-testing section.  The code won't stop, and will skip
   the "bad" configuration.  These problems seem to arise from the pseudo
   wave functions acquiring nodes at small radii, which can happen for the
   fully-non-local radial Schroedinger equation.  My present solution for
   this is to only start counting nodes for r>0.5 a_B (in the routine vkboutwf).
   [Since then: such a node enters at r = 0 as the energy changes, without
   an eigenvalue passing, and a count from a fixed radius jumps where it
   crosses that radius. The pseudo-atom solvers with projectors (lschvkbb,
   lschvkbbe) count the eigenvalues below the energy instead (outwf_index
   in outwf_m), which no such node changes.]

8) The shell scripts run.sh run_nr.sh etc. find the executables in
   build/bin of the CMake build of the tree they are in.  If you move the
   executables elsewhere, set ONCVPSP_BIN_DIR to their directory.

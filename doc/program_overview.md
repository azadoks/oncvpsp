## OVERVIEW OF ONCVPSP

References to "the paper" are to oncvpsp6.pdf, which is reproduced in this
directory (now published as Phys. Rev. B 88, 085117 (2013)).

The calculation is the subroutine oncvpsp_main (oncvpsp_main.F90).  Each
executable (oncvpsp.x, oncvpspnr.x, oncvpspr.x and their _ae.x and _ps.x
stages) is a short program calling it with the relativistic treatment,
irel, and the stage.  It has the following major sections:

1) Data readin: Oncvpsp reads an input file, say *32_Ge.dat*, from its standard
   input.  While this file may contain comment lines with the initial character
   '#', its format is fixed, and input data must appear in the proper order,
   line by line.  Examples in tests/data as well as *32_Ge_annotated.dat* in
   this directory should clarify this.  The program will stop, identifying
   the first bad line encountered in the input file if it does not find what
   it expects.  After a successful readin, the data is scanned by the
   routine *check_data* for violations of bounds, etc., and offending data will
   be identified prior to halting.  There are, quite deliberately, no defaults.
   Everything must be thought about and specified (but there is not that much
   input data).

2) The routine sratom (relatom for the fully-relativistic code, both
   called through aeatom) will then run to perform a self-consistent scalar-
   relativistic calculation for the reference configuration of the atom,
   as described in the input.  The screening (Hartree plus exchange-
   correlation) potential is mixed, and the SCF is converged when no
   eigenvalue moves and the root mean square change of the screening
   potential over the electrons (ae_residual, as in RadialAtomicSolvers.jl)
   is below 1e-8 Ha.  The radial solvers (lschfb, ldiracfb) take two more
   steps from RadialAtomicSolvers.jl.  They start the inward integration
   where the WKB decay from the turning point reaches machine precision
   (inward_start); the start at 10 times the turning point used before cut
   the tail of a shallow 3d or 4f state, which made the energy first order
   in the mesh.  And an occupied state with no bound solution in an early
   SCF potential (the 4f of Yb) gets its box state, zero at the end of the
   mesh, and the SCF goes on; it cannot converge with one.  The scalar-
   relativistic SCF with a GGA also follows RadialAtomicSolvers.jl (its note
   SR_SHOOTING_VPRIME.md), and converges on meshes down to amesh 1.00075
   and rr1 1e-8/Z, where it diverged for light atoms.  Inside
   r* ~ Z/(2c^2) the relativistic mass is set by the potential, and with a
   GGA that closes a loop potential -> orbital exponent -> density gradient
   -> potential, whose modes grow with the mesh points per decade.  lschfb,
   and lschfs for the scattering states, integrate the Koelling-Harmon
   equation in its first-order flux form (sr_flux_m), which needs no dV/dr;
   they start from the series at the origin of that
   equation for -Z/r plus a constant (with the exponent
   sqrt(l(l+1)+1-(Z/c)^2)), whose logarithmic derivative at rr1 removes
   what was a 1e-5 Ha error at the default rr1 for Na.  That series converges
   only inside r = Z/(2c^2), for every l, so points beyond a quarter of that
   (rr1 of the default mesh for H) are reached from there by RK4 in ln r on
   the same model equation.  lschfb converges
   each eigenvalue to 1e-12 relative (1e-10 left 4e-7 Ha in the 1s of U).
   The mixing (ae_mix) is Anderson over the last 8 iterations, with the
   residual filtered in ln r inside 10 r* (KerkerMixing).  The built-in PBE
   (excggc) takes its potential in the divergence form of the Libxc path,
   from the first derivative of the density only; the second derivative it
   used was round-off noise at the first points of a fine mesh.  The total
   energy and its kinetic, electron-nucleus,
   Hartree and exchange-correlation terms (ae_energies, as NIST defines them,
   the kinetic energy from the eigenvalue sum) are printed with the
   reference-configuration results.  I strongly believe that neutral ground-state
   reference configurations provide the best pseudopotentials for solids,
   and oncvpsp incorporates a number of features to facilitate their use.
   After some bookkeeping, the main program echoes the input data in a format
   similar to a typical input file, with the addition of the reference-
   configuration eigenvalues to the configuration description.  All output
   printing is to the standard output.

3) Next the main program loops through angular momenta from 0 to lmax to
   calculate one or two optimized pseudo wave functions.  For each angular
   momentum, there are four major sub-sections:

  A) The all-electron wave function is calculated for the lowest occupied
     state declared as valence or for a barrier-confined "well state" for a
     specified positive energy corresponding to a scattering state
     (subroutines lschfb and wellstate).  The barrier potential is given
     in Eq.(17) of the paper.  To find a bound state at the desired
     energy ep, wellstate starts with a large r_b = 8*r_c and searches
     for a suitable v_inf within the range [ep, ep + 1 Ha].  The barrier
     starts at well_onset*r_c (TOML input, at least 1, default 1.006**5,
     where the mesh-dependent onset of before sat on the default mesh);
     it is flat to every order there, so 1 is the natural limit.

     With well_tail = "free" (TOML input; default "well") the well states
     do not depend on the well (free_tail_m): the state is the atom's
     solution at exactly ep, integrated outward (scatterstate) and
     normalized to 1 inside r_c, and beyond r_c it is the free solution
     at ep under an erfc window, t = u w, w = erfc((r - R0)/L)/2,
     R0 = r_c + 6L, with L from qcut and ep so that the window's share of
     the residual energy above qcut is about 1e-13; a negative ep (with
     ep_as_given) gives a growing tail under the same window. wellstate
     then only resolves the energy (its replacements of ep) and builds no
     well, except in the fully-relativistic code, where the j = l-1/2
     state takes the eigenvalue of the j = l+1/2 well as its energy.

  B) If two projectors using the norm-conserving Vanderbilt method are
     specified, the all-electron wave function with one additional node is
     calculated for either a second bound valence state (when the first
     is a shallow core) or for a second "well state" at a specified energy
     above the lowest state at the current angular momentum.

  C) One or two pseudo wave functions are calculated using the RRKJ
     convergence optimization approach.  Both will satisfy continuity
     and norm-conserving conditions with their corresponding all-electron
     wave functions.  The second will have one node and will also satisfy
     the off-diagonal Vanderbilt generalized norm-conserving condition,
     which is dependent of the first pseudo wave function as well as the
     overlap of the two all-electron wave functions.  Details of the
     routine run_optimize are expanded below. Diagnostic information is
     printed.

  D) The semi-local pseudopotentials are calculated.  More diagnostics
     and the convergence profile (Energy error per electron vs. cutoff)
     are printed.

4) The routine *run_vkb* then is then called.  It constructs the local
   potential from a polynomial extrapolation of the all-electron potential
   to zero (the recommended option, especially to take advantage of the
   Vanderbilt two-projector method), or uses one of the semi-local
   potentials, which may be satisfactory when using a single Kleinman=
   Bylander projector.  It then loops over angular momenta, calculating
   the "raw" projectors (chi in Eq.(18) of the paper) and the B matrix
   (Eq.(21)).  For two projectors, any residual asymmetry  of B is
   removed (see the discussion following Eq.(24)), it is diagonalized,
   and the final diagonal and normalized separable non-local potential
   of Eq.(24) is produced.  A small amount of diagnostic information is
   printed.  As of release 2.0, orthonormal projectors are subsequently
   constructed, and overlap information is omitted.  The "OPTION"
   variable in *run_vkb* can be reset to restore the old non-orthogonal
   output, although only for non- and scalar-relativistic cases.

5) The valence pseudocharge density is calculated based on the separable
   potential, and then that of any "extra" valence state (ie., when cores
   are included), using subroutine lschvkbb.  If specified, a monotonic
   polynomial model core charge is calculated which joins smoothly to
   the all-electron core charge at the radius at which this charge density
   equals a specified fraction of the valence pseudocharge density
   (subroutine modcore).  The Hartree and exchange-correlation potentials
   of the pseudocharge are computed (subroutine vout)and the pseudopotentials
   are "unscreened."

   Beyond the largest core radius rcmax the screened potentials are the
   all-electron potential, screened by the all-electron valence charge.
   The pseudo valence charge differs from it there (by about 1e-4 of it
   for a scalar- or fully relativistic reference), so unscreening with it
   would leave the Hartree and exchange-correlation response to that
   difference in the local potential: a long-ranged term that is not
   ionic, and that plane-wave codes see through the G=0 term. Instead
   (subroutine ionic_screening in ionic_tail.f90, after the ionic tail of
   PseudoPotentialGenerators.jl) the screened potentials are switched,
   from rcmax to where the frozen-core ionic potential

       V_ion = -Z/r + V_H[n_c] + V_xc[n_c + n_v] - V_xc[n_v + rhomod]

   is -zion/r to 1e-10 (or to 2 rcmax if that is further out), to
   V_ion + V_Hxc[rho + rhomod], with the C-infinity switch
   f(1-x)/(f(1-x)+f(x)), f(t) = exp(-1/t). The unscreened potentials are
   then V_ion beyond the window, which keeps the core penetration and the
   core-valence exchange-correlation. The pseudo valence states are
   iterated with the switched local potential to a fixed point (1e-10),
   so that the reference configuration stays a fixed point of the
   pseudo-atom (PS_ref - PS_tst of test configuration 0 is round-off).
   The model core is built before, and kept. The out file gives the
   window and the number of steps:

       Ionic tail switched over r =  2.7030 to  5.4060 a.u., converged in  5 steps

   and, after the ghost tests, how closely the local potential is
   -zion/r (subroutine coulomb_tail):

       Coulomb tail of the local potential
         alpha_z = 4 pi int r**2 (V_loc + zion/r) dr  2.2381270191E+01 Ha bohr**3
         alpha_z change, table cut at r =  6.00:  3.9343E-07
         alpha_z change, table cut at r = 10.00:  1.7178E-07
         |V_loc + zion/r| < 1e-10 beyond r =    6.8728
         max |V_loc - V_ion| beyond rcmax  8.1246E-06

   alpha_z is the G=0 term of the energy over the whole log mesh (QE's
   alpha = int r (r V_loc + 2 zion) dr in Ry is alpha_z/(2 pi)).

6) The routine *run_diag* is next called to test the semi-local and separable
   pseudopotentials.  For bound states, the normalized bound-state
   eigenvalues and pseudo wave functions are calculated.  The eigenvalues
   and amplitudes and first derivatives at `r_c` are compared to the all-
   electron results (as ratios).  For positive-energy states, the all-
   electron log derivative is computed at *r_c*, and the "eigenvalue" of
   the pseudo wave function which matches this is calculated.  Value and
   slope ratios (equal by construction in this case) are printed.
   For Vanderbilt two-projector calculations, these calculations are
   repeated at the energies used to construct the second projector.
   "Extra" valence states are treated as bound states here.

7) The routine *run_config* performs all-electron and VKB pseudopotential
   self-consistent atom calculations on several configurations.
   Configuration 0 is always the reference configuration, and minimal
   errors here are primarily a consistency check.  Up to 4 additional
   valence configurations may be specified in the input data.
   Eigenvalues are compared to all-electron results, and energy differences
   relative to the reference configuration indicate the excitation energy
   error of the pseudopotential.  The VKB calculations are performed with
   the specified number of projectors and non-linear core correction.
   This routine duplicates a capability of the OPIUM code, but its utility
   is not clear (see paragraph 3 of Sec. V of the paper).

8) The routine *run_plot* generates output that will be automatically
   plotted by gnuplot when running the main shell script.  These include the
   unscreened semi-local potentials, with the local potential added in
   the more general case of a polynomial interpolation.  The pseudo-valence
   charge, core charge (up to a cutoff), and model core charge if calculated
   are plotted next.  Then, all-electron and pseudo wave functions are
   compared, looping from l=0 to l = lmax.  In the shallow core with occupied
   valence case, both wave functions are computed, the single-node pseudo-
   valence wave function being calculated with the specified number of VKB
   projectors.  For two projectors, the first- and second-projector pseudo
   wave functions are printed for comparison.  The cutoff energy convergence
   predictions based on the residual kinetic energy are output in form
   suitable for plotting as of version 2.0.  The *run_plot* output is also
   printed on the standard output below the heading "DATA FOR PLOTTING."

9) The routine *run_phsft* compares the log derivatives of the all-electron
   and pseudo wave functions (bound or scattering) at the radius rxpsh, or
   by default at rc(l) times amesh**2 (rc of the local potential for l
   beyond lmax), the second mesh point beyond rc when rc is one; the
   solutions are interpolated there, so the radius need not be a mesh
   point.  Since simple plots of log derivatives contain multiple divergences, and
   are difficult to interpret visually, I plot atan(r * ((d psi(r)/dr)/psi(r)))
   at this radius, which is akin but not equal to a scattering phase shift.
   These are plotted using the specified number of vkb projectors.  If the
   local potential is that of a single l, no projectors are used for this l.
   These results are printed to the standard output, and plotted after
   the wave functions plots for the corresponding l.  A smoothed downward step
   ~pi in these plots at negative energies indicates a bound state.  The energy
   range of these plots is specified in the input data, and starting the
   range at an energy somewhat below the lowest relevant core state energy
   serves as a ghost test.  A step in the pseudopotential plot which is not
   matched by one in the all-electron plot indicates a ghost.  The
   pseudopotential results are shifted by an appropriate multiple of
   pi to align the plots at the first reference energy, where the two
   solutions agree beyond rc.  The phase is that of the point
   (psi, r d psi/dr), which is defined only up to pi, as a solution's sign
   is arbitrary; *phsft_m* follows twice it, which is sign-invariant, from
   each energy to the next, halving the step until no step turns it by more
   than pi/2.  The curves are then continuous functions of the energy, the
   radius and the potentials, with no corrections for jumps (see
   ghosts.md).

10) The routine gnuscript produces an input script for gnuplot, which is
    extracted and appropriately handled for plotting by the main shell
    script scripts/run.sh.

    The routine *diagout* then writes the block "Begin DIAGNOSTICS" to
    "End DIAGNOSTICS", for programs that read the output: the cutoffs at
    which each projector's convergence curve reaches every error level
    (1e-2 Ha down by half decades to 1e-7, past the tables' 1e-5; a level
    the curve does not reach on its q grid has no line; the curve, here
    and in the tables and plots, is the residual to infinity, the q
    integral to q_infinity plus its tail, while the optimization
    minimizes the residual to q_infinity), the
    curve's tail (`ekin_tail l ikap iproj t`, the residual kinetic
    energy from q_infinity to three times that, up to 1.4e-7 Ha), the
    r-space kinetic energy of each wave less the Fourier one to three
    times q_infinity (`ekin_check l ikap iproj t`: by Parseval the
    energy above that, below 1e-10 Ha, plus the errors of the q grid's
    quadrature, up to 8e-9 Ha, and of the r-space integral, up to 5e-9
    Ha at the default mesh), and the log
    derivatives of *run_phsft*
    with their radius, and the margins of the ghost tests, to 17
    significant digits, one tagged value per line (`ecut l ikap iproj err
    ecut`, `ld_r l ikap r`, `ld l ikap e ae ps`; ikap is 1 for j=l+1/2 or a
    scalar-relativistic l, 2 for j=l-1/2). The ghost margins are, per row
    of the bound-state test, the energy of the nonlocal Schroedinger
    equation and the basis eigenvalue matched to it (`ghost l ikap n e_ref
    e_basis`; a basis eigenvalue below min(e_ref - 1e-4, e_ref (1 + 1e-4))
    is a GHOST(-)), and per channel the basis state of least <r>/rc more
    than 0.05 Ha above the reference (`ghost_r l ikap e r`, a GHOST(+)
    below 1). States from 0.025 Ha above the reference enter with their
    <r>/rc divided by a weight rising smoothly from 0 there to 1 at
    0.05 Ha, in the margin and in the GHOST(+) test alike, so that the
    least does not jump as a state crosses 0.05 Ha. The ghost tests'
    basis is the states of the local potential (with its switch to the
    ionic tail beyond rcmax, as in every other solver) inside a hard
    wall at exactly 3 rc, each coupled to the projectors with
    a weight that falls smoothly (C-infinity) from 1 at 1.25 times the
    kinetic energy of the largest sbf of the optimization to 0 at 1.5
    times that, so these values are smooth functions of the input.
    With budget_cutoff = "weighted" (icmod 6, 7 with fcfact <= 0), the
    block also gives the model core's budget, the cutoff of the
    occupation-weighted mean of the convergence curves (`ecut_w err
    ecut`), and the weights (`ecut_w_weight l ikap iproj w`, see
    doc/core_correction.md). It also gives the conditioning of the
    projectors per channel: the eigenvalues of the symmetrized B matrix
    (`vkb_b l ikap i b`), their condition number (`vkb_cond l ikap c`)
    and the orthonormal projectors' coefficients (`vkb_e l ikap i e`, the
    eigenvalues of the nonlocal kernel, unlike the projectors themselves
    independent of the basis); the test configurations' energy
    differences from the reference configuration, all-electron, pseudo
    and their difference (`test_config icfg de_ae de_ps err`, or
    `test_config_failed icfg`), and their eigenvalues (`test_eig icfg n
    l ikap f eae eps`); and, with a model core, its extent: the AE and
    model core charges at rcmax (`core_rc r ae model`), the radii beyond
    rcmax where each falls to 1e-2 and 1e-4 of that (`core_extent frac
    r_ae r_model`), and the model's largest and smallest value on (0,
    rcmax] with their radii (`core_model_max r rho`, `core_model_min r
    rho`; a stationary point refined by Brent's method on the
    derivative of its piece, so a hollow core's maximum moves smoothly).
    Last, the Fourier convergence of the local potential (*vloc_cutoff*):
    the wave cutoffs at which a plane-wave code's local energy of the
    reference configuration's pseudo valence density n, summed over the
    density's sphere |G| < 2 sqrt(2 Ecut), is converged to each level
    from 1e-2 Ha per electron down by half decades to 1e-5 (`ecut_vloc
    err ecut`). The truncation is bounded by (Cauchy-Schwarz)
    [int_Q q^2 n(q)^2 dq int_Q q^2 V(q)^2 dq]^(1/2) / (2 pi^2 N), a
    smooth function of the input (the tail of |n(q) V(q)|, up to 35 %
    lower, has kinks at the zeros of n V); a spike of V_loc near the
    origin shows in it. The transforms are taken on a uniform r grid to
    q = 60 bohr^-1, and `vloc_check e` gives their local energy per
    electron less the r-space one (3e-14 to 1.5e-9 Ha in the references).

11) The routine linout creates a linear radial mesh, whose size and spacing
    are specified in the input data.  All relevant potentials and projectors
    are interpolated onto this mesh, which is printed to the output file in
    a format that can be read by ABINIT as its psp format 8.  The linear
    mesh has advantages of economy for pseudopotentials at small radii, and
    of greater accuracy for the Fourier transforms that will by performed by
    ABINIT at larger radii.  In addition, format 8 accepts numerical
    projectors and their coefficients, which allows the multi-projector
    flexibility that ONCVPSP achieves to be used.  The script extract.sh
    transfers this part of the single ONCVPSP output file to an ABINIT
    pseudopotential file in a specified directory.  The routine upfout
    uses the same linear mesh but produces a PWSCF-conpatible pseudopotential
    file in the UPF format.  Which one is called depends on the "psfile"
    input variable; psfile = none calls neither, so no pseudopotential
    file is written.  While all the psp files presently linked to PWSCF use
    a log mesh, the linear mesh works fine and produces results in excellent
    agreement with ABINIT results.  With upf_mesh > 0 (the optional third
    value of the output grid line), upfout writes the UPF file on the log
    mesh instead: every istr-th point of it out to the radius of the linear
    mesh, istr the smallest stride leaving at most upf_mesh points, with the
    tables as generated at those points rather than interpolated (as
    PseudoPotentialGenerators.jl does). QE stops on a mesh of more than its
    ndmx = 3500 points, so 3500 is the natural value. The default, 0, keeps
    the linear mesh, the one the reference tests compare.


The auxiliary program *oncvpspnr.x* runs oncvpsp_main with irel=0, which
does non-relativistic all-electron calculations.  The only
reason for including this feature is to support the discussion following
Eq.(24) of the paper. The script *scripts/run_nr.sh* calls it.

The subroutine *run_optimize* executes the math in Sec. II of the paper for
each l, and has the following steps:

1) The routine *wf_rc_der* calculates the value and 4 derivatives of the
   input all-electron wave function at *r_c*.

2) The routine qroots selects the wave vectors *q_i* for the xi^B basis,
   Eq.(3) according to the principles discussed following Eq.(10).

3) The routine *sbf_basis_con* creates the xi^O basis using Eq.(4).

4) The routine *Const_basis* executes Eqs.(6-10), producing the xi^N basis
   and the *phi_0* component of the pseudo wave function.  Diagnostics
   confirming that the constraints are satisfied are printed.

5) The routine eresid calculates the matrix elements in Eq.(11) of the
   E^r "operator" defined in Eqs.(1-2).  The q integration is performed
   inward from a large value ("infinity"), saving snapshots of the matrix
   elements at regular intervals as well as the values at the specific
   specified *q_c* cutoff.  These snapshots allow the convergence profile
   of the optimized pseudo wave function to be calculated very economically.
   This is the single most time-consuming step in the code, and setting
   the integration increments dq and dr to larger values will save time
   at some expense in accuracy.  These are set in the first executable
   statements of *run_optimize*.  With the new sbf8 routine of release 2.1.1
   you shouldn't need to  bother.

6) The routine optimize diagonalizes the E^r_ij matrix, calculates the
   xi_R basis functions, and the coefficients in Eq.(14).  It then
   executes the minimization procedure discussed in the paragraph
   surrounding Eqs.(15-16), giving the optimum x_i in Eq.(13).
   It finds the coefficients of the optimized pseudo wave function
   in the xi^N basis and uses the "snapshot" values of the terms in
   Eq.(11) to compute the convergence profile.  Finally, it finds
   the coefficients in the xi^O and xi^B basis sets.

7) The routine pspot calculates the semi-local potential and the
   (epsilon_i - T)|ph_i> component of the |chi_i> projector, Eq.(18).
   It also calculates the total kinetic energy from the Laplacian to
   be compared with the Fourier-transform result (the q_c = 0 end
   of the "snapshot" convergence profile).

At this point, run_optimize exits if only one projector is wanted.
If a second projector is wanted, it repeats step 2 above, computing
derivatives of the second all-electron wave function.  The existing
xi^B basis appears to work well for the second projector, so step 3
(choosing a set of q_i) is skipped.  The constraint vector and matrix,
Eqs.(6-7) is then incremented following the discussion surrounding
Eq.(25), and steps 4-7 are completed to provide an optimized second
projector.

*run_optimize* prints some additional diagnostic comparisons following
several of the subroutine calls, and prints the plane-wave energy
cutoffs corresponding to kinetic-energy-per-electron errors of 10^-2 to 10^-5 Ha.

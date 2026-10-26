## Changes from ONCVPSP 4.0.1

This page lists how the `smooooth` branch differs from the released
ONCVPSP 4.0.1. It is kept up to date: a change that moves an output, or
changes what an input means, adds or edits an entry here in the same
commit. Each entry says what 4.0.1 did, what this code does, and the
commits that did it, whose messages hold the measurements.
`doc/smoothness_plan.md` cites the same changes by item (R-, Q-, M-, O-,
S-, X-) and records the decisions (D-1 to D-7). `doc/numerics_audit.md`
cites them as C-1 and so on.

The changes fall into two families:

- **Methods** change what is computed: how the physics is approximated,
  which assumptions are made, and how the pseudopotential is built. A
  methods change gives a different pseudopotential even when it is
  computed exactly.
- **Numerics** change how it is computed: the solvers, quadratures,
  interpolants, root finders, iterations and linear algebra, and the
  tools around them. A numerics change moves the outputs only by the
  error of the old algorithm, and the outputs converge to the same answer
  as the mesh is refined and the tolerances are tightened.

Some entries could go in either family. The rule used here is whether a
perfect computer would still give a different answer. Using `rc` as
given rather than snapping it to the mesh is numerics. Replacing the
confining well of a well state with the free solution is a method.

Where an entry has an input that selects it, the 4.0.1 behaviour is the
default unless the entry says otherwise. Inputs are written as their
TOML keys. The text input has the same values as optional fields (see
`doc/32_Ge_annotated.dat` and `doc/users_guide.md`).

### Methods

#### Smoothness in the continuous inputs

The overarching goal of the branch: every output is a smooth function of
every continuous input (core radii, `qcut`, projector energies, `debl`,
`rcloc`, the model-core parameters, `rxpsh`). The aim is an output that
an optimizer can differentiate by finite differences or follow along a
path. In 4.0.1 many outputs were staircases or had kinks in the inputs.
The causes were radii snapped to mesh points, bisections stopped at loose
tolerances, discrete searches, and stencils across kinks. Most of the
numerics entries below serve this goal. They are tested by the
`smooth`-labelled scans (`scripts/smooth_scan.py`), which fit each output
over a narrow window of one input and require the residual to sit at
round-off (about 1e-10 relative or below).

What stays discrete on purpose is listed in `doc/smoothness_plan.md`
under "Truly discrete (kept)". The list includes the principal quantum
number of a state, whether a bound state or a ghost exists, the
model-core fill type, and the basin a search lands in. Where such an
event happens, the margin to it is printed (the DIAGNOSTICS block)
rather than hidden.

#### The tail of the local potential: the frozen-core ionic potential

- **4.0.1.** Beyond the largest core radius, the screened potential is
  the all-electron one, screened by the all-electron valence charge. It
  was unscreened with the pseudo valence charge, which is about 1e-4
  short of the all-electron charge there in scalar- and fully
  relativistic runs. So V_loc kept the Hartree-xc response to that
  deficit. Further out, the unscreened potentials were set to -zion/r at
  the mesh points where the density had underflowed to 0, an index that
  moved with the inputs.
- **Now.** Beyond rcmax the screened potentials are switched to
  V_ion + V_Hxc[rho + rhomod], with the frozen-core ionic potential
  V_ion = -Z/r + V_H[n_c] + V_xc[n_c + n_v] - V_xc[n_v + rhomod].
  The switch is a C-infinity step from rcmax to where V_ion reaches
  -zion/r to 1e-10. The pseudo valence states are iterated with it to a
  fixed point (`src/ionic_tail_m.f90`). The method follows
  PseudoPotentialGenerators.jl. Over the outer half of the log mesh
  (fixed radii far beyond any rc) the potentials are switched to
  -zion/r. The reference configuration's pseudo-atom then reproduces
  the all-electron one to round-off again (73_Ta: 8e-5 -> 5e-9 Ha).
  A "Coulomb tail of the local potential" block reports alpha_z and the
  radius where V_loc reaches -zion/r.
- Commits: 5e706ad (replaces the METAPSP smootherstep of 8b5cddd),
  25de811, a4f44c5.

#### Model core charges

- **icmod 6, budgeted core (new).** The model core is the all-electron
  core charge wherever the plane-wave cutoff can afford it. Inside the
  smallest radius r_s whose model needs no more than the cutoff budget
  for an xc-energy error below `rcfact`, it is a fill. The fill is the
  "minimum tail" polynomial sum_k c_k (r/r_s)^2k, matched C2 to the core
  charge, whose free coefficients minimize the density's high-q weight
  above the budget's q. Where that fill is not monotone, a Martins form
  exp(a + b r^2 + c r^4) is used. The budget is `fcfact` (Ha) when
  positive. When `fcfact` <= 0 it is taken from the pseudo-waves'
  convergence: by default the largest Ecut of the tables at `rcfact` per
  electron over the ground-state waves (those of the occupied states and
  each empty channel's first; over every wave with `budget_cutoff =
  "max_all"`), or, with `budget_cutoff = "weighted"`, the
  occupation-weighted mean (smooth where two waves' crossings pass each
  other, which the maximum is not). The mean runs to the end of the
  longest curve, each shorter one held at its last value: the s curves
  end near 200 Ha, below where a hard p wave at a small rc converges.
  An empty channel's first wave weighs
  2(2l+1), or its `ecut_weight` ([vkb_projectors]; 0 leaves it out of
  the maximum too). With `fcfact` < 0 the budget is -fcfact times that
  cutoff. `slack` (0 < slack <= 1) makes the model fit slack times the
  budget. The method is ported from the BudgetedCore of
  PseudoPotentialGenerators.jl. Commits: 8fef166, dc1883b, f098dc9,
  9b90162, 59d0ba8, d0a99ca, 9ed4066 (the ground-state waves), 270c1e6
  (ecut_weight).
- **icmod 7, hardness-matched core (new).** Of the radii that fit the
  budget, the one where the model core best reproduces the all-electron
  d2Exc/dn_i dn_j of the bound pseudo states. This is the
  HardnessMatchedCore of PseudoPotentialGenerators.jl. Commit: 9b90162.
- **Teter cores (icmod 2-4).** The Teter function is now evaluated
  analytically everywhere, and is tapered to 0 over x = 3 to 4 with a
  C-infinity switch. The bare function's x^-10 tail is much wider than a
  core charge, and it moved the end of the ionic switch from 6 to 23
  bohr in 79_Au_lxc. Commit: 3b812b9 (M-5).
- **icmod 4 search grid** is an input (`fcfact`/`rcfact` ranges and
  steps). The defaults are 4.0.1's hard-coded values. Commit: feb98a3.

#### Well states (unbound projector references)

- **Onset of the well.** In 4.0.1 the confining well rose from
  rc*amesh^5, so the well states, and through them the kernels, changed
  at first order in the mesh spacing. The onset is now rc*`well_onset`
  (default 1.006^5, where it sat on the default mesh). The well is
  multiplied by exp(-0.01/x), which makes it C-infinity at its onset and
  keeps the confined states those of before. Commits: e3a587f (R-12),
  5a16005.
- **`ep_as_given` (D-1, off by default).** 4.0.1 overrides the input
  `ep` of a well state: a bound state whose outermost peak is inside
  0.75 rc takes its own energy, one inside 1.5 rc forces ep = 0.05, and
  a negative ep becomes 0.25. With `ep_as_given = true`, ep is always
  used. A bound state at or above it is an input error, and one below it
  is raised to it by a wall as close as that takes. Commit: e3a587f.
- **`well_tail = "free"` (D-7, off by default; the optimizer sets it).**
  The kernels of a well state depended on its tail beyond rc, through
  the residual kinetic energy, and so on the arbitrary well: its onset,
  height and width, and the node count it confines. A 7 % jump in 14_Si
  l=1 came from a node entering. With the free tail the well state is
  the atom's own solution at exactly ep, normalized to 1 inside rc.
  Beyond rc it is cut off by an erfc window whose width comes from qcut
  and the state's own oscillation or decay. The window is smooth across
  ep = 0, so a growing tail for ep < 0 is allowed. Nothing then depends
  on a well (outputs identical for `well_onset` 1 to 2.5). The fully
  relativistic code keeps the j = l+1/2 well's eigenvalue for j = l-1/2.
  Commits: 7323a34, d603e40. The window's bound holds while kappa L
  <= 12 (kappa the tail's growth rate, L the window's width): beyond
  that the windowed tail of an ep < 0 state peaks at exp(6 kappa L +
  (kappa L)**2/4), and a low qcut made it overflow silently (E_resid
  6e171 Ha for an I d well state, a NaN in the optimization for N s).
  Such an input is now an error, like qcut not above the state's wave
  vector, and the well's share is undefined there. Commit: e781800.
  That bound is not enough: below it, E_resid is computed from the
  transform of a windowed tail whose peak grows as exp(6 kappa L), and
  its round-off grows with the peak squared, from kappa L ~ 1.5
  (oncvpsp-opt rows: I d at -1.61 Ha, kappa L 3.6, E_resid 2.5e3 Ha; N s
  at -0.40, 6.1, 4.5e9; Cs p at -0.49, 7.8, 5e20; about 1 % of its
  completed rows). So the free tail is now for ep >= 0 only (D-7, the
  user's decision): a well state below 0 with `well_tail = "free"` is an
  input error, which names the state and suggests a higher ep or
  `well_tail = "confined"` or `"well"` (the input's ep; the j = l-1/2
  state of the fully relativistic code keeps the well's eigenvalue,
  D-7c). The test inputs' free-tail states are all above 0, so nothing
  printed changes; the smoothness row of 14_Si's ep[2] with the free
  tail now scans 0 to 0.04 Ha instead of across 0.
- **`well_tail = "confined"` (D-8, off by default).** 4.0.1's well rises
  to ep + 0.5 Ha. As ep rises its wall steepens until no state of the
  asked n exists, and the state moves to n + 1, a jump in its tail and
  the kernels (14_Si l=0 at debl 0.84: 1.5e-2). The confined wall,
  x^3 exp(-0.01/x) Ha, rises without bound, and the state keeps its n,
  so every ep above the atom's own state of that n has one well and
  nothing switches. The n of each l's first well state is the optional
  `well_n` (default: the n the well states start from); a higher n gives
  a gentler wall, whose tail adds less to E_resid (the well's-share line
  says how much). The states are normalized inside rc, as with the free
  tail. Commit: 4286074.
- **`well_norm` (off by default).** 4.0.1 normalized every well state
  over all space, which for a confined state is the well's choice: a
  wider well spreads the state and lowers its E_resid. `well_norm`
  ("space" or "rc") sets the normalization of the states the well
  confines; by default it is that of `well_tail` (all space for "well",
  inside rc for "free" and "confined"), and "space" with "free" is an
  error, as the free tail has no finite norm. The kernels do not depend
  on it (6e-13 in 14_Si); the residual energies, convergence tables and
  cutoffs of the well states do (49 % at 1e-5 in 14_Si's d), so cutoffs
  of different `well_tail` modes compare at the same `well_norm`.

#### Pseudo-wave optimization

- **The residual kinetic energy at the exact qcut.** 4.0.1 evaluated
  E_resid at the q mesh point nearest qcut, rounding qcut by up to
  0.008 bohr^-1. It is now at qcut itself (see Numerics for how). This
  changes the kernels by about 2e-3. Commit: 6c310b9 (Q-1).
- **`wave_regularization` w (default 0).** Where qcut*rc is large,
  E_resid(qcut) is flat to round-off in the null space. The sign of the
  null-space coefficient, which decides between a nodeless wave and one
  with a node, was then set by the quadrature's rounding. With w > 0 the
  optimization minimizes E_resid(qcut) + w E_resid(qcut/4). With w = 0 a
  flat E_resid is reported with a warning. Commit: 9bb82b7.

#### Relativistic correction of the projectors

- **The taper (D-4).** 4.0.1 subtracted f v_rel phi from chi inside the
  window 0.95 rc to rc, with v_rel = Delta/2u formed by dividing by the
  all-electron u. The code needed a threshold and a neighbour average at
  nodes of u. It now subtracts f Delta[u_AE]/2, where
  Delta = u'' - [l(l+1)/r^2 + 2(V - eps)]u, from the first-order system.
  It never divides by u and, in the fully relativistic code, holds the
  spin-orbit term too. The two differ by O((0.05 rc)^ncon). The method
  follows PseudoPotentialGenerators.jl. Commit: e3a587f (R-11).
- **Semi-local local potential (lloc < 4).** The local potential took
  the semi-local potential of channel lloc untapered, which in
  relativistic runs jumps by vr(rc) at rc. The other channels' projectors
  inherited the jump. It now gets the same taper as the projectors.
  Commit: d1393c5.

#### Exchange-correlation and relativity inputs

- `iexc = 6`: Slater exchange with VWN5 correlation (Libxc -001007).
- `iexc = 7`: VWN with the MacDonald-Vosko relativistic exchange (Libxc
  -532007, `lda_x_rel`). These are the functionals of the NIST atomic
  reference data. `iexc = 5` is refused: it is r2SCAN in METAPSP.
- `c_light`: the speed of light is an input (default 137.036; NIST uses
  137.0359895).
- The non-relativistic solvers are truly non-relativistic. In 4.0.1 they
  used fss = 1e-12, a speed of light of 1e6 (4e-5 Ha in U), and the
  scalar-relativistic starting exponent for l > 0.
- Commits: 34a4981, 522b58a, 88db006.

#### Diagnostics: what is measured

These additions change what a run reports, not the pseudopotential. All
are in the `Begin DIAGNOSTICS ... End DIAGNOSTICS` block at full
precision.

- **Convergence curves to infinity.** The residual kinetic energy above
  q_infinity, the wave's own E_resid from q_infinity to 3 q_infinity, is
  added to the convergence curve. The Ecut tables, plots, model-core
  budget and every level use the corrected curve. Levels run down to
  1e-7 Ha, with `ekin_tail` and `ekin_check` (the r-space kinetic energy
  less the Fourier one). The optimization itself keeps E_resid to
  q_infinity. Commits: 9889e33, 0a7971f, 1a6da96, 2e4ab37.
- **Local-potential convergence (`ecut_vloc`).** The cutoff at which a
  plane-wave code's local energy of the reference valence density is
  converged, from a Cauchy-Schwarz bound. A spike in V_loc near the
  origin is invisible to the projectors' cutoffs. Commit: e7655cd.
- **Projector conditioning** (`vkb_b`, `vkb_cond`, `vkb_e`), test
  configurations, and the model core's extent. Commit: 9889e33.
- **The ghost test, continuous.** The hard wall sits at 3 rc exactly,
  not at the mesh point below it. Barrier states enter the basis with a
  C-infinity weight in energy, not at a sharp cutoff. The basis is built
  on the same switched V_loc as the reference energies. The margins are
  printed (`ghost`, `ghost_r`). Commits: 99e97fc, 123954b.
- **Log derivatives** are at rxpsh exactly. The default radius is
  rc*amesh^2, which was 4.0.1's mesh node whenever rc was one. The phase
  is followed by continuation in its doubled angle, which removes 4.0.1's
  pi and 2 pi corrections, and is aligned at ep itself. `oncvpsp_ld.x`
  gives them at any radii and energies from a saved run. Commits:
  ab51d65, e6b38ed, e0c94b0.
- **The confined atom (new), `oncvpsp_confine.x`.** From the
  all-electron file and the `--ps-file` of a run, the band edges of a
  Wigner-Seitz sphere of radius rs, all-electron against pseudopotential
  at frozen potentials. The bottom is where r phi'/phi = 0 at rs, the
  centre where it is -l-1, and the top where phi = 0. Each edge is the root
  of a monotonic phase (atan(r phi'/phi) less pi per node inside rs,
  `phsft_node_phase`), independent of the scan step. The run also gives the
  cutoff each pseudo state needs in the sphere, from its variational energy
  in the spherical Bessel functions of the same boundary condition
  (`ws_basis_m`). This covers every channel at the energies of the valence
  bands and replaces the band cutoffs removed earlier, which measured an
  invented interstitial tail. Si d at rs = 1.25 rc: 1e-5 Ha at 8 to 28 Ha
  from the bottom to the top of the 3d band (15.6 Ha from the old tail).
  The program also compresses the atom: the self-consistent all-electron
  atom and pseudo-atom of the reference configuration in a smooth wall,
  1/2 k (r - r0)^2 exp(-0.5/(r - r0)) beyond r0, the largest core radius
  by default, at a few strengths k. It gives their total-energy changes
  from the free atom and their valence eigenvalues. For 14_Si at k = 0.1
  Ha/bohr^2, the AE change is 55.280 mHa and the PS change 55.248 mHa.
  `aeatom` and `psatom` take the wall as an argument, which is zero in
  the run, and the pseudopotential file (format version 2) carries the
  unscreened local potential, the valence and model-core densities and
  iexc. Nothing in the run's output changes. An edge's partner is the
  pseudopotential edge of the same phase, and the output names what does
  not pair: the band of a core state (`ws_core`, such as Nd 4f), a
  pseudopotential edge with no all-electron partner (`ws_extra`, a ghost
  in the window, such as the localized s ghosts of I and Pb near 4.4 Ha
  that the run's GHOST(+) test reports),
  and an edge of another phase nearer than the partner
  (`ws_paired_off`). A last `confine_summary` line gives the largest
  |e_ps - e_ae| and cutoffs at 1e-4 and 1e-5 Ha of the edges below 0,
  the counts of those lines and of failed walls, the largest
  |de_ps - de_ae| of the walls, and the cutoffs of the top edges alone,
  free of the boundary's tail (in the optimizer's Li fronts they still
  need 45-51 Ry where the tables give 11-13, from the empty p channel the
  budget leaves out). The pseudopotential's phase counts the
  poles of its r phi'/phi at rs rather than its nodes, as its wave
  function can gain a single node inside the projectors (Pb p, I s and Nd
  f in the optimizer's runs, where it paired every edge above it one
  band off). The poles are the eigenvalues of the sphere with phi(rs) = 0
  (`ws_eigenvalues`), each refined to the root of u(rs). The two phases
  are aligned by the multiple of pi they differ by about the first
  projector's reference energy, where the pseudopotential matches the
  all-electron log derivative by construction (on the whole scan only
  for a channel without projectors). The most frequent multiple of the
  whole scan had paired every band with the next one down in a window
  wide enough to take in a core band and the energies of several Ha
  where the pseudopotential no longer follows the atom (the PseudoDojo
  inputs' -12 Ha: max |e_ps - e_ae| 7.3 Ha for Cs, 3.6e-4 Ha now).
  Each wall's atoms start from the wall before. The free all-electron
  atom starts from the data file's potential (`aeatom`'s `warm`), and
  the pseudo-atom ramps its wall in over 10 iterations rather than 50.
  Each edge of a sphere shares the quadrature of its channel
  (`ws_grid_init`), and the scan of a channel takes the phases of every
  radius from one solution at each energy (`phsft_node_phases`; a root
  still integrates only to its own radius). On the test inputs this cuts
  the instructions by 45-52 % (Bi 26.6 to 13.2 G), and the output
  changes at round-off: the walls' energies by up to 1.3e-10 Ha, the
  edges by 2e-14 Ha, and the cutoffs by up to 1e-7 relative. Commits: 7041579 (the sphere), 3eb250c (the
  walls), 5c91c31 (the unpaired edges and the summary).
- **Well-state diagnostics** (`run_diag`) take the state's index from
  the pseudo-atom at epa rather than the all-electron node count, which
  can differ from it (see the eigenvalue index under Numerics). Commit:
  d159196.
- **The well's share of the cutoffs (new).** 4.0.1 did not say how much
  of a well state's residual kinetic energy came from its arbitrary
  well. For each well state confined by its well, the output now gives
  E_resid at qcut of the optimized wave with the well's tail and with the
  free tail at the same energy (`well_share`), and a WARNING where they
  differ by more than 10 % (84 of the 358 well states of the test
  inputs, 74_W's up to 94 %) or where qcut is not above the state's wave
  vector. Nothing else changes. Commit: 74b1c54.

### Numerics

#### Radii, mesh and discretization

- **Core radii as given (R-1 to R-3).** 4.0.1 truncated rc to 6
  decimals and moved it to the next mesh point, so every output was a
  staircase in rc. Radii are now used as given. Everything built inside
  rc works at rc itself: projectors as pieces split at 0.95 rc, rc(5)
  and rc, the polynomial local potential matched at rc(5), integrals per
  segment, and the well and ionic windows at continuous radii. On inputs
  whose rc was already on 4.0.1's rounding the outputs are unchanged.
  Commits: 54c5267, e3a587f (R-6 to R-13).
- **Log mesh as input.** `amesh`, `rr1` and `rmax` are optional inputs.
  The defaults are 4.0.1's mesh. Commit: 8490a59.
- **Matching data at rc from the radial equation (R-4).** The value and
  four derivatives of u/r at rc came from 7-point stencils, which were
  off by up to 1e-4 in the 4th derivative of well states. They now come
  from a Taylor recursion of the first-order system the solver
  integrates (`ode_rc_der`). Commits: b47bc4d, 79b1cd8, 693c78e.
- **Log-mesh operations (`logmesh_m`).** These are values, derivatives,
  integrals and roots at any radius. The interpolant is Floater-Hormann
  (d = 6) in ln r, windowed by a C-infinity bump. Derivatives come from
  11-point stencil fields, cumulative integrals from 8th-order cells with
  compensated summation, and roots from Brent in ln r. `ugrid_m` gives
  the same on uniform grids. Commits: a1b2a06, b47bc4d, 6c310b9. The
  bump's smooth step takes one exponential instead of two, and each
  node's weight one division instead of two, for 2 to 3% of a run; this
  moves outputs by up to 1e-9 relative where they amplify round-off
  (c60094f).
- **Piecewise functions (`pieces_m`).** A function with kinks is held as
  pieces, each a smooth extension over the whole mesh, so no stencil or
  interpolant spans a kink. Commit: e3a587f.

#### Radial solvers

- **Scalar-relativistic flux form.** `lschfb` and `lschfs` integrate the
  Koelling-Harmon equation in its first-order flux form, which needs no
  dV/dr. Finite differences of a GGA potential near the nucleus made
  4.0.1's SCF diverge on fine meshes. The outward start is the series at
  the origin with the KH exponent, carried out by RK4 in ln r over a
  range independent of the energy. The method follows
  RadialAtomicSolvers.jl. Commits: 77c1355, 4352a49, 37d4f5c.
- **Dirac outward start.** The outward start is now the series of the
  Dirac equation at the origin, with PECECE steps and an inward
  re-integration of the first points. 4.0.1 started from the leading
  term only, with starting derivatives of the wrong sign of the energy
  term. Commits: 0b1c655, 34befb5 (C-4).
- **Inward start at practical infinity.** 4.0.1 started the inward
  integration at 10 times the classical turning point and set u = 0
  beyond it. For shallow localized states (Tm 4f, Cu 3d) that truncated
  the norm at about exp(-5). The start is now where the WKB decay
  reaches machine precision. Commits: 9cc675a, 60d4ad3.
- **Well states past the mesh's end (C-19).** 4.0.1 solved a well state
  on the log mesh to rmax, so a wide well (a small `ep`), whose state
  reaches rmax, was a box at rmax: the state, its tail in E_resid and the
  kernels depended on rmax (14_Si's d at 0.02 Ha: kernel 2e-2 between
  rmax 45 and 60). The mesh is now continued past rmax, for that state
  only, until its decay beyond the well reaches machine precision. The
  well rises only 0.5 Ha above the state, so most well states use it
  (the mesh to 56 to 70 bohr); the kernels of the test inputs move by
  at most 3e-8, 74_W's by 1.5e-6. Commit: 0a37862.
- **Box state.** When an early SCF potential has no bound state for an
  occupied 4f (the f-collapse instability: Yb, Tb, Dy, Eu), the solver
  returns the box state with the right node count instead of failing.
  Commit: 9cc675a.
- **Continuous matching point (S-8).** Bound states are matched at the
  continuous classical turning point, with the energy correction
  averaged over a C-infinity window. 4.0.1 used the turning point's mesh
  index, so the eigenvalue jumped by about 1e-13 Ha whenever that index
  moved. Commit: 8eae4dd.
- **Converged to round-off (S-2 to S-5).** Every eigen-solver iterates
  until the correction is at a few ulps and returns the energy and
  function of the last pass. 4.0.1 stopped on a small correction and
  returned the energy from before it, so the result depended on the
  start by up to the tolerance (1e-8 to 1e-12). Commits: 60d4ad3,
  8eae4dd.
- **Pseudo-atom per segment (R-10).** The outward integration with
  projectors runs on its own grid in each segment between kinks of the
  potential or the projectors. Each segment is restarted from the one
  before and the result is interpolated back to the nodes. Projector
  overlaps end at rc instead of at the matching point. Commits: e809a29,
  dc24b82.
- **Eigenvalue index instead of node count.** With projectors, a node
  can enter from r = 0 without an eigenvalue passing. 4.0.1's count
  beyond Hamann's small-r cutoff, max(0.1 rsrc, 0.05) bohr, then jumped.
  `lschvkbb` and `lschvkbbe` now bisect on the number of eigenvalues
  below E (`outwf_index`, from the inertia of the Schur complement of
  the projectors, by Haynsworth's theorem). Commits: d159196, bd574c9.
  Where every energy tried has more eigenvalues below it than the state's
  index allows, a state lies below the bracket (a bound ghost of the
  projectors): `lschvkbb` returns ierr 200 instead of 100, and the
  pseudo-atom's checks say so, where 4.0.1 ran out of iterations
  silently.
- **Normalization by the Floater-Hormann integral.** The bound and
  scattering solvers (`lschfb`, `lschfs`, `ldiracfb`, `ldiracfs`,
  `lschvkbb`, `lschvkbbe`, `lschvkbs`, `lschpsbar`) and `renorm_r`
  normalize with `quad_m`'s `log_norm`: the power law u(1)² rr(1)/(2γ+1)
  below rr(1) and `lm_integral` above it, from the solution's points
  alone, as `fpovlp`'s overlaps. 4.0.1 summed al r u² over the nodes with
  a 3-point end correction and a power law below rr(1)/√amesh that
  assumed u = r^γ. The ghost test's <r> takes the same integral on the
  barrier's mesh. The kernels of the test inputs move by at most 3e-8,
  V_loc by 2e-9, the AE energy terms by 6e-8 Ha (their total by 1e-10),
  and the ghost test's <r>/rc by up to 1.3e-4. On 14_Si's meshes from
  amesh 1.006 to 1.00075, <r>/rc is closer to its converged value than
  before at every mesh.
- Performance: the Adams formulas inlined, all projectors' outward
  solutions in one loop, energy-independent terms hoisted out of
  `lschfb`. Commits: 13c5b30, 9e59c13, 77d9cb0, c411754.

#### Self-consistency

- **All-electron SCF.** 4.0.1 stopped when no eigenvalue moved, which
  could happen at once from a poor start (H on a fine mesh). The SCF now
  mixes the screening potential (Anderson over 8 iterations, with a
  Kerker-type filter in ln r inside 10 r* for scalar-relativistic GGA).
  It converges on the screening residual (1e-11 Ha) and the eigenvalue
  changes, accepting a residual stalled at its round-off floor with a
  warning. A failed reference SCF stops the run. The method follows
  RadialAtomicSolvers.jl. Commits: 35f1150, 77c1355, e0ac9fe, a9e131c,
  b3a7c7d, fa0b95e (C-2).
- **All-electron energy.** The total energy is now computed as RAS and
  NIST compute it, with the kinetic energy as the eigenvalue sum less
  the potential energy in the input potential. The terms are printed.
  The head of each energy integral below rr(1) follows the density's
  power law, not a constant density. Commits: 35f1150, 260d2af (C-3).
- **Pseudo-atom SCF.** It converges on the eigenvalues and residual to
  1e-12 Ha, only at the target occupations, with Anderson mixing. From
  the 50th iteration on, the occupations are exactly the targets. 4.0.1
  could stop during the occupation ramp, or never converge where the
  ramp did not round to the target. Commits: 60d4ad3, fa0b95e (C-1).
  During the ramp's 50 iterations each bound state stops at a relative
  energy correction of 1e-4, not at round-off (`lschvkbb_tol`): those
  iterations only lead to the target's start. This moves the test
  configurations' energies by about 1e-13 Ha, within the SCF's
  tolerance, and saves 4 to 13% of a run. Commit: 3b60ff4.
  Like the all-electron SCF, it accepts a residual stalled at its
  round-off floor, with the same warning (`ae_converged`). Li with an s
  ghost almost bound stalled at 1.4-3e-12 Ha from the first iteration
  and failed after 100, in the test configurations and the confined
  atom's free atom. No test input reaches that floor, so their outputs
  are unchanged.
- **Repeated all-electron SCFs reused.** The all-electron atom depends
  only on its states and occupations. Test configuration 0 repeats the
  reference at both levels, and a configuration that only removes
  charge repeats its intermediate level at its excited one. 4.0.1 solved
  each again; a run now reuses the earlier solution, to the same bits,
  with its warning of a stop on the residual's round-off floor. This
  saves 2 to 5 of the 9 SCFs of the test inputs. Commit: 410344c.
- **Ionic-tail fixed point** to 1e-13 Ha. Commit: 25de811 (S-6).
- Tested against the NIST atomic reference data (LDA, ScRLDA, RLDA,
  H to U and cations) at NIST's stated accuracy (`ctest -L nist`;
  expected failures listed in `scripts/nist_test.py`). Commits: 4f7b5b7,
  0fadb26.

#### Bessel basis and the pseudo-wave optimization

- **Residual kinetic energy (`eresid`, Q-1 to Q-4, C-10).** Inside rc,
  the transform is the exact one of the Bessel sum, from Lommel's closed
  form. The all-electron tail is on a linear mesh continuous in rc. The
  q grid is anchored at qmax with a continuous spacing. The integrals
  accumulate from qmax by an 8th-order rule with Kahan summation, and
  the value at the exact qcut comes from a Gauss-Legendre partial cell.
  qmax is a C-infinity maximum. Ecut crossings are roots of the
  interpolant of ln E, not log-linear interpolation. Commits: 4c0ea5f
  (the sin/cos rotation, 1.3-2x faster), 6c310b9, 23d54d2, 8e777ba
  (M-13). The curve's tail above qmax reuses the curve's tail mesh and
  function, which it had rebuilt at a spacing equal up to round-off. This
  moves its cutoffs by about 1e-15 relative (410344c). The tail's
  transforms on the q grid are chirp-z transforms (`czt_m`, Bluestein's
  convolution by FFT): for l <= 3, r j_l(qr) is a sum of
  q^-k r^(1-k) sin or cos(qr), so each of those l+1 sums over the tail
  mesh comes at every node at once, in O(L log L) instead of O(nr nq).
  Their round-off is 1e-16 of the transform's largest value, below the
  4e-16 of the sums of sines they replace. At the lowest nodes, where
  those terms cancel, and for l > 3, the tail is still summed directly.
  E_resid takes a fifth of the time it did; kernels move by at most
  2e-11 and cutoffs by about 1e-10 relative.
- **Basis roots (`qroots`).** The roots are found on the pole-free
  g(q) = ulgd j_l(q rc) - q j_l'(q rc) and Newton-polished. 4.0.1 tested
  the log-derivative mismatch and rejected large sign changes as poles,
  which for a large ulgd rejected every root. Commits: feb98a3, 1d5bfc0
  (Q-5), 40c219f.
- **Basis linear algebra (R-5, Q-8).** The overlap matrix is formed in
  closed form. The constraint rows are scaled before the solve. The
  basis's linear algebra (overlap, orthonormalization, constraints, null
  space and inside transforms) runs in quad precision. Hamann's basis is
  nearly dependent (overlap eigenvalues down to 5e-9), and double
  precision lost up to 5e-10 in the kernels. Commits: bc00d9b, 170bc99.
- **Spherical Bessel functions (Q-10).** j_0 to j_3 switch to the power
  series below x = 0.5 to 3.5, where the closed forms cancel. At
  x = 1e-3 the old j_3 was off by 2.6e5 times its value. `sbf8` sets
  every order and starts its recursion consistently. Commit: 0971ea2.
- **Roots polished.** The norm constraint of `optimize` and the
  well-state width (Illinois regula falsi to a few ulps) are polished to
  round-off instead of stopped at a bisection tolerance. Commit:
  1d5bfc0.
- **V_l inside rc** uses T u = q^2 u/2 for each Bessel term (Bessel's
  equation), not l(l+1)u/2r^2 - u''/2, which cancelled near the origin.
  Commit: bbae0de (C-11).
- **r-space kinetic energy** is split at rc. The origin term of 4.0.1 is
  corrected: it is ro^(2l+3)/(2l+3). Commits: 5e3d9e0, 0a7971f.

#### Model core numerics

- **Crossover and matching radii** are continuous roots of
  ln rhoc - ln(fcfact rhotps), not the last mesh point below. The core
  charge's derivatives there come from `lm_derivs`. icmod 1's constant
  is in closed form, and icmod 2's scale is found by Newton's method.
  Commit: 3b812b9 (M-1 to M-4).
- **icmod 4** refines the Nelder-Mead optimum by Newton's method along
  the Hessian's eigenvectors, so the result does not depend on where the
  simplex stopped. Commit: 8974b5c (M-6).
- **icmod 6/7.** r_s is a Brent root, and icmod 7's optimum a Brent root
  of the slope. The fill is in a Legendre basis. The error sphere is
  tapered. The quadrature below r_s runs on its own mesh ending at r_s.
  Iteration limits are errors. Commits: ca09107 (M-8 to M-12), ca1ba1e,
  86f8e58 (C-6).
- **d2Exc stencil.** 4.0.1's five-point stencil evaluated -2h twice,
  giving d2Exc 4/3 too large, second order only. Commits: 333c55e,
  a672619.

#### Output

- **Output meshes from the pieces (O-1).** 4.0.1 put the log-mesh arrays
  on the linear mesh with `dpnint`, an 8-point Lagrange interpolant
  across rc. Projectors rang beyond rc (up to 9e-4) and the kernels were
  1e-8 off smooth. Each output point is now taken from its own piece,
  and the projectors are exactly 0 beyond rc. The model core's columns
  are taken from their pieces at the joint. Commits: 9796df1, 6f87629
  (M-14).
- **Fully relativistic psp8 (O-2, D-2).** 4.0.1 dropped scalar and
  spin-orbit terms below 1e-3 of the first, so a term could drop out as
  an input moved. Every term is now written. Commit: 870835b.
- **UPF.** Every number has 17 significant digits; 4.0.1 wrote 11, and
  the linear mesh 4 decimals. `upf_mesh` writes the file on the log mesh
  with no interpolation. Exponents past 99 are written as 0. Commits:
  5469c9f, 6278f71, 8d7fc0a.
- **DIAGNOSTICS block**: tagged values at 17 digits, for scripts and the
  smoothness scans. Commit: e188eee.

#### Corrections to 4.0.1

Bugs that moved outputs, beyond those above:

- Fully relativistic all-electron valence density: the j = l+1/2 part of
  each l > 0 state was counted twice (80_Hg's integrated to 18 for
  zion = 12). This also changed the icmod 1-4 fits. Commit: b700909.
- Fully relativistic phase shifts were aligned at the input ep (often
  0.0) instead of the projector energy. Commit: e5cdf20.
- `psatom` never reported non-convergence (`it == 400` in a 100-step
  loop). `psatom_r` lost a failed state's charge. Commits: 8c7ae75,
  a672619.
- `run_config` set the intermediate occupations of one stale index.
  Commits: 3dfb268, a672619.
- `run_optimize` stored up to 9 constraints in `cons(6)`. `eresid`
  never filled qcut's matrices when qcut was within 4 dq of q_infinity.
  Commit: 3dfb268.
- `excwig` left vxc unset below a density of 1e-15. The `vkbphsft`
  alignment was unset outside the scan. The `run_vkb` division by 0 for
  an l with no projectors. Commit: 3dfb268.
- psp8, UPF and gnuplot writers: an empty l header in fully
  relativistic psp8, the "NOX+HL" label for iexc = 2 (QE dropped
  exchange), i1 formats for 10 or more projectors, and wrong plot
  columns. Commit: 2ca45c3.
- `lschpsbar`'s inward start, `wellstate`'s unconverged trial states,
  and `modcore3`'s unset start. Commit: 5e3d9e0.
- Inputs that overran arrays, `depsh` = 0, NaN and Inf reals, and inputs
  read from a pipe. Commits: f7c2e11, 17f5b9b, 0dfa6cf.

#### Tools

- **Inputs.** The TOML input (`-t`, from upstream master) is checked and
  has every key the text input has. Every echo of the input (output,
  psp8, UPF, PSML) reruns to the same output. `psfile = none` writes no
  file. Commits: e083200, e42dd70, 8722505, 6d39ae9.
- **Executables.** `oncvpsp*_ae.x` solve the all-electron atom once and
  `oncvpsp*_ps.x` reuse it (the optimizer's loop). `oncvpsp_ld.x` gives
  log derivatives from saved runs, `oncvpsp_confine.x` the confined atom,
  and `--ps-file` saves the pseudopotential data. Every run prints its commit and build. Commits:
  1ca8f73, 37bcac2, 48b6d6a, e6b38ed, 1f11258.
- **One code for nr, sr and fr.** The `_r` copies were merged into one
  driver, `oncvpsp_main(irel, istage)`. The physics is unchanged.
  Commits: 0b301d7 to 3ca0b78.
- **Build.** CMake only (the make build was dropped), with dependency
  fetching, presets and a `study` preset (-O3 -march=native
  -ffast-math). xmlf90 comes from its repository. PSML goes into the
  output instead of a file. Commits: abf71da, 394f1fc, b61415b, 739a92b,
  8b2a016.
- **Tests.** Outputs are compared with pyoncvpsp quantity by quantity,
  and pseudopotentials by their nonlocal kernels
  sum_ij |beta_i> D_ij <beta_j| and V_loc, never by the projectors,
  which are not unique. Tolerances sit at the noise of the compiler
  flags. The suite was cut to the fewest inputs that keep its coverage.
  There are smoothness scans, NIST tests and coverage tooling. See
  `doc/developer_setup.md` and `doc/test_suite_speed.md`. Commits:
  5a8f039, 10369e0, e02bcf3, 350fe06, 43be171, e2f4665.
- **Code standards.** Intents on every argument, generic intrinsics,
  documented routines, a constants module, pre-commit hooks and
  fortitude. See `doc/coding_standards.md`.

Upstream master since 4.0.1 (inherited here) adds the CMake build, Libxc
5-7, the TOML input, PSML output, PW92 in UPF and the `rxpsh` input.

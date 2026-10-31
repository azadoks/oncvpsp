!
! Copyright (c) 2026 by Austin Zadoks
!
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!> The kinds of failure that a routine's ifail reports. No library routine
!> stops: each failure writes its message and returns one of these, and
!> every caller returns the ifail it got, unchanged, up to the program,
!> which stops with exit status 1 for any of them (oncvpsp_main). A caller
!> from another language reads the kind from ifail; the message says which
!> failure of that kind it was.
!>
!> The site of a failure sets ifail=failed(kind,msg,ierr): the kind in its
!> last two digits, which message (msg_*, one per site) above them, and the
!> solver status it stops on (ierr, if any) above that, so that a binding
!> tells the failures apart without their text, and every caller still
!> only tests ifail/=0 and passes it on. fail_context splits it into a
!> failure_t, with where the run was (the stage, l, ikap) that the stages
!> record before each call. The values of msg_* are kept: a new site takes
!> the next value, and a removed one's value is not reused.
!>
!> A warning does not stop the run: its site writes its message and adds
!> a record to the warnings_t passed down to it (warn), which wrn_* names
!> as msg_* names the failures, with the same rule for their values.
!>
!> ierr, where a routine has one, is a status that its caller interprets
!> (no root in a bracket, a box state, ...), documented with the routine;
!> it is not a failure until the caller makes it one. The statuses that
!> pass between modules, those of the radial solvers and the
!> self-consistent atoms, are named here (ierr_*).

 module failure_m
 use, intrinsic :: iso_c_binding, only: c_int
 implicit none
 private
 public :: fail_none, fail_input, fail_convergence, fail_not_found, &
&          fail_construction, fail_numerical, fail_io, fail_internal, &
&          fail_name
 public :: failed, fail_kind, fail_msg, fail_ierr, fail_context, failure_t
 public :: stage_none, stage_input, stage_ae, stage_channel, stage_vkb, &
&          stage_valence, stage_model_core, stage_screening, stage_diag, &
&          stage_ghosts, stage_unscreen, stage_configs, stage_plot, &
&          stage_phsft, stage_tests, stage_output, stage_confine
 public :: ierr_box_state, ierr_iterations, ierr_not_converged, &
&          ierr_ps_scf, ierr_ghost
 public :: ierr_lxc_unknown, ierr_lxc_no_vxc, ierr_lxc_family
 public :: warnings_t, warn, warn_merge, mxwarn

!> no failure
 integer, parameter :: fail_none=0
!> the input is invalid, inconsistent or asks for something not built in
!> (check_data, the readers, a functional without Libxc)
 integer, parameter :: fail_input=1
!> an iteration reached its limit (a root search, a self-consistent loop,
!> a bound state, a Newton polish)
 integer, parameter :: fail_convergence=2
!> something the construction needs does not exist for these parameters:
!> a bound state, a root in a bracket, a crossing
 integer, parameter :: fail_not_found=3
!> the pseudopotential cannot be constructed as asked: a node in a pseudo
!> wave function, a ghost, a norm that cannot be kept, a parameter out of
!> its range for this atom
 integer, parameter :: fail_construction=4
!> a numerical operation failed: a singular or non-positive matrix, a
!> LAPACK routine's info, a value that is not finite
 integer, parameter :: fail_numerical=5
!> a file cannot be opened, read or written, or its contents do not match
 integer, parameter :: fail_io=6
!> arguments that the calling code should never pass: a bug, not the input
 integer, parameter :: fail_internal=7

! The statuses (ierr) of the radial solvers and of the self-consistent
! atoms, which their callers interpret; the values are those of 4.0.1,
! which some messages print. 0 is success.

!> no bound state: the box state stands in for it (bs_solve, lschfb,
!> ldiracfb); the self-consistent atoms accept it during their iteration
 integer, parameter :: ierr_box_state=10
!> the energy iteration of a bound-state solver did not converge (bs_solve,
!> lschfb, ldiracfb, lschvkbbe, lschpsbar)
 integer, parameter :: ierr_iterations=60
!> lschvkbb's energy iteration, or sratom's and relatom's self-consistency,
!> did not converge
 integer, parameter :: ierr_not_converged=100
!> psatom's self-consistency did not converge
 integer, parameter :: ierr_ps_scf=101
!> lschvkbb: every energy tried had more states of l below it than the
!> node count expects, so a state it does not expect (a bound ghost of the
!> projectors) lies below emin
 integer, parameter :: ierr_ghost=200

! The statuses of lxc_init and exc_libxc_check, which check_data reports

!> Libxc does not know the functional
 integer, parameter :: ierr_lxc_unknown=1
!> the functional gives no potential
 integer, parameter :: ierr_lxc_no_vxc=2
!> the functional is neither an LDA nor a GGA (meta-GGA, hybrid, ...)
 integer, parameter :: ierr_lxc_family=3

!> ifail is kind+msg_scale*(msg+msg_limit*ierr)
 integer, parameter :: msg_scale=100
!> the messages' values are below this
 integer, parameter :: msg_limit=1000
!> the largest ierr that ifail holds; failed drops a larger one
 integer, parameter :: ierr_limit=huge(1)/(msg_scale*msg_limit)-1

! Where a run was, the stages of oncvpsp_run in order (failure_t%stage)

!> no stage recorded
 integer, parameter :: stage_none=0
!> the input completed and checked (complete_input)
 integer, parameter :: stage_input=1
!> the all-electron atom and its configurations (solve_ae)
 integer, parameter :: stage_ae=2
!> the core radii and the channel l's optimized pseudo wave functions
!> (optimize_channel)
 integer, parameter :: stage_channel=3
!> the local potential and the projectors (run_vkb)
 integer, parameter :: stage_vkb=4
!> the pseudo valence states and densities (valence_densities, sr_so_r)
 integer, parameter :: stage_valence=5
!> the model core charge (modcore*, modcore_budget)
 integer, parameter :: stage_model_core=6
!> the screening of the pseudo atom (vout, ionic_screening)
 integer, parameter :: stage_screening=7
!> the diagnostics of the channel's states (run_diag)
 integer, parameter :: stage_diag=8
!> the ghost test, of ikap (run_ghosts)
 integer, parameter :: stage_ghosts=9
!> the unscreened potentials (coulomb_tail)
 integer, parameter :: stage_unscreen=10
!> the test configurations (run_config_ae, run_config_ps)
 integer, parameter :: stage_configs=11
!> the plotted wave functions (run_plot)
 integer, parameter :: stage_plot=12
!> the log derivatives (run_phsft)
 integer, parameter :: stage_phsft=13
!> the cutoff tests (vloc_cutoff, diagout)
 integer, parameter :: stage_tests=14
!> the pseudopotential data file (write_ps_file)
 integer, parameter :: stage_output=15
!> the confined atom of a run (confine_atom)
 integer, parameter :: stage_confine=16

!> a failure, for a binding: what (kind, msg, ierr) and where (stage, l,
!> ikap, iproj); 0 where it does not apply, l -1
 type, bind(c) :: failure_t
!> the kind of failure (fail_*)
   integer(c_int) :: kind
!> which message (msg_*), 0 if none
   integer(c_int) :: msg
!> the solver status the failure stopped on (ierr_*), 0 if none
   integer(c_int) :: ierr
!> the stage of the run (stage_*)
   integer(c_int) :: stage
!> the angular momentum of the channel, -1 if none
   integer(c_int) :: ll
!> the spin-orbit index (1 or 2), 0 if none
   integer(c_int) :: ikap
!> the projector, 0 if none
   integer(c_int) :: iproj
 end type failure_t

!> the warnings records kept
 integer, parameter :: mxwarn=64

!> the warnings of a calculation, in the order written (warn)
 type :: warnings_t
!> their number, all of them
   integer :: nwarn=0
!> the first mxwarn: the warning (wrn_*), the stage (stage_*), l (-1 if
!> none), ikap (0 if none), and the configuration or the state's n that
!> it is of, as its wrn_* says (0 if none)
   integer :: rec(5,mxwarn)=0
!> the stage of the next ones, which the stages set
   integer :: stage=0
!> their l where the site does not give it, which the stages set (-1 for
!> none)
   integer :: ll=-1
!> their ikap where the site does not give it (0 for none)
   integer :: ikap=0
 end type warnings_t

! The messages of the failure sites (failed's msg), each with its routine,
! file and message

!> write_head (ae_file_m.f90): write_ae_file: ERROR cannot open
 integer, parameter, public :: msg_ae_file_open_write=1
!> read_head (ae_file_m.f90): read_ae_file: ERROR file was written with a
!> different relativistic treatment
 integer, parameter, public :: msg_ae_file_irel=2
!> read_head (ae_file_m.f90): read_ae_file: ERROR atsym, z, nc, nv or iexc
!> differ from the input
 integer, parameter, public :: msg_ae_file_atom=3
!> read_head (ae_file_m.f90): read_ae_file: ERROR reference configuration
!> differs from the input
 integer, parameter, public :: msg_ae_file_reference=4
!> read_head (ae_file_m.f90): read_ae_file: ERROR radial mesh size differs
 integer, parameter, public :: msg_ae_file_mmax=5
!> read_head (ae_file_m.f90): read_ae_file: ERROR radial mesh differs from the
!> one this program builds
 integer, parameter, public :: msg_ae_file_mesh=6
!> read_head (ae_file_m.f90): read_ae_file: ERROR c_light differs from the
!> input
 integer, parameter, public :: msg_ae_file_c_light=7
!> open_ae (ae_file_m.f90): read_ae_file: ERROR cannot open
 integer, parameter, public :: msg_ae_file_open_read=8
!> open_ae (ae_file_m.f90): read_ae_file: ERROR not an all-electron data file:
 integer, parameter, public :: msg_ae_file_not_ae=9
!> open_ae (ae_file_m.f90): read_ae_file: ERROR file was written on a machine
!> with the other byte order
 integer, parameter, public :: msg_ae_file_byte_order=10
!> open_ae (ae_file_m.f90): read_ae_file: ERROR file format version but this
!> program reads version
 integer, parameter, public :: msg_ae_file_version=11
!> check_found (ae_file_m.f90): read_ae_file: ERROR test configuration is not
!> in the all-electron data file
 integer, parameter, public :: msg_ae_file_test_config=12
!> ae_bound_projectors (ae_projectors_m.f90): oncvpsp-387: lschfb convergence
!> ERROR n,l,iter=
 integer, parameter, public :: msg_ae_proj_bound_conv=13
!> ae_well_projectors (ae_projectors_m.f90): oncvpsp: ERROR debl = for l= ERROR
!> not allowed with 2 or more scattering states program will stop
 integer, parameter, public :: msg_ae_proj_debl_scattering=14
!> ae_well_projectors (ae_projectors_m.f90): ae_well_projectors: ERROR well_n
!> for l= is , below the first well state's n,
 integer, parameter, public :: msg_ae_proj_well_n=15
!> ae_well_projectors (ae_projectors_m.f90): oncvpsp_r 474: ldiracfb
!> convergence ERROR, n= l= kap= epa=
 integer, parameter, public :: msg_ae_proj_well_conv=16
!> free_tail_state (ae_projectors_m.f90): free_tail_state: ERROR well_tail =
!> "free" needs ep >= 0: the well state of l= , n= at ep= Ha has a growing
!> tail; raise its ep (debl), or use well_tail = "confined" or "well"
 integer, parameter, public :: msg_free_tail_ep=17
!> ae_delta_window (ae_projectors_m.f90): ae_delta_window: ERROR node of the
!> local channel l= all-electron wave near rc, at r=
 integer, parameter, public :: msg_delta_window_node=18
!> ecut_weights (ae_projectors_m.f90): ecut_weights: ERROR ecut_weight for l=
!> is for a channel with no occupied state; this one has one
 integer, parameter, public :: msg_ecut_weight_unoccupied=19
!> load_pair (ae_ps_pair_m.f90): : ERROR the files differ in relativistic
!> treatment, atom, mesh size, c_light or functional
 integer, parameter, public :: msg_pair_differ_iexc=20
!> load_pair (ae_ps_pair_m.f90): : ERROR the files differ in relativistic
!> treatment, atom, mesh size or c_light
 integer, parameter, public :: msg_pair_differ=21
!> load_pair (ae_ps_pair_m.f90): : ERROR the files' radial meshes differ
 integer, parameter, public :: msg_pair_mesh=22
!> bs_box_energy (boundsolve_m.f90): bs_box_energy: ERROR no energy with more
!> than nodes up to
 integer, parameter, public :: msg_box_energy_nodes=23
!> bs_box_energy (boundsolve_m.f90): bs_box_energy: ERROR bisection not
!> converged in 200 steps, e=
 integer, parameter, public :: msg_box_energy_bisection=24
!> check_data (check_data_m.f90): ERROR: test_data found ERROR; stopping
 integer, parameter, public :: msg_check_data=25
!> confine_walls (confine_wall_m.f90): confine_walls: ERROR the free atom
!> failed
 integer, parameter, public :: msg_confine_free_atom=26
!> add_edge (confine_ws_m.f90): add_edge: ERROR too many edges in the scan,
 integer, parameter, public :: msg_confine_edges=27
!> const_basis (const_basis_m.f90): const_basis: ERROR the constraints are
!> linearly dependent, singular value ratio
 integer, parameter, public :: msg_const_basis_dependent=28
!> fall_to (diagout_m.f90): diagout: ERROR the core charge's bracket lost its
!> sign
 integer, parameter, public :: msg_diagout_bracket=29
!> origin_series (dirac_m.f90): dirac_outward: ERROR origin series not
!> converged, l= , kap= , ee= , r=
 integer, parameter, public :: msg_dirac_origin_series=30
!> eresid (eresid_m.f90): eresid: ERROR q grid below 0, nq= , dq= , qmax=
 integer, parameter, public :: msg_eresid_q_grid=31
!> eresid (eresid_m.f90): eresid: ERROR tail too short, nr=
 integer, parameter, public :: msg_eresid_tail=32
!> at_qout (eresid_m.f90): eresid: ERROR q off the grid, q= qmax=
 integer, parameter, public :: msg_eresid_q_off_grid=33
!> exc_libxc_init (exc_libxc_m.f90): exc_libxc_init: ERROR Libxc cannot
!> evaluate functional
 integer, parameter, public :: msg_libxc_init=34
!> exc_libxc (exc_libxc_stub_m.f90): exc_libxc_stub: ERROR iexc = requires
!> libxc
 integer, parameter, public :: msg_libxc_stub_exc=35
!> exc_libxc_fxc (exc_libxc_stub_m.f90): exc_libxc_fxc: ERROR without libxc,
!> iexc =
 integer, parameter, public :: msg_libxc_stub_fxc=36
!> exc_libxc_info (exc_libxc_stub_m.f90): exc_libxc_stub: ERROR iexc = requires
!> libxc
 integer, parameter, public :: msg_libxc_stub_info=37
!> free_tail_mesh (free_tail_m.f90): free_tail_mesh: ERROR wave function
!> vanishes too close to rc: rc= , last nonzero point
 integer, parameter, public :: msg_free_tail_mesh=38
!> free_tail_width (free_tail_m.f90): free_tail_width: ERROR qcut= is not above
!> the wave vector of the well state: raise qcut or lower its energy
 integer, parameter, public :: msg_free_tail_qcut=39
!> free_tail_width (free_tail_m.f90): free_tail_width: ERROR the free tail of
!> the well state at Ha with qcut= grows too fast for its window (kappa L = >
!> 12): raise qcut or the state's energy
 integer, parameter, public :: msg_free_tail_growth=40
!> read_error (input_text_m.f90): Read ERROR, input data file line
 integer, parameter, public :: msg_input_text_read=41
!> input_error (input_toml_m.f90): TOML input ERROR:
 integer, parameter, public :: msg_input_toml=42
!> ionic_valence (ionic_tail_m.f90): ionic_valence: ERROR lschvkbb n,l,ierr=
 integer, parameter, public :: msg_ionic_valence_construction=43
!> ionic_valence (ionic_tail_m.f90): ionic_valence: ERROR lschvkbb n,l,ierr=
 integer, parameter, public :: msg_ionic_valence_conv=44
!> ionic_screening (ionic_tail_m.f90): ionic_screening: ERROR the ionic tail
!> did not converge in steps, change =
 integer, parameter, public :: msg_ionic_screening_conv=45
!> ldiracfb (ldiracfb_m.f90): ldiracfb: ERROR ll = > nn =
 integer, parameter, public :: msg_ldiracfb_ll=46
!> ldiracfb (ldiracfb_m.f90): ldiracfb: ERROR kap = ll =
 integer, parameter, public :: msg_ldiracfb_kap=47
!> ldiracfb (ldiracfb_m.f90): ldiracfb: ERROR zz =
 integer, parameter, public :: msg_ldiracfb_zz=48
!> ldiracfs (ldiracfs_m.f90): ldiracfs: ERROR kap = ll =
 integer, parameter, public :: msg_ldiracfs_kap=49
!> ldiracfs (ldiracfs_m.f90): ldiracfs: ERROR zz =
 integer, parameter, public :: msg_ldiracfs_zz=50
!> lm_root (logmesh_m.f90): lm_root: ERROR not converged in 200 steps, r= , f=
 integer, parameter, public :: msg_lm_root_conv=51
!> lschfb (lschfb_m.f90): lschfb: ERROR ll = > nn =
 integer, parameter, public :: msg_lschfb_ll=52
!> write_reference_results (main_setup_m.F90): : ERROR all-electron reference
!> atom not converged
 integer, parameter, public :: msg_ae_reference_conv=53
!> core_radius_points (main_setup_m.F90): core_radius_points: ERROR rc inside
!> the mesh start for l= : rc=
 integer, parameter, public :: msg_core_radius_rc=54
!> match_point (match_m.f90): match_point: ERROR matching window at mch= leaves
!> the mesh, mmax=
 integer, parameter, public :: msg_match_point_window=55
!> match_point (match_m.f90): match_point: ERROR no sign change below mch=
 integer, parameter, public :: msg_match_point_sign=56
!> modcore (modcore1_m.f90): rhomod: ERROR ircc (core-valence charge crossover)
!> not found
 integer, parameter, public :: msg_modcore_ircc=57
!> modcore2 (modcore2_m.f90): modcore2: ERROR core charge never exceeds valence
!> pseudocharge
 integer, parameter, public :: msg_modcore2_core_below=58
!> modcore2 (modcore2_m.f90): modcore2: ERROR ircc (core-valence charge
!> crossover) not found
 integer, parameter, public :: msg_modcore2_ircc=59
!> modcore2 (modcore2_m.f90): modcore2: ERROR no Teter scale for the slope
 integer, parameter, public :: msg_modcore2_teter=60
!> modcore3 (modcore3_m.f90): modcore3: ERROR ircc (core-valence charge
!> crossover) not found
 integer, parameter, public :: msg_modcore3_ircc=61
!> modcore3 (modcore3_m.f90): modcore3: ERROR Nelder-Mead not converged in 100
!> steps
 integer, parameter, public :: msg_modcore3_nelder_mead=62
!> newton_refine (modcore3_m.f90): newton_refine: ERROR not converged in 30
!> steps, last step=
 integer, parameter, public :: msg_newton_refine_conv=63
!> budget_init (modcore_budget_m.f90): modcore_budget: ERROR ecut beyond Ha
 integer, parameter, public :: msg_budget_ecut=64
!> budget_init (modcore_budget_m.f90): modcore_budget: ERROR no core charge
 integer, parameter, public :: msg_budget_no_core=65
!> budget_init (modcore_budget_m.f90): modcore_budget: ERROR r_s up to bohr
!> needs more nodes of the transform
 integer, parameter, public :: msg_budget_nodes=66
!> budget_tables (modcore_budget_m.f90): modcore_budget: ERROR no core charge
 integer, parameter, public :: msg_budget_tables_no_core=67
!> no_model (modcore_budget_m.f90): modcore_budget: ERROR no model core inside
!> rc = needs less than Ha for an error of Ha
 integer, parameter, public :: msg_budget_no_model=68
!> set_model (modcore_budget_m.f90): modcore_budget: ERROR r_s beyond the table
 integer, parameter, public :: msg_budget_rs_table=69
!> tail_monotone (modcore_budget_m.f90): modcore_budget: ERROR dhseqr info =
 integer, parameter, public :: msg_budget_dhseqr=70
!> min_tail (modcore_budget_m.f90): modcore_budget: ERROR dgglse info =
 integer, parameter, public :: msg_budget_dgglse=71
!> zeroin (modcore_budget_m.f90): modcore_budget: ERROR no root between
 integer, parameter, public :: msg_budget_no_root=72
!> zeroin (modcore_budget_m.f90): modcore_budget: ERROR zeroin not converged in
!> 200 steps, x= , f=
 integer, parameter, public :: msg_budget_zeroin_conv=73
!> modcore_budget_table (modcore_budget_m.f90): modcore_budget_table: ERROR no
!> pseudo-wave reaches Ha per electron in its table
 integer, parameter, public :: msg_budget_table_unreached=74
!> modcore_budget_table (modcore_budget_m.f90): modcore_budget_table: ERROR
!> pseudo-wave has a weight and no convergence curve
 integer, parameter, public :: msg_budget_table_weight_curve=75
!> modcore_budget_table (modcore_budget_m.f90): modcore_budget_table: ERROR no
!> weighted curves
 integer, parameter, public :: msg_budget_table_no_weights=76
!> modcore_budget_table (modcore_budget_m.f90): modcore_budget_table: ERROR the
!> weighted mean does not reach Ha per electron
 integer, parameter, public :: msg_budget_table_mean_unreached=77
!> modcore_budget_table (modcore_budget_m.f90): modcore_budget_table: ERROR the
!> bracket lost its sign change
 integer, parameter, public :: msg_budget_table_bracket=78
!> modcore_budget_table (modcore_budget_m.f90): modcore_budget_table: ERROR the
!> weighted mean is below Ha per electron where the curves start
 integer, parameter, public :: msg_budget_table_mean_below=79
!> poly_core_fit (modcore_m.f90): poly_core_fit: ERROR singular match, dgesv
!> info =
 integer, parameter, public :: msg_poly_core_dgesv=80
!> poly_core_fit (modcore_m.f90): poly_core_fit: ERROR dhseqr info =
 integer, parameter, public :: msg_poly_core_dhseqr=81
!> teter_scale (modcore_m.f90): teter_scale: ERROR not converged
 integer, parameter, public :: msg_teter_scale_conv=82
!> ode_rc_der (ode_rc_der_m.f90): ode_rc_der: ERROR uu and up disagree for l=
!> at rc= uu' from up and from stencils on uu:
 integer, parameter, public :: msg_ode_rc_der_disagree=83
!> complete_input (oncvpsp_run_m.f90): oncvpsp_run: ERROR well_norm =
!> "space" (1) is not defined for well_tail = "free": the free tail has no
!> finite norm
 integer, parameter, public :: msg_run_well_norm_free=84
!> complete_input (oncvpsp_run_m.f90): oncvpsp_run: ERROR well_n is for
!> well_tail = "confined" only
 integer, parameter, public :: msg_run_well_n_confined=85
!> run_ae (oncvpsp_run_m.f90): oncvpsp_run: ERROR aeatom failed for the
!> reference configuration, ierr= , it= , res=
 integer, parameter, public :: msg_run_aeatom=86
!> optimize (optimize_m.f90): optimize: residual-energy eigenvalue ERROR, info=
 integer, parameter, public :: msg_optimize_dsyev_residual=87
!> optimize (optimize_m.f90): optimize: residual-energy eigenvalue ERROR, info=
 integer, parameter, public :: msg_optimize_dsyev_regularized=88
!> solve_norm_constraint (optimize_m.f90): optimize: ERROR ps0norm > uunorm,
!> code will stop
 integer, parameter, public :: msg_norm_constraint_ps0norm=89
!> solve_norm_constraint (optimize_m.f90): optimize: ERROR interval-halving
!> search not converged
 integer, parameter, public :: msg_norm_constraint_halving=90
!> solve_norm_constraint (optimize_m.f90): optimize: ERROR Newton polish not
!> converged, last step
 integer, parameter, public :: msg_norm_constraint_newton=91
!> om_local_values (outmesh_m.f90): om_local_values: ERROR the pieces end
!> inside the mesh
 integer, parameter, public :: msg_om_local_pieces=92
!> check_range (outmesh_m.f90): check_range: ERROR output mesh to bohr, beyond
!> the log mesh,
 integer, parameter, public :: msg_outmesh_beyond=93
!> check_range (outmesh_m.f90): check_range: ERROR more than one output point
!> below the log mesh
 integer, parameter, public :: msg_outmesh_below=94
!> origin (outmesh_m.f90): origin: ERROR log mesh too short for the value at
!> the origin
 integer, parameter, public :: msg_outmesh_origin=95
!> outwf_setup (outwf_m.F90): outwf_setup: ERROR the potential's pieces end
!> inside the mesh
 integer, parameter, public :: msg_outwf_setup_pieces=96
!> breakpoints (outwf_m.F90): outwf_setup: ERROR the sources' segments differ
 integer, parameter, public :: msg_outwf_setup_segments=97
!> breakpoints (outwf_m.F90): outwf_setup: ERROR too many breakpoints
 integer, parameter, public :: msg_outwf_setup_breakpoints=98
!> first_grid (outwf_m.F90): outwf_setup: ERROR breakpoint below the 4th node,
!> r=
 integer, parameter, public :: msg_outwf_setup_first_node=99
!> first_grid (outwf_m.F90): outwf_setup: ERROR breakpoint too close to the
!> mesh end, r=
 integer, parameter, public :: msg_outwf_setup_first_end=100
!> later_grid (outwf_m.F90): outwf_setup: ERROR breakpoint too close to the
!> mesh end, r=
 integer, parameter, public :: msg_outwf_setup_later_end=101
!> later_grid (outwf_m.F90): outwf_setup: ERROR breakpoint too close to the
!> mesh end, r=
 integer, parameter, public :: msg_outwf_setup_later_stencil=102
!> outwf_index (outwf_m.F90): outwf_index: ERROR R inside the sources, R, rsrc=
 integer, parameter, public :: msg_outwf_index_sources=103
!> outwf_index (outwf_m.F90): outwf_index: dsyev ERROR, stopping info =
 integer, parameter, public :: msg_outwf_index_dsyev=104
!> outwf_solve (outwf_m.F90): outwf_solve: dgesv ERROR, stopping info =
 integer, parameter, public :: msg_outwf_solve_dgesv=105
!> common_init (phsft_m.f90): common_init: ERROR the radius is too close to the
!> end of the mesh rq= rr(mmax)=
 integer, parameter, public :: msg_phsft_radius_end=106
!> follow (phsft_m.f90): follow: ERROR the phase turns by more than pi/2
!> between and
 integer, parameter, public :: msg_phsft_phase_turn=107
!> pw_alloc (pieces_m.f90): pw_alloc: ERROR number of segments
 integer, parameter, public :: msg_pw_alloc_segments=108
!> pw_alloc (pieces_m.f90): pw_alloc: ERROR bounds out of order at segment
 integer, parameter, public :: msg_pw_alloc_order=109
!> pw_product_integral (pieces_m.f90): pw_product_integral: ERROR the segments
!> differ
 integer, parameter, public :: msg_pw_product_segments=110
!> write_ps_file (ps_file_m.f90): write_ps_file: ERROR cannot open
 integer, parameter, public :: msg_ps_file_open_write=111
!> read_ps_file (ps_file_m.f90): read_ps_file: ERROR cannot open
 integer, parameter, public :: msg_ps_file_open_read=112
!> read_ps_file (ps_file_m.f90): read_ps_file: ERROR not a pseudopotential data
!> file:
 integer, parameter, public :: msg_ps_file_not_ps=113
!> read_ps_file (ps_file_m.f90): read_ps_file: ERROR file was written on a
!> machine with the other byte order
 integer, parameter, public :: msg_ps_file_byte_order=114
!> read_ps_file (ps_file_m.f90): read_ps_file: ERROR file format version but
!> this program reads version
 integer, parameter, public :: msg_ps_file_version=115
!> read_pw (ps_file_m.f90): read_ps_file: ERROR number of segments
 integer, parameter, public :: msg_ps_file_segments=116
!> valence_densities (ps_valence_m.f90): oncvpsp-461: lschfb convergence ERROR
!> n,l,iter=
 integer, parameter, public :: msg_valence_bound_conv=117
!> valence_densities (ps_valence_m.f90): oncvpsp: lschvkbb ERROR
 integer, parameter, public :: msg_valence_vkb_construction=118
!> valence_densities (ps_valence_m.f90): oncvpsp: lschvkbb ERROR
 integer, parameter, public :: msg_valence_vkb_conv=119
!> add_zero_r (psmlout_m.f90): add_zero_r: ERROR f0 is not one longer than ff
 integer, parameter, public :: msg_add_zero_r_size=120
!> add_zero_r (psmlout_m.f90): add_zero_r: ERROR log mesh too short for
!> lm_origin
 integer, parameter, public :: msg_add_zero_r_mesh=121
!> pspot (pspot_m.f90): ERROR pspot: try changing psp parameters for this l
 integer, parameter, public :: msg_pspot_node=122
!> qroots (qroots_m.f90): qroots: ERROR - failed to find nroot roots
 integer, parameter, public :: msg_qroots_not_found=123
!> polish (qroots_m.f90): qroots: ERROR Newton polish not converged, ll= , q= ,
!> last step
 integer, parameter, public :: msg_qroots_newton=124
!> quad_gauss_legendre (quad_m.f90): quad_gauss_legendre: ERROR node not
!> converged, nn= , node=
 integer, parameter, public :: msg_gauss_legendre_node=125
!> relatom (relatom_m.f90): relatom: ldiracfb convergence ERROR
!> n,l,kap,iter,ierr=
 integer, parameter, public :: msg_relatom_conv=126
!> renorm_r (renorm_r_m.f90): renorm_r: ERROR wave function is zero, ll= , kap=
 integer, parameter, public :: msg_renorm_r_zero=127
!> root_brent (root_m.f90): root_brent: ERROR not converged in 200 steps, x= ,
!> f=
 integer, parameter, public :: msg_root_brent_conv=128
!> root_fmin (root_m.f90): root_fmin: ERROR not converged in 100 steps, x= ,
!> bracket width=
 integer, parameter, public :: msg_root_fmin_conv=129
!> run_diag (run_diag_m.f90): : ERROR matching radius beyond the mesh, rm=
 integer, parameter, public :: msg_diag_match_radius=130
!> run_diag (run_diag_m.f90): run_diag_r: ldiracfb convergence ERROR
!> n,l,,kap,ierr=
 integer, parameter, public :: msg_diag_dirac_conv=131
!> run_diag (run_diag_m.f90): run_diag: lschfb convergence ERROR n,l,ierr=
 integer, parameter, public :: msg_diag_lschfb_conv=132
!> run_ghosts (run_ghosts_m.f90): run_ghosts: hamiltonian matrix eigenvalue
!> ERROR, info=
 integer, parameter, public :: msg_ghosts_dsyev=133
!> barrier_mesh (run_ghosts_m.f90): run_ghosts: ERROR barrier radius beyond
!> radial mesh, mbar= , mmax=
 integer, parameter, public :: msg_ghosts_barrier_mesh=134
!> run_optimize (run_optimize_m.f90): run_optimize: ERROR q grid longer than
!> the buffer of the convergence curve, nq= , ncurve= (ncurve in oncvpsp_main;
!> the grid grows with qcut*rc)
 integer, parameter, public :: msg_run_optimize_q_grid=135
!> run_phsft (run_phsft_m.f90): : ERROR rxpsh beyond radial mesh , rxpsh= ,
!> rr(mmax)=
 integer, parameter, public :: msg_phsft_rxpsh=136
!> write_waves (run_plot_m.f90): run_plot_r 409: ldiracfb convergence ERROR, n=
!> l= kap- ierr=
 integer, parameter, public :: msg_plot_dirac_conv=137
!> write_waves (run_plot_m.f90): run_plot: lschfb convergence ERROR n,l,ierr=
 integer, parameter, public :: msg_plot_lschfb_conv=138
!> local_potential (run_vkb_m.f90): run_vkb: ERROR node of the local channel's
!> pseudo wave within reach of rc, l= r=
 integer, parameter, public :: msg_vkb_local_node=139
!> channel_vkb (run_vkb_m.f90): run_vkb: B matrix eigenvalue ERROR, info=
 integer, parameter, public :: msg_vkb_b_dsyev=140
!> channel_vkb (run_vkb_m.f90): run_vkb: sovl matrix eigenvalue ERROR, info=
 integer, parameter, public :: msg_vkb_sovl_dsyev=141
!> channel_vkb (run_vkb_m.f90): run_vkb: B^(-1)* matrix eigenvalue ERROR, info=
 integer, parameter, public :: msg_vkb_binv_dsyev=142
!> sbf_basis_con (sbf_basis_con_m.f90): sbfbasis: negative eigenvalue of
!> overlap matrix
 integer, parameter, public :: msg_sbf_basis_overlap=143
!> sbf_rc_der (sbf_m.f90): sbf_rc_der: argument ERROR, llin = >5
 integer, parameter, public :: msg_sbf_rc_der_llin=144
!> sbfq_rc_der (sbf_q_m.f90): sbfq_rc_der: argument ERROR, llin = >4
 integer, parameter, public :: msg_sbfq_rc_der_llin=145
!> series_error (sr_flux_m.f90): sr_outward: ERROR origin series not converged,
!> l= , ee= , r=
 integer, parameter, public :: msg_sr_origin_series=146
!> sr_so_r (sr_so_r_m.f90): sr_so_r: S matrix eigenvalue ERROR, info=
 integer, parameter, public :: msg_sr_so_s_dsyev=147
!> sr_so_r (sr_so_r_m.f90): sr_so_r: A* matrix eigenvalue ERROR, info=
 integer, parameter, public :: msg_sr_so_sr_dsyev=148
!> sr_so_r (sr_so_r_m.f90): sr_so_r: A* matrix eigenvalue ERROR, info=
 integer, parameter, public :: msg_sr_so_so_dsyev=149
!> sratom (sratom_m.f90): sratom123: lschfb convergence ERROR n,l,iter=
 integer, parameter, public :: msg_sratom_conv=150
!> ug_weights (ugrid_m.f90): ug_weights: ERROR fewer nodes than the rule, nn=
 integer, parameter, public :: msg_ug_weights_nodes=151
!> ug_cumint_down (ugrid_m.f90): ug_cumint_down: ERROR fewer nodes than the
!> rule, nn=
 integer, parameter, public :: msg_ug_cumint_nodes=152
!> ug_interp (ugrid_m.f90): ug_interp: ERROR position off the grid, tt= nn=
 integer, parameter, public :: msg_ug_interp_off_grid=153
!> ug_cross_down (ugrid_m.f90): ug_cross_down: ERROR the bracket lost its sign
!> change
 integer, parameter, public :: msg_ug_cross_bracket=154
!> upf_functional (upfout_m.f90): upfout: ERROR iexc = is presently unsupported
!> for UPF output
 integer, parameter, public :: msg_upf_iexc=155
!> vploc_fit (vploc_m.f90): vploc_fit: ERROR lpopt
 integer, parameter, public :: msg_vploc_fit_lpopt=156
!> wm_short (well_mesh_m.f90): wm_short: ERROR the well state at Ha, l= does
!> not decay within bohr
 integer, parameter, public :: msg_well_decay=157
!> wellstate (wellstate_m.f90): wellstate: ERROR E(rwell) jumps across ep, n=
!> l= E-ep=
 integer, parameter, public :: msg_wellstate_jump=158
!> wellstate (wellstate_m.f90): wellstate: ERROR lschfb convergence error, n=
!> l=
 integer, parameter, public :: msg_wellstate_conv=159
!> wellstate (wellstate_m.f90): wellstate: ERROR with the confined well, ep
!> must be above the atom's own state of this n: raise ep or lower well_n
 integer, parameter, public :: msg_wellstate_confined_ep=160
!> target_energy (wellstate_m.f90): wellstate: ERROR a bound state exists for
!> n= l= at energy wellstate: ERROR above the requested ep=
 integer, parameter, public :: msg_wellstate_bound_below=161
!> ws_roots (ws_basis_m.f90): ws_roots: ERROR root_brent failed, ierr=
 integer, parameter, public :: msg_ws_roots_brent=162
!> ws_convergence (ws_basis_m.f90): ws_convergence: ERROR dsyev info
 integer, parameter, public :: msg_ws_convergence_dsyev=163
!> ws_eigenvalues (ws_basis_m.f90): ws_eigenvalues: ERROR dsyev info
 integer, parameter, public :: msg_ws_eigenvalues_dsyev=164
!> exc_dispatch (xc_dispatch_m.f90): : ERROR bad input iexc =
 integer, parameter, public :: msg_xc_dispatch_iexc=165
!> open_output (oncvpsp_c_m.F90): open_output: ERROR cannot open the text
!> output
 integer, parameter, public :: msg_c_output=166
!> read_toml (oncvpsp_c_m.F90): read_toml: ERROR TOML input needs a build
!> with TOML
 integer, parameter, public :: msg_c_no_toml=167
!> lookup (oncvpsp_c_m.F90): lookup: ERROR the handle has no results of the
!> field
 integer, parameter, public :: msg_c_no_results=168
!> lookup (oncvpsp_c_m.F90): lookup: ERROR no field
 integer, parameter, public :: msg_c_field=169
!> check_size (oncvpsp_c_m.F90): check_size: ERROR the field is of the other
!> kind
 integer, parameter, public :: msg_c_kind=170
!> check_size (oncvpsp_c_m.F90): check_size: ERROR the field has values, the
!> buffer
 integer, parameter, public :: msg_c_size=171
!> oncv_write_files (oncvpsp_c_m.F90): oncv_write_files: ERROR the format
!> is psp8, upf or psml, not
 integer, parameter, public :: msg_c_format=172
!> oncv_write_files (oncvpsp_c_m.F90): oncv_write_files: ERROR cannot open
 integer, parameter, public :: msg_c_file=173
!> write_psp_file (psp_files_m.f90): write_psp_file: ERROR no file format
 integer, parameter, public :: msg_psp_file_kind=174
!> check_ae (ae_state_m.f90): check_ae: ERROR the all-electron atom is of
!> another irel
 integer, parameter, public :: msg_ae_given_irel=175
!> check_ae (ae_state_m.f90): check_ae: ERROR the all-electron atom is of
!> another atom, reference configuration or functional
 integer, parameter, public :: msg_ae_given_atom=176
!> check_ae (ae_state_m.f90): check_ae: ERROR the all-electron atom is on
!> another mesh
 integer, parameter, public :: msg_ae_given_mesh=177
!> check_ae (ae_state_m.f90): check_ae: ERROR the all-electron atom has
!> another speed of light
 integer, parameter, public :: msg_ae_given_c_light=178
!> run_channel (oncvpsp_run_m.f90): run_channel: ERROR no channel l=
 integer, parameter, public :: msg_run_channel_l=179
!> source (oncvpsp_c_m.F90): source: ERROR the handle holds no
!> all-electron atom
 integer, parameter, public :: msg_c_no_ae=180
!> oncv_coupled (oncvpsp_c_m.F90): oncv_coupled: ERROR the handle holds no
!> stages
 integer, parameter, public :: msg_c_no_stages=181
!> confine_atom (confine_atom_m.f90): confine_atom: ERROR the radius
!> ... is not beyond the largest core radius
 integer, parameter, public :: msg_confine_radius=182
!> oncv_confine (oncvpsp_c_m.F90): oncv_confine: ERROR the handle holds no
!> run
 integer, parameter, public :: msg_c_no_run=183
!> oncv_confine (oncvpsp_c_m.F90): oncv_confine: ERROR the channels are 0
!> to 3, the walls positive and the window emin < emax with de > 0
 integer, parameter, public :: msg_c_confine_args=184
!> run_channel (oncvpsp_run_m.f90): run_channel: ERROR the mode is 0, 1
!> or 2
 integer, parameter, public :: msg_run_channel_mode=185
!> run_stages (oncvpsp_run_m.f90): run_stages: ERROR the channel l= is
!> not of this input or atom
 integer, parameter, public :: msg_run_stages_channel=186
!> run_channel (oncvpsp_run_m.f90): run_channel: ERROR the local potential
!> is not of this input or atom
 integer, parameter, public :: msg_run_channel_local=187
!> run_channel (oncvpsp_run_m.f90): run_channel: ERROR the tests need mode
!> 1; the window is emin < emax, de > 0
 integer, parameter, public :: msg_run_channel_tests=188
!> oncv_channel_local_toml (oncvpsp_c_m.F90): local_of: ERROR the handle
!> holds no local potential
 integer, parameter, public :: msg_c_no_local=189
!> oncv_assemble_toml (oncvpsp_c_m.F90): oncv_assemble_toml: ERROR a
!> handle holds no channel, two have one l, or an l has none
 integer, parameter, public :: msg_c_assemble=190
!> modcore_budget (modcore_budget_m.f90): modcore_budget: ERROR r_s= bohr
!> outside [r_lo, r_top] = [ , ]
 integer, parameter, public :: msg_budget_rs_range=191
!> modcore_budget (modcore_budget_m.f90): modcore_budget: ERROR the fill is
!> not positive at r_s= bohr
 integer, parameter, public :: msg_budget_fill_positive=192
!> budget_measures (modcore_budget_m.f90): modcore_budget: ERROR the error
!> share s= has its cutoff beyond Ha
 integer, parameter, public :: msg_budget_share_range=193

! The warnings (warn's id), each with its routine, file and message, and
! what the record's index is

!> ae_floor_warning (ae_scf_m.f90): ae_converged: WARNING scf stopped on
!> the round-off floor of the residual
 integer, parameter, public :: wrn_ae_floor=1
!> sratom (sratom_m.f90): sratom: WARNING failed to converge
 integer, parameter, public :: wrn_sratom_conv=2
!> relatom (relatom_m.f90): relatom: WARNING failed to converge
 integer, parameter, public :: wrn_relatom_conv=3
!> config_ae_warning (config_ae_m.f90): run_config: WARNING for AE atom,
!> no output for configuration, self-consistency failed to converge; the
!> configuration
 integer, parameter, public :: wrn_config_ae_scf=4
!> config_ae_warning (config_ae_m.f90): runconfig: WARNING lschfb
!> convergence ERROR; the configuration
 integer, parameter, public :: wrn_config_ae_lschfb=5
!> config_ae_warning (config_ae_m.f90): runconfig: WARNING ldiracfb
!> convergence ERROR; the configuration
 integer, parameter, public :: wrn_config_ae_ldiracfb=6
!> run_config_ps (run_config_m.f90): run_config: WARNING for fully
!> non-local PS atom, no output for configuration; the configuration
 integer, parameter, public :: wrn_config_ps=7
!> psatom (psatom_m.f90): psatom: WARNING lschvkbb convergence error; n
 integer, parameter, public :: wrn_psatom_state=8
!> psatom (psatom_m.f90): psatom: potential failed to converge
 integer, parameter, public :: wrn_psatom_scf=9
!> optimize (optimize_m.f90): optimize: WARNING E_resid flat to round-off
!> at qcut
 integer, parameter, public :: wrn_optimize_flat=10
!> well_share (run_optimize_m.f90): WARNING the well's tail changes it
!> by; the projector
 integer, parameter, public :: wrn_well_tail=11
!> well_share (run_optimize_m.f90): WARNING qcut is not above this well
!> state's wave vector; the projector
 integer, parameter, public :: wrn_well_qcut=12
!> well_share (run_optimize_m.f90): WARNING this well state's free tail
!> grows too fast; the projector
 integer, parameter, public :: wrn_well_free_tail=13
!> target_energy (wellstate_m.f90): WARNING wellstate: localized bound
!> state found; n
 integer, parameter, public :: wrn_well_localized=14
!> target_energy (wellstate_m.f90): WARNING wellstate: moderately
!> localized bound state found; n
 integer, parameter, public :: wrn_well_moderate=15
!> target_energy (wellstate_m.f90): WARNING wellstate: negative energy
!> specified; n
 integer, parameter, public :: wrn_well_negative=16
!> run_diag (run_diag_m.f90): run_diag: lschvkbb ERROR; the projector
 integer, parameter, public :: wrn_diag_lschvkbb=17
!> run_diag (run_diag_m.f90): run_diag: lschvkbbe ERROR; the projector
 integer, parameter, public :: wrn_diag_lschvkbbe=18
!> bound_tests (run_ghosts_m.f90): run_ghosts: lschvkbb ERROR; n
 integer, parameter, public :: wrn_ghost_lschvkbb=19
!> bound_tests (run_ghosts_m.f90): WARNING - GHOST(-)
 integer, parameter, public :: wrn_ghost_minus=20
!> run_ghosts (run_ghosts_m.f90): WARNING - GHOST(+)
 integer, parameter, public :: wrn_ghost_plus=21
!> barrier_basis (run_ghosts_m.f90): run_ghosts: WARNING basis cut there
 integer, parameter, public :: wrn_ghost_basis_cut=22
!> barrier_basis (run_ghosts_m.f90): run_ghosts: WARNING basis full below
!> the cutoff
 integer, parameter, public :: wrn_ghost_basis_full=23
!> write_waves (run_plot_m.f90): write_waves: WARNING lschvkbb did not
!> converge; the projector
 integer, parameter, public :: wrn_plot_lschvkbb=24
!> sr_so_r (sr_so_r_m.f90): sr_so_r: WARNING directions of the
!> projectors' overlap at its round-off dropped
 integer, parameter, public :: wrn_sr_so_drop=25
!> budget_hardness7 (modcore_budget_m.f90): modcore_budget: WARNING the
!> radius of the least mismatch does not fit the budget
 integer, parameter, public :: wrn_budget_misfit=26
!> minimum_radius (modcore_budget_m.f90): modcore_budget: WARNING no
!> minimum of the mismatch bracketed
 integer, parameter, public :: wrn_budget_no_bracket=27
!> zeroin (modcore_budget_m.f90): modcore_budget: WARNING root at a jump
!> of the function
 integer, parameter, public :: wrn_budget_jump=28
!> modcore_budget (modcore_budget_m.f90): modcore_budget: WARNING icmod 6
!> and 7 are deprecated: icmod 8 and 9 take r_s
 integer, parameter, public :: wrn_budget_deprecated=29

 contains

 pure function failed(kfail,msg,ierr) result(ifail)

!! the ifail of a failure: its kind, its message and the solver status it
!! stops on, if any (ifail of the module's description)
!!
!! @param kfail  the kind of failure
!! @param msg  the message (msg_*)
!! @param ierr  the solver status, if the failure stops on one

!Input variables
 integer, intent(in) :: kfail,msg
 integer, intent(in), optional :: ierr

!Output variables
 integer :: ifail

 ifail=kfail+msg_scale*msg
 if(present(ierr)) then
   if(ierr>0 .and. ierr<=ierr_limit) ifail=ifail+msg_scale*msg_limit*ierr
 end if
 return
 end function failed

 pure function fail_kind(ifail) result(kfail)

!! the kind of failure of an ifail (fail_*)
!!
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

!Input variables
 integer, intent(in) :: ifail

!Output variables
 integer :: kfail

 kfail=mod(ifail,msg_scale)
 return
 end function fail_kind

 pure function fail_msg(ifail) result(msg)

!! the message of an ifail (msg_*), 0 if it has none
!!
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

!Input variables
 integer, intent(in) :: ifail

!Output variables
 integer :: msg

 msg=mod(ifail/msg_scale,msg_limit)
 return
 end function fail_msg

 pure function fail_ierr(ifail) result(ierr)

!! the solver status an ifail stopped on, 0 if none
!!
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

!Input variables
 integer, intent(in) :: ifail

!Output variables
 integer :: ierr

 ierr=ifail/(msg_scale*msg_limit)
 return
 end function fail_ierr

 pure function fail_context(ifail,at) result(failure)

!! a failure, for a binding: ifail split into its kind, message and
!! status, where the run was (at) beside them
!!
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
!! @param at  where the run was (failure_m): its stage, l and ikap

!Input variables
 integer, intent(in) :: ifail
 type(failure_t), intent(in) :: at

!Output variables
 type(failure_t) :: failure

 failure=at
 failure%kind=fail_kind(ifail)
 failure%msg=fail_msg(ifail)
 failure%ierr=fail_ierr(ifail)
 return
 end function fail_context

 pure subroutine warn(warns,id,ll,ikap,idx)

!! adds a warning to warns, with warns's stage, and its l and ikap as
!! given or else warns's
!!
!! @param warns  the warnings, to which it adds its own (failure_m)
!! @param id  the warning (wrn_*)
!! @param ll  angular momentum, if not warns's
!! @param ikap  its spin-orbit index, if not warns's
!! @param idx  the configuration, n or projector it is of (wrn_*), if any

!Input variables
 integer, intent(in) :: id
 integer, intent(in), optional :: ll,ikap,idx

!Input/output variables
 type(warnings_t), intent(inout) :: warns

 warns%nwarn=warns%nwarn+1
 if(warns%nwarn>mxwarn) return
 warns%rec(:,warns%nwarn)=[id,warns%stage,warns%ll,warns%ikap,0]
 if(present(ll)) warns%rec(3,warns%nwarn)=ll
 if(present(ikap)) warns%rec(4,warns%nwarn)=ikap
 if(present(idx)) warns%rec(5,warns%nwarn)=idx
 return
 end subroutine warn

 pure subroutine warn_merge(warns,more)

!! adds the warnings more to warns, after its own
!!
!! @param warns  the warnings, to which it adds its own (failure_m)
!! @param more  those to add

!Input variables
 type(warnings_t), intent(in) :: more

!Input/output variables
 type(warnings_t), intent(inout) :: warns

!Local variables
 integer :: ii

 do ii=1,min(more%nwarn,mxwarn)
   warns%nwarn=warns%nwarn+1
   if(warns%nwarn<=mxwarn) warns%rec(:,warns%nwarn)=more%rec(:,ii)
 end do
 warns%nwarn=warns%nwarn+max(more%nwarn-mxwarn,0)
 return
 end subroutine warn_merge

 pure function fail_name(kfail) result(name)

!! the name of a kind of failure, for messages and bindings
!!
!! @param kfail  the kind of failure

!Input variables
 integer, intent(in) :: kfail

!Output variables
 character(len=12) :: name

 select case(kfail)
 case(fail_none)
   name='none'
 case(fail_input)
   name='input'
 case(fail_convergence)
   name='convergence'
 case(fail_not_found)
   name='not_found'
 case(fail_construction)
   name='construction'
 case(fail_numerical)
   name='numerical'
 case(fail_io)
   name='io'
 case(fail_internal)
   name='internal'
 case default
   name='unknown'
 end select
 return
 end function fail_name

 end module failure_m

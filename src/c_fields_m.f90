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
!> The fields of the results that the C interface copies out
!> (oncvpsp_c_m): for each its kind, its name, its axes, its units and what
!> it is, which oncv_field_name and oncv_field_info give a binding.

 module c_fields_m
 implicit none
 private
 public :: field_info_t, fields

!> a field of the results
 type :: field_info_t
!> its kind (r real, i integer, a logical as 0 or 1), a blank and its name,
!> the component of results_t (pot.*, tests.*: pseudo_m, test_results_m),
!> of the all-electron atom (ae.*: ae_state_t), of a channel (channel.*:
!> channel_t) or of the confined atom's summary (confine.*:
!> confine_summary_t)
   character(len=32) :: name
!> the names of its axes, blank-separated, the first index first: r the
!> radial mesh, iproj the projector, l, ikap (1 for j=l+1/2 and nkap 1, 2
!> for j=l-1/2), istate the state, level the error level, iq the q grid,
!> ie the energy, icfg the test configuration, irow its row, ighost,
!> deriv the derivative, iv the valence state, iwarn the warning, icrad
!> the basin radius, ipow the power; an axis of different quantities as
!> their names in brackets, [ecut,err]
   character(len=48) :: axes
!> its units, blank for none, as the axes for an axis of quantities
   character(len=24) :: units
!> what it is, one sentence
   character(len=160) :: doc
 end type field_info_t

!> the run's pseudopotential and its files' data
 type(field_info_t), parameter :: fields_run(*)=[ &
&  field_info_t('r pot.zz','','', &
&    'atomic number'), &
&  field_info_t('r pot.zion','','e', &
&    'ion charge'), &
&  field_info_t('r pot.c_light','','', &
&    'speed of light of the run (Hartree atomic units)'), &
&  field_info_t('r pot.rcmax','','bohr', &
&    'largest core radius'), &
&  field_info_t('r pot.rr','r','bohr', &
&    'log radial mesh'), &
&  field_info_t('i pot.nproj','l','', &
&    'number of projectors of each l'), &
&  field_info_t('r pot.evkb','iproj l ikap','Ha', &
&    'the projectors'' coefficients: the nonlocal operator of l and '// &
&    'ikap is sum_i vkb_i(r) evkb_i vkb_i(r'')'), &
&  field_info_t('r pot.vkb','r iproj l ikap','bohr^-1/2', &
&    'the projectors r*beta on the mesh, normalized to integral '// &
&    'vkb**2 dr = 1'), &
&  field_info_t('r pot.esr','iproj l','Ha', &
&    'the scalar-relativistic projectors'' coefficients of the fully '// &
&    'relativistic potential, 2*mxprj per l'), &
&  field_info_t('r pot.eso','iproj l','Ha', &
&    'the spin-orbit projectors'' coefficients, as pot.esr'), &
&  field_info_t('r pot.vpuns','r l','Ha', &
&    'the unscreened semi-local potentials, l = 0..mxl+1; column '// &
&    'lloc+1 is the local potential (lloc 4 the polynomial)'), &
&  field_info_t('r pot.rho','r','e/bohr^3', &
&    'the valence pseudo-density of the reference configuration, '// &
&    '4*pi times the density: integral rho r**2 dr is the charge'), &
&  field_info_t('r pot.rhomod','r deriv','e/bohr^3', &
&    'the model core charge, 4*pi times the density, and its first '// &
&    'four radial derivatives'), &
&  field_info_t('r epstot','','Ha', &
&    'the pseudo-atom''s total energy'), &
&  field_info_t('r ea','istate ikap','Ha', &
&    'the all-electron eigenvalues of the reference configuration, '// &
&    'the UPF''s orbital energies'), &
&  field_info_t('r rhotae','r','e/bohr^3', &
&    'the all-electron valence density, 4*pi times the density'), &
&  field_info_t('r rhoc','r','e/bohr^3', &
&    'the core density, 4*pi times the density, that of the psp8 '// &
&    'file'), &
&  field_info_t('r rcore','','bohr', &
&    'the PSML file''s radius beyond which the model core is the core'), &
&  field_info_t('r vpsml','r l ikap','Ha', &
&    'the PSML file''s semi-local potentials'), &
&  field_info_t('r uupsa','r ikap iv','bohr^-1/2', &
&    'the pseudo-atomic orbitals r*R of the valence states, the '// &
&    'UPF''s')]

!> the run's tests
 type(field_info_t), parameter :: fields_tests(*)=[ &
&  field_info_t('i tests.logder.ne','','', &
&    'number of energies of the log derivatives, 0 with skip_ld'), &
&  field_info_t('r tests.logder.ee','ie','Ha', &
&    'their energies, from epsh2 down by depsh'), &
&  field_info_t('r tests.logder.rl','l','bohr', &
&    'the radius of each l the log derivatives are taken at'), &
&  field_info_t('r tests.logder.ae','ie l ikap','rad', &
&    'the all-electron log derivative at rl, atan(r psi''/psi)'), &
&  field_info_t('r tests.logder.ps','ie l ikap','rad', &
&    'the pseudo log derivative at rl, atan(r psi''/psi)'), &
&  field_info_t('i tests.ghosts.nrow','l ikap','', &
&    'number of bound states of the nonlocal Schroedinger equation'), &
&  field_info_t('i tests.ghosts.nn','istate l ikap','', &
&    'their n'), &
&  field_info_t('r tests.ghosts.ee','[e_bound,e_basis] istate l ikap','Ha', &
&    'each bound state''s energy and that of the barrier-basis '// &
&    'eigenstate matched to it; e_basis below e_bound is a GHOST(-)'), &
&  field_info_t('r tests.ghosts.least','[energy,r_rc] l ikap','[Ha,1]', &
&    'of the basis eigenstates more than 0.05 Ha above the '// &
&    'reference, the one of least <r>/rc: its energy and <r>/rc, a '// &
&    'GHOST(+) below 1; 0 if none'), &
&  field_info_t('i tests.configs.ncfg','','', &
&    'number of test configurations, the reference first'), &
&  field_info_t('i tests.configs.nrow','icfg','', &
&    'the rows of each one''s eigenvalue table, -1 if either atom '// &
&    'failed'), &
&  field_info_t('i tests.configs.state','[n,l,ikap] irow icfg','', &
&    'n, l and ikap of each row'), &
&  field_info_t('r tests.configs.ee','[f,e_ae,e_ps] irow icfg','[e,Ha,Ha]', &
&    'the occupation, the all-electron and the pseudo eigenvalue of '// &
&    'each row'), &
&  field_info_t('r tests.configs.de','[de_ae,de_ps,err] icfg','Ha', &
&    'the all-electron and pseudo total-energy differences from the '// &
&    'reference and the excitation error de_ae - de_ps'), &
&  field_info_t('r tests.convergence.cross', &
&    '[ecut,err] level iproj l ikap','Ha', &
&    '(Ecut, err) where the residual kinetic energy per electron '// &
&    'crosses each level, 1e-2 down by half decades; Ecut 0 for a '// &
&    'level not reached'), &
&  field_info_t('r tests.convergence.curve','[ecut,err] iq iproj l ikap','Ha', &
&    'the convergence curves, (Ecut, residual kinetic energy per '// &
&    'electron) on the q grid, 0 beyond their ends'), &
&  field_info_t('r tests.convergence.tail','[tail,check] iproj l ikap','Ha', &
&    'the residual from q_infinity to three times it, and the '// &
&    'r-space less the Fourier kinetic energy to three times '// &
&    'q_infinity'), &
&  field_info_t('r tests.convergence.opt', &
&    '[eresid,reg,e1,f1] iproj l ikap','[Ha,Ha,Ha,1]', &
&    'the residual at qcut, the regularization''s term, and the '// &
&    'smallest eigenvalue e1 and |f1| of the minimization at qcut'), &
&  field_info_t('r tests.convergence.wshare', &
&    '[well,free,flag] iproj l ikap','[Ha,Ha,1]', &
&    'of a well state, the residual at qcut with the well''s tail and'// &
&    ' with the free tail; flag 1 for these, -1 where the free tail '// &
&    'is undefined, 0 otherwise'), &
&  field_info_t('i tests.weighted_cut.used','','', &
&    '1 if the model core''s budget took the weighted cutoff'), &
&  field_info_t('r tests.weighted_cut.errw','','Ha', &
&    'its error per electron'), &
&  field_info_t('r tests.weighted_cut.ecutw','','Ha', &
&    'the weighted cutoff, the weighted mean of the pseudo-waves'' '// &
&    'convergence curves'), &
&  field_info_t('r tests.weighted_cut.weight','iproj l ikap','', &
&    'the weight of each pseudo-wave in the mean'), &
&  field_info_t('r tests.vloc_cut.cut','[ecut,err] level','Ha', &
&    '(Ecut, err) of the local energy of the reference pseudo '// &
&    'valence density per electron, 1e-2 down by half decades; Ecut '// &
&    '0 for a level not reached'), &
&  field_info_t('r tests.vloc_cut.check','','Ha', &
&    'the local energy per electron from the q grid less that from r'// &
&    ' space, the error of the transforms'), &
&  field_info_t('r tests.bbeig','iproj l ikap','', &
&    'the eigenvalues of the projectors'' symmetrized B matrices, '// &
&    'their condition'), &
&  field_info_t('i tests.ghosts.nminus','l ikap','', &
&    'number of GHOST(-) states, basis eigenvalues below a bound '// &
&    'state (or the vacuum level)'), &
&  field_info_t('r tests.ghosts.minus','ighost l ikap','Ha', &
&    'their energies, the first 20'), &
&  field_info_t('i tests.ghosts.nplus','l ikap','', &
&    'number of GHOST(+) states, basis states above the reference '// &
&    'with weighted <r>/rc below 1'), &
&  field_info_t('r tests.ghosts.plus','[energy,r_rc] ighost l ikap','[Ha,1]', &
&    'their energies and weighted <r>/rc, the first 20'), &
&  field_info_t('r tests.ghosts.ecut','l ikap','Ha', &
&    'the energy the basis goes to, the cutoff of the test'), &
&  field_info_t('i tests.ghosts.basis','l ikap','', &
&    'the basis: 0 complete to its cutoff, 1 cut at a barrier state '// &
&    'that did not converge, 2 full below the cutoff')]

!> the run's model core and warnings
 type(field_info_t), parameter :: fields_core(*)=[ &
&  field_info_t('r core.teter', &
&    '[fcfact,rcfact,rmatch,rhocmatch]','[1,1,bohr,e/bohr^3]', &
&    'icmod 3, 4: the Teter function''s amplitude and scale '// &
&    'prefactors, the core-valence crossover radius and the core '// &
&    'charge there'), &
&  field_info_t('r core.budget.ecut','','Ha', &
&    'icmod 6 to 9: the cutoff budget'), &
&  field_info_t('r core.budget.error','','Ha', &
&    'the exchange-correlation energy error that defines the needed '// &
&    'cutoff'), &
&  field_info_t('r core.budget.slack','','', &
&    'the share of the budget the model must fit'), &
&  field_info_t('r core.budget.ecutw','','Ha', &
&    'the pseudo-waves'' cutoff the budget is factor times when '// &
&    'fcfact <= 0, else 0'), &
&  field_info_t('r core.budget.factor','','', &
&    'that multiple, 0 when fcfact > 0'), &
&  field_info_t('r core.budget.r_s','','bohr', &
&    'the model radius r_s'), &
&  field_info_t('r core.budget.rc_min','','bohr', &
&    'the smallest rc'), &
&  field_info_t('r core.budget.r_core','','bohr', &
&    'the last radius with core charge when it is below rc_min and '// &
&    'so bounds r_s, else 0'), &
&  field_info_t('i core.budget.fill','','', &
&    'the fill inside r_s: 1 the minimum tail, 2 Martins'''), &
&  field_info_t('r core.budget.cpoly','ipow','', &
&    'the minimum tail''s coefficients c0..c8 of (r/r_s)**(2k)'), &
&  field_info_t('r core.budget.martins','[a,b,c]','[1,bohr^-2,bohr^-4]', &
&    'Martins'' fill exp(a + b r**2 + c r**4)'), &
&  field_info_t('r core.budget.ecut_needed','','Ha', &
&    'the cutoff the model core needs'), &
&  field_info_t('r core.budget.r_s6','','bohr', &
&    'icmod 7: the radius of icmod 6, the smallest that fits'), &
&  field_info_t('r core.budget.mismatch','','Ha', &
&    'icmod 7: the rms hardness mismatch of the model'), &
&  field_info_t('r core.budget.crad','icrad','bohr', &
&    'icmod 7: the radii of the basin'), &
&  field_info_t('i core.budget.cfill','icrad','', &
&    'their fills, 0 for none'), &
&  field_info_t('r core.budget.cmis','icrad','Ha', &
&    'their hardness mismatches'), &
&  field_info_t('i warnings.nwarn','','', &
&    'number of the run''s warnings'), &
&  field_info_t('i warnings.rec','[warning,stage,l,ikap,index] iwarn','', &
&    'the first 64 warnings, in the order written: the warning '// &
&    '(wrn_*), the stage, l (-1 none), ikap (0 none) and the '// &
&    'configuration or n it is of')]

!> the all-electron atom
 type(field_info_t), parameter :: fields_ae(*)=[ &
&  field_info_t('i ae.irel','','', &
&    '0 non-relativistic, 1 scalar-relativistic, 2 fully '// &
&    'relativistic'), &
&  field_info_t('i ae.nkap','','', &
&    'number of kappas per l: 1, or 2 for irel 2'), &
&  field_info_t('r ae.zz','','', &
&    'atomic number'), &
&  field_info_t('i ae.nc','','', &
&    'number of core states'), &
&  field_info_t('i ae.nv','','', &
&    'number of valence states'), &
&  field_info_t('i ae.iexc','','', &
&    'the functional'), &
&  field_info_t('i ae.na','istate','', &
&    'principal quantum number of each state, core then valence'), &
&  field_info_t('i ae.la','istate','', &
&    'angular momentum of each state'), &
&  field_info_t('r ae.fa','istate','e', &
&    'occupation of each state'), &
&  field_info_t('r ae.c_light','','', &
&    'speed of light (Hartree atomic units)'), &
&  field_info_t('i ae.mmax','','', &
&    'size of the log radial mesh'), &
&  field_info_t('r ae.al','','', &
&    'log of the mesh ratio'), &
&  field_info_t('r ae.rr','r','bohr', &
&    'log radial mesh'), &
&  field_info_t('i ae.it','','', &
&    'the self-consistent iterations'), &
&  field_info_t('r ae.etot','','Ha', &
&    'total energy'), &
&  field_info_t('r ae.eterm','[kinetic,nuclear,hartree,xc]','Ha', &
&    'the total energy''s terms: kinetic, electron-nucleus, Hartree, '// &
&    'exchange-correlation'), &
&  field_info_t('r ae.res','','Ha', &
&    'the residual of the screening potential'), &
&  field_info_t('r ae.ea','istate ikap','Ha', &
&    'eigenvalues of each state per kappa; for irel 2 the l=0 '// &
&    'states'' second column equals the first'), &
&  field_info_t('r ae.rhoc','r','e/bohr^3', &
&    'the core density, 4*pi times the density'), &
&  field_info_t('r ae.rho','r','e/bohr^3', &
&    'the valence density, 4*pi times the density'), &
&  field_info_t('r ae.vfull','r','Ha', &
&    'the all-electron potential'), &
&  field_info_t('i ae.warnings.nwarn','','', &
&    'number of the all-electron atom''s warnings'), &
&  field_info_t('i ae.warnings.rec','[warning,stage,l,ikap,index] iwarn','', &
&    'the first 64 warnings, as warnings.rec')]

!> a channel and the confined atom
 type(field_info_t), parameter :: fields_stage(*)=[ &
&  field_info_t('i channel.l','','', &
&    'the channel''s l'), &
&  field_info_t('i channel.mkap','','', &
&    'number of kappas: 1, or 2 for irel 2 and l > 0'), &
&  field_info_t('i channel.nproj','','', &
&    'number of projectors'), &
&  field_info_t('i channel.irc','','', &
&    'the last log mesh point at or inside rc'), &
&  field_info_t('i channel.npa','iproj','', &
&    'the principal quantum number of each projector''s state'), &
&  field_info_t('i channel.wellprj','iproj','', &
&    '1 for a projector that is a well state'), &
&  field_info_t('r channel.qmsbf','','bohr^-1', &
&    'the q of the last basis function'), &
&  field_info_t('r channel.epa','iproj ikap','Ha', &
&    'the projectors'' energies'), &
&  field_info_t('r channel.pswf','r iproj ikap','bohr^-1/2', &
&    'the optimized pseudo-waves r*R'), &
&  field_info_t('r channel.vp','r ikap','Ha', &
&    'the screened semi-local potential'), &
&  field_info_t('r channel.cvgplt','[ecut,err] level iproj ikap','Ha', &
&    'the channel''s tests.convergence.cross'), &
&  field_info_t('r channel.ecurve','[ecut,err] iq iproj ikap','Ha', &
&    'the channel''s tests.convergence.curve'), &
&  field_info_t('r channel.ektail','[tail,check] iproj ikap','Ha', &
&    'the channel''s tests.convergence.tail'), &
&  field_info_t('r channel.optdg', &
&    '[eresid,reg,e1,f1] iproj ikap','[Ha,Ha,Ha,1]', &
&    'the channel''s tests.convergence.opt'), &
&  field_info_t('r channel.wshare','[well,free,flag] iproj ikap','[Ha,Ha,1]', &
&    'the channel''s tests.convergence.wshare'), &
&  field_info_t('r channel.wcurve','iproj ikap','', &
&    'the channel''s tests.weighted_cut.weight'), &
&  field_info_t('i channel.warnings.nwarn','','', &
&    'number of the channel''s warnings'), &
&  field_info_t('i channel.warnings.rec', &
&    '[warning,stage,l,ikap,index] iwarn','', &
&    'the first 64 warnings, as warnings.rec'), &
&  field_info_t('r channel.vloc','r','Ha', &
&    'mode projectors: the screened local potential'), &
&  field_info_t('r channel.evkb','iproj ikap','Ha', &
&    'mode projectors: the projectors'' coefficients, as pot.evkb'), &
&  field_info_t('r channel.vkb','r iproj ikap','bohr^-1/2', &
&    'mode projectors: the projectors, as pot.vkb'), &
&  field_info_t('r channel.bbeig','iproj ikap','', &
&    'mode projectors: the channel''s tests.bbeig'), &
&  field_info_t('i confine.nmiss','','', &
&    'number of ws_missing lines, -1 for a part left out'), &
&  field_info_t('i confine.nextra','','', &
&    'number of ws_extra lines'), &
&  field_info_t('i confine.noff','','', &
&    'number of ws_paired_off lines'), &
&  field_info_t('i confine.nwfail','','', &
&    'number of wall_failed lines'), &
&  field_info_t('r confine.max_de','','Ha', &
&    'the largest |e_ps - e_ae| of the edges below 0'), &
&  field_info_t('r confine.cut','[all4,all5,top4,top5]','Ha', &
&    'the cutoffs at 1e-4 and 1e-5 Ha of all the edges below 0, then'// &
&    ' of the top edges alone'), &
&  field_info_t('r confine.dde','','Ha', &
&    'the largest |de_ps - de_ae| of the walls'), &
&  field_info_t('i confine.inf','[all4,all5,top4,top5]','', &
&    '1 where an edge never reaches the tolerance of confine.cut')]

!> those added since: the local potential alone and a channel's tests
!> (api 6), the passes of the budget's search (api 7)
 type(field_info_t), parameter :: fields_tiles(*)=[ &
&  field_info_t('i local.lloc','','', &
&    'l of the local potential, 4 for the polynomial'), &
&  field_info_t('r local.rcloc','','bohr', &
&    'its matching radius'), &
&  field_info_t('r local.vloc','r','Ha', &
&    'the screened local potential on the mesh, before a whole run''s'// &
&    ' ionic tail beyond rcmax'), &
&  field_info_t('i channel.logder.ne','','', &
&    'number of energies of the channel''s log derivatives, 0 without'// &
&    ' its tests'), &
&  field_info_t('r channel.logder.ee','ie','Ha', &
&    'their energies, from the window''s top down'), &
&  field_info_t('r channel.logder.rl','','bohr', &
&    'their radius, as the run''s'), &
&  field_info_t('r channel.logder.ae','ie ikap','rad', &
&    'the all-electron log derivative at rl, atan(r psi''/psi)'), &
&  field_info_t('r channel.logder.ps','ie ikap','rad', &
&    'the pseudo log derivative at rl, on the local potential before'// &
&    ' the ionic tail'), &
&  field_info_t('r channel.logder.rmax','','bohr', &
&    'the largest radius of the local potential they read'), &
&  field_info_t('r channel.logder.rcmax','','bohr', &
&    'the input''s largest core radius, where a whole run''s ionic '// &
&    'tail starts'), &
&  field_info_t('i channel.logder.exact','','', &
&    '1 if they are the input''s whole run''s bit for bit (rmax at or '// &
&    'inside rcmax)'), &
&  field_info_t('i channel.ghosts.nrow','ikap','', &
&    'the channel''s tests.ghosts.nrow, on the local potential before'// &
&    ' the ionic tail'), &
&  field_info_t('i channel.ghosts.nn','istate ikap','', &
&    'the channel''s tests.ghosts.nn'), &
&  field_info_t('r channel.ghosts.ee','[e_bound,e_basis] istate ikap','Ha', &
&    'the channel''s tests.ghosts.ee'), &
&  field_info_t('r channel.ghosts.least','[energy,r_rc] ikap','[Ha,1]', &
&    'the channel''s tests.ghosts.least'), &
&  field_info_t('i channel.ghosts.nminus','ikap','', &
&    'the channel''s tests.ghosts.nminus'), &
&  field_info_t('r channel.ghosts.minus','ighost ikap','Ha', &
&    'the channel''s tests.ghosts.minus'), &
&  field_info_t('i channel.ghosts.nplus','ikap','', &
&    'the channel''s tests.ghosts.nplus'), &
&  field_info_t('r channel.ghosts.plus','[energy,r_rc] ighost ikap','[Ha,1]', &
&    'the channel''s tests.ghosts.plus'), &
&  field_info_t('r channel.ghosts.ecut','ikap','Ha', &
&    'the channel''s tests.ghosts.ecut'), &
&  field_info_t('i channel.ghosts.basis','ikap','', &
&    'the channel''s tests.ghosts.basis'), &
&  field_info_t('i core.budget.passes','','', &
&    'icmod 6, 7: the passes of the search, more than 1 when the '// &
&    'error peaks beyond the first cutoffs (0 for icmod 8, 9)'), &
&  field_info_t('r core.budget.qext','','bohr^-1', &
&    'the last density cutoff the radius was found with'), &
&  field_info_t('r core.budget.pnorm','','Ha', &
&    'icmod 6 to 9: the 16-mean of the xc energy error over the '// &
&    'density cutoffs qwin, ((1/(q_hi-q_b)) int |e|**16 dq)**(1/16), '// &
&    'smooth in the inputs'), &
&  field_info_t('r core.budget.emax','','Ha', &
&    'the largest xc energy error over qwin, C0 in the inputs'), &
&  field_info_t('r core.budget.qwin','[q_b,q_hi]','bohr^-1', &
&    'the window of pnorm and emax'), &
&  field_info_t('r core.budget.share','ishare','', &
&    'the shares s of the cutoff budget of the signed errors'), &
&  field_info_t('r core.budget.err_q','ishare','bohr^-1', &
&    'their density cutoffs, 2 sqrt(2 s ecut)'), &
&  field_info_t('r core.budget.err_signed','ishare','Ha', &
&    'the signed xc energy errors there'), &
&  field_info_t('r core.budget.margin','','', &
&    'the fill''s monotonicity margin, the least of '// &
&    '-(d rho/du)/rho(r_s) over u = (r/r_s)**2 in [0, 1]')]

!> the fields, in the order oncv_field_name numbers them (a new group
!> last, so that the numbers of the others stay)
 type(field_info_t), parameter :: fields(*)=[fields_run,fields_tests, &
&  fields_core,fields_ae,fields_stage,fields_tiles]

 end module c_fields_m

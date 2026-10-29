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
!> The C interface of the library (include/oncvpsp.h): a run handle that
!> holds one calculation, its input given as TOML text, its results and its
!> failure. The handle is the only state: two handles are two independent
!> calculations, in one thread or in several. Every function returns 0 or
!> the ifail of the failure (failure_m); a binding raises from the kind and
!> the message of the record (oncv_failure). Beside whole runs, a handle
!> holds a stage of one (oncvpsp_run_m): the all-electron atom, from which
!> runs and channels start; one channel; or a run's stages, from which the
!> model core and the tests are made again with another model core.

 module oncvpsp_c_m
 use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_char, c_size_t, &
&                 c_double, c_associated, c_loc, c_f_pointer
 use input_m, only: input_t, config_t
 use results_m, only: results_t
 use failure_m, only: failure_t
 use ae_state_m, only: ae_state_t
 use channel_m, only: channel_t
 use run_vkb_m, only: local_t, projectors_t
 use confine_summary_m, only: confine_summary_t
 use oncvpsp_run_m, only: stages_t
 use constants_m, only: dp, mxlin, mxcnf
 implicit none
 private
 public :: oncv_api_version, oncv_new, oncv_free, oncv_set_output, &
&          oncv_run_toml, oncv_failure, oncv_dims, dims_t, oncv_field_name, &
&          oncv_shape, oncv_get_double, oncv_get_int, oncv_ae_toml, &
&          oncv_run_ae_toml, oncv_channel_toml, oncv_coupled, &
&          oncv_set_output_buffer, oncv_output_text, oncv_confine

!> the version of the interface: a change to a function or a struct
!> raises it
 integer(c_int), parameter :: api_version=3

!> 1 if c_double is the kind of the reals, else a division by zero that
!> stops the compilation
 integer, parameter :: same_double=1/merge(1,0,c_double==dp)

!> the largest rank of a field
 integer, parameter :: mxrank=5

!> the fields of the results that oncv_shape and oncv_get_* copy out, each
!> its kind (r real, i integer, a logical as 0 or 1) and its name, the
!> component of results_t (pot.*, tests.*: pseudo_m, test_results_m), of
!> the all-electron atom (ae.*: ae_state_t), of a channel (channel.*:
!> channel_t) or of the confined atom's summary (confine.*:
!> confine_summary_t)
 character(len=32), parameter :: fields(*)=[character(len=32) :: &
&  'r pot.zz', 'r pot.zion', 'r pot.c_light', 'r pot.rcmax', 'r pot.rr', &
&  'i pot.nproj', 'r pot.evkb', 'r pot.vkb', 'r pot.esr', 'r pot.eso', &
&  'r pot.vpuns', 'r pot.rho', 'r pot.rhomod', 'r epstot', 'r ea', 'r rhotae', &
&  'r rhoc', 'r rcore', 'r vpsml', 'r uupsa', 'i tests.logder.ne', &
&  'r tests.logder.ee', 'r tests.logder.rl', 'r tests.logder.ae', &
&  'r tests.logder.ps', 'i tests.ghosts.nrow', 'i tests.ghosts.nn', &
&  'r tests.ghosts.ee', 'r tests.ghosts.least', 'i tests.configs.ncfg', &
&  'i tests.configs.nrow', 'i tests.configs.state', 'r tests.configs.ee', &
&  'r tests.configs.de', 'r tests.convergence.cross', &
&  'r tests.convergence.curve', 'r tests.convergence.tail', &
&  'r tests.convergence.opt', 'r tests.convergence.wshare', &
&  'i tests.weighted_cut.used', 'r tests.weighted_cut.errw', &
&  'r tests.weighted_cut.ecutw', 'r tests.weighted_cut.weight', &
&  'r tests.vloc_cut.cut', 'r tests.vloc_cut.check', 'r tests.bbeig', &
&  'i tests.ghosts.nminus', 'r tests.ghosts.minus', 'i tests.ghosts.nplus', &
&  'r tests.ghosts.plus', 'r tests.ghosts.ecut', 'i tests.ghosts.basis', &
&  'r core.teter', 'r core.budget.ecut', 'r core.budget.error', &
&  'r core.budget.slack', 'r core.budget.ecutw', 'r core.budget.factor', &
&  'r core.budget.r_s', 'r core.budget.rc_min', 'r core.budget.r_core', &
&  'i core.budget.fill', 'r core.budget.cpoly', 'r core.budget.martins', &
&  'r core.budget.ecut_needed', 'r core.budget.r_s6', &
&  'r core.budget.mismatch', 'r core.budget.crad', 'i core.budget.cfill', &
&  'r core.budget.cmis', 'i warnings.nwarn', 'i warnings.rec', &
&  'i ae.irel', 'i ae.nkap', 'r ae.zz', 'i ae.nc', 'i ae.nv', 'i ae.iexc', &
&  'i ae.na', 'i ae.la', 'r ae.fa', 'r ae.c_light', 'i ae.mmax', 'r ae.al', &
&  'r ae.rr', 'i ae.it', 'r ae.etot', 'r ae.eterm', 'r ae.res', 'r ae.ea', &
&  'r ae.rhoc', 'r ae.rho', 'r ae.vfull', 'i ae.warnings.nwarn', &
&  'i ae.warnings.rec', 'i channel.l', 'i channel.mkap', 'i channel.nproj', &
&  'i channel.irc', 'i channel.npa', 'i channel.wellprj', &
&  'r channel.qmsbf', 'r channel.epa', 'r channel.pswf', 'r channel.vp', &
&  'r channel.cvgplt', 'r channel.ecurve', 'r channel.ektail', &
&  'r channel.optdg', 'r channel.wshare', 'r channel.wcurve', &
&  'i channel.warnings.nwarn', 'i channel.warnings.rec', 'r channel.vloc', &
&  'r channel.evkb', 'r channel.vkb', 'r channel.bbeig', &
&  'i confine.nmiss', 'i confine.nextra', 'i confine.noff', &
&  'i confine.nwfail', 'r confine.max_de', 'r confine.cut', 'r confine.dde', &
&  'i confine.inf']

!> a field copied out of the results: its values, flat, and its shape
 type :: field_t
!> whether the name is a field
   logical :: found=.false.
!> whether its values are integers (iv), else reals (rv)
   logical :: isint=.false.
!> its rank, 0 for a scalar
   integer :: rank=0
!> its shape, shp(1:rank)
   integer :: shp(mxrank)=0
!> its real values, in Fortran order
   real(dp), allocatable :: rv(:)
!> its integer values, in Fortran order
   integer, allocatable :: iv(:)
 end type field_t

!> the sizes of a run's results (oncv_dims), for the caller's arrays
 type, bind(c) :: dims_t
!> the version of the interface that filled it (api_version)
   integer(c_int) :: version
!> 0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
   integer(c_int) :: irel
!> number of kappas per l: 1, or 2 for irel=2
   integer(c_int) :: nkap
!> the functional (iexc)
   integer(c_int) :: iexc
!> size of the log radial mesh
   integer(c_int) :: mmax
!> the largest l of the projectors
   integer(c_int) :: lmax
!> the l of the local potential, 4 for the polynomial one
   integer(c_int) :: lloc
!> the projectors of each l, l=0..mxlin
   integer(c_int) :: nproj(mxlin+1)
!> points of the linear output mesh
   integer(c_int) :: nrl
 end type dims_t

!> a calculation: its settings, input, results and failure, or a stage of
!> one
 type :: handle_t
!> the file of the text output
   character(len=:), allocatable :: out_path
!> whether the text output goes to text rather than to out_path
   logical :: buffered=.false.
!> the text output of the last run, with buffered
   character(len=:), allocatable :: text
!> the input values, completed by the run
   type(input_t) :: inp
!> the results
   type(results_t) :: psp
!> the reference and test configurations (complete_input)
   type(config_t) :: cnf(mxcnf)
!> the input as the stages left it, from which oncv_coupled starts
   type(input_t) :: inp_stages
!> the all-electron atom of oncv_ae_toml
   type(ae_state_t) :: ae
!> the stages of a run before the model core (run_stages)
   type(stages_t) :: st
!> the channel of oncv_channel_toml
   type(channel_t) :: chan
!> its local potential, with its projectors
   type(local_t) :: lp
!> its projectors (oncv_channel_toml with mode 1)
   type(projectors_t) :: prj
!> where the run was when it failed
   type(failure_t) :: at=failure_t(0,0,0,0,-1,0,0)
!> the run's ifail, 0 before a run
   integer :: ifail=0
!> the run's relativistic treatment
   integer :: irel=0
!> whether a run finished, so that the results are there
   logical :: done=.false.
!> whether it holds an all-electron atom (oncv_ae_toml)
   logical :: has_ae=.false.
!> whether it holds a run's stages (st, inp_stages and cnf)
   logical :: has_stages=.false.
!> whether it holds a channel (oncv_channel_toml)
   logical :: has_chan=.false.
!> whether the channel has its local potential and projectors
   logical :: has_vkb=.false.
!> the summary of oncv_confine
   type(confine_summary_t) :: sm
!> whether it holds the summary of a confined atom (oncv_confine)
   logical :: has_confine=.false.
 end type handle_t

!> a field's values and shape from an array or scalar of the results
 interface put
   module procedure put_r0, put_r1, put_r2, put_r3, put_r4, put_r5, &
&                   put_i0, put_i1, put_i2, put_i3, put_l0, put_l1
 end interface put

 contains

 function oncv_api_version() result(version) bind(c,name='oncv_api_version')

!! the version of the interface (api_version)

!Output variables
 integer(c_int) :: version

 version=api_version
 return
 end function oncv_api_version

 function oncv_new() result(handle) bind(c,name='oncv_new')

!! a new run handle, its text output going nowhere (/dev/null); free it
!! with oncv_free

!Output variables
 type(c_ptr) :: handle

!Local variables
 type(handle_t), pointer :: hh

 allocate(hh)
 hh%out_path='/dev/null'
 handle=c_loc(hh)
 return
 end function oncv_new

 subroutine oncv_free(handle) bind(c,name='oncv_free')

!! frees a run handle and all it holds; a null handle is ignored
!!
!! @param handle  the run handle

!Input variables
 type(c_ptr), value :: handle

!Local variables
 type(handle_t), pointer :: hh

 if(.not. c_associated(handle)) return
 call c_f_pointer(handle,hh)
 deallocate(hh)
 return
 end subroutine oncv_free

 function oncv_set_output(handle,path,npath) result(ifail) &
&         bind(c,name='oncv_set_output')

!! the file the runs of the handle write their text output to, replaced
!! at each run (/dev/null for none)
!!
!! @param handle  the run handle
!! @param path  the file name, npath characters, not terminated
!! @param npath  its length

!Input variables
 type(c_ptr), value :: handle
 character(kind=c_char), intent(in) :: path(npath)
 integer(c_size_t), value :: npath

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh

 ifail=0
 call c_f_pointer(handle,hh)
 hh%out_path=from_c(path,npath)
 return
 end function oncv_set_output

 function oncv_set_output_buffer(handle,on) result(ifail) &
&         bind(c,name='oncv_set_output_buffer')

!! with on non-zero, the runs of the handle write their text output to a
!! buffer the handle keeps, replaced at each run and read with
!! oncv_output_text, rather than to the file of oncv_set_output; with on
!! 0, to that file again
!!
!! @param handle  the run handle
!! @param on  non-zero for the buffer

!Input variables
 type(c_ptr), value :: handle
 integer(c_int), value :: on

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh

 ifail=0
 call c_f_pointer(handle,hh)
 hh%buffered=on/=0
 return
 end function oncv_set_output_buffer

 function oncv_output_text(handle,text,ntext,nout) result(ifail) &
&         bind(c,name='oncv_output_text')

!! the text output of the handle's last run, or stage, from its buffer
!! (oncv_set_output_buffer): its first ntext bytes into text, not
!! terminated, and its whole length in nout (0 without a buffer), so that a
!! call with ntext 0 gives the size to allocate
!!
!! @param handle  the run handle
!! @param text  the caller's buffer, ntext characters
!! @param ntext  its length
!! @param nout  the length of the text output

!Input variables
 type(c_ptr), value :: handle
 integer(c_size_t), value :: ntext

!Output variables
 character(kind=c_char), intent(out) :: text(ntext)
 integer(c_size_t), intent(out) :: nout
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh
 integer(c_size_t) :: ii

 ifail=0
 call c_f_pointer(handle,hh)
 nout=0
 if(.not. allocated(hh%text)) return
 nout=len(hh%text,kind=c_size_t)
 do ii=1,min(ntext,nout)
   text(ii)=hh%text(ii:ii)
 end do
 return
 end function oncv_output_text

 function oncv_run_toml(handle,irel,text,ntext,skip_ld) result(ifail) &
&         bind(c,name='oncv_run_toml')

!! the pseudopotential and its tests of the TOML input text, with the
!! relativistic treatment irel, into the handle, whose earlier results it
!! replaces (oncvpsp_run)
!!
!! @param handle  the run handle
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param text  the TOML input, ntext characters, not terminated
!! @param ntext  its length
!! @param skip_ld  non-zero to leave out the log derivatives

!Input variables
 type(c_ptr), value :: handle
 integer(c_int), value :: irel,skip_ld
 character(kind=c_char), intent(in) :: text(ntext)
 integer(c_size_t), value :: ntext

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh
 integer :: iout,jfail

 call c_f_pointer(handle,hh)
 call reset(hh,irel)
 call open_output(hh,iout,ifail)
 if(ifail/=0) return
 call read_toml(hh,text,ntext,iout,jfail)
 if(jfail==0) call run_all(hh,skip_ld/=0,iout,jfail)
 call close_output(hh,iout)
 hh%ifail=jfail
 ifail=jfail
 return
 end function oncv_run_toml

 function oncv_ae_toml(handle,irel,text,ntext) result(ifail) &
&         bind(c,name='oncv_ae_toml')

!! the all-electron atom of the TOML input text, with the relativistic
!! treatment irel, into the handle, whose earlier contents it replaces
!! (run_ae): the stage that oncv_run_ae_toml and oncv_channel_toml start
!! from, its fields ae.*
!!
!! @param handle  the run handle
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param text  the TOML input, ntext characters, not terminated
!! @param ntext  its length

 use oncvpsp_run_m, only: run_ae

!Input variables
 type(c_ptr), value :: handle
 integer(c_int), value :: irel
 character(kind=c_char), intent(in) :: text(ntext)
 integer(c_size_t), value :: ntext

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh
 integer :: iout,jfail

 call c_f_pointer(handle,hh)
 call reset(hh,irel)
 call open_output(hh,iout,ifail)
 if(ifail/=0) return
 call read_toml(hh,text,ntext,iout,jfail)
 if(jfail==0) then
   call run_ae(irel,0,hh%inp,'',hh%cnf,hh%ae,hh%at,iout,jfail)
 end if
 call close_output(hh,iout)
 hh%has_ae=jfail==0
 hh%ifail=jfail
 ifail=jfail
 return
 end function oncv_ae_toml

 function oncv_run_ae_toml(handle,ae_handle,text,ntext,skip_ld) &
&         result(ifail) bind(c,name='oncv_run_ae_toml')

!! the pseudopotential and its tests of the TOML input text into the
!! handle, as oncv_run_toml, from the all-electron atom of ae_handle
!! (oncv_ae_toml, or a run's) in place of solving it: the input must be of
!! that atom, its irel, reference configuration, functional and mesh
!! (check_ae). The handle may be ae_handle
!!
!! @param handle  the run handle
!! @param ae_handle  the handle of the all-electron atom
!! @param text  the TOML input, ntext characters, not terminated
!! @param ntext  its length
!! @param skip_ld  non-zero to leave out the log derivatives

!Input variables
 type(c_ptr), value :: handle,ae_handle
 integer(c_int), value :: skip_ld
 character(kind=c_char), intent(in) :: text(ntext)
 integer(c_size_t), value :: ntext

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh,src
 type(ae_state_t) :: ae
 integer :: iout,jfail

 call c_f_pointer(handle,hh)
 call c_f_pointer(ae_handle,src)
 call source(src,ae,jfail)
 call reset(hh,ae%irel)
 if(jfail/=0) then
   hh%ifail=jfail
   ifail=jfail
   return
 end if
 call open_output(hh,iout,ifail)
 if(ifail/=0) return
 call read_toml(hh,text,ntext,iout,jfail)
 if(jfail==0) call run_all(hh,skip_ld/=0,iout,jfail,ae)
 call close_output(hh,iout)
 hh%ifail=jfail
 ifail=jfail
 return
 end function oncv_run_ae_toml

 function oncv_channel_toml(handle,ae_handle,text,ntext,ll,mode) &
&         result(ifail) bind(c,name='oncv_channel_toml')

!! the channel of angular momentum ll alone of the TOML input text into
!! the handle, whose earlier contents it replaces, from the all-electron
!! atom of ae_handle (run_channel): its fields channel.*. The input must be
!! of that atom (check_ae). A first pseudo-wave with a node is the failure
!! of pspot (oncv_failure); the handle may be ae_handle. With mode 1, also
!! the local potential and the channel's projectors (channel.vloc, evkb,
!! vkb, bbeig), whose failures are at the stage of the projectors; with
!! mode 2, only the failures up to the node test of each kappa's first
!! projector (pspot), the handle then holding no fields channel.*
!!
!! @param handle  the run handle
!! @param ae_handle  the handle of the all-electron atom
!! @param text  the TOML input, ntext characters, not terminated
!! @param ntext  its length
!! @param ll  angular momentum
!! @param mode  0 the channel, 1 also the local potential and projectors,
!!              2 the node test of its first projector only

 use oncvpsp_run_m, only: run_channel

!Input variables
 type(c_ptr), value :: handle,ae_handle
 integer(c_int), value :: ll,mode
 character(kind=c_char), intent(in) :: text(ntext)
 integer(c_size_t), value :: ntext

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh,src
 type(ae_state_t) :: ae
 integer :: iout,jfail

 call c_f_pointer(handle,hh)
 call c_f_pointer(ae_handle,src)
 call source(src,ae,jfail)
 call reset(hh,ae%irel)
 if(jfail/=0) then
   hh%ifail=jfail
   ifail=jfail
   return
 end if
 call open_output(hh,iout,ifail)
 if(ifail/=0) return
 call read_toml(hh,text,ntext,iout,jfail)
 if(jfail==0) then
   call run_channel(ae,hh%inp,ll,int(mode),hh%chan,hh%lp,hh%prj,hh%at, &
&                   iout,jfail)
 end if
 call close_output(hh,iout)
 hh%has_chan=jfail==0 .and. mode/=2
 hh%has_vkb=jfail==0 .and. mode==1
 hh%ifail=jfail
 ifail=jfail
 return
 end function oncv_channel_toml

 function oncv_coupled(handle,base,icmod,fcfact,rcfact,skip_ld) &
&         result(ifail) bind(c,name='oncv_coupled')

!! the model core, the screening and the tests again into the handle
!! (run_coupled), from the stages of the run of base (oncv_run_toml,
!! oncv_run_ae_toml or oncv_coupled, whose model core may have failed)
!! with the model core icmod, fcfact and rcfact in place of the input's:
!! the handle then holds those stages and a run's results, as oncv_run_toml
!! with that model core gives them. The handle may be base
!!
!! @param handle  the run handle
!! @param base  the handle of the run whose stages it starts from
!! @param icmod  the model core's kind (fcfact 0 for icmod 0)
!! @param fcfact  its amplitude factor
!! @param rcfact  its radius factor
!! @param skip_ld  non-zero to leave out the log derivatives

 use oncvpsp_run_m, only: run_coupled
 use check_data_m, only: check_data
 use failure_m, only: fail_input, failed, msg_c_no_stages
 use, intrinsic :: iso_fortran_env, only: error_unit

!Input variables
 type(c_ptr), value :: handle,base
 integer(c_int), value :: icmod,skip_ld
 real(c_double), value :: fcfact,rcfact

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh,bb
 type(input_t) :: inp
 real(dp) :: fcgrid(3),rcgrid(3)
 integer :: iout,irel,jfail

 call c_f_pointer(handle,hh)
 call c_f_pointer(base,bb)
 if(.not. bb%has_stages) then
   write(error_unit,'(/a)') 'oncv_coupled: ERROR the handle holds no stages'
   call reset(hh,0)
   hh%ifail=failed(fail_input,msg_c_no_stages)
   ifail=hh%ifail
   return
 end if
 irel=bb%irel
 if(.not. c_associated(handle,base)) then
   hh%st=bb%st
   hh%cnf(:)=bb%cnf(:)
   hh%inp_stages=bb%inp_stages
 end if
 call reset(hh,irel)
 hh%has_stages=.true.
 call open_output(hh,iout,ifail)
 if(ifail/=0) return
 inp=hh%inp_stages
 inp%icmod=icmod
 inp%fcfact=fcfact
 inp%rcfact=rcfact
 if(icmod==0) inp%fcfact=0.0d0
 fcgrid=[inp%fcfact_min,inp%fcfact_max,inp%fcfact_step]
 rcgrid=[inp%rcfact_min,inp%rcfact_max,inp%rcfact_step]
 call check_data(inp,fcgrid,rcgrid,iout,jfail)
 if(jfail==0) then
   call run_coupled(irel,0,inp,hh%cnf,hh%st,'',skip_ld/=0,hh%psp,hh%at, &
&                   iout,jfail)
 end if
 call close_output(hh,iout)
 hh%inp=inp
 hh%done=jfail==0
 hh%ifail=jfail
 ifail=jfail
 return
 end function oncv_coupled

 function oncv_confine(handle,run,radii,nradii,channels,nchannels,energies, &
&         walls,nwalls,r0,ecutmax) result(ifail) bind(c,name='oncv_confine')

!! the confined atom of the run of the handle run (oncv_run_toml,
!! oncv_run_ae_toml or oncv_coupled) into the handle, whose earlier
!! contents it replaces, as oncvpsp_confine.x makes it of the run's
!! all-electron and pseudopotential data files (confine_atom): its lines
!! in the text output and its summary, the fields confine.*. The handle
!! may be run
!!
!! @param handle  the handle
!! @param run  the handle of the run
!! @param radii  the radii of the sphere (bohr), nradii of them
!! @param nradii  their number: 0 for the default, less for no sphere
!! @param channels  the l of the sphere, nchannels of them
!! @param nchannels  their number, 0 or less for all four
!! @param energies  emin, emax and the step of the scan (Ha), or null for
!!                  the input's window of the log derivatives
!! @param walls  the strengths of the wall (Ha/bohr**2), nwalls of them
!! @param nwalls  their number: less than 0 for the default, 0 for none
!! @param r0  the onset of the wall (bohr), 0 or less for the largest
!!            core radius
!! @param ecutmax  the largest cutoff of the basis (Ha), 0 or less for 200

 use confine_atom_m, only: confine_atom, kdef
 use ps_file_m, only: ps_file_t, ps_file_of
 use failure_m, only: fail_input, failed, msg_c_no_run, msg_c_confine_args, &
&                    stage_confine
 use constants_m, only: mxstate, mxkap
 use, intrinsic :: iso_fortran_env, only: error_unit

!Input variables
 type(c_ptr), value :: handle,run,radii,channels,energies,walls
 integer(c_int), value :: nradii,nchannels,nwalls
 real(c_double), value :: r0,ecutmax

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh,rh
 type(ps_file_t) :: pf
 type(confine_summary_t) :: sm
 real(c_double), pointer :: rq(:),ee(:),kw(:)
 integer(c_int), pointer :: ll(:)
 real(dp), allocatable :: rr(:),vfull(:),rs(:),kwall(:)
 real(dp) :: c_light,zz,ew(3),ecut
 real(dp) :: fa(mxstate),ea(mxstate,mxkap)
 integer :: na(mxstate),la(mxstate)
 integer :: irel,mmax,nc,nv,iexc,ii,iout,jfail
 logical :: chan(0:3),bad

 call c_f_pointer(handle,hh)
 call c_f_pointer(run,rh)
 if(.not. (rh%done .and. rh%has_stages)) then
   write(error_unit,'(/a)') 'oncv_confine: ERROR the handle holds no run'
   call reset(hh,0)
   hh%ifail=failed(fail_input,msg_c_no_run)
   ifail=hh%ifail
   return
 end if

!the options, the defaults those of oncvpsp_confine.x
 bad=.false.
 allocate(rs(max(nradii,0)))
 if(nradii>0) then
   call c_f_pointer(radii,rq,[nradii])
   rs(:)=rq(:)
 end if
 chan(:)=nchannels<=0
 if(nchannels>0) then
   call c_f_pointer(channels,ll,[nchannels])
   do ii=1,nchannels
     if(ll(ii)<0 .or. ll(ii)>3) then
       bad=.true.
     else
       chan(ll(ii))=.true.
     end if
   end do
 end if
 ew(:)=0.0d0
 if(c_associated(energies)) then
   call c_f_pointer(energies,ee,[3])
   ew(:)=ee(:)
   bad=bad .or. ew(1)>=ew(2) .or. ew(3)<=0.0d0
 end if
 if(nwalls<0) then
   kwall=kdef
 else
   allocate(kwall(nwalls))
   if(nwalls>0) then
     call c_f_pointer(walls,kw,[nwalls])
     kwall(:)=kw(:)
     bad=bad .or. any(kwall<=0.0d0)
   end if
 end if
 if(bad) then
   write(error_unit,'(/a)') 'oncv_confine: ERROR the channels are 0 to 3, '// &
&       'the walls positive and the window emin < emax with de > 0'
   call reset(hh,0)
   hh%ifail=failed(fail_input,msg_c_confine_args)
   ifail=hh%ifail
   return
 end if
 ecut=200.0d0
 if(ecutmax>0.0d0) ecut=ecutmax

!what oncvpsp_confine.x reads of the run's two files
 call ps_file_of(rh%psp%pot,rh%inp%rc,rh%inp%epsh1,rh%inp%epsh2, &
&                rh%inp%depsh,rh%inp%rxpsh,rh%st%epa,pf)
 irel=rh%irel
 associate(ae=>rh%st%ae)
   c_light=ae%c_light
   zz=ae%zz
   mmax=ae%mmax
   rr=ae%rr(1:mmax)
   vfull=ae%vfull(1:mmax)
   nc=ae%nc
   nv=ae%nv
   iexc=ae%iexc
   na(:)=ae%na(:)
   la(:)=ae%la(:)
   fa(:)=ae%fa(:)
   ea(:,:)=ae%ea(:,:)
 end associate

 call reset(hh,irel)
 call open_output(hh,iout,ifail)
 if(ifail/=0) return
 call confine_atom(pf,irel,c_light,zz,mmax,rr,vfull,nc,nv,iexc,na,la,fa,ea, &
&                  chan,nradii>=0,size(rs),rs,ew,ecut,.false.,size(kwall), &
&                  kwall,r0,sm,iout,jfail)
 call close_output(hh,iout)
 hh%sm=sm
 hh%has_confine=jfail==0
 if(.not. hh%has_confine) hh%at%stage=stage_confine
 hh%ifail=jfail
 ifail=jfail
 return
 end function oncv_confine

 subroutine reset(hh,irel)

!! the handle emptied for a new run or stage of the relativistic treatment
!! irel: no results, stages, atom or channel, and no failure
!!
!! @param hh  the handle
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic

 use failure_m, only: stage_input

!Input variables
 integer, intent(in) :: irel

!Input/output variables
 type(handle_t), intent(inout) :: hh

 hh%done=.false.
 hh%has_ae=.false.
 hh%has_stages=.false.
 hh%has_chan=.false.
 hh%has_vkb=.false.
 hh%has_confine=.false.
 hh%text=''
 hh%ifail=0
 hh%irel=irel
 hh%at=failure_t(0,0,0,stage_input,-1,0,0)
 return
 end subroutine reset

 subroutine open_output(hh,iout,ifail)

!! the handle's text output opened, replaced: its file (oncv_set_output),
!! or a scratch file for its buffer (oncv_set_output_buffer)
!!
!! @param hh  the handle
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written, also the handle's

 use failure_m, only: fail_io, failed, msg_c_output
 use, intrinsic :: iso_fortran_env, only: error_unit

!Input/output variables
 type(handle_t), intent(inout) :: hh

!Output variables
 integer, intent(out) :: iout,ifail

!Local variables
 integer :: ios

 ifail=0
 if(hh%buffered) then
   open(newunit=iout,status='scratch',action='readwrite',iostat=ios)
 else
   open(newunit=iout,file=hh%out_path,action='write',status='replace', &
&       iostat=ios)
 end if
 if(ios/=0) then
   write(error_unit,'(/2a)') &
&       'open_output: ERROR cannot open the text output ',hh%out_path
   ifail=failed(fail_io,msg_c_output)
   hh%ifail=ifail
   return
 end if
 return
 end subroutine open_output

 subroutine close_output(hh,iout)

!! the handle's text output closed, its scratch file read into its buffer
!! first when it has one
!!
!! @param hh  the handle
!! @param iout  unit of the text output

!Input variables
 integer, intent(in) :: iout

!Input/output variables
 type(handle_t), intent(inout) :: hh

!Local variables
 character(len=:), allocatable :: line,text
 integer :: ios,nn,nt

 if(hh%buffered) then
   rewind(iout)
   allocate(character(len=65536) :: text)
   nt=0
   do
     call read_line(iout,line,ios)
     if(ios/=0) exit
     nn=len(line)+1
     do while(nt+nn>len(text))
       text=text//repeat(' ',len(text))
     end do
     text(nt+1:nt+nn)=line//new_line('a')
     nt=nt+nn
   end do
   hh%text=text(1:nt)
 end if
 close(iout)
 return
 end subroutine close_output

 subroutine read_toml(hh,text,ntext,iout,ifail)

!! the TOML input text into the handle's input
!!
!! @param hh  the handle
!! @param text  the TOML input, ntext characters, not terminated
!! @param ntext  its length
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

#ifdef WITH_TOML
 use input_toml_m, only: read_input_toml_string
#else
 use failure_m, only: fail_input, failed, msg_c_no_toml
#endif

!Input variables
 character(kind=c_char), intent(in) :: text(ntext)
 integer(c_size_t), intent(in) :: ntext
 integer, intent(in) :: iout

!Input/output variables
 type(handle_t), intent(inout) :: hh

!Output variables
 integer, intent(out) :: ifail

#ifdef WITH_TOML
 call read_input_toml_string(from_c(text,ntext),hh%inp,iout,ifail)
#else
 write(iout,'(/a)') 'read_toml: ERROR TOML input needs a build with TOML'
 ifail=failed(fail_input,msg_c_no_toml)
#endif
 return
 end subroutine read_toml

 subroutine run_all(hh,skip_ld,iout,ifail,ae_in)

!! the stages of the handle's input (run_stages), kept for oncv_coupled,
!! then the model core and the tests (run_coupled): oncvpsp_run, with the
!! all-electron atom ae_in if given
!!
!! @param hh  the handle
!! @param skip_ld  whether to leave out the log derivatives
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
!! @param ae_in  the all-electron atom solved before, if given

 use oncvpsp_run_m, only: run_stages, run_coupled

!Input variables
 logical, intent(in) :: skip_ld
 integer, intent(in) :: iout
 type(ae_state_t), intent(in), optional :: ae_in

!Input/output variables
 type(handle_t), intent(inout) :: hh

!Output variables
 integer, intent(out) :: ifail

 call run_stages(hh%irel,0,hh%inp,'',hh%cnf,hh%st,hh%at,iout,ifail,ae_in)
 if(ifail/=0) return
 hh%has_stages=.true.
 hh%inp_stages=hh%inp
 call run_coupled(hh%irel,0,hh%inp,hh%cnf,hh%st,'',skip_ld,hh%psp,hh%at, &
&                 iout,ifail)
 hh%done=ifail==0
 return
 end subroutine run_all

 subroutine source(src,ae,ifail)

!! a copy of the all-electron atom that the handle src holds, its own
!! (oncv_ae_toml) or its run's
!!
!! @param src  the handle
!! @param ae  the all-electron atom
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_input, failed, msg_c_no_ae
 use, intrinsic :: iso_fortran_env, only: error_unit

!Input variables
 type(handle_t), intent(in) :: src

!Output variables
 type(ae_state_t), intent(out) :: ae
 integer, intent(out) :: ifail

 ifail=0
 if(src%has_stages) then
   ae=src%st%ae
 else if(src%has_ae) then
   ae=src%ae
 else
   write(error_unit,'(/a)') &
&       'source: ERROR the handle holds no all-electron atom'
   ifail=failed(fail_input,msg_c_no_ae)
 end if
 return
 end subroutine source

 function oncv_write_files(handle,format,nformat,path,npath) result(ifail) &
&         bind(c,name='oncv_write_files')

!! writes the pseudopotential of the handle's run as the file format
!! ("psp8", "upf" or "psml") to path: the lines the executables append to
!! the text output between its Begin line and END_PSP (write_psp_file),
!! without the blank line before END_PSP
!!
!! @param handle  the run handle
!! @param format  the format, nformat characters, not terminated
!! @param nformat  its length
!! @param path  the file name, npath characters, not terminated
!! @param npath  its length

 use psp_files_m, only: write_psp_file
 use failure_m, only: fail_input, fail_io, failed, msg_c_no_results, &
&                    msg_c_format, msg_c_file
 use, intrinsic :: iso_fortran_env, only: error_unit

!Input variables
 type(c_ptr), value :: handle
 character(kind=c_char), intent(in) :: format(nformat),path(npath)
 integer(c_size_t), value :: nformat,npath

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh
 character(len=:), allocatable :: kind,file,line,held
 integer :: iunit,iscr,ios,jfail
 logical :: inside,hold

 call c_f_pointer(handle,hh)
 if(.not. hh%done) then
   write(error_unit,'(/a)') 'oncv_write_files: ERROR the handle has no results'
   ifail=failed(fail_input,msg_c_no_results)
   return
 end if
 kind=from_c(format,nformat)
 if(kind/='psp8' .and. kind/='upf' .and. kind/='psml') then
   write(error_unit,'(/2a)') &
&       'oncv_write_files: ERROR the format is psp8, upf or psml, not ',kind
   ifail=failed(fail_input,msg_c_format)
   return
 end if
 file=from_c(path,npath)
 open(newunit=iunit,file=file,action='write',status='replace',iostat=ios)
 if(ios/=0) then
   write(error_unit,'(/2a)') 'oncv_write_files: ERROR cannot open ',file
   ifail=failed(fail_io,msg_c_file)
   return
 end if
 open(newunit=iscr,status='scratch',action='readwrite')
 call write_psp_file(kind,hh%inp,hh%psp,iscr,jfail)
 if(jfail==0) then
   rewind(iscr)
   inside=.false.
   hold=.false.
   do
     call read_line(iscr,line,ios)
     if(ios/=0) exit
     if(.not. inside) then
       inside=index(line,'Begin ')==1
       cycle
     end if
     if(line=='END_PSP') exit
     if(hold) write(iunit,'(a)') held
     held=line
     hold=.true.
   end do
   if(hold .and. len(held)>0) write(iunit,'(a)') held
 end if
 close(iscr)
 close(iunit)
 ifail=jfail
 return
 end function oncv_write_files

 subroutine read_line(iunit,line,ios)

!! the next line of iunit, whole, its trailing blanks kept
!!
!! @param iunit  the unit
!! @param line  the line
!! @param ios  0, or the iostat at the end of the file

!Input variables
 integer, intent(in) :: iunit

!Output variables
 character(len=:), allocatable, intent(out) :: line
 integer, intent(out) :: ios

!Local variables
 character(len=256) :: buf
 integer :: nread

 line=''
 do
   read(iunit,'(a)',advance='no',size=nread,iostat=ios) buf
   line=line//buf(1:nread)
   if(is_iostat_eor(ios)) then
     ios=0
     return
   end if
   if(ios/=0) return
 end do
 end subroutine read_line

 function oncv_failure(handle,failure) result(ifail) &
&         bind(c,name='oncv_failure')

!! the failure of the handle's last run (fail_context): its kind, message
!! and status, and where the run was; kind 0 if it did not fail
!!
!! @param handle  the run handle
!! @param failure  the record

 use failure_m, only: fail_context

!Input variables
 type(c_ptr), value :: handle

!Output variables
 type(failure_t), intent(out) :: failure
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh

 call c_f_pointer(handle,hh)
 failure=fail_context(hh%ifail,hh%at)
 ifail=0
 return
 end function oncv_failure

 function oncv_dims(handle,dims) result(ifail) bind(c,name='oncv_dims')

!! the sizes of the handle's results, for the caller's arrays; without a
!! finished run, those of its all-electron atom (irel, nkap, iexc, mmax)
!! if it holds one, all 0 but the version else
!!
!! @param handle  the run handle
!! @param dims  the sizes

!Input variables
 type(c_ptr), value :: handle

!Output variables
 type(dims_t), intent(out) :: dims
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh

 call c_f_pointer(handle,hh)
 dims=dims_t(api_version,0,0,0,0,0,0,0,0)
 ifail=0
 if(.not. hh%done) then
   if(hh%has_stages) then
     call ae_dims(hh%st%ae,dims)
   else if(hh%has_ae) then
     call ae_dims(hh%ae,dims)
   end if
   return
 end if
 associate(pot=>hh%psp%pot)
   dims%irel=pot%irel
   dims%nkap=pot%nkap
   dims%iexc=pot%iexc
   dims%mmax=pot%mmax
   dims%lmax=pot%lmax
   dims%lloc=pot%lloc
   dims%nproj(:)=pot%nproj(:)
 end associate
 dims%nrl=hh%psp%nrl
 return
 end function oncv_dims

 subroutine ae_dims(ae,dims)

!! the sizes of an all-electron atom
!!
!! @param ae  the all-electron atom
!! @param dims  the sizes, irel, nkap, iexc and mmax set

!Input variables
 type(ae_state_t), intent(in) :: ae

!Input/output variables
 type(dims_t), intent(inout) :: dims

 dims%irel=ae%irel
 dims%nkap=ae%nkap
 dims%iexc=ae%iexc
 dims%mmax=ae%mmax
 return
 end subroutine ae_dims

 function oncv_field_name(ifield,name,nname) result(nout) &
&         bind(c,name='oncv_field_name')

!! the field ifield (1, 2, ...) of the results that oncv_shape and
!! oncv_get_* copy out: its kind (r real, i integer), a blank and its name,
!! into name, nname bytes at most, not terminated; the number of bytes
!! written, 0 past the last field
!!
!! @param ifield  the field
!! @param name  its kind and name
!! @param nname  the size of name

!Input variables
 integer(c_int), value :: ifield
 integer(c_size_t), value :: nname

!Input/output variables
 character(kind=c_char), intent(inout) :: name(nname)

!Output variables
 integer(c_int) :: nout

!Local variables
 integer :: ii

 nout=0
 if(ifield<1 .or. ifield>size(fields)) return
 nout=int(min(int(len_trim(fields(ifield)),c_size_t),nname),c_int)
 do ii=1,nout
   name(ii)=fields(ifield)(ii:ii)
 end do
 return
 end function oncv_field_name

 function oncv_shape(handle,name,nname,rank,shp) result(ifail) &
&         bind(c,name='oncv_shape')

!! the rank and shape of a field of the handle's results (oncv_field_name),
!! in Fortran order (the first index fastest); a field the run did not
!! compute has the size 0
!!
!! @param handle  the run handle
!! @param name  the field's name, nname characters, not terminated
!! @param nname  its length
!! @param rank  its rank, 0 for a scalar
!! @param shp  its shape, shp(1:rank)

!Input variables
 type(c_ptr), value :: handle
 integer(c_size_t), value :: nname
 character(kind=c_char), intent(in) :: name(nname)

!Output variables
 integer(c_int), intent(out) :: rank
 integer(c_int), intent(out) :: shp(mxrank)
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh
 type(field_t) :: ff

 rank=0
 shp(:)=0
 call c_f_pointer(handle,hh)
 call lookup(hh,from_c(name,nname),ff,ifail)
 if(ifail/=0) return
 rank=ff%rank
 shp(:)=ff%shp(:)
 return
 end function oncv_shape

 function oncv_get_double(handle,name,nname,values,nvalues) result(ifail) &
&         bind(c,name='oncv_get_double')

!! the values of a real field of the handle's results, in Fortran order,
!! into values, whose size nvalues must be the field's (oncv_shape)
!!
!! @param handle  the run handle
!! @param name  the field's name, nname characters, not terminated
!! @param nname  its length
!! @param values  the values
!! @param nvalues  their number

!Input variables
 type(c_ptr), value :: handle
 integer(c_size_t), value :: nname,nvalues
 character(kind=c_char), intent(in) :: name(nname)

!Input/output variables
 real(c_double), intent(inout) :: values(nvalues)

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh
 type(field_t) :: ff

 call c_f_pointer(handle,hh)
 call lookup(hh,from_c(name,nname),ff,ifail)
 if(ifail/=0) return
 call check_size(from_c(name,nname),ff,.false.,nvalues,ifail)
 if(ifail/=0) return
 values(:)=ff%rv(:)
 return
 end function oncv_get_double

 function oncv_get_int(handle,name,nname,values,nvalues) result(ifail) &
&         bind(c,name='oncv_get_int')

!! the values of an integer field of the handle's results (a logical as 0
!! or 1), in Fortran order, into values, whose size nvalues must be the
!! field's (oncv_shape)
!!
!! @param handle  the run handle
!! @param name  the field's name, nname characters, not terminated
!! @param nname  its length
!! @param values  the values
!! @param nvalues  their number

!Input variables
 type(c_ptr), value :: handle
 integer(c_size_t), value :: nname,nvalues
 character(kind=c_char), intent(in) :: name(nname)

!Input/output variables
 integer(c_int), intent(inout) :: values(nvalues)

!Output variables
 integer(c_int) :: ifail

!Local variables
 type(handle_t), pointer :: hh
 type(field_t) :: ff

 call c_f_pointer(handle,hh)
 call lookup(hh,from_c(name,nname),ff,ifail)
 if(ifail/=0) return
 call check_size(from_c(name,nname),ff,.true.,nvalues,ifail)
 if(ifail/=0) return
 values(:)=ff%iv(:)
 return
 end function oncv_get_int

 subroutine check_size(name,ff,isint,nvalues,ifail)

!! whether a caller's buffer of nvalues values of the kind isint fits the
!! field ff
!!
!! @param name  the field's name
!! @param ff  the field
!! @param isint  whether the buffer holds integers
!! @param nvalues  its size
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_input, failed, msg_c_kind, msg_c_size
 use, intrinsic :: iso_fortran_env, only: error_unit

!Input variables
 character(len=*), intent(in) :: name
 type(field_t), intent(in) :: ff
 logical, intent(in) :: isint
 integer(c_size_t), intent(in) :: nvalues

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: nff

 ifail=0
 if(ff%isint.neqv.isint) then
   write(error_unit,'(/3a)') 'check_size: ERROR the field ',trim(name), &
&       ' is of the other kind (oncv_get_double, oncv_get_int)'
   ifail=failed(fail_input,msg_c_kind)
   return
 end if
 nff=product(ff%shp(1:ff%rank))
 if(nvalues/=int(nff,c_size_t)) then
   write(error_unit,'(/3a,i0,a,i0)') 'check_size: ERROR the field ', &
&       trim(name),' has ',nff,' values, the buffer ',nvalues
   ifail=failed(fail_input,msg_c_size)
   return
 end if
 return
 end subroutine check_size

 subroutine lookup(hh,name,ff,ifail)

!! the field name of the handle's results, copied (fields): of its
!! all-electron atom (ae.*: its own or its run's), its channel (channel.*),
!! its confined atom (confine.*) or its run (the others)
!!
!! @param hh  the handle
!! @param name  the field's name
!! @param ff  the field
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_input, failed, msg_c_field, msg_c_no_results
 use, intrinsic :: iso_fortran_env, only: error_unit

!Input variables
 type(handle_t), intent(in) :: hh
 character(len=*), intent(in) :: name

!Output variables
 type(field_t), intent(out) :: ff
 integer, intent(out) :: ifail

 ifail=0
 if(index(name,'ae.')==1) then
   if(hh%has_stages) then
     call lookup_ae(hh%st%ae,name,ff)
   else if(hh%has_ae) then
     call lookup_ae(hh%ae,name,ff)
   end if
 else if(index(name,'channel.')==1) then
   if(hh%has_chan) call lookup_channel(hh,name,ff)
 else if(index(name,'confine.')==1) then
   if(hh%has_confine) call lookup_confine(hh%sm,name,ff)
 else if(hh%done) then
   call lookup_run(hh,name,ff)
 end if
 if(ff%found) return
 if(any(fields(:)(3:)==name)) then
   write(error_unit,'(/2a)') 'lookup: ERROR the handle has no results '// &
&       'of the field ',name
   ifail=failed(fail_input,msg_c_no_results)
 else
   write(error_unit,'(/3a)') 'lookup: ERROR no field ',name, &
&       ' (oncv_field_name)'
   ifail=failed(fail_input,msg_c_field)
 end if
 return
 end subroutine lookup

 subroutine lookup_run(hh,name,ff)

!! the field name of the results of the handle's run, copied, not found if
!! it is none of theirs
!!
!! @param hh  the handle
!! @param name  the field's name
!! @param ff  the field

 use failure_m, only: mxwarn

!Input variables
 type(handle_t), intent(in) :: hh
 character(len=*), intent(in) :: name

!Output variables
 type(field_t), intent(out) :: ff

 associate(pot=>hh%psp%pot,tt=>hh%psp%tests,bb=>hh%psp%core%budget)
   select case(name)
   case('pot.zz')
     call put(ff,pot%zz)
   case('pot.zion')
     call put(ff,pot%zion)
   case('pot.c_light')
     call put(ff,pot%c_light)
   case('pot.rcmax')
     call put(ff,pot%rcmax)
   case('pot.rr')
     if(allocated(pot%rr)) then
       call put(ff,pot%rr)
     else
       call put_none(ff,1)
     end if
   case('pot.nproj')
     call put(ff,pot%nproj)
   case('pot.evkb')
     if(allocated(pot%evkb)) then
       call put(ff,pot%evkb)
     else
       call put_none(ff,3)
     end if
   case('pot.vkb')
     if(allocated(pot%vkb)) then
       call put(ff,pot%vkb)
     else
       call put_none(ff,4)
     end if
   case('pot.esr')
     if(allocated(pot%esr)) then
       call put(ff,pot%esr)
     else
       call put_none(ff,2)
     end if
   case('pot.eso')
     if(allocated(pot%eso)) then
       call put(ff,pot%eso)
     else
       call put_none(ff,2)
     end if
   case('pot.vpuns')
     if(allocated(pot%vpuns)) then
       call put(ff,pot%vpuns)
     else
       call put_none(ff,2)
     end if
   case('pot.rho')
     if(allocated(pot%rho)) then
       call put(ff,pot%rho)
     else
       call put_none(ff,1)
     end if
   case('pot.rhomod')
     if(allocated(pot%rhomod)) then
       call put(ff,pot%rhomod)
     else
       call put_none(ff,2)
     end if
   case('epstot')
     call put(ff,hh%psp%epstot)
   case('ea')
     call put(ff,hh%psp%ea)
   case('rhotae')
     if(allocated(hh%psp%rhotae)) then
       call put(ff,hh%psp%rhotae)
     else
       call put_none(ff,1)
     end if
   case('rhoc')
     if(allocated(hh%psp%rhoc)) then
       call put(ff,hh%psp%rhoc)
     else
       call put_none(ff,1)
     end if
   case('rcore')
     call put(ff,hh%psp%rcore)
   case('vpsml')
     if(allocated(hh%psp%vpsml)) then
       call put(ff,hh%psp%vpsml)
     else
       call put_none(ff,3)
     end if
   case('uupsa')
     if(allocated(hh%psp%uupsa)) then
       call put(ff,hh%psp%uupsa)
     else
       call put_none(ff,3)
     end if
   case('tests.logder.ne')
     call put(ff,tt%logder%ne)
   case('tests.logder.ee')
     if(allocated(tt%logder%ee)) then
       call put(ff,tt%logder%ee)
     else
       call put_none(ff,1)
     end if
   case('tests.logder.rl')
     call put(ff,tt%logder%rl)
   case('tests.logder.ae')
     if(allocated(tt%logder%ae)) then
       call put(ff,tt%logder%ae)
     else
       call put_none(ff,3)
     end if
   case('tests.logder.ps')
     if(allocated(tt%logder%ps)) then
       call put(ff,tt%logder%ps)
     else
       call put_none(ff,3)
     end if
   case('tests.ghosts.nrow')
     call put(ff,tt%ghosts%nrow)
   case('tests.ghosts.nn')
     call put(ff,tt%ghosts%nn)
   case('tests.ghosts.ee')
     call put(ff,tt%ghosts%ee)
   case('tests.ghosts.least')
     call put(ff,tt%ghosts%least)
   case('tests.configs.ncfg')
     call put(ff,tt%configs%ncfg)
   case('tests.configs.nrow')
     call put(ff,tt%configs%nrow)
   case('tests.configs.state')
     call put(ff,tt%configs%state)
   case('tests.configs.ee')
     call put(ff,tt%configs%ee)
   case('tests.configs.de')
     call put(ff,tt%configs%de)
   case('tests.convergence.cross')
     if(allocated(tt%convergence%cross)) then
       call put(ff,tt%convergence%cross)
     else
       call put_none(ff,5)
     end if
   case('tests.convergence.curve')
     if(allocated(tt%convergence%curve)) then
       call put(ff,tt%convergence%curve)
     else
       call put_none(ff,5)
     end if
   case('tests.convergence.tail')
     if(allocated(tt%convergence%tail)) then
       call put(ff,tt%convergence%tail)
     else
       call put_none(ff,4)
     end if
   case('tests.convergence.opt')
     if(allocated(tt%convergence%opt)) then
       call put(ff,tt%convergence%opt)
     else
       call put_none(ff,4)
     end if
   case('tests.convergence.wshare')
     if(allocated(tt%convergence%wshare)) then
       call put(ff,tt%convergence%wshare)
     else
       call put_none(ff,4)
     end if
   case('tests.weighted_cut.used')
     call put(ff,tt%weighted_cut%used)
   case('tests.weighted_cut.errw')
     call put(ff,tt%weighted_cut%errw)
   case('tests.weighted_cut.ecutw')
     call put(ff,tt%weighted_cut%ecutw)
   case('tests.weighted_cut.weight')
     if(allocated(tt%weighted_cut%weight)) then
       call put(ff,tt%weighted_cut%weight)
     else
       call put_none(ff,3)
     end if
   case('tests.vloc_cut.cut')
     call put(ff,tt%vloc_cut%cut)
   case('tests.vloc_cut.check')
     call put(ff,tt%vloc_cut%check)
   case('tests.bbeig')
     if(allocated(tt%bbeig)) then
       call put(ff,tt%bbeig)
     else
       call put_none(ff,3)
     end if
   case('tests.ghosts.nminus')
     call put(ff,tt%ghosts%nminus)
   case('tests.ghosts.minus')
     call put(ff,tt%ghosts%minus)
   case('tests.ghosts.nplus')
     call put(ff,tt%ghosts%nplus)
   case('tests.ghosts.plus')
     call put(ff,tt%ghosts%plus)
   case('tests.ghosts.ecut')
     call put(ff,tt%ghosts%ecut)
   case('tests.ghosts.basis')
     call put(ff,tt%ghosts%basis)
   case('core.teter')
     call put(ff,hh%psp%core%teter)
   case('core.budget.ecut')
     call put(ff,bb%ecut)
   case('core.budget.error')
     call put(ff,bb%error)
   case('core.budget.slack')
     call put(ff,bb%slack)
   case('core.budget.ecutw')
     call put(ff,bb%ecutw)
   case('core.budget.factor')
     call put(ff,bb%factor)
   case('core.budget.r_s')
     call put(ff,bb%r_s)
   case('core.budget.rc_min')
     call put(ff,bb%rc_min)
   case('core.budget.r_core')
     call put(ff,bb%r_core)
   case('core.budget.fill')
     call put(ff,bb%fill)
   case('core.budget.cpoly')
     call put(ff,bb%cpoly)
   case('core.budget.martins')
     call put(ff,bb%martins)
   case('core.budget.ecut_needed')
     call put(ff,bb%ecut_needed)
   case('core.budget.r_s6')
     call put(ff,bb%r_s6)
   case('core.budget.mismatch')
     call put(ff,bb%mismatch)
   case('core.budget.crad')
     if(allocated(bb%crad)) then
       call put(ff,bb%crad)
     else
       call put_none(ff,1)
     end if
   case('core.budget.cfill')
     if(allocated(bb%cfill)) then
       call put(ff,bb%cfill)
     else
       call put(ff,[integer ::])
     end if
   case('warnings.nwarn')
     call put(ff,hh%psp%warnings%nwarn)
   case('warnings.rec')
     call put(ff,hh%psp%warnings%rec(:,1:min(hh%psp%warnings%nwarn,mxwarn)))
   case('core.budget.cmis')
     if(allocated(bb%cmis)) then
       call put(ff,bb%cmis)
     else
       call put_none(ff,1)
     end if
   case default
     continue
   end select
 end associate
 return
 end subroutine lookup_run

 subroutine lookup_ae(ae,name,ff)

!! the field name of an all-electron atom (ae.*), copied, not found if it
!! is none of its
!!
!! @param ae  the all-electron atom
!! @param name  the field's name
!! @param ff  the field

 use failure_m, only: mxwarn

!Input variables
 type(ae_state_t), intent(in) :: ae
 character(len=*), intent(in) :: name

!Output variables
 type(field_t), intent(out) :: ff

!Local variables
 integer :: ncv

 ncv=ae%nc+ae%nv
 select case(name)
 case('ae.irel')
   call put(ff,ae%irel)
 case('ae.nkap')
   call put(ff,ae%nkap)
 case('ae.zz')
   call put(ff,ae%zz)
 case('ae.nc')
   call put(ff,ae%nc)
 case('ae.nv')
   call put(ff,ae%nv)
 case('ae.iexc')
   call put(ff,ae%iexc)
 case('ae.na')
   call put(ff,ae%na(1:ncv))
 case('ae.la')
   call put(ff,ae%la(1:ncv))
 case('ae.fa')
   call put(ff,ae%fa(1:ncv))
 case('ae.c_light')
   call put(ff,ae%c_light)
 case('ae.mmax')
   call put(ff,ae%mmax)
 case('ae.al')
   call put(ff,ae%al)
 case('ae.rr')
   call put(ff,ae%rr)
 case('ae.it')
   call put(ff,ae%it)
 case('ae.etot')
   call put(ff,ae%etot)
 case('ae.eterm')
   call put(ff,ae%eterm)
 case('ae.res')
   call put(ff,ae%res)
 case('ae.ea')
   call put(ff,ae%ea(1:ncv,1:ae%nkap))
 case('ae.rhoc')
   call put(ff,ae%rhoc)
 case('ae.rho')
   call put(ff,ae%rho)
 case('ae.vfull')
   call put(ff,ae%vfull)
 case('ae.warnings.nwarn')
   call put(ff,ae%warn%nwarn)
 case('ae.warnings.rec')
   call put(ff,ae%warn%rec(:,1:min(ae%warn%nwarn,mxwarn)))
 case default
   continue
 end select
 return
 end subroutine lookup_ae

 subroutine lookup_channel(hh,name,ff)

!! the field name of a channel (channel.*), copied, not found if it is none
!! of its; the arrays per projector and kappa as the run's (mxprj, nkap),
!! and the convergence curves to the end of the longest; the local
!! potential and projectors of size 0 without them, and the projectors of
!! the local channel too
!!
!! @param hh  the handle, which holds a channel
!! @param name  the field's name
!! @param ff  the field

 use failure_m, only: mxwarn

!Input variables
 type(handle_t), target, intent(in) :: hh
 character(len=*), intent(in) :: name

!Output variables
 type(field_t), intent(out) :: ff

!Local variables
 type(channel_t), pointer :: ch
 integer :: ii,nkap
 logical :: has_prj

 ch=>hh%chan
 nkap=size(ch%vp,2)
 has_prj=hh%has_vkb .and. ch%l1/=hh%inp%lloc+1
 select case(name)
 case('channel.l')
   call put(ff,ch%l1-1)
 case('channel.mkap')
   call put(ff,ch%mkap)
 case('channel.nproj')
   call put(ff,ch%nproj)
 case('channel.irc')
   call put(ff,ch%irc)
 case('channel.npa')
   call put(ff,ch%npa)
 case('channel.wellprj')
   call put(ff,ch%wellprj)
 case('channel.qmsbf')
   call put(ff,ch%qmsbf)
 case('channel.epa')
   call put(ff,ch%epa(:,1:size(ch%vp,2)))
 case('channel.pswf')
   call put(ff,ch%pswf)
 case('channel.vp')
   call put(ff,ch%vp)
 case('channel.cvgplt')
   call put(ff,ch%cvgplt)
 case('channel.ecurve')
   nkap=0
   do ii=1,size(ch%ecurve,2)
     if(any(ch%ecurve(1,ii,:,:)>0.0d0)) nkap=ii
   end do
   call put(ff,ch%ecurve(:,1:nkap,:,:))
 case('channel.ektail')
   call put(ff,ch%ektail)
 case('channel.optdg')
   call put(ff,ch%optdg)
 case('channel.wshare')
   call put(ff,ch%wshare)
 case('channel.wcurve')
   call put(ff,ch%wcurve)
 case('channel.warnings.nwarn')
   call put(ff,ch%warn%nwarn)
 case('channel.warnings.rec')
   call put(ff,ch%warn%rec(:,1:min(ch%warn%nwarn,mxwarn)))
 case('channel.vloc')
   if(hh%has_vkb) then
     call put(ff,hh%lp%vloc)
   else
     call put_none(ff,1)
   end if
 case('channel.evkb')
   if(has_prj) then
     call put(ff,hh%prj%evkb(:,1:nkap))
   else
     call put_none(ff,2)
   end if
 case('channel.vkb')
   if(has_prj) then
     call put(ff,hh%prj%vkb)
   else
     call put_none(ff,3)
   end if
 case('channel.bbeig')
   if(has_prj) then
     call put(ff,hh%prj%bbeig(:,1:nkap))
   else
     call put_none(ff,2)
   end if
 case default
   continue
 end select
 return
 end subroutine lookup_channel

 subroutine lookup_confine(sm,name,ff)

!! the field name of the summary of a confined atom (confine.*), copied,
!! not found if it is none of its; -1 for a part left out
!!
!! @param sm  the summary
!! @param name  the field's name
!! @param ff  the field

!Input variables
 type(confine_summary_t), intent(in) :: sm
 character(len=*), intent(in) :: name

!Output variables
 type(field_t), intent(out) :: ff

 select case(name)
 case('confine.nmiss')
   call put(ff,sm%nmiss)
 case('confine.nextra')
   call put(ff,sm%nextra)
 case('confine.noff')
   call put(ff,sm%noff)
 case('confine.nwfail')
   call put(ff,sm%nwfail)
 case('confine.max_de')
   call put(ff,sm%smax_de)
 case('confine.cut')
   call put(ff,sm%scut)
 case('confine.dde')
   call put(ff,sm%sdde)
 case('confine.inf')
   call put(ff,sm%sinf)
 case default
   continue
 end select
 return
 end subroutine lookup_confine

 subroutine put_none(ff,rank)

!! a field the run did not compute: rank rank, size 0
!!
!! @param ff  the field
!! @param rank  its rank

!Input variables
 integer, intent(in) :: rank

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%rank=rank
 allocate(ff%rv(0))
 return
 end subroutine put_none

 subroutine put_r0(ff,xx)

!! a real field of rank 0
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 real(dp), intent(in) :: xx

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%rank=0
 ff%rv=[xx]
 return
 end subroutine put_r0

 subroutine put_r1(ff,xx)

!! a real field of rank 1
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 real(dp), intent(in) :: xx(:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%rank=1
 ff%shp(1:1)=shape(xx)
 ff%rv=reshape(xx,[size(xx)])
 return
 end subroutine put_r1

 subroutine put_r2(ff,xx)

!! a real field of rank 2
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 real(dp), intent(in) :: xx(:,:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%rank=2
 ff%shp(1:2)=shape(xx)
 ff%rv=reshape(xx,[size(xx)])
 return
 end subroutine put_r2

 subroutine put_r3(ff,xx)

!! a real field of rank 3
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 real(dp), intent(in) :: xx(:,:,:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%rank=3
 ff%shp(1:3)=shape(xx)
 ff%rv=reshape(xx,[size(xx)])
 return
 end subroutine put_r3

 subroutine put_r4(ff,xx)

!! a real field of rank 4
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 real(dp), intent(in) :: xx(:,:,:,:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%rank=4
 ff%shp(1:4)=shape(xx)
 ff%rv=reshape(xx,[size(xx)])
 return
 end subroutine put_r4

 subroutine put_r5(ff,xx)

!! a real field of rank 5
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 real(dp), intent(in) :: xx(:,:,:,:,:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%rank=5
 ff%shp(1:5)=shape(xx)
 ff%rv=reshape(xx,[size(xx)])
 return
 end subroutine put_r5

 subroutine put_i0(ff,xx)

!! an integer field of rank 0
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 integer, intent(in) :: xx

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%isint=.true.
 ff%rank=0
 ff%iv=[xx]
 return
 end subroutine put_i0

 subroutine put_i1(ff,xx)

!! an integer field of rank 1
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 integer, intent(in) :: xx(:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%isint=.true.
 ff%rank=1
 ff%shp(1:1)=shape(xx)
 ff%iv=reshape(xx,[size(xx)])
 return
 end subroutine put_i1

 subroutine put_i2(ff,xx)

!! an integer field of rank 2
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 integer, intent(in) :: xx(:,:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%isint=.true.
 ff%rank=2
 ff%shp(1:2)=shape(xx)
 ff%iv=reshape(xx,[size(xx)])
 return
 end subroutine put_i2

 subroutine put_i3(ff,xx)

!! an integer field of rank 3
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 integer, intent(in) :: xx(:,:,:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%isint=.true.
 ff%rank=3
 ff%shp(1:3)=shape(xx)
 ff%iv=reshape(xx,[size(xx)])
 return
 end subroutine put_i3

 subroutine put_l0(ff,xx)

!! a logical field, as the integer 0 or 1
!!
!! @param ff  the field
!! @param xx  its value

!Input variables
 logical, intent(in) :: xx

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%isint=.true.
 ff%iv=[merge(1,0,xx)]
 return
 end subroutine put_l0

 subroutine put_l1(ff,xx)

!! a logical field of rank 1, as the integers 0 or 1
!!
!! @param ff  the field
!! @param xx  its values

!Input variables
 logical, intent(in) :: xx(:)

!Output variables
 type(field_t), intent(out) :: ff

 ff%found=.true.
 ff%isint=.true.
 ff%rank=1
 ff%shp(1)=size(xx)
 ff%iv=merge(1,0,xx)
 return
 end subroutine put_l1

 function from_c(chars,nchars) result(string)

!! a Fortran string of nchars C characters
!!
!! @param chars  the characters
!! @param nchars  their number

!Input variables
 character(kind=c_char), intent(in) :: chars(nchars)
 integer(c_size_t), intent(in) :: nchars

!Output variables
 character(len=:), allocatable :: string

!Local variables
 integer(c_size_t) :: ii

 allocate(character(len=nchars) :: string)
 do ii=1,nchars
   string(ii:ii)=chars(ii)
 end do
 return
 end function from_c

 end module oncvpsp_c_m

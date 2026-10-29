!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
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
!
!> The calculation of a pseudopotential, from its input to its tests and
!> the results (results_m)
 module oncvpsp_run_m
 use constants_m, only: dp, mxlin
 use ae_state_m, only: ae_state_t
 use pieces_m, only: piecewise_t
 use failure_m, only: warnings_t
 implicit none
 private
 public :: oncvpsp_run, complete_input, run_stages, run_coupled, stages_t, &
&          run_ae, run_channel

!> the results of the stages before the model core (oncvpsp_run), from
!> which run_coupled makes the model core, the screening and the tests;
!> run_coupled does not change them, so it can be called again with
!> another model core
 type :: stages_t
!> the all-electron atom, with the test configurations' all-electron
!> halves read with its file
   type(ae_state_t) :: ae
!> number of kappas per l: 1, or 2 for irel=2
   integer :: nkap=0
!> ion charge
   real(dp) :: zion=0.0d0
!> largest core radius (bohr)
   real(dp) :: rcmax=0.0d0
!> points of the output mesh
   integer :: nrl=0
!> the q of each channel's last basis function
   real(dp) :: qmsbf(mxlin+1)=0.0d0
!> the core radii as the echoes give them
   real(dp) :: rc_in(mxlin+1)=0.0d0
!> fcfact as the echoes give it
   real(dp) :: fcfact_in=0.0d0
!> the principal quantum number of each projector's state
   integer, allocatable :: npa(:,:)
!> the projectors' energies (Ha)
   real(dp), allocatable :: epa(:,:,:)
!> whether each projector is a well state
   logical, allocatable :: wellprj(:,:)
!> the screened semi-local and local potentials (Ha)
   real(dp), allocatable :: vp(:,:,:)
!> the projectors on the mesh
   real(dp), allocatable :: vkb(:,:,:,:)
!> their coefficients (Ha)
   real(dp), allocatable :: evkb(:,:,:)
!> the projectors as pieces
   type(piecewise_t), allocatable :: vkbpw(:,:,:)
!> the screened local potential as pieces
   type(piecewise_t) :: vlocpw
!> the cutoffs of the error levels (run_optimize)
   real(dp), allocatable :: cvgplt(:,:,:,:,:)
!> the convergence curves (run_optimize)
   real(dp), allocatable :: ecurve(:,:,:,:,:)
!> the curves' tails and checks (run_optimize)
   real(dp), allocatable :: ektail(:,:,:,:)
!> the optimization's diagnostics (run_optimize)
   real(dp), allocatable :: optdg(:,:,:,:)
!> the waves' shares (run_optimize)
   real(dp), allocatable :: wshare(:,:,:,:)
!> the weights of the pseudo-waves in the weighted cutoff
   real(dp), allocatable :: wcurve(:,:,:)
!> the eigenvalues of the B matrices (run_vkb)
   real(dp), allocatable :: bbeig(:,:,:)
!> the pseudo valence charge density (valence_densities)
   real(dp), allocatable :: rho(:)
!> the all-electron valence charges, state by state
   real(dp), allocatable :: rhoae(:,:)
!> the pseudo valence charges, state by state
   real(dp), allocatable :: rhops(:,:)
!> the all-electron valence charge density
   real(dp), allocatable :: rhotae(:)
!> the pseudo-atomic orbitals
   real(dp), allocatable :: uupsa(:,:,:)
!> the pseudo valence eigenvalue sum (Ha)
   real(dp) :: eeig=0.0d0
!> valence charge
   real(dp) :: zval=0.0d0
!> fully relativistic: the scalar-relativistic projectors' coefficients
   real(dp), allocatable :: esr(:,:)
!> fully relativistic: the spin-orbit projectors' coefficients
   real(dp), allocatable :: eso(:,:)
!> fully relativistic: the scalar-relativistic projectors as pieces
   type(piecewise_t), allocatable :: vsrpw(:,:)
!> fully relativistic: the spin-orbit projectors as pieces
   type(piecewise_t), allocatable :: vsopw(:,:)
!> the warnings of the stages
   type(warnings_t) :: warn
 end type stages_t

 contains

 subroutine oncvpsp_run(irel,istage,inp,ae_filename,ps_filename,skip_ld, &
&                       psp,at,iout,ifail)
!! the pseudopotential of the input inp, and its tests: the stages before
!! the model core (run_stages), then the model core and the tests
!! (run_coupled). The text output goes to iout, and the pseudopotential
!! and the tests to psp. inp is completed (configurations, core radii) as
!! the writers echo it. With istage=1 only the all-electron atom, written
!! to ae_filename; with istage=2 that file in place of it. A ps_filename,
!! if not blank, receives the pseudopotential for oncvpsp_ld.x; skip_ld
!! leaves out the log derivatives. On a failure, at says where the run
!! was (fail_context).
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param istage  0 the whole calculation, 1 the all-electron stage only, 2 the
!!                pseudopotential stage only
!! @param inp  the input values (input_m), as read
!! @param ae_filename  the all-electron data file of istage 1 and 2
!! @param ps_filename  the pseudopotential file for oncvpsp_ld.x, or blank
!! @param skip_ld  whether to leave out the log derivatives
!! @param psp  the results (results_m)
!! @param at  where the run was (failure_m): its stage, l and ikap
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use input_m, only: input_t, config_t
 use results_m, only: results_t
 use constants_m, only: mxcnf
 use failure_m, only: failure_t

!Input variables
 integer, intent(in) :: irel,istage
 character(len=*), intent(in) :: ae_filename,ps_filename
 logical, intent(in) :: skip_ld
 integer, intent(in) :: iout

!Input/output variables
 type(input_t), intent(inout) :: inp

!Output variables
 type(results_t), intent(out) :: psp
 type(failure_t), intent(out) :: at
 integer, intent(out) :: ifail

!Local variables
 type(config_t) :: cnf(mxcnf)
 type(stages_t) :: st

 call run_stages(irel,istage,inp,ae_filename,cnf,st,at,iout,ifail)
 if(ifail/=0 .or. istage==1) return
 call run_coupled(irel,istage,inp,cnf,st,ps_filename,skip_ld,psp,at,iout, &
&                 ifail)
 return
 end subroutine oncvpsp_run

 subroutine run_stages(irel,istage,inp,ae_filename,cnf,st,at,iout,ifail, &
&                      ae_in)

!! the stages of a run before the model core, from the input inp: the
!! input completed (complete_input) and echoed, the all-electron atom
!! (solve_ae), the channels (optimize_channel), the local potential and
!! the projectors (run_vkb), and the valence densities, in st for
!! run_coupled. With istage=1 only the all-electron atom, written to
!! ae_filename, and st is left incomplete; with istage=2 that file in place
!! of it; with ae_in (istage 0), that atom in place of it (given_ae)
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param istage  0 the whole calculation, 1 the all-electron stage only, 2 the
!!                pseudopotential stage only
!! @param inp  the input values (input_m), as read on input, completed on
!!             output
!! @param ae_filename  the all-electron data file of istage 1 and 2
!! @param cnf  the reference and test configurations (input_m)
!! @param st  the stages' results
!! @param at  where the run was (failure_m): its stage, l and ikap
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
!! @param ae_in  the all-electron atom solved before, if given (solve_ae)
 use failure_m, only: failure_t, fail_none, stage_input, &
&                    stage_channel, stage_vkb, stage_valence, stage_configs, &
&                    warn_merge
 use config_ae_m, only: run_config_ae, ae_config_t
 use ae_file_m, only: write_ae_file
 use ae_state_m, only: ae_state_t
 use channel_m, only: channel_t, optimize_channel, ncurve
 use sr_so_r_m, only: sr_so_r
 use run_vkb_m, only: run_vkb
 use ps_valence_m, only: valence_densities
 use constants_m, only: nerrdg
 use input_echo_m, only: write_atom, write_lmax, write_icmod, write_epsh, &
&                        write_mesh, fx, vkb_head, vkb_line
 use main_setup_m, only: main_names, write_reference_results, core_radius_points
 use constants_m, only: mxprj, mxcnf
 use input_m, only: input_t, config_t
 use constants_m, only: dp, mxl, mxlin, mxstate, mxkap
 use pieces_m, only: piecewise_t
 use adams_m, only: adams_t
!Input variables
 integer, intent(in) :: irel,istage
 character(len=*), intent(in) :: ae_filename
 integer, intent(in) :: iout
 type(ae_state_t), intent(in), optional :: ae_in

!Input/output variables
 type(input_t), intent(inout) :: inp

!Output variables
 type(config_t), intent(out) :: cnf(mxcnf)
 type(stages_t), intent(out) :: st
 type(failure_t), intent(out) :: at
 integer, intent(out) :: ifail

!Local variables
!the warnings, those of solve_ae and each channel's and their own
 type(warnings_t) :: warns
 integer :: ii,ierr,it
 integer :: jj,l1
 integer :: mmax,nlim,nrl
 integer :: nkap
 integer,allocatable :: npa(:,:)
!
 integer :: irc(mxlin+1)

 real(dp) :: al
 real(dp) :: eeig
 real(dp) :: fcfact_in
 real(dp) :: fcgrid(3),rcgrid(3)
 real(dp) :: rcmax
 real(dp) :: log_mesh(4)
 type(adams_t) :: adams
 real(dp) :: zion,zval,etot,res
!> the speed of light of the run, from the input
 real(dp) :: c_light
 real(dp) :: eterm(4)
!
 real(dp) :: ea(mxstate,mxkap)
 real(dp) :: qmsbf(mxlin+1),rc_in(mxlin+1)

 real(dp), allocatable :: evkb(:,:,:),cvgplt(:,:,:,:,:)
 real(dp), allocatable :: ecurve(:,:,:,:,:)
!the r-space less the Fourier kinetic energy of each pseudo-wave
!(run_optimize)
 real(dp), allocatable :: ektail(:,:,:,:),optdg(:,:,:,:),wshare(:,:,:,:)
!the eigenvalues of the B matrices (run_vkb)
 real(dp), allocatable :: bbeig(:,:,:)
 real(dp) :: ecutw
!weights of the pseudo-waves' curves in the weighted cutoff (budget_cutoff)
 real(dp), allocatable :: wcurve(:,:,:)
 real(dp), allocatable :: rr(:)
 real(dp), allocatable :: rho(:),rhoc(:)
 real(dp), allocatable :: uu(:,:),up(:,:)
 logical, allocatable :: wellprj(:,:)
 real(dp), allocatable :: vp(:,:,:),vfull(:),vkb(:,:,:,:),pswf(:,:,:,:)
!the valence charges, state by state, and the all-electron density
 real(dp), allocatable :: rhoae(:,:),rhops(:,:),rhotae(:)
 real(dp), allocatable :: uupsa(:,:,:) !pseudo-atomic orbitals array
 real(dp), allocatable :: epa(:,:,:)
 real(dp), allocatable :: dlta(:,:,:,:),vrloc(:,:),pswfx(:,:,:,:)
 type(piecewise_t), allocatable :: vkbpw(:,:,:)
 type(piecewise_t) :: vlocpw
 real(dp), allocatable :: chix(:,:,:,:)
!fully-relativistic only: scalar-relativistic and spin-orbit projectors
 real(dp), allocatable :: vsr(:,:,:),vso(:,:,:),esr(:,:),eso(:,:)
 type(piecewise_t), allocatable :: vsrpw(:,:),vsopw(:,:)

!all-electron results for the test configurations (see run_config_ae)
 type(ae_config_t) :: aec(mxcnf)
!the all-electron atom (ae_state_m) and a channel (channel_m)
 type(ae_state_t) :: ae
 type(channel_t) :: chan

!names in the messages: the equation solved and the program
 character(len=19) :: treatment
 character(len=16) :: pname
 character(len=16) :: exe

 ifail=0
 at=failure_t(fail_none,0,0,stage_input,-1,0,0)
 call main_names(irel,istage,nkap,treatment,pname,exe)
 call run_ae(irel,istage,inp,ae_filename,cnf,ae,at,iout,ifail,ae_in)
 if(ifail/=0) return

!icmod=4 coarse-search grids: min, max, step
 fcgrid=[inp%fcfact_min,inp%fcfact_max,inp%fcfact_step]
 rcgrid=[inp%rcfact_min,inp%rcfact_max,inp%rcfact_step]
 warns=ae%warn
 c_light=ae%c_light
 al=ae%al
 mmax=ae%mmax
 rr=ae%rr
 log_mesh(:)=ae%log_mesh(:)
 adams=ae%adams

!PWSCF wants an even number of mesh points (not made so when fully
!relativistic, see doc/fr_merge_plan.md): one more, or one fewer where
!one more would pass the log mesh
 nrl=int((inp%rlmax/inp%drl)-0.5d0)+1
 if(mod(nrl,2)/=0 .and. irel/=2) then
   if(inp%drl*dble(nrl)<=rr(mmax)) then
     nrl=nrl+1
   else
     nrl=nrl-1
   end if
 end if

!calculate zion for output
 zion=inp%zz
 do ii=1,inp%nc
   zion=zion-inp%fa(ii)
 end do


 allocate(rho(mmax),rhoc(mmax))
 allocate(uu(mmax,2),up(mmax,2),uupsa(mmax,nkap,inp%nv))
 allocate(evkb(mxprj,mxl+1,nkap))
!the projectors beyond nproj and of lloc stay 0
 evkb(:,:,:)=0.0d0
 allocate(cvgplt(2,nerrdg,mxprj,mxl+1,nkap))
 allocate(ecurve(2,ncurve,mxprj,mxl+1,nkap))
 allocate(ektail(2,mxprj,mxl+1,nkap),bbeig(mxprj,mxl+1,nkap))
 allocate(optdg(4,mxprj,mxl+1,nkap),wshare(3,mxprj,mxl+1,nkap))
 allocate(wcurve(mxprj,mxl+1,nkap))
 allocate(vp(mmax,mxl+2,nkap),vfull(mmax),vkb(mmax,mxprj,mxl+1,nkap))
 allocate(pswf(mmax,mxprj,mxl+1,nkap))
 allocate(rhoae(mmax,inp%nv),rhops(mmax,inp%nv),rhotae(mmax))
 allocate(npa(mxprj,mxlin+1))
 allocate(wellprj(mxprj,mxlin+1))
 wellprj(:,:)=.false.
 allocate(epa(mxprj,mxlin+1,nkap))
 allocate(vsr(mmax,2*mxprj,mxl+1),vso(mmax,2*mxprj,mxl+1))
 allocate(esr(2*mxprj,mxl+1),eso(2*mxprj,mxl+1))
!the results hold them for every irel (sr_so_r sets them for irel=2)
 esr(:,:)=0.0d0
 eso(:,:)=0.0d0
 allocate(dlta(mmax,mxprj,mxlin+1,nkap),vrloc(mmax,nkap))
 allocate(vkbpw(mxprj,mxl+1,nkap))
 allocate(vsrpw(2*mxprj,mxl+1),vsopw(2*mxprj,mxl+1))
 allocate(pswfx(mmax,mxprj,mxl+1,nkap),chix(mmax,mxprj,mxl+1,nkap))

 dlta(:,:,:,:)=0.0d0
 vrloc(:,:)=0.0d0
 pswfx(:,:,:,:)=0.0d0
 chix(:,:,:,:)=0.0d0
 vp(:,:,:)=0.0d0
 vkb(:,:,:,:)=0.0d0
 epa(:,:,:)=0.0d0

 rhoc(:)=ae%rhoc(:)
 rho(:)=ae%rho(:)
 vfull(:)=ae%vfull(:)
 ea(:,:)=ae%ea(:,:)
 it=ae%it
 etot=ae%etot
 eterm(:)=ae%eterm(:)
 res=ae%res
 aec(:)=ae%aec(:)
!
!

 rc_in(:)=inp%rc(:)
!input value, for the echoes (modcore3 resets fcfact for icmod=4)
 fcfact_in=inp%fcfact

! output printing (echos input data, with all-electron eigenvalues added)

 call write_atom(inp%atsym,inp%zz,inp%nc,inp%nv,inp%iexc,c_light,inp%psfile, &
&                iout)
 write(iout,'(a/a)') '#','#   n    l    f        energy (Ha)'
 do ii=1,inp%nc+inp%nv
   if(irel/=2) then
     write(iout,'(2i5,f8.2,1pe18.7)') inp%na(ii),inp%la(ii),inp%fa(ii),ea(ii,1)
   else if(inp%la(ii)==0) then
     write(iout,'(2i5,f8.2,1p,e18.7)') inp%na(ii),inp%la(ii),inp%fa(ii),ea(ii,1)
   else
     write(iout,'(2i5,f8.2,1p,2d18.7)') inp%na(ii),inp%la(ii),inp%fa(ii), &
&          ea(ii, &
&          1),ea(ii,2)
   end if
 end do

 call write_lmax(inp%lmax,inp%ep_as_given,inp%well_onset,inp%well_tail, &
&                inp%wave_regularization,inp%well_norm,iout)
 write(iout,'(a/a)') '#','#   l,   rc,      ep,       ncon, nbas, qcut'
 do l1=1,inp%lmax+1
   write(iout,'(i5,2a,2i5,a)') l1-1,fx(inp%rc(l1),10,5),fx(inp%ep(l1,1),10,5), &
&        inp%ncon(l1),inp%nbas(l1),fx(inp%qcut(l1),10,5)
 end do

 write(iout,'(a/a/a,a)') '#','# LOCAL POTENTIAL','# lloc, lpopt,  rc(5),', &
&      '   dvloc0'
 write(iout,'(2i5,3a)') inp%lloc,inp%lpopt,fx(inp%rc(5),10,5),'   ', &
&      fx(inp%dvloc0,10,5)

 write(iout,'(a/a/a)') '#','# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', &
&      trim(vkb_head(inp%lmax,inp%well_n,inp%ecut_weight))
 do l1=1,inp%lmax+1
   write(iout,'(a)') trim(vkb_line(inp%lmax,l1,inp%nproj,inp%debl,inp%well_n, &
&        inp%ecut_weight))
 end do

 call write_icmod(inp%icmod,inp%fcfact,inp%rcfact,inp%slack,inp%budget_cutoff, &
&                 fcgrid,rcgrid,iout)

 call write_epsh(inp%epsh1,inp%epsh2,inp%depsh,inp%rxpsh,iout)

 if(inp%upf_mesh==0) then
   write(iout,'(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl'
   write(iout,'(2a)') fx(inp%rlmax,8,4),fx(inp%drl,8,4)
 else
   write(iout,'(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl, upf_mesh'
   write(iout,'(2a,i8)') fx(inp%rlmax,8,4),fx(inp%drl,8,4),inp%upf_mesh
 end if

 write(iout,'(a/a/a)') '#','# TEST CONFIGURATIONS','# ncnf'
 write(iout,'(i5)') inp%ncnf
 write(iout,'(a/a)') '# nvcnf','#   n    l    f'
 do jj=2,inp%ncnf+1
   write(iout,'(i5)') cnf(jj)%nv
   do ii=inp%nc+1,inp%nc+cnf(jj)%nv
     write(iout,'(2i5,f8.2)') cnf(jj)%na(ii),cnf(jj)%la(ii), &
&                             cnf(jj)%fa(ii)
   end do
   write(iout,'(a)') '#'
 end do

 call write_mesh(log_mesh,iout)

 call write_reference_results(pname,it,etot,res,eterm,iout,ifail)
 if(ifail/=0) return

 if(istage==1) then
!all-electron halves of the test configurations; their pseudo-atom halves
!and all output about them come from the matching _ps.x
   warns%stage=stage_configs
   do jj=1,inp%ncnf+1
     call run_config_ae(ae%cache,jj,cnf,inp%nc,rr, &
&                       inp%zz, &
&                       mmax,adams,inp%iexc,ea,irel,c_light,nkap,aec(jj), &
&                       warns,iout, &
&                       ifail)
     if(ifail/=0) return
   end do

   call write_ae_file(ae_filename,irel,c_light,nkap,inp%atsym,inp%zz,inp%nc, &
&                     inp%nv,inp%iexc,inp%na,inp%la,inp%fa,mmax,rr,it,etot, &
&                     eterm,res,ea,rhoc,rho,vfull,inp%ncnf,cnf,aec,iout, &
&                     ifail)
   if(ifail/=0) return

   write(iout,'(/2a)') 'All-electron data written to ',trim(ae_filename)
   return
 end if

!the last log mesh point at or inside each rc (R-2)
 at%stage=stage_channel
 call core_radius_points(inp%lmax,inp%lloc,inp%rc,rr,mmax,irc,rcmax,iout, &
&                        ifail)
 if(ifail/=0) return

!
 cvgplt(:,:,:,:,:)=0.0d0
 ecurve(:,:,:,:,:)=0.0d0
 ecutw=0.0d0
 ektail(:,:,:,:)=0.0d0
 optdg(:,:,:,:)=0.0d0
 wshare(:,:,:,:)=0.0d0
 wcurve(:,:,:)=0.0d0
!
! loop to construct pseudopotentials for all angular momenta
!
 write(iout,'(/a/a)') &
&      'Begin loop to  construct optimized pseudo wave functions',&
&      'and semi-local pseudopoentials for all angular momenta'

!temporarily set this to 1 so that the pseudo wave function needed for the
!local potential will be generated.  Reset after run_vkb.
 npa(:,:)=0
 inp%nproj(inp%lloc+1)=1
 do l1=1,inp%lmax+1
   at%ll=l1-1
   call optimize_channel(ae,inp,l1,inp%nproj(l1),irc(l1),chan,iout,ifail)
   if(ifail/=0) return
   call warn_merge(warns,chan%warn)
   irc(l1)=chan%irc
   qmsbf(l1)=chan%qmsbf
   npa(:,l1)=chan%npa(:)
   wellprj(:,l1)=chan%wellprj(:)
   epa(:,l1,1:nkap)=chan%epa(:,1:nkap)
   dlta(:,:,l1,:)=chan%dlta(:,:,:)
   if(l1==inp%lloc+1) vrloc(:,:)=chan%vrloc(:,:)
   pswf(:,:,l1,:)=chan%pswf(:,:,:)
   pswfx(:,:,l1,:)=chan%pswfx(:,:,:)
   chix(:,:,l1,:)=chan%chix(:,:,:)
   vp(:,l1,:)=chan%vp(:,:)
   vkb(:,:,l1,:)=chan%vkb(:,:,:)
   cvgplt(:,:,:,l1,:)=chan%cvgplt(:,:,:,:)
   ecurve(:,:,:,l1,:)=chan%ecurve(:,:,:,:)
   ektail(:,:,l1,:)=chan%ektail(:,:,:)
   optdg(:,:,l1,:)=chan%optdg(:,:,:)
   wshare(:,:,l1,:)=chan%wshare(:,:,:)
   wcurve(:,l1,:)=chan%wcurve(:,:)
 end do !l1

! construct Vanderbilt / Kleinman-Bylander projectors

 write(iout,'(/a,a)') 'Construct Vanderbilt / Kleinmman-Bylander projectors'

 at%stage=stage_vkb
 at%ll=-1
 call run_vkb(inp%lmax,inp%lloc,inp%lpopt,inp%dvloc0,inp%rc,irc,inp%nproj,rr, &
&             mmax,pswf,pswfx, &
&             chix,vfull,vp,evkb,vkb,nlim,dlta,vrloc,irel,nkap,vkbpw, &
&             vlocpw,bbeig,iout,ifail)
 if(ifail/=0) return

!restore this to its proper value
 inp%nproj(inp%lloc+1)=0


! accumulate charge and eigenvalues
! pseudo wave functions are calculated with VKB projectors for
! maximum consistency of unscreening
! get all-electron and pseudopotential valence-state by valence-state
! charge densities

 at%stage=stage_valence
 warns%stage=stage_valence
 call valence_densities(irel,c_light,nkap,inp%nc,inp%nv,inp%na,inp%la,inp%fa, &
&                       ea,it,inp%lloc,inp%nproj,rr,inp%zz,vfull,vp,vlocpw, &
&                       vkbpw,evkb,mmax,adams,uu,up,uupsa,eeig,zval,rho,rhops, &
&                       rhoae, &
&                       rhotae,iout,ifail)
 if(ifail/=0) return

 if(irel==2) then
! decomposition into scalar-relativistic and spin-orbit projectors

   call sr_so_r(inp%lmax,inp%nproj,rr,mmax,evkb,vkbpw, &
&               vsr,esr,vso,eso,vsrpw,vsopw,warns,iout,ifail)
   if(ifail/=0) return
 end if

!the earlier stages' results, from which the model core and the tests
!follow (run_coupled)
 st%ae=ae
 st%ae%aec(:)=aec(:)
 st%nkap=nkap
 st%zion=zion
 st%rcmax=rcmax
 st%nrl=nrl
 st%qmsbf(:)=qmsbf(:)
 st%npa=npa
 st%epa=epa
 st%wellprj=wellprj
 st%vp=vp
 st%vkb=vkb
 st%evkb=evkb
 st%vkbpw=vkbpw
 st%vlocpw=vlocpw
 st%cvgplt=cvgplt
 st%ecurve=ecurve
 st%ektail=ektail
 st%optdg=optdg
 st%wshare=wshare
 st%wcurve=wcurve
 st%bbeig=bbeig
 st%rc_in(:)=rc_in(:)
 st%fcfact_in=fcfact_in
 st%rho=rho
 st%rhoae=rhoae
 st%rhops=rhops
 st%rhotae=rhotae
 st%uupsa=uupsa
 st%eeig=eeig
 st%zval=zval
 st%esr=esr
 st%eso=eso
 st%vsrpw=vsrpw
 st%vsopw=vsopw
 st%warn=warns

 return
 end subroutine run_stages

 subroutine run_coupled(irel,istage,inp,cnf,st,ps_filename,skip_ld,psp,at, &
&                       iout,ifail)

!! the model core, the screening with the ionic tail, the diagnostics, the
!! ghost tests, the test configurations, the log derivatives and the
!! results psp, from the earlier stages' results st, which it does not
!! change: called again with another icmod, fcfact or rcfact in inp, it
!! gives what a whole run with them gives
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param istage  0 the whole calculation, 1 the all-electron stage only, 2 the
!!                pseudopotential stage only
!! @param inp  the input values (input_m), completed (complete_input)
!! @param cnf  the reference and test configurations (input_m)
!! @param st  the earlier stages' results
!! @param ps_filename  the pseudopotential file for oncvpsp_ld.x, or blank
!! @param skip_ld  whether to leave out the log derivatives
!! @param psp  the results (results_m)
!! @param at  where the run was (failure_m): its stage, l and ikap
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use config_ae_m, only: run_config_ae, ae_config_t
 use failure_m, only: failure_t, fail_none, stage_model_core, &
&                    stage_screening, stage_diag, stage_ghosts, &
&                    stage_unscreen, stage_configs, stage_plot, stage_phsft, &
&                    stage_tests, stage_output
 use ps_file_m, only: write_ps_file
 use vout_m, only: vout
 use ionic_tail_m, only: coulomb_tail, ionic_screening
 use modcore1_m, only: modcore
 use modcore2_m, only: modcore2
 use modcore3_m, only: modcore3
 use modcore_budget_m, only: modcore_budget, modcore_budget_table, &
&                            modcore_probes
 use modcore_m, only: core_pieces_smooth
 use vloc_cutoff_m, only: vloc_cutoff
 use diagout_m, only: diagout
 use gnu_script_m, only: gnu_script
 use run_config_m, only: run_config_ps
 use run_diag_m, only: run_diag
 use run_ghosts_m, only: run_ghosts
 use run_phsft_m, only: run_phsft
 use run_plot_m, only: run_plot
 use ps_valence_m, only: unscreen_potentials, coulomb_switch
 use channel_m, only: ncurve
 use pieces_m, only: pw_free
 use input_m, only: input_t, config_t
 use results_m, only: results_t
 use adams_m, only: adams_t
 use constants_m, only: mxl, mxlin, mxstate, mxkap, mxcnf, nerrlv, mxprj, &
&                       nrhomod

!Input variables
 integer, intent(in) :: irel,istage
 type(config_t), intent(in) :: cnf(mxcnf)
 type(stages_t), intent(in) :: st
 character(len=*), intent(in) :: ps_filename
 logical, intent(in) :: skip_ld
 integer, intent(in) :: iout

!Input/output variables
 type(input_t), intent(inout) :: inp

!Output variables
 type(results_t), intent(out) :: psp
 type(failure_t), intent(out) :: at
 integer, intent(out) :: ifail

!Local variables
!the warnings, the stages' and their own
 type(warnings_t) :: warns
!the earlier stages' results, copied, as the steps below change some
 type(ae_state_t) :: ae
 type(ae_config_t) :: aec(mxcnf)
 type(adams_t) :: adams
 integer :: mmax,nkap,nrl
 real(dp) :: c_light,zion,rcmax,etot,eeig,zval,fcfact_in,ecutw
 real(dp) :: fcgrid(3),rcgrid(3),log_mesh(4),rc_in(mxlin+1),qmsbf(mxlin+1)
 real(dp) :: ea(mxstate,mxkap)
 integer, allocatable :: npa(:,:)
 logical, allocatable :: wellprj(:,:)
 real(dp), allocatable :: rr(:),rhoc(:),vfull(:),rho(:),rhot(:),vo(:),vxc(:)
 real(dp), allocatable :: vpuns(:,:),vion(:),vpsml(:,:,:)
 real(dp), allocatable :: epa(:,:,:),vp(:,:,:),vkb(:,:,:,:),evkb(:,:,:)
 real(dp), allocatable :: cvgplt(:,:,:,:,:),ecurve(:,:,:,:,:)
 real(dp), allocatable :: ektail(:,:,:,:),optdg(:,:,:,:),wshare(:,:,:,:)
 real(dp), allocatable :: wcurve(:,:,:),bbeig(:,:,:)
 real(dp), allocatable :: rhomod(:,:),rhoae(:,:),rhops(:,:),rhotae(:)
 real(dp), allocatable :: uupsa(:,:,:),esr(:,:),eso(:,:)
!probe densities and their l for the hardness of icmod=7
 real(dp), allocatable :: probes(:,:)
 integer, allocatable :: lprobe(:)
 integer :: nprobe
 type(piecewise_t), allocatable :: vkbpw(:,:,:),vsrpw(:,:),vsopw(:,:)
 type(piecewise_t) :: vlocpw
 type(piecewise_t) :: corepw(nrhomod)
 integer :: ii,jj,l1,npsh,nq
 real(dp) :: eeel,eexc,epstot,rcore,vlchk
 real(dp) :: cvloc(2,nerrlv),rphs(4)
!the test configurations (run_config_ps): rows of eigenvalues (-1 if one
!failed), their n, l, ikap and f, eae, eps, and the energy differences
 integer :: ntcf(mxcnf),itcf(3,2*mxstate,mxcnf)
 real(dp) :: etcf(3,2*mxstate,mxcnf),dtcf(3,mxcnf)
 real(dp), allocatable :: pshae(:,:,:),pshps(:,:,:)

 ifail=0
 at=failure_t(fail_none,0,0,stage_model_core,-1,0,0)
 warns=st%warn
 warns%stage=stage_model_core
 ae=st%ae
 aec(:)=st%ae%aec(:)
 adams=ae%adams
 mmax=ae%mmax
 c_light=ae%c_light
 rr=ae%rr
 rhoc=ae%rhoc
 vfull=ae%vfull
 ea(:,:)=ae%ea(:,:)
 etot=ae%etot
 log_mesh(:)=ae%log_mesh(:)
 nkap=st%nkap
 zion=st%zion
 rcmax=st%rcmax
 nrl=st%nrl
 qmsbf(:)=st%qmsbf(:)
 rc_in(:)=st%rc_in(:)
 fcfact_in=st%fcfact_in
 npa=st%npa
 epa=st%epa
 wellprj=st%wellprj
 vp=st%vp
 vkb=st%vkb
 evkb=st%evkb
 vkbpw=st%vkbpw
 vlocpw=st%vlocpw
 cvgplt=st%cvgplt
 ecurve=st%ecurve
 ektail=st%ektail
 optdg=st%optdg
 wshare=st%wshare
 wcurve=st%wcurve
 bbeig=st%bbeig
 rho=st%rho
 rhoae=st%rhoae
 rhops=st%rhops
 rhotae=st%rhotae
 uupsa=st%uupsa
 eeig=st%eeig
 zval=st%zval
 esr=st%esr
 eso=st%eso
 vsrpw=st%vsrpw
 vsopw=st%vsopw
 ecutw=0.0d0
 allocate(rhot(mmax),vo(mmax),vxc(mmax),vpuns(mmax,mxl+2))
!the columns past max(lmax,lloc)+1 stay 0
 vpuns(:,:)=0.0d0

!icmod=4 coarse-search grids: min, max, step
 fcgrid=[inp%fcfact_min,inp%fcfact_max,inp%fcfact_step]
 rcgrid=[inp%rcfact_min,inp%rcfact_max,inp%rcfact_step]

 allocate(rhomod(mmax,nrhomod))

 rhomod(:,:)=0.0d0

! construct model core charge based on monotonic polynomial fit
! or Teter function fit

 if(inp%icmod==1) then
   call modcore(inp%icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod,corepw, &
&               inp%fcfact,inp%rcfact,rcmax,mmax,adams,rr,inp%nc,inp%nv, &
&               inp%la, &
&               zion,inp%iexc,c_light,iout,ifail)
   if(ifail/=0) return

 else if(inp%icmod==2) then
   call modcore2(inp%icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod,corepw, &
&               inp%fcfact,inp%rcfact,mmax,adams,rr,inp%nc,inp%nv,inp%la,zion, &
&               inp%iexc,c_light,iout,ifail)
   if(ifail/=0) return

 else if(inp%icmod==6 .or. inp%icmod==7) then
   if(inp%fcfact<=0.0d0) then
     call modcore_budget_table(ecurve,ncurve,mxprj*4*nkap,wcurve, &
&                              inp%budget_cutoff==1,inp%budget_cutoff==2, &
&                              inp%rcfact,ecutw,iout,ifail)
     if(ifail/=0) return
   end if
!the probe densities of the hardness of icmod=7
   allocate(probes(mmax,mxprj*(inp%lmax+1)),lprobe(mxprj*(inp%lmax+1)))
   nprobe=0
   if(inp%icmod==7) then
     call modcore_probes(inp%lmax,inp%lloc,inp%nproj,epa,wellprj,vp,vlocpw, &
&                        vkbpw,evkb, &
&                        rr,mmax,adams, &
&                        nkap,probes,lprobe,nprobe,iout,ifail)
     if(ifail/=0) return
   end if
   call modcore_budget(inp%icmod,rho,rhoc,rhomod,corepw,inp%fcfact,ecutw, &
&               inp%budget_cutoff==1,inp%rcfact, &
&               inp%slack, &
&               minval(inp%rc(1:inp%lmax+1)),mmax,adams,rr,inp%iexc,c_light, &
&               zion, &
&               probes(:,1:nprobe),lprobe(1:nprobe),nprobe,psp%core%budget, &
&               warns,iout,ifail)
   if(ifail/=0) return
   deallocate(probes,lprobe)

 else if(inp%icmod>=3) then
   call modcore3(inp%icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod,corepw, &
&               inp%fcfact,inp%rcfact,fcgrid,rcgrid,mmax,adams,rr,inp%nc, &
&               inp%nv, &
&               inp%la,zion,inp%iexc,c_light,psp%core%teter,iout,ifail)
   if(ifail/=0) return

 end if
! the pieces of no model
 if(corepw(1)%nseg==0) then
   call core_pieces_smooth(rr,rhomod,mmax,corepw,ifail=ifail)
   if(ifail/=0) return
 end if

! screening potential for pseudocharge

 at%stage=stage_screening
 warns%stage=stage_screening
 call vout(.true.,rho,rhomod(1,1),vo,vxc,zval,eeel,eexc,rr,mmax,adams, &
&          inp%iexc, &
&          c_light,ifail)
 if(ifail/=0) return

! switch the screened potentials to the ionic tail beyond rcmax, iterating
! the pseudo valence states with it to a fixed point (after the model core,
! which V_ion uses)

 allocate(vion(mmax))
 call ionic_screening(nkap,inp%nc,inp%nv,inp%la,inp%fa,ea,inp%lmax,inp%lloc, &
&                     inp%nproj,rr,vp,vlocpw,vkbpw,evkb,mmax,adams,inp%zz, &
&                     zion,zval,inp%iexc,c_light,rcmax,rhoc,rhotae,rhomod(1, &
&                     1),rho,rhops,uupsa,eeig,vo,vxc,eeel,eexc,vion,iout, &
&                     ifail)
 if(ifail/=0) return

! total energy output

 epstot= eeig + eexc - 0.5d0*eeel
 write(iout,'(/a,f12.6/)') 'Pseudoatom total energy', epstot

 at%stage=stage_diag
 warns%stage=stage_diag
 at%ikap=1
 call run_diag(inp%lmax,npa,epa,wellprj,inp%lloc,inp%rc, &
&               vlocpw,vkbpw,evkb,inp%nproj,rr,vfull,vp,inp%zz,mmax,adams, &
&               irel,c_light,nkap, &
&               .false.,vsrpw,esr,vsopw,eso,warns,iout,ifail)
 if(ifail/=0) return

 if(irel==2) then
   at%ikap=2
   call run_diag(inp%lmax,npa,epa,wellprj,inp%lloc,inp%rc, &
&                 vlocpw,vkbpw,evkb,inp%nproj,rr,vfull,vp,inp%zz,mmax,adams, &
&                 irel,c_light, &
&                 nkap,.true.,vsrpw,esr,vsopw,eso,warns,iout,ifail)
   if(ifail/=0) return

! save full info for psml, switched to -zion/r far out as vpuns is
   allocate(vpsml(mmax,mxl+2,mxkap))
   vpsml(:,:,:)=0.0d0
   do l1=1,max(inp%lmax+1,inp%lloc+1)
     vpsml(:,l1,1) = vp(:,l1,1) - vo(:)
     vpsml(:,l1,2) = vp(:,l1,2) - vo(:)
     call coulomb_switch(rr,mmax,zion,vpsml(:,l1,1))
     call coulomb_switch(rr,mmax,zion,vpsml(:,l1,2))
   end do

! ghost testing

   write(iout,'(/a)') 'Ghost test for J = L + 1/2'
 end if

 at%stage=stage_ghosts
 warns%stage=stage_ghosts
 at%ikap=1
 call run_ghosts(inp%lmax,inp%la,ea(1,1),inp%nc,inp%nv,inp%lloc,inp%rc,qmsbf, &
&                vlocpw,vkbpw(:,:,1),evkb(1,1,1),inp%nproj,rr,vp,mmax,adams, &
&                1,psp%tests%ghosts,warns,iout,ifail)
 if(ifail/=0) return

 if(irel==2) then
   write(iout,'(/a)') 'Ghost test for J = L - 1/2'
   vkb(:,:,1,2)=vkb(:,:,1,1)
   evkb(:,1,2)=evkb(:,1,1)
   do jj=1,inp%nproj(1)
     vkbpw(jj,1,2)=vkbpw(jj,1,1)
   end do
   at%ikap=2
   call run_ghosts(inp%lmax,inp%la,ea(1,2),inp%nc,inp%nv,inp%lloc,inp%rc, &
&                  qmsbf, &
&                  vlocpw,vkbpw(:,:,2),evkb(1,1,2),inp%nproj,rr,vp,mmax,adams, &
&                  2,psp%tests%ghosts,warns,iout,ifail)
   if(ifail/=0) return
   vkb(:,:,1,2)=0.0d0
   evkb(:,1,2)=0.0d0
   do jj=1,inp%nproj(1)
     call pw_free(vkbpw(jj,1,2))
   end do
 end if

! unscreen semi-local potentials, taking the scalar-relativistic average
! when fully relativistic; -zion/r far out
 call unscreen_potentials(irel,inp%lmax,inp%lloc,vp,vo,zion,rr,mmax,nkap,vpuns)

! how closely the local potential is -zion/r beyond rcmax
 at%stage=stage_unscreen
 warns%stage=stage_unscreen
 at%ikap=0
 call coulomb_tail(vpuns(1,inp%lloc+1),vion,rr,mmax,zion,rcmax,iout,ifail)
 if(ifail/=0) return


! loop over reference plus test atom configurations
 rhot(:)=rho(:)
 ntcf(:)=0
 itcf(:,:,:)=0
 etcf(:,:,:)=0.0d0
 dtcf(:,:)=0.0d0
 at%stage=stage_configs
 warns%stage=stage_configs
 do jj=1,inp%ncnf+1

   write(iout,'(/a,i2)') 'Test configuration',jj-1

! charge density is initialized to that of reference configuration

   rhot(:)=rho(:)

   if(istage/=2) then
     call run_config_ae(ae%cache,jj,cnf,inp%nc,rr, &
&                       inp%zz, &
&                       mmax,adams,inp%iexc,ea,irel,c_light,nkap,aec(jj), &
&                       warns,iout, &
&                       ifail)
     if(ifail/=0) return
   end if

   call run_config_ps(jj,aec(jj),inp%nc,rhot,rhomod,rr, &
&                     rcmax,mmax,adams,inp%iexc,etot,epstot,inp%nproj,vpuns, &
&                     inp%lloc,vlocpw,vkbpw,evkb,irel,c_light,nkap,ntcf(jj), &
&                     itcf(:,:,jj),etcf(:,:,jj),dtcf(:,jj),warns,iout,ifail)
   if(ifail/=0) return

 end do !jj

 at%stage=stage_plot
 warns%stage=stage_plot
 call run_plot(inp%lmax,npa,epa,wellprj,inp%lloc,inp%rc, &
&               vkb,vlocpw,vkbpw,evkb,inp%nproj,rr,vfull,vp,vpuns,inp%zz, &
&               mmax,adams, &
&               inp%drl,nrl, &
&               rho,rhoc,rhomod,irel,c_light,nkap,cvgplt,warns,iout,ifail)
 if(ifail/=0) return




!the energies of the log derivatives, from epsh2 down by depsh; none with
!--no-log-derivatives
 if(skip_ld) then
   npsh=0
   rphs(:)=0.0d0
   allocate(pshae(npsh,4,nkap),pshps(npsh,4,nkap))
 else
   npsh=int(((inp%epsh2-inp%epsh1)/inp%depsh)-0.5d0)+1
   allocate(pshae(npsh,4,nkap),pshps(npsh,4,nkap))
   at%stage=stage_phsft
   warns%stage=stage_phsft
   call run_phsft(inp%lmax,inp%lloc,inp%nproj,epa,npsh,inp%epsh2,inp%depsh, &
&                 inp%rxpsh,vlocpw,vkbpw,evkb,rr,vfull,vp,inp%zz,mmax,adams, &
&                 inp%rc,irel,c_light,nkap,rphs,pshae,pshps,iout,ifail)
   if(ifail/=0) return
 end if

 call gnu_script(epa,evkb,inp%lmax,inp%lloc,inp%nproj,nkap,iout)

 at%stage=stage_tests
 warns%stage=stage_tests
 call vloc_cutoff(rr,mmax,rho,vpuns(1,inp%lloc+1),zion,cvloc,vlchk,ifail)
 if(ifail/=0) return

!the tests' results
 psp%tests%logder%ne=npsh
 allocate(psp%tests%logder%ee(npsh))
 do ii=1,npsh
   psp%tests%logder%ee(ii)=inp%epsh2-inp%depsh*real(ii-1,dp)
 end do
 psp%tests%logder%rl(:)=rphs(:)
 call move_alloc(pshae,psp%tests%logder%ae)
 call move_alloc(pshps,psp%tests%logder%ps)
 psp%tests%configs%ncfg=inp%ncnf+1
 psp%tests%configs%nrow(:)=ntcf(:)
 psp%tests%configs%state(:,:,:)=itcf(:,:,:)
 psp%tests%configs%ee(:,:,:)=etcf(:,:,:)
 psp%tests%configs%de(:,:)=dtcf(:,:)
 call move_alloc(cvgplt,psp%tests%convergence%cross)
!the curves to the end of the longest
 nq=0
 do ii=1,ncurve
   if(any(ecurve(1,ii,:,:,:)>0.0d0)) nq=ii
 end do
 psp%tests%convergence%curve=ecurve(:,1:nq,:,:,:)
 call move_alloc(ektail,psp%tests%convergence%tail)
 call move_alloc(optdg,psp%tests%convergence%opt)
 call move_alloc(wshare,psp%tests%convergence%wshare)
 psp%tests%weighted_cut%used=(inp%icmod==6 .or. inp%icmod==7) &
&                         .and. inp%fcfact<=0.0d0 .and. inp%budget_cutoff==1
 psp%tests%weighted_cut%errw=inp%rcfact
 psp%tests%weighted_cut%ecutw=ecutw
 call move_alloc(wcurve,psp%tests%weighted_cut%weight)
 psp%tests%vloc_cut%cut(:,:)=cvloc(:,:)
 psp%tests%vloc_cut%check=vlchk

 call move_alloc(bbeig,psp%tests%bbeig)

 call diagout(inp%lmax,inp%nproj,nkap,psp%tests,inp%icmod,rr,mmax, &
&             rcmax,rhoc,corepw,inp%lloc,evkb,iout,ifail)
 if(ifail/=0) return

! Write info for PSML format
!
!
! rcore: the radius beyond which the model core equals the all-electron
! core, its joint (doc/smoothness_plan.md, M-14); with no joint, the
! first point beyond which they are equal
!
 if(corepw(1)%nseg==2) then
   rcore=corepw(1)%bound(1)
 else
   rcore=rr(mmax)
   do ii=mmax,1,-1
     if(rhomod(ii,1)/=rhoc(ii)) exit
     rcore=rr(ii)
   end do
 end if
!the fully-relativistic semi-local potentials are vp-vo for each kappa,
!without the corrections of vpuns (see doc/fr_merge_plan.md)
 if(irel/=2) then
   allocate(vpsml(mmax,mxl+2,1))
   vpsml(:,:,1)=vpuns(:,:)
 end if

!the pseudopotential
 psp%pot%atsym=inp%atsym
 psp%pot%zz=inp%zz
 psp%pot%zion=zion
 psp%pot%c_light=c_light
 psp%pot%irel=irel
 psp%pot%nkap=nkap
 psp%pot%iexc=inp%iexc
 psp%pot%mmax=mmax
 call move_alloc(rr,psp%pot%rr)
 psp%pot%lmax=inp%lmax
 psp%pot%lloc=inp%lloc
 psp%pot%nproj(:)=inp%nproj(:)
 psp%pot%rcmax=rcmax
 call move_alloc(evkb,psp%pot%evkb)
 call move_alloc(vkb,psp%pot%vkb)
 call move_alloc(vkbpw,psp%pot%vkbpw)
 call move_alloc(esr,psp%pot%esr)
 call move_alloc(eso,psp%pot%eso)
 call move_alloc(vsrpw,psp%pot%vsrpw)
 call move_alloc(vsopw,psp%pot%vsopw)
 psp%pot%vlocpw=vlocpw
 psp%pot%vloc=vp(:,inp%lloc+1,1:nkap)
 call move_alloc(vpuns,psp%pot%vpuns)
 call move_alloc(rho,psp%pot%rho)
 call move_alloc(rhomod,psp%pot%rhomod)
 psp%pot%corepw(:)=corepw(:)
!the rest the writers take
 psp%nrl=nrl
 psp%epstot=epstot
 psp%ea(:,:)=ea(:,:)
 call move_alloc(rhotae,psp%rhotae)
 call move_alloc(rhoc,psp%rhoc)
 psp%rcore=rcore
 call move_alloc(vpsml,psp%vpsml)
 call move_alloc(uupsa,psp%uupsa)
 psp%echo%rc_in(:)=rc_in(:)
 psp%echo%fcfact_in=fcfact_in
 psp%echo%fcgrid(:)=fcgrid(:)
 psp%echo%rcgrid(:)=rcgrid(:)
 psp%echo%log_mesh(:)=log_mesh(:)
 psp%warnings=warns

!the pseudopotential for oncvpsp_ld.x and oncvpsp_confine.x
 if(len_trim(ps_filename)/=0) then
   at%stage=stage_output
   warns%stage=stage_output
   call write_ps_file(trim(ps_filename),psp%pot,inp%rc,inp%epsh1, &
&                     inp%epsh2,inp%depsh,inp%rxpsh,epa,vp,iout,ifail)
   if(ifail/=0) return
 end if
 return
 end subroutine run_coupled

 subroutine run_ae(irel,istage,inp,ae_filename,cnf,ae,at,iout,ifail,ae_in)

!! the input inp completed (complete_input) and its all-electron atom
!! (solve_ae), or ae_in in its place (given_ae); an aeatom that does not
!! converge is a failure
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param istage  0 the whole calculation, 1 the all-electron stage only, 2 the
!!                pseudopotential stage only
!! @param inp  the input values (input_m), as read on input, completed on
!!             output
!! @param ae_filename  the all-electron data file of istage 1 and 2
!! @param cnf  the reference and test configurations (input_m)
!! @param ae  the all-electron atom (ae_state_m)
!! @param at  where the run was (failure_m): its stage, l and ikap
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
!! @param ae_in  the all-electron atom solved before, if given (solve_ae)

 use ae_state_m, only: ae_state_t, solve_ae, given_ae
 use input_m, only: input_t, config_t
 use constants_m, only: mxcnf
 use failure_m, only: failure_t, fail_none, fail_convergence, failed, &
&                    msg_run_aeatom, stage_input, stage_ae

!Input variables
 integer, intent(in) :: irel,istage
 character(len=*), intent(in) :: ae_filename
 integer, intent(in) :: iout
 type(ae_state_t), intent(in), optional :: ae_in

!Input/output variables
 type(input_t), intent(inout) :: inp

!Output variables
 type(config_t), intent(out) :: cnf(mxcnf)
 type(ae_state_t), intent(out) :: ae
 type(failure_t), intent(out) :: at
 integer, intent(out) :: ifail

!Local variables
 integer :: ierr,nkap

 at=failure_t(fail_none,0,0,stage_input,-1,0,0)
 call complete_input(inp,cnf,iout,ifail)
 if(ifail/=0) return

!the all-electron atom: its mesh, and the reference configuration solved
!or read
 at%stage=stage_ae
 nkap=1
 if(irel==2) nkap=2
 ierr=0
 if(present(ae_in)) then
   call given_ae(ae_in,irel,inp,ae,iout,ifail)
 else
   call solve_ae(irel,nkap,istage,inp,cnf,ae_filename,ae,ierr,iout,ifail)
 end if
 if(ifail/=0) return
 if(ierr/=0) then
   write(iout,'(/a,i0,a,i0,a,es12.4)') &
&        'oncvpsp_run: ERROR aeatom failed for the reference '// &
&        'configuration'// &
&        ', ierr=',ierr,', it=',ae%it,', res=',ae%res
   ifail=failed(fail_convergence,msg_run_aeatom,ierr)
   return
 end if
 return
 end subroutine run_ae

 subroutine run_channel(ae,inp,ll,mode,ch,lp,prj,at,iout,ifail)

!! the channel of angular momentum ll alone from the input inp and the
!! all-electron atom ae, solved before (solve_ae) and checked against it
!! (check_ae): the input completed (complete_input), the core radii's mesh
!! points, and the channel (optimize_channel) with the input's projectors,
!! or one for the local channel. A first pseudo-wave with a node is
!! pspot's failure. With mode 1, also the local potential
!! (local_potential, for lloc<4 from the local channel, optimized for it
!! when ll is another) and the channel's projectors (channel_vkb, none for
!! the local channel), the same as the whole run's; with mode 2, the
!! channel only up to the node test of each kappa's first projector, its
!! results unset
!!
!! @param ae  the all-electron atom (ae_state_m)
!! @param inp  the input values (input_m), as read on input, completed on
!!             output
!! @param ll  angular momentum
!! @param mode  0 the channel, 1 also the local potential and projectors,
!!              2 the node test of its first projector only
!! @param ch  the channel (channel_m)
!! @param lp  the local potential (run_vkb_m), with mode 1
!! @param prj  the channel's projectors (run_vkb_m), with mode 1
!! @param at  where the run was (failure_m): its stage, l and ikap
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use ae_state_m, only: ae_state_t, check_ae
 use channel_m, only: channel_t, optimize_channel
 use run_vkb_m, only: local_t, projectors_t, local_potential, channel_vkb
 use main_setup_m, only: core_radius_points
 use input_m, only: input_t, config_t
 use constants_m, only: dp, mxcnf, mxlin, mxprj
 use failure_m, only: failure_t, fail_none, fail_input, failed, &
&                    msg_run_channel_l, msg_run_channel_mode, stage_input, &
&                    stage_ae, stage_channel, stage_vkb

!Input variables
 type(ae_state_t), intent(in) :: ae
 integer, intent(in) :: ll,mode
 integer, intent(in) :: iout

!Input/output variables
 type(input_t), intent(inout) :: inp

!Output variables
 type(channel_t), intent(out) :: ch
 type(local_t), intent(out) :: lp
 type(projectors_t), intent(out) :: prj
 type(failure_t), intent(out) :: at
 integer, intent(out) :: ifail

!Local variables
 type(config_t) :: cnf(mxcnf)
 type(channel_t) :: chl
 integer :: irc(mxlin+1)
 integer :: nprj,nkap
 real(dp) :: rcmax

 at=failure_t(fail_none,0,0,stage_input,-1,0,0)
 call complete_input(inp,cnf,iout,ifail)
 if(ifail/=0) return
 if(ll<0 .or. ll>inp%lmax) then
   write(iout,'(/a,i0,a,i0)') 'run_channel: ERROR no channel l=',ll, &
&       ', lmax=',inp%lmax
   ifail=failed(fail_input,msg_run_channel_l)
   return
 end if
 if(mode<0 .or. mode>2) then
   write(iout,'(/a,i0)') 'run_channel: ERROR the mode is 0, 1 or 2, not ',mode
   ifail=failed(fail_input,msg_run_channel_mode)
   return
 end if
 at%stage=stage_ae
 call check_ae(ae,ae%irel,inp,iout,ifail)
 if(ifail/=0) return
 at%stage=stage_channel
 call core_radius_points(inp%lmax,inp%lloc,inp%rc,ae%rr,ae%mmax,irc,rcmax, &
&                        iout,ifail)
 if(ifail/=0) return
 nprj=inp%nproj(ll+1)
 if(ll==inp%lloc) nprj=1
 at%ll=ll
 call optimize_channel(ae,inp,ll+1,nprj,irc(ll+1),ch,iout,ifail,mode==2)
 if(ifail/=0 .or. mode/=1) return

! the local potential, from the local channel's first wave for lloc<4
 nkap=size(ch%vp,2)
 if(inp%lloc<4 .and. ll/=inp%lloc) then
   at%ll=inp%lloc
   call optimize_channel(ae,inp,inp%lloc+1,1,irc(inp%lloc+1),chl,iout,ifail)
   if(ifail/=0) return
 else
   chl=ch
 end if
 if(inp%lloc<4) irc(inp%lloc+1)=chl%irc
 irc(ll+1)=ch%irc
 at%stage=stage_vkb
 at%ll=-1
 call local_potential(inp%lloc,inp%lpopt,inp%dvloc0,inp%rc,irc,ae%rr, &
&                     ae%mmax,ae%vfull,nkap,chl%pswfx,chl%chix,chl%vrloc,lp, &
&                     iout,ifail)
 if(ifail/=0) return

! the channel's projectors
 allocate(prj%vkb(ae%mmax,mxprj,nkap))
 prj%vkb(:,:,:)=ch%vkb(:,:,:)
 if(ll==inp%lloc) return
 at%ll=ll
 call channel_vkb(ll+1,inp%lloc,lp,inp%rc(ll+1),irc(ll+1),inp%nproj(ll+1), &
&                 ae%rr,ae%mmax,ae%vfull,ch%pswfx,ch%chix,ch%dlta,ae%irel, &
&                 nkap,prj%evkb(:,1:nkap),prj%vkb,prj%vkbpw(:,1:nkap), &
&                 prj%bbeig(:,1:nkap),iout,ifail)
 return
 end subroutine run_channel

 subroutine complete_input(inp,cnf,iout,ifail)

!! the input as the stages take it: the well's options checked; the
!! configurations, the reference first and the core states in each;
!! fcfact 0 without a model core; the values checked (check_data); a core
!! radius of 0 the largest of the others (doc/smoothness_plan.md, R-1); and
!! no projectors for lloc
!!
!! @param inp  the input values (input_m), as read on input
!! @param cnf  the reference and test configurations (input_m)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use input_m, only: input_t, config_t
 use constants_m, only: dp, mxcnf
 use check_data_m, only: check_data
 use failure_m, only: fail_input, failed, msg_run_well_n_confined, &
&                    msg_run_well_norm_free

!Input variables
 integer, intent(in) :: iout

!Input/output variables
 type(input_t), intent(inout) :: inp

!Output variables
 type(config_t), intent(out) :: cnf(mxcnf)
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: fcgrid(3),rcgrid(3),rcmax
 integer :: jj,l1

 ifail=0
!the free tail has no finite norm over all space
 if(inp%well_tail==1 .and. inp%well_norm==1) then
   write(iout,'(/a)') &
&        'oncvpsp_run: ERROR well_norm = "space" (1) is not defined '// &
&       'for well_tail = "free": the free tail has no finite norm'
   ifail=failed(fail_input,msg_run_well_norm_free)
   return
 end if
 if(inp%well_tail/=2 .and. any(inp%well_n(:)>0)) then
   write(iout,'(/a)') &
&       'oncvpsp_run: ERROR well_n is for well_tail = "confined" only'
   ifail=failed(fail_input,msg_run_well_n_confined)
   return
 end if

!the configurations, the reference first, the core states in each
 cnf(:)=inp%cnf(:)
 cnf(1)%nv=inp%nv
 cnf(1)%na(:)=inp%na(:)
 cnf(1)%la(:)=inp%la(:)
 cnf(1)%fa(:)=inp%fa(:)
 do jj=2,inp%ncnf+1
   cnf(jj)%na(1:inp%nc)=inp%na(1:inp%nc)
   cnf(jj)%la(1:inp%nc)=inp%la(1:inp%nc)
   cnf(jj)%fa(1:inp%nc)=inp%fa(1:inp%nc)
 end do

 if(inp%icmod==0) then
   inp%fcfact=0.0d0
 end if

!icmod=4 coarse-search grids: min, max, step
 fcgrid=[inp%fcfact_min,inp%fcfact_max,inp%fcfact_step]
 rcgrid=[inp%rcfact_min,inp%rcfact_max,inp%rcfact_step]

 call check_data(inp,fcgrid,rcgrid,iout,ifail)
 if(ifail/=0) return

! the core radii are used as given (doc/smoothness_plan.md, R-1)
 rcmax=0.0d0
 do l1=1,inp%lmax+1
   rcmax=max(rcmax,inp%rc(l1))
 end do
 do l1=1,inp%lmax+1
   if(inp%rc(l1)==0.0d0) then
     inp%rc(l1)=rcmax
   end if
 end do

 inp%nproj(inp%lloc+1)=0
 return
 end subroutine complete_input

 end module oncvpsp_run_m

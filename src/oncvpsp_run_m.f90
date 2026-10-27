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

 use constants_m, only: dp, mxl, mxlin, mxstate, mxkap, nrhomod
 use pieces_m, only: piecewise_t
 use results_m, only: results_t
 implicit none
 private
 public :: oncvpsp_run


 contains

 subroutine oncvpsp_run(irel,istage,inp,ae_filename,ps_filename,skip_ld, &
&                       psp,iout,ifail)
!! the pseudopotential of the input inp, and its tests: the text output
!! to iout, and the pseudopotential and the tests in psp. inp is completed
!! (configurations, core radii) as the writers echo it. With istage=1
!! only the all-electron atom, written to ae_filename; with istage=2 that
!! file in place of it. A ps_filename, if not blank, receives the
!! pseudopotential for oncvpsp_ld.x; skip_ld leaves out the log
!! derivatives.
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param istage  0 the whole calculation, 1 the all-electron stage only, 2 the
!!                pseudopotential stage only
!! @param inp  the input, as read
!! @param ae_filename  the all-electron data file of istage 1 and 2
!! @param ps_filename  the pseudopotential file for oncvpsp_ld.x, or blank
!! @param skip_ld  whether to leave out the log derivatives
!! @param psp  the results (results_m)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written
 use config_ae_m, only: run_config_ae, aeatom_cached, ae_config_t, &
&                       ae_cache_t
 use ae_file_m, only: write_ae_file, read_ae_file
 use ps_file_m, only: write_ps_file
 use exc_libxc_m, only: exc_libxc_info
 use vout_m, only: vout
 use sr_so_r_m, only: sr_so_r
 use ionic_tail_m, only: coulomb_tail, ionic_screening
 use modcore1_m, only: modcore
 use modcore2_m, only: modcore2
 use modcore3_m, only: modcore3
 use modcore_budget_m, only: modcore_budget, modcore_budget_table, &
&                            modcore_probes
 use run_optimize_m, only: run_optimize
 use vloc_cutoff_m, only: vloc_cutoff
 use check_data_m, only: check_data
 use diagout_m, only: diagout
 use gnu_script_m, only: gnu_script
 use run_config_m, only: run_config_ps
 use run_diag_m, only: run_diag
 use run_ghosts_m, only: run_ghosts
 use run_phsft_m, only: run_phsft
 use run_plot_m, only: run_plot
 use run_vkb_m, only: run_vkb
 use kappa_m, only: kap_count, kap_of, j_share
 use modcore_m, only: core_pieces_smooth
 use well_mesh_m, only: well_mesh_t
 use ps_valence_m, only: valence_densities, unscreen_potentials
 use ae_projectors_m, only: ae_bound_projectors, ae_well_projectors, &
&                           ae_delta_window, ae_overlap, ecut_weights
 use constants_m, only: nerrdg, nerrlv
 use input_echo_m, only: write_atom, write_lmax, write_icmod, write_epsh, &
&                        write_mesh, fx, vkb_head, vkb_line
 use relativity_m, only: c_light_of
 use main_setup_m, only: main_names, &
&                        make_log_mesh, write_reference_results, &
&                        core_radius_points
 use constants_m, only: mxprj, mxcnf
 use pieces_m, only: pw_free
 use input_m, only: input_t
 implicit none

!Input variables
 integer, intent(in) :: irel,istage
 character(len=*), intent(in) :: ae_filename,ps_filename
 logical, intent(in) :: skip_ld
 integer, intent(in) :: iout

!Input/output variables
 type(input_t), intent(inout) :: inp

!Output variables
 type(results_t), intent(out) :: psp
 integer, intent(out) :: ifail




!Local variables
 integer :: ii,ierr,ios,it,nq
 integer :: jj,ll,l1,inline
 integer :: mmax,nlim,nrl
 integer :: iprj,nkap,nnwell,npsh
 integer :: ikap,kap,mkap
 integer,allocatable :: npa(:,:)
!
 integer :: irc(mxlin+1)
! the ghost tests' margins for the diagnostics block (run_ghosts), per kap
 integer :: ngh(mxl+1,mxkap),ghn(mxstate,mxl+1,mxkap)
 real(dp) :: ghe(2,mxstate,mxl+1,mxkap),ghp(2,mxl+1,mxkap)

 real(dp) :: al,eeel
 real(dp) :: eeig,eexc
 real(dp) :: fcfact0,soscale
 real(dp) :: fcgrid(3),rcgrid(3)
 real(dp) :: rcmax
 real(dp) :: log_mesh(3)
 real(dp) :: zion,zval,etot,res
!> the speed of light of the run, from the input
 real(dp) :: c_light
 real(dp) :: cvloc(2,nerrlv),vlchk
 real(dp) :: eterm(4),rphs(4)
 real(dp) :: fj,rcore
!
 real(dp) :: ea(mxstate,mxkap)
 real(dp) :: qmsbf(mxlin+1),rc0(mxlin+1)
 real(dp) :: rpk(mxstate,mxkap)
 real(dp) :: epstot

 real(dp), allocatable :: evkb(:,:,:),cvgplt(:,:,:,:,:),qq(:,:)
 real(dp), allocatable :: ecurve(:,:,:,:,:)
!the r-space less the Fourier kinetic energy of each pseudo-wave
!(run_optimize)
 real(dp), allocatable :: ektail(:,:,:,:),optdg(:,:,:,:),wshare(:,:,:,:)
!the eigenvalues of the B matrices (run_vkb)
 real(dp), allocatable :: bbeig(:,:,:)
!the test configurations (run_config_ps): rows of eigenvalues (-1 if one
!failed), their n, l, ikap and f, eae, eps, and the energy differences
 integer :: ntcf(mxcnf),itcf(3,2*mxstate,mxcnf)
 real(dp) :: etcf(3,2*mxstate,mxcnf),dtcf(3,mxcnf)
 real(dp), allocatable :: pshae(:,:,:),pshps(:,:,:)
 real(dp) :: ecutw
!weights of the pseudo-waves' curves in the weighted cutoff (budget_cutoff)
 real(dp), allocatable :: wcurve(:,:,:)
! the points of a convergence curve: run_optimize's q grid,
! 1 + 24 rc max(2 q of the last basis function, 2 qcut, 20), ~2400 for
! nbas=14 and less than 3000 for qcut up to 15 at rc up to 4; a buffer,
! with room for qcut rc up to 160
 integer, parameter :: ncurve=8000
 real(dp), allocatable :: rr(:)
 real(dp), allocatable :: rho(:),rhoc(:),rhot(:)
 real(dp), allocatable :: uu(:,:),up(:,:),vwell(:,:)
!> whether each projector of an l is a well state, confined by its well,
!> and of each l
 logical :: inwell(mxprj)
!> each projector's confined well state on its continued mesh (well_mesh_m)
 type(well_mesh_t) :: wmesh(mxprj)
 logical, allocatable :: wellprj(:,:)
 real(dp), allocatable :: vp(:,:,:),vfull(:),vkb(:,:,:,:),pswf(:,:,:,:)
 real(dp), allocatable :: vpuns(:,:)
!frozen-core ionic potential of the tail of the local potential
 real(dp), allocatable :: vion(:)
 real(dp), allocatable :: vo(:),vxc(:)
 real(dp), allocatable :: rhomod(:,:),rhoae(:,:),rhops(:,:),rhotae(:)
 real(dp), allocatable :: uupsa(:,:,:) !pseudo-atomic orbitals array
!probe densities and their l for the hardness of icmod=7
 real(dp), allocatable :: probes(:,:)
 integer, allocatable :: lprobe(:)
 integer :: nprobe
 real(dp), allocatable :: epa(:,:,:)
 real(dp), allocatable :: uua(:,:,:),upa(:,:,:)
 real(dp), allocatable :: dlta(:,:,:,:),vrloc(:,:),pswfx(:,:,:,:)
 type(piecewise_t), allocatable :: vkbpw(:,:,:)
 type(piecewise_t) :: vlocpw
 type(piecewise_t) :: corepw(nrhomod)
 logical :: well_tail_free
!> the cutoff of the pseudo-waves for icmod 6, 7 (budget_cutoff): 0 the
!> largest of the ground-state waves, 1 their weighted mean, 2 the largest
!> of all
!> the tail of the well states (0 "well", 1 "free", 2 "confined"), and the
!> n of each l's first well state with "confined" (0 for the default)
!> the normalization of the well states: 0 that of well_tail, 1 all
!> space, 2 inside rc; and whether they are normalized inside rc
 logical :: well_norm_rc
!> the weight of the first pseudo-wave of each l with no occupied state in
!> the cutoff of the pseudo-waves (ecut_weight; negative for 2(2l+1))
 real(dp), allocatable :: chix(:,:,:,:)
!fully-relativistic only: scalar-relativistic and spin-orbit projectors, and
!the semi-local potentials for psmlout_r
 real(dp), allocatable :: vsr(:,:,:),vso(:,:,:),esr(:,:),eso(:,:)
 type(piecewise_t), allocatable :: vsrpw(:,:),vsopw(:,:)
 real(dp), allocatable :: vpsml(:,:,:)

!all-electron results for the test configurations (see run_config_ae)
 type(ae_config_t) :: aec(mxcnf)
!the all-electron atoms of this run, for the test configurations
 type(ae_cache_t) :: cache


!names in the messages: the equation solved and the program
 character(len=19) :: treatment
 character(len=9) :: pname
 character(len=16) :: exe

 ifail=0
 call main_names(irel,istage,nkap,treatment,pname,exe)
 well_tail_free=inp%well_tail==1
!the free tail has no finite norm over all space
 if(inp%well_tail==1 .and. inp%well_norm==1) then
   write(iout,'(/a)') &
&        'oncvpsp: ERROR well_norm = "space" (1) is not defined '// &
&       'for well_tail = "free": the free tail has no finite norm'
   ifail=1
   return
 end if
 well_norm_rc=inp%well_norm==2 .or. (inp%well_norm==0 .and. inp%well_tail/=0)
 if(inp%well_tail/=2 .and. any(inp%well_n(:)>0)) then
   write(iout,'(/a)') 'oncvpsp: ERROR well_n is for well_tail = "confined" only'
   ifail=1
   return
 end if

 inp%cnf(1)%nv=inp%nv
 do ii=1,inp%nc+inp%nv
   inp%cnf(1)%n(ii)=inp%na(ii)
   inp%cnf(1)%l(ii)=inp%la(ii)
   inp%cnf(1)%f(ii)=inp%fa(ii)
 end do

 do jj=2,inp%ncnf+1
   do ii=1,inp%nc
     inp%cnf(jj)%n(ii)=inp%na(ii)
     inp%cnf(jj)%l(ii)=inp%la(ii)
     inp%cnf(jj)%f(ii)=inp%fa(ii)
   end do
 end do

 if(inp%icmod==0) then
   inp%fcfact=0.0d0
 end if

!icmod=4 coarse-search grids: min, max, step
 fcgrid=[inp%fcfact_min,inp%fcfact_max,inp%fcfact_step]
 rcgrid=[inp%rcfact_min,inp%rcfact_max,inp%rcfact_step]

 call check_data(inp,fcgrid,rcgrid,iout,ifail)
 if(ifail/=0) return
 c_light=c_light_of(inp%c_light)

 nrl=int((inp%rlmax/inp%drl)-0.5d0)+1

!PWSCF wants an even number of mesh pointe (not made so when fully
!relativistic, see doc/fr_merge_plan.md)
 if(mod(nrl,2)/=0 .and. irel/=2) nrl=nrl+1

 call make_log_mesh(irel,inp%zz,inp%amesh,inp%rr1,inp%rmax,al,mmax,rr,log_mesh)

!calculate zion for output
 zion=inp%zz
 do ii=1,inp%nc
   zion=zion-inp%fa(ii)
 end do


 allocate(rho(mmax),rhoc(mmax),rhot(mmax))
 allocate(uu(mmax,2),up(mmax,2),uupsa(mmax,nkap,inp%nv))
 allocate(vwell(mmax,mxprj))
 allocate(evkb(mxprj,mxl+1,nkap),qq(mxprj,mxprj))
!the projectors beyond nproj and of lloc stay 0
 evkb(:,:,:)=0.0d0
 allocate(cvgplt(2,nerrdg,mxprj,mxl+1,nkap))
 allocate(ecurve(2,ncurve,mxprj,mxl+1,nkap))
 allocate(ektail(2,mxprj,mxl+1,nkap),bbeig(mxprj,mxl+1,nkap))
 allocate(optdg(4,mxprj,mxl+1,nkap),wshare(3,mxprj,mxl+1,nkap))
 allocate(wcurve(mxprj,mxl+1,nkap))
 allocate(vp(mmax,mxl+2,nkap),vfull(mmax),vkb(mmax,mxprj,mxl+1,nkap))
 allocate(pswf(mmax,mxprj,mxl+1,nkap))
 allocate(vpuns(mmax,mxl+2))
 allocate(vo(mmax),vxc(mmax))
 allocate(rhoae(mmax,inp%nv),rhops(mmax,inp%nv),rhotae(mmax))
 allocate(npa(mxprj,mxlin+1))
 allocate(wellprj(mxprj,mxlin+1))
 wellprj(:,:)=.false.
 allocate(epa(mxprj,mxlin+1,nkap))
 allocate(uua(mmax,mxprj,nkap),upa(mmax,mxprj,nkap))
 allocate(vsr(mmax,2*mxprj,mxl+1),vso(mmax,2*mxprj,mxl+1))
 allocate(esr(2*mxprj,mxl+1),eso(2*mxprj,mxl+1))
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

!
! full potential atom solution
!
 if(inp%iexc<0) call exc_libxc_info(inp%iexc,iout)
 if(istage==2) then
   call read_ae_file(ae_filename,irel,c_light,nkap,inp%atsym,inp%zz,inp%nc, &
&                    inp%nv,inp%iexc,inp%na,inp%la,inp%fa,mmax,rr,inp%ncnf, &
&                    inp%cnf,it,etot,eterm,res,ea,rhoc,rho,vfull,aec,iout, &
&                    ifail)
   if(ifail/=0) return
 else
   ierr=0
   call aeatom_cached(cache,irel,c_light,nkap,inp%na,inp%la,ea,inp%fa,rpk, &
&                     inp%nc,inp%nc+inp%nv,it,rhoc,rho,rr,vfull,inp%zz,mmax, &
&                     inp%iexc,etot,eterm,res,ierr,iout,ifail)
   if(ifail/=0) return
   if(ierr/=0) then
     write(iout,'(/a,i0,a,i0,a,es12.4)') &
&          'oncvpsp: ERROR aeatom failed for the reference configuration'// &
&          ', ierr=',ierr,', it=',it,', res=',res
     ifail=1
     return
   end if
 end if
!
!

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
 rc0(:)=inp%rc(:)
!input value, for the echoes (modcore3 resets fcfact for icmod=4)
 fcfact0=inp%fcfact

! output printing (echos input data, with all-electron eigenvalues added)

 call write_atom(inp%atsym,inp%zz,inp%nc,inp%nv,inp%iexc,c_light,inp%psfile, &
&                iout)
 write(iout,'(a/a)') '#','#   n    l    f        energy (Ha)'
 do ii=1,inp%nc+inp%nv
   if(irel/=2) then
     write(iout,'(2i5,f8.2,1pe18.7)') inp%na(ii),inp%la(ii),inp%fa(ii),ea(ii,1)
   else if(inp%la(ii)==0) then
     write(iout,'(2i5,f8.2,1p,d18.7)') inp%na(ii),inp%la(ii),inp%fa(ii),ea(ii,1)
     ea(ii,2)=ea(ii,1)
   else
     write(iout,'(2i5,f8.2,1p,2d18.7)') inp%na(ii),inp%la(ii),inp%fa(ii), &
&          ea(ii, &
&          1),ea(ii,2)
   end if
 end do

 call write_lmax(inp%lmax,inp%ep_as_given,inp%well_onset,inp%well_tail, &
&                inp%wreg,inp%well_norm,iout)
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
&      trim(vkb_head(inp%lmax,inp%well_n,inp%ecut_wt))
 do l1=1,inp%lmax+1
   write(iout,'(a)') trim(vkb_line(inp%lmax,l1,inp%nproj,inp%debl,inp%well_n, &
&        inp%ecut_wt))
 end do

 call write_icmod(inp%icmod,inp%fcfact,inp%rcfact,inp%slack,inp%budget_cut, &
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
   write(iout,'(i5)') inp%cnf(jj)%nv
   do ii=inp%nc+1,inp%nc+inp%cnf(jj)%nv
     write(iout,'(2i5,f8.2)') inp%cnf(jj)%n(ii),inp%cnf(jj)%l(ii), &
&                             inp%cnf(jj)%f(ii)
   end do
   write(iout,'(a)') '#'
 end do

 call write_mesh(log_mesh,iout)

 call write_reference_results(pname,it,etot,res,eterm,iout)

 if(istage==1) then
!all-electron halves of the test configurations; their pseudo-atom halves
!and all output about them come from the matching _ps.x
   do jj=1,inp%ncnf+1
     call run_config_ae(cache,jj,inp%cnf,inp%nc,rr, &
&                       inp%zz, &
&                       mmax,inp%iexc,ea,irel,c_light,nkap,aec(jj),iout,ifail)
     if(ifail/=0) return
   end do

   call write_ae_file(ae_filename,irel,c_light,nkap,inp%atsym,inp%zz,inp%nc, &
&                     inp%nv,inp%iexc,inp%na,inp%la,inp%fa,mmax,rr,it,etot, &
&                     eterm,res,ea,rhoc,rho,vfull,inp%ncnf,inp%cnf,aec,iout, &
&                     ifail)
   if(ifail/=0) return

   write(iout,'(/2a)') 'All-electron data written to ',trim(ae_filename)
   return
 end if

!the last log mesh point at or inside each rc (R-2)
 call core_radius_points(inp%lmax,inp%lloc,inp%rc,rr,mmax,irc,rcmax)

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
   vwell(:,:)=0.0d0
   do jj=1,mxprj
     wmesh(jj)%mx=0
   end do
   ll=l1-1
   mkap=kap_count(nkap,ll)
! loop on J = ll +/- 1/2 when nkap=2; fj is the weight of each J
   do ikap=1,mkap
     kap=kap_of(ll,ikap)
     fj=1.0d0
     if(nkap==2) fj=j_share(ll,ikap,1.0d0)

     uu(:,:)=0.0d0; qq(:,:)=0.0d0
     if(ikap==1) inwell(:)=.false.

     call ae_bound_projectors(irel,c_light,ll,kap,ikap,inp%nc,inp%nv,inp%na, &
&                             inp%la,ea,inp%nproj(l1),it,rr,inp%zz,vfull,mmax, &
&                             uu,up,iprj,epa(:,l1,ikap),npa(:,l1),uua(:,:, &
&                             ikap),upa(:,:,ikap),iout,ifail)
     if(ifail/=0) return
     call ae_well_projectors(irel,c_light,ll,kap,ikap,inp%nproj(l1),inp%ep(l1, &
&                            ikap),inp%debl(l1),inp%well_n(l1), &
&                            inp%ep_as_given,inp%well_onset,inp%well_tail, &
&                            well_norm_rc,inp%rc(l1),irc(l1),inp%qcut(l1),rr, &
&                            inp%zz,vfull,mmax,uu,up,nnwell,iprj,epa(:,l1, &
&                            ikap),npa(:,l1),uua(:,:,ikap),upa(:,:,ikap), &
&                            vwell,inwell,wmesh,iout,ifail)
     if(ifail/=0) return
     call ae_delta_window(irel,c_light,ll,kap,inp%nproj(l1),epa(:,l1,ikap),rr, &
&                         vfull,uua(:,:,ikap),upa(:,:,ikap),mmax,al, &
&                         inp%rc(l1),irc(l1),l1==inp%lloc+1,dlta(:,:,l1,ikap), &
&                         vrloc(:,ikap),iout,ifail)
     if(ifail/=0) return
     call ae_overlap(irel,c_light,ll,kap,inp%nproj(l1),uua(:,:,ikap),rr, &
&                    inp%zz,mmax,inp%rc(l1),qq)
     call run_optimize(epa(1,l1,ikap),ll,kap,irel,c_light,mmax,mxprj,rr, &
&                      uua(1,1,ikap),upa(1,1,ikap),qq, &
&                      well_tail_free .and. inwell(:),wmesh,inp%rc(l1), &
&                      irc(l1), &
&                      inp%qcut(l1),qmsbf(l1),inp%ncon(l1),inp%nbas(l1), &
&                      inp%nproj(l1), &
&                      pswf(1,1,l1,ikap),vp(1,l1,ikap),vkb(1,1,l1,ikap), &
&                      vfull,cvgplt(1,1,1,l1,ikap),ncurve, &
&                      ecurve(1,1,1,l1,ikap), &
&                      pswfx(1,1,l1,ikap),chix(1,1,l1,ikap), &
&                      ektail(1,1,l1,ikap),inp%wreg,optdg(1,1,l1,ikap), &
&                      wshare(1,1,l1,ikap),iout,ifail)
     if(ifail/=0) return

!the weights of its pseudo-waves in the cutoff of the pseudo-waves
     call ecut_weights(ll,inp%nc,inp%nv,inp%la,inp%fa,fj,inp%nproj(l1), &
&                      inp%ecut_wt(l1), &
&                      wcurve(:,l1,ikap),iout,ifail)
     if(ifail/=0) return

   end do !ikap
   wellprj(:,l1)=inwell(:)
 end do !l1

! construct Vanderbilt / Kleinman-Bylander projectors

 write(iout,'(/a,a)') 'Construct Vanderbilt / Kleinmman-Bylander projectors'

 call run_vkb(inp%lmax,inp%lloc,inp%lpopt,inp%dvloc0,inp%rc,irc,inp%nproj,rr, &
&             mmax,mxprj,pswf,pswfx, &
&             chix,vfull,vp,evkb,vkb,nlim,dlta,vrloc,irel,nkap,vkbpw, &
&             vlocpw,bbeig,iout,ifail)
 if(ifail/=0) return

!restore this to its proper value
 inp%nproj(inp%lloc+1)=0

 deallocate(uua,upa)

! accumulate charge and eigenvalues
! pseudo wave functions are calculated with VKB projectors for
! maximum consistency of unscreening
! get all-electron and pseudopotential valence-state by valence-state
! charge densities

 call valence_densities(irel,c_light,nkap,inp%nc,inp%nv,inp%na,inp%la,inp%fa, &
&                       ea,it,inp%lloc,inp%nproj,rr,inp%zz,vfull,vp,vlocpw, &
&                       vkbpw,evkb,mmax,uu,up,uupsa,eeig,zval,rho,rhops,rhoae, &
&                       rhotae,iout,ifail)
 if(ifail/=0) return

 if(irel==2) then
! decomposition into scalar-relativistic and spin-orbit projectors

   call sr_so_r(inp%lmax,inp%nproj,rr,mmax,mxprj,evkb,vkbpw, &
&               vsr,esr,vso,eso,vsrpw,vsopw,iout,ifail)
   if(ifail/=0) return
 end if

 allocate(rhomod(mmax,nrhomod))

 rhomod(:,:)=0.0d0

! construct model core charge based on monotonic polynomial fit
! or Teter function fit

 if(inp%icmod==1) then
   call modcore(inp%icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod,corepw, &
&               inp%fcfact,inp%rcfact,rcmax,mmax,rr,inp%nc,inp%nv,inp%la, &
&               zion,inp%iexc,c_light,iout,ifail)
   if(ifail/=0) return

 else if(inp%icmod==2) then
   call modcore2(inp%icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
&               inp%fcfact,inp%rcfact,mmax,rr,inp%nc,inp%nv,inp%la,zion, &
&               inp%iexc,c_light,iout,ifail)
   if(ifail/=0) return

 else if(inp%icmod==6 .or. inp%icmod==7) then
   if(inp%fcfact<=0.0d0) then
     call modcore_budget_table(ecurve,ncurve,mxprj*4*nkap,wcurve, &
&                              inp%budget_cut==1,inp%budget_cut==2, &
&                              inp%rcfact,ecutw,iout,ifail)
     if(ifail/=0) return
   end if
!the probe densities of the hardness of icmod=7
   allocate(probes(mmax,mxprj*(inp%lmax+1)),lprobe(mxprj*(inp%lmax+1)))
   nprobe=0
   if(inp%icmod==7) then
     call modcore_probes(inp%lmax,inp%lloc,inp%nproj,epa,wellprj,vp,vlocpw, &
&                        vkbpw,evkb, &
&                        rr,mmax,mxprj, &
&                        nkap,probes,lprobe,nprobe,iout,ifail)
     if(ifail/=0) return
   end if
   call modcore_budget(inp%icmod,rho,rhoc,rhomod,corepw,inp%fcfact,ecutw, &
&               inp%budget_cut==1,inp%rcfact, &
&               inp%slack, &
&               minval(inp%rc(1:inp%lmax+1)),mmax,rr,inp%iexc,c_light,zion, &
&               probes(:,1:nprobe),lprobe(1:nprobe),nprobe,iout,ifail)
   if(ifail/=0) return
   deallocate(probes,lprobe)

 else if(inp%icmod>=3) then
   call modcore3(inp%icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
&               inp%fcfact,inp%rcfact,fcgrid,rcgrid,mmax,rr,inp%nc,inp%nv, &
&               inp%la,zion,inp%iexc,c_light,iout,ifail)
   if(ifail/=0) return

 end if
! the pieces of a model with no joint (the Teter function, or none)
 if(corepw(1)%nseg==0) call core_pieces_smooth(rr,rhomod,mmax,corepw)

! screening potential for pseudocharge

 call vout(1,rho,rhomod(1,1),vo,vxc,zval,eeel,eexc,rr,mmax,inp%iexc,c_light)

! switch the screened potentials to the ionic tail beyond rcmax, iterating
! the pseudo valence states with it to a fixed point (after the model core,
! which V_ion uses)

 allocate(vion(mmax))
 call ionic_screening(nkap,inp%nc,inp%nv,inp%la,inp%fa,ea,inp%lmax,inp%lloc, &
&                     inp%nproj,rr,vp,vlocpw,vkbpw,evkb,mmax,mxprj,inp%zz, &
&                     zion,zval,inp%iexc,c_light,rcmax,rhoc,rhotae,rhomod(1, &
&                     1),rho,rhops,uupsa,eeig,vo,vxc,eeel,eexc,vion,iout, &
&                     ifail)
 if(ifail/=0) return

! total energy output

 epstot= eeig + eexc - 0.5d0*eeel
 write(iout,'(/a,f12.6/)') 'Pseudoatom total energy', epstot

 call run_diag(inp%lmax,npa,epa,wellprj,inp%lloc,irc, &
&               vlocpw,vkbpw,evkb,inp%nproj,rr,vfull,vp,inp%zz,mmax,mxprj, &
&               irel,c_light,nkap, &
&               0,vsrpw,esr,vsopw,eso,iout,ifail)
 if(ifail/=0) return

 if(irel==2) then
   call run_diag(inp%lmax,npa,epa,wellprj,inp%lloc,irc, &
&                 vlocpw,vkbpw,evkb,inp%nproj,rr,vfull,vp,inp%zz,mmax,mxprj, &
&                 irel,c_light, &
&                 nkap,1,vsrpw,esr,vsopw,eso,iout,ifail)
   if(ifail/=0) return

! save full info for psml
   allocate(vpsml(mmax,mxl+2,mxkap))
   do l1=1,max(inp%lmax+1,inp%lloc+1)
     vpsml(:,l1,1) = vp(:,l1,1) - vo(:)
     vpsml(:,l1,2) = vp(:,l1,2) - vo(:)
   end do

! ghost testing

   write(iout,'(/a)') 'Ghost test for J = L + 1/2'
 end if

 ngh(:,:)=0
 ghn(:,:,:)=0
 ghe(:,:,:,:)=0.0d0
 ghp(:,:,:)=0.0d0
 call run_ghosts(inp%lmax,inp%la,ea(1,1),inp%nc,inp%nv,inp%lloc,inp%rc,qmsbf, &
&                vlocpw,vkbpw(:,:,1),evkb(1,1,1),inp%nproj,rr,vp,mmax,mxprj, &
&                ngh(:,1),ghn(:,:,1),ghe(:,:,:,1),ghp(:,:,1),iout,ifail)
 if(ifail/=0) return

 if(irel==2) then
   write(iout,'(/a)') 'Ghost test for J = L - 1/2'
   vkb(:,:,1,2)=vkb(:,:,1,1)
   evkb(:,1,2)=evkb(:,1,1)
   do jj=1,inp%nproj(1)
     vkbpw(jj,1,2)=vkbpw(jj,1,1)
   end do
   call run_ghosts(inp%lmax,inp%la,ea(1,2),inp%nc,inp%nv,inp%lloc,inp%rc, &
&                  qmsbf, &
&                  vlocpw,vkbpw(:,:,2),evkb(1,1,2),inp%nproj,rr,vp,mmax,mxprj, &
&                  ngh(:,2),ghn(:,:,2),ghe(:,:,:,2),ghp(:,:,2),iout,ifail)
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
 call coulomb_tail(vpuns(1,inp%lloc+1),vion,rr,mmax,zion,rcmax,iout)


! loop over reference plus test atom configurations
 rhot(:)=rho(:)
 ntcf(:)=0
 itcf(:,:,:)=0
 etcf(:,:,:)=0.0d0
 dtcf(:,:)=0.0d0
 do jj=1,inp%ncnf+1

   write(iout,'(/a,i2)') 'Test configuration',jj-1

! charge density is initialized to that of reference configuration

   rhot(:)=rho(:)

   if(istage/=2) then
     call run_config_ae(cache,jj,inp%cnf,inp%nc,rr, &
&                       inp%zz, &
&                       mmax,inp%iexc,ea,irel,c_light,nkap,aec(jj),iout,ifail)
     if(ifail/=0) return
   end if

   call run_config_ps(jj,aec(jj),inp%nc,rhot,rhomod,rr, &
&                     rcmax,mmax,mxprj,inp%iexc,etot,epstot,inp%nproj,vpuns, &
&                     inp%lloc,vlocpw,vkbpw,evkb,irel,c_light,nkap,ntcf(jj), &
&                     itcf(:,:,jj),etcf(:,:,jj),dtcf(:,jj),iout,ifail)
   if(ifail/=0) return

 end do !jj

 call run_plot(inp%lmax,npa,epa,wellprj,inp%lloc,irc, &
&               vkb,vlocpw,vkbpw,evkb,inp%nproj,rr,vfull,vp,vpuns,inp%zz, &
&               mmax,mxprj, &
&               inp%drl,nrl, &
&               rho,rhoc,rhomod,irel,c_light,nkap,cvgplt,iout,ifail)
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
   call run_phsft(inp%lmax,inp%lloc,inp%nproj,epa,npsh,inp%epsh2,inp%depsh, &
&                 inp%rxpsh,vlocpw,vkbpw,evkb,rr,vfull,vp,inp%zz,mmax,mxprj, &
&                 inp%rc,irel,c_light,nkap,rphs,pshae,pshps,iout,ifail)
   if(ifail/=0) return
 end if

 call gnu_script(epa,evkb,inp%lmax,inp%lloc,mxprj,inp%nproj,nkap,iout)

 call vloc_cutoff(rr,mmax,rho,vpuns(1,inp%lloc+1),zion,cvloc,vlchk)

!the tests' results
 psp%tests%logder%ne=npsh
 allocate(psp%tests%logder%e(npsh))
 do ii=1,npsh
   psp%tests%logder%e(ii)=inp%epsh2-inp%depsh*real(ii-1,dp)
 end do
 psp%tests%logder%r(:)=rphs(:)
 call move_alloc(pshae,psp%tests%logder%ae)
 call move_alloc(pshps,psp%tests%logder%ps)
 psp%tests%ghosts%nrow(:,:)=ngh(:,:)
 psp%tests%ghosts%n(:,:,:)=ghn(:,:,:)
 psp%tests%ghosts%e(:,:,:,:)=ghe(:,:,:,:)
 psp%tests%ghosts%least(:,:,:)=ghp(:,:,:)
 psp%tests%configs%ncfg=inp%ncnf+1
 psp%tests%configs%nrow(:)=ntcf(:)
 psp%tests%configs%state(:,:,:)=itcf(:,:,:)
 psp%tests%configs%e(:,:,:)=etcf(:,:,:)
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
&                         .and. inp%fcfact<=0.0d0 .and. inp%budget_cut==1
 psp%tests%weighted_cut%errw=inp%rcfact
 psp%tests%weighted_cut%ecutw=ecutw
 call move_alloc(wcurve,psp%tests%weighted_cut%weight)
 psp%tests%vloc_cut%cut(:,:)=cvloc(:,:)
 psp%tests%vloc_cut%check=vlchk

 call move_alloc(bbeig,psp%tests%bbeig)

 call diagout(inp%lmax,inp%nproj,mxprj,nkap,psp%tests,inp%icmod,rr,mmax, &
&             rcmax,rhoc,corepw,inp%lloc,evkb,iout)

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
 psp%echo%rc0(:)=rc0(:)
 psp%echo%fcfact0=fcfact0
 psp%echo%fcgrid(:)=fcgrid(:)
 psp%echo%rcgrid(:)=rcgrid(:)
 psp%echo%log_mesh(:)=log_mesh(:)

!the pseudopotential for oncvpsp_ld.x and oncvpsp_confine.x
 if(len_trim(ps_filename)/=0) then
   call write_ps_file(trim(ps_filename),psp%pot,inp%rc,inp%epsh1, &
&                     inp%epsh2,inp%depsh,inp%rxpsh,epa,vp,iout,ifail)
   if(ifail/=0) return
 end if
 return
 end subroutine oncvpsp_run

 end module oncvpsp_run_m

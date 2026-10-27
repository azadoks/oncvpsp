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

 module oncvpsp_main_m

 implicit none
 private
 public :: oncvpsp_main

 contains

 subroutine oncvpsp_main(irel,istage)
!! Creates and tests optimized norm-conserving Vanderbilt or Kleinman-Bylander
!! pseudopotentials based on D. R. Hamann, Phys. Rev. B 88, 085117 (2013)
!! and references therein.
!!
!!   D. R. Hamann
!!   Mat-Sim Research LLC
!!   P.O. Box 742
!!   Murray Hill, NJ 07974
!!   USA
!!
!!   Developed from original "gncpp" code of March 8,1987
!!
!!   Output format for ABINIT pspcod=8 and upf format for quantumespresso
!!
!!   The whole of each executable: oncvpsp.x, oncvpspnr.x and oncvpspr.x and
!!   their _ae.x and _ps.x stages are short programs calling this with their
!!   irel and istage.  With istage=1 (the _ae.x), it solves only the
!!   all-electron atom and the all-electron halves of the test
!!   configurations, and writes them to the file given by --ae-file.  With
!!   istage=2 (the _ps.x), it reads that file instead of solving the
!!   all-electron atom and otherwise runs as with istage=0, with the same
!!   output.  Both read the same input.  See ae_file_m.  Except in the
!!   all-electron stage, --ps-file writes the pseudopotential for
!!   oncvpsp_ld.x, which takes the log derivatives from it and the
!!   all-electron file at any radii and energies (ps_file_m), and
!!   --no-log-derivatives leaves them out of this run.
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param istage  0 the whole calculation, 1 the all-electron stage only, 2 the
!!                pseudopotential stage only
 use, intrinsic :: iso_fortran_env, only: stdin => input_unit
 use psmlout_m, only: psmlout
 use input_text_m, only: read_input_text
 use config_ae_m, only: run_config_ae, aeatom_cached
 use ae_file_m, only: write_ae_file, read_ae_file
 use ps_file_m, only: write_ps_file
 use exc_libxc_m, only: exc_libxc
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
 use linout_m, only: linout
 use run_config_m, only: run_config_ps
 use run_diag_m, only: run_diag
 use run_ghosts_m, only: run_ghosts
 use run_phsft_m, only: run_phsft
 use run_plot_m, only: run_plot
 use run_vkb_m, only: run_vkb
 use upfout_m, only: upfout
#if (defined WITH_TOML)
 use input_toml_m, only: read_input_toml
#endif
 use constants_m, only: dp, mxprj, mxstate, mxcnf, mxl, mxlin, mxkap, nrhomod
 use kappa_m, only: kap_count, kap_of, j_share
 use pieces_m, only: piecewise_t, pw_free
 use modcore_m, only: core_pieces_smooth
 use well_mesh_m, only: well_mesh_t
 use ps_valence_m, only: valence_densities, unscreen_potentials
 use ae_projectors_m, only: ae_bound_projectors, ae_well_projectors, &
&                           ae_delta_window, ae_overlap, ecut_weights
 use constants_m, only: nlev, nerrdg, nerrlv, well_onset_default
 use input_echo_m, only: write_atom, write_lmax, write_icmod, write_epsh, &
&                        write_mesh, fx, vkb_head, vkb_line
 use relativity_m, only: set_relativity
 use main_setup_m, only: main_names, parse_command_line, write_banner, &
&                        make_log_mesh, write_reference_results, &
&                        core_radius_points, INPUT_STDIN, INPUT_TEXT, &
&                        INPUT_TOML
 use input_m, only: input_t
 implicit none
 integer, parameter :: iout=6

 type(input_t) :: inp

!Input variables
 integer, intent(in) :: irel,istage

!Local variables
 integer :: ifail
 integer :: ii,ierr,ios,it
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
 real(dp) :: cvloc(2,nerrlv),vlchk
 real(dp) :: eterm(4),rphs(4)
 real(dp) :: fj,rcore
!
 real(dp) :: ea(mxstate,mxkap)
 real(dp) :: qmsbf(mxlin+1),rc0(mxlin+1)
 real(dp) :: rpk(mxstate,mxkap)
 real(dp) :: epstot
 real(dp), parameter :: eps=1.0d-8

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
 integer :: istcnf(mxcnf),ierrcnf(mxcnf),itcnf(mxcnf),kerrcnf(mxcnf)
 integer :: nvtcnf(mxcnf),natcnf(mxstate,mxcnf),latcnf(mxstate,mxcnf)
 real(dp) :: fatcnf(mxstate,nlev,mxcnf),etstcnf(mxcnf)
 real(dp), allocatable :: eatcnf(:,:,:,:)


!names in the messages: the equation solved and the program
 character(len=19) :: treatment
 character(len=9) :: pname
 character(len=16) :: exe

 integer :: input_mode
 integer :: unit
 character(len=1024) :: input_filename
 character(len=1024) :: ae_filename
 character(len=1024) :: ps_filename
 logical :: skip_ld

 call main_names(irel,istage,nkap,treatment,pname,exe)
 allocate(eatcnf(mxstate,nkap,nlev,mxcnf))
 call parse_command_line(istage,treatment,pname,exe,input_mode, &
&                        input_filename,ae_filename,ps_filename,skip_ld)
 call write_banner(irel,treatment,exe,iout)

 if(istage/=0 .and. len_trim(ae_filename) == 0) then
   write(iout,'(/4a)') trim(pname),': ERROR ',trim(exe),' requires --ae-file'
   error stop 1
 end if

 ios = 0
 inp%ep_as_given=.false.
 inp%well_onset=well_onset_default
 inp%well_tail=0
 inp%well_n(:)=0
 inp%well_norm=0
 inp%ecut_wt(:)=-1.0d0
 inp%wreg=0.0d0
 select case(input_mode)
 case(INPUT_STDIN)
   unit=stdin
 case(INPUT_TEXT,INPUT_TOML)
   open(newunit=unit, file=input_filename, status='old', action='read', &
&       iostat=ios)
   if(ios /= 0) then
     write(iout, '(/3a)') trim(pname),': ERROR cannot open input file, file=', &
&       trim(input_filename)
     error stop 1
   end if
 case default
   write(iout,'(/2a,i0)') trim(pname), &
&        ': ERROR unknown input mode, input_mode=', &
&        input_mode
   error stop 1
 end select
 if(input_mode==INPUT_TOML) then
#if (defined WITH_TOML)
   call read_input_toml(unit, inp, iout, ifail)
   if(ifail/=0) error stop 1
#endif
 else
   call read_input_text(unit, inline, inp, iout, ifail)
   if(ifail/=0) error stop 1
 end if
 if(input_mode/=INPUT_STDIN) close(unit)
 well_tail_free=inp%well_tail==1
!the free tail has no finite norm over all space
 if(inp%well_tail==1 .and. inp%well_norm==1) then
   write(iout,'(/a)') &
&        'oncvpsp: ERROR well_norm = "space" (1) is not defined '// &
&       'for well_tail = "free": the free tail has no finite norm'
   error stop 1
 end if
 well_norm_rc=inp%well_norm==2 .or. (inp%well_norm==0 .and. inp%well_tail/=0)
 if(inp%well_tail/=2 .and. any(inp%well_n(:)>0)) then
   write(iout,'(/a)') 'oncvpsp: ERROR well_n is for well_tail = "confined" only'
   error stop 1
 end if

 inp%nvcnf(1)=inp%nv
 do ii=1,inp%nc+inp%nv
   inp%nacnf(ii,1)=inp%na(ii)
   inp%lacnf(ii,1)=inp%la(ii)
   inp%facnf(ii,1)=inp%fa(ii)
 end do

 do jj=2,inp%ncnf+1
   do ii=1,inp%nc
     inp%nacnf(ii,jj)=inp%na(ii)
     inp%lacnf(ii,jj)=inp%la(ii)
     inp%facnf(ii,jj)=inp%fa(ii)
   end do
 end do

 if(inp%icmod==0) then
   inp%fcfact=0.0d0
 end if

!icmod=4 coarse-search grids: min, max, step
 fcgrid=[inp%fcfact_min,inp%fcfact_max,inp%fcfact_step]
 rcgrid=[inp%rcfact_min,inp%rcfact_max,inp%rcfact_step]

 call check_data(inp,fcgrid,rcgrid,iout,ifail)
 if(ifail/=0) error stop 1
 call set_relativity(inp%c_light)

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
 if(istage==2) then
   call read_ae_file(ae_filename,irel,nkap,inp%atsym,inp%zz,inp%nc,inp%nv, &
&                    inp%iexc,inp%na,inp%la,inp%fa, &
&                    mmax,rr,inp%ncnf,inp%nvcnf,inp%nacnf,inp%lacnf,inp%facnf, &
&                    it,etot,eterm,res,ea,rhoc,rho,vfull, &
&                    istcnf,ierrcnf,itcnf, &
&                    kerrcnf,nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf, &
&                    iout,ifail)
   if(ifail/=0) error stop 1
!aeatom prints the libxc functional description on its first exc_libxc call
   if(inp%iexc<0) then
     block
       real(dp) :: vxct(mmax),exct(mmax)
       call exc_libxc(inp%iexc,al,rho,vxct,exct,rr,mmax)
     end block
   end if
 else
   ierr=0
   call aeatom_cached(irel,nkap,inp%na,inp%la,ea,inp%fa,rpk,inp%nc, &
&                     inp%nc+inp%nv,it,rhoc,rho, &
&                     rr,vfull,inp%zz,mmax,inp%iexc,etot,eterm,res,ierr, &
&                     iout,ifail)
   if(ifail/=0) error stop 1
   if(ierr/=0) then
     write(iout,'(/a,i0,a,i0,a,es12.4)') &
&          'oncvpsp: ERROR aeatom failed for the reference configuration'// &
&          ', ierr=',ierr,', it=',it,', res=',res
     error stop 1
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

 call write_atom(inp%atsym,inp%zz,inp%nc,inp%nv,inp%iexc,inp%psfile,iout)
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
   write(iout,'(i5)') inp%nvcnf(jj)
   do ii=inp%nc+1,inp%nc+inp%nvcnf(jj)
     write(iout,'(2i5,f8.2)') inp%nacnf(ii,jj),inp%lacnf(ii,jj),inp%facnf(ii,jj)
   end do
   write(iout,'(a)') '#'
 end do

 call write_mesh(log_mesh,iout)

 call write_reference_results(pname,it,etot,res,eterm,iout)

 if(istage==1) then
!all-electron halves of the test configurations; their pseudo-atom halves
!and all output about them come from the matching _ps.x
   do jj=1,inp%ncnf+1
     call run_config_ae(jj,inp%nacnf,inp%lacnf,inp%facnf,inp%nc,inp%nvcnf,rr, &
&                       inp%zz, &
&                       mmax,inp%iexc,ea,irel,nkap,istcnf,ierrcnf,itcnf, &
&                       kerrcnf,nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf, &
&                       iout,ifail)
     if(ifail/=0) error stop 1
   end do

   call write_ae_file(ae_filename,irel,nkap,inp%atsym,inp%zz,inp%nc,inp%nv, &
&                     inp%iexc,inp%na,inp%la,inp%fa, &
&                     mmax,rr,it,etot,eterm,res,ea,rhoc,rho,vfull, &
&                     inp%ncnf,inp%nvcnf,inp%nacnf,inp%lacnf,inp%facnf, &
&                     istcnf,ierrcnf,itcnf, &
&                     kerrcnf,nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf, &
&                     iout,ifail)
   if(ifail/=0) error stop 1

   write(iout,'(/2a)') 'All-electron data written to ',trim(ae_filename)
   stop
 end if

!the last log mesh point at or inside each rc (R-2)
 call core_radius_points(inp%lmax,inp%lloc,inp%rc,rr,mmax,irc,rcmax)

!
 cvgplt(:,:,:,:,:)=0.0d0
 ecurve(:,:,:,:,:)=0.0d0
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

     call ae_bound_projectors(irel,ll,kap,ikap,inp%nc,inp%nv,inp%na,inp%la, &
&                             ea,inp%nproj(l1),it, &
&                             rr,inp%zz,vfull,mmax,uu,up,iprj,epa(:,l1,ikap), &
&                             npa(:,l1),uua(:,:,ikap),upa(:,:,ikap),iout,ifail)
     if(ifail/=0) error stop 1
     call ae_well_projectors(irel,ll,kap,ikap,inp%nproj(l1),inp%ep(l1,ikap), &
&                            inp%debl(l1),inp%well_n(l1),inp%ep_as_given, &
&                            inp%well_onset, &
&                            inp%well_tail,well_norm_rc,inp%rc(l1),irc(l1), &
&                            inp%qcut(l1),rr,inp%zz,vfull,mmax,uu,up,nnwell, &
&                            iprj, &
&                            epa(:,l1,ikap),npa(:,l1),uua(:,:,ikap), &
&                            upa(:,:,ikap),vwell,inwell,wmesh,iout,ifail)
     if(ifail/=0) error stop 1
     call ae_delta_window(irel,ll,kap,inp%nproj(l1),epa(:,l1,ikap),rr,vfull, &
&                         uua(:,:,ikap),upa(:,:,ikap),mmax,al,inp%rc(l1), &
&                         irc(l1),l1==inp%lloc+1,dlta(:,:,l1,ikap), &
&                         vrloc(:,ikap),iout,ifail)
     if(ifail/=0) error stop 1
     call ae_overlap(irel,ll,kap,inp%nproj(l1),uua(:,:,ikap),rr,inp%zz,mmax, &
&                    inp%rc(l1), &
&                    qq)
     call run_optimize(epa(1,l1,ikap),ll,kap,irel,mmax,mxprj,rr, &
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
     if(ifail/=0) error stop 1

!the weights of its pseudo-waves in the cutoff of the pseudo-waves
     call ecut_weights(ll,inp%nc,inp%nv,inp%la,inp%fa,fj,inp%nproj(l1), &
&                      inp%ecut_wt(l1), &
&                      wcurve(:,l1,ikap),iout,ifail)
     if(ifail/=0) error stop 1

   end do !ikap
   wellprj(:,l1)=inwell(:)
 end do !l1

! construct Vanderbilt / Kleinman-Bylander projectors

 write(iout,'(/a,a)') 'Construct Vanderbilt / Kleinmman-Bylander projectors'

 call run_vkb(inp%lmax,inp%lloc,inp%lpopt,inp%dvloc0,inp%rc,irc,inp%nproj,rr, &
&             mmax,mxprj,pswf,pswfx, &
&             chix,vfull,vp,evkb,vkb,nlim,dlta,vrloc,irel,nkap,vkbpw, &
&             vlocpw,bbeig,iout,ifail)
 if(ifail/=0) error stop 1

!restore this to its proper value
 inp%nproj(inp%lloc+1)=0

 deallocate(uua,upa)

! accumulate charge and eigenvalues
! pseudo wave functions are calculated with VKB projectors for
! maximum consistency of unscreening
! get all-electron and pseudopotential valence-state by valence-state
! charge densities

 call valence_densities(irel,nkap,inp%nc,inp%nv,inp%na,inp%la,inp%fa,ea,it, &
&                       inp%lloc,inp%nproj,rr,inp%zz, &
&                       vfull,vp,vlocpw,vkbpw,evkb,mmax,uu,up,uupsa,eeig, &
&                       zval,rho,rhops,rhoae,rhotae,iout,ifail)
 if(ifail/=0) error stop 1

 if(irel==2) then
! decomposition into scalar-relativistic and spin-orbit projectors

   call sr_so_r(inp%lmax,inp%nproj,rr,mmax,mxprj,evkb,vkbpw, &
&               vsr,esr,vso,eso,vsrpw,vsopw,iout,ifail)
   if(ifail/=0) error stop 1
 end if

 allocate(rhomod(mmax,nrhomod))

 rhomod(:,:)=0.0d0

! construct model core charge based on monotonic polynomial fit
! or Teter function fit

 if(inp%icmod==1) then
   call modcore(inp%icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod,corepw, &
&               inp%fcfact,inp%rcfact,rcmax,mmax,rr,inp%nc,inp%nv,inp%la, &
&               zion,inp%iexc,iout,ifail)
   if(ifail/=0) error stop 1

 else if(inp%icmod==2) then
   call modcore2(inp%icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
&               inp%fcfact,inp%rcfact,mmax,rr,inp%nc,inp%nv,inp%la,zion, &
&               inp%iexc,iout,ifail)
   if(ifail/=0) error stop 1

 else if(inp%icmod==6 .or. inp%icmod==7) then
   ecutw=0.0d0
   if(inp%fcfact<=0.0d0) then
     call modcore_budget_table(ecurve,ncurve,mxprj*4*nkap,wcurve, &
&                              inp%budget_cut==1,inp%budget_cut==2, &
&                              inp%rcfact,ecutw,iout,ifail)
     if(ifail/=0) error stop 1
   end if
!the probe densities of the hardness of icmod=7
   allocate(probes(mmax,mxprj*(inp%lmax+1)),lprobe(mxprj*(inp%lmax+1)))
   nprobe=0
   if(inp%icmod==7) then
     call modcore_probes(inp%lmax,inp%lloc,inp%nproj,epa,wellprj,vp,vlocpw, &
&                        vkbpw,evkb, &
&                        rr,mmax,mxprj, &
&                        nkap,probes,lprobe,nprobe,iout,ifail)
     if(ifail/=0) error stop 1
   end if
   call modcore_budget(inp%icmod,rho,rhoc,rhomod,corepw,inp%fcfact,ecutw, &
&               inp%budget_cut==1,inp%rcfact, &
&               inp%slack, &
&               minval(inp%rc(1:inp%lmax+1)),mmax,rr,inp%iexc,zion, &
&               probes(:,1:nprobe),lprobe(1:nprobe),nprobe,iout,ifail)
   if(ifail/=0) error stop 1
   deallocate(probes,lprobe)

 else if(inp%icmod>=3) then
   call modcore3(inp%icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
&               inp%fcfact,inp%rcfact,fcgrid,rcgrid,mmax,rr,inp%nc,inp%nv, &
&               inp%la,zion,inp%iexc,iout,ifail)
   if(ifail/=0) error stop 1

 end if
! the pieces of a model with no joint (the Teter function, or none)
 if(corepw(1)%nseg==0) call core_pieces_smooth(rr,rhomod,mmax,corepw)

! screening potential for pseudocharge

 call vout(1,rho,rhomod(1,1),vo,vxc,zval,eeel,eexc,rr,mmax,inp%iexc)

! switch the screened potentials to the ionic tail beyond rcmax, iterating
! the pseudo valence states with it to a fixed point (after the model core,
! which V_ion uses)

 allocate(vion(mmax))
 call ionic_screening(nkap,inp%nc,inp%nv,inp%la,inp%fa,ea,inp%lmax,inp%lloc, &
&                     inp%nproj,rr,vp,vlocpw, &
&                     vkbpw,evkb,mmax,mxprj,inp%zz,zion,zval,inp%iexc,rcmax, &
&                     rhoc, &
&                     rhotae,rhomod(1,1),rho,rhops,uupsa,eeig,vo,vxc, &
&                     eeel,eexc,vion,iout,ifail)
 if(ifail/=0) error stop 1

! total energy output

 epstot= eeig + eexc - 0.5d0*eeel
 write(iout,'(/a,f12.6/)') 'Pseudoatom total energy', epstot

 call run_diag(inp%lmax,npa,epa,wellprj,inp%lloc,irc, &
&               vlocpw,vkbpw,evkb,inp%nproj,rr,vfull,vp,inp%zz,mmax,mxprj, &
&               irel,nkap, &
&               0,vsrpw,esr,vsopw,eso,iout,ifail)
 if(ifail/=0) error stop 1

 if(irel==2) then
   call run_diag(inp%lmax,npa,epa,wellprj,inp%lloc,irc, &
&                 vlocpw,vkbpw,evkb,inp%nproj,rr,vfull,vp,inp%zz,mmax,mxprj, &
&                 irel, &
&                 nkap,1,vsrpw,esr,vsopw,eso,iout,ifail)
   if(ifail/=0) error stop 1

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
 if(ifail/=0) error stop 1

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
   if(ifail/=0) error stop 1
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
 do jj=1,inp%ncnf+1

   write(iout,'(/a,i2)') 'Test configuration',jj-1

! charge density is initialized to that of reference configuration

   rhot(:)=rho(:)

   if(istage/=2) then
     call run_config_ae(jj,inp%nacnf,inp%lacnf,inp%facnf,inp%nc,inp%nvcnf,rr, &
&                       inp%zz, &
&                       mmax,inp%iexc,ea,irel,nkap,istcnf,ierrcnf,itcnf, &
&                       kerrcnf,nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf, &
&                       iout,ifail)
     if(ifail/=0) error stop 1
   end if

   call run_config_ps(jj,istcnf,ierrcnf,itcnf,kerrcnf,nvtcnf, &
&                     natcnf,latcnf,fatcnf,eatcnf,etstcnf,inp%nc,rhot,rhomod, &
&                     rr, &
&                     rcmax,mmax,mxprj,inp%iexc,etot,epstot,inp%nproj,vpuns, &
&                     inp%lloc,vlocpw,vkbpw,evkb,irel,nkap,ntcf(jj), &
&                     itcf(:,:,jj),etcf(:,:,jj),dtcf(:,jj),iout,ifail)
   if(ifail/=0) error stop 1

 end do !jj

 call run_plot(inp%lmax,npa,epa,wellprj,inp%lloc,irc, &
&               vkb,vlocpw,vkbpw,evkb,inp%nproj,rr,vfull,vp,vpuns,inp%zz, &
&               mmax,mxprj, &
&               inp%drl,nrl, &
&               rho,rhoc,rhomod,irel,nkap,cvgplt,iout,ifail)
 if(ifail/=0) error stop 1



!the pseudopotential for oncvpsp_ld.x and oncvpsp_confine.x
 if(len_trim(ps_filename)/=0) then
   call write_ps_file(trim(ps_filename),irel,nkap,inp%atsym,inp%zz,mmax,rr, &
&                     inp%lmax, &
&                     inp%lloc,inp%nproj,inp%rc,inp%epsh1,inp%epsh2, &
&                     inp%depsh,inp%rxpsh,epa,vp,vlocpw, &
&                     evkb,vkbpw,inp%iexc,rcmax,vpuns(:,inp%lloc+1),rho, &
&                     rhomod(:,1),iout,ifail)
   if(ifail/=0) error stop 1
 end if

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
&                 inp%rxpsh,vlocpw,vkbpw, &
&                 evkb, &
&                 rr,vfull,vp,inp%zz,mmax,mxprj,inp%rc,irel,nkap,rphs,pshae, &
&                 pshps,iout,ifail)
   if(ifail/=0) error stop 1
 end if

 call gnu_script(epa,evkb,inp%lmax,inp%lloc,mxprj,inp%nproj,nkap,iout)

 call vloc_cutoff(rr,mmax,rho,vpuns(1,inp%lloc+1),zion,cvloc,vlchk)

 call diagout(inp%lmax,inp%nproj,mxprj,nkap,cvgplt,npsh,inp%epsh2,inp%depsh, &
&             rphs,pshae, &
&             pshps,ngh,ghn,ghe,ghp,(inp%icmod==6 .or. inp%icmod==7) &
&             .and. inp%fcfact<=0.0d0 .and. inp%budget_cut==1,inp%rcfact, &
&             ecutw,wcurve, &
&             ektail,inp%icmod,rr,mmax,rcmax,rhoc,corepw,inp%lloc,evkb,bbeig, &
&             inp%ncnf+1,ntcf,itcf,etcf,dtcf,cvloc,vlchk,optdg,wshare,iout)
 deallocate(pshae,pshps)

 if(trim(inp%psfile)=='psp8' .or. trim(inp%psfile)=='both') then
   soscale=1.0d0
   call linout(inp,vkbpw,evkb,vsrpw,esr,vsopw,eso,rr,vpuns,vlocpw,rho,corepw, &
&              rhotae,rhoc,zion,mmax,mxprj,nrl,soscale,rc0,fcfact0,fcgrid, &
&              rcgrid,log_mesh,nkap,iout,ifail)
   if(ifail/=0) error stop 1
 end if

 if(trim(inp%psfile)=='upf' .or. trim(inp%psfile)=='both') then
   call upfout(inp,vkbpw,evkb,rr,vpuns,vlocpw,rho,corepw,zion,mmax,mxprj,nrl, &
&              epstot,rc0,fcfact0,fcgrid,rcgrid,uupsa,ea,log_mesh,nkap, &
&              iout,ifail)
   if(ifail/=0) error stop 1
 end if


 if(trim(inp%psfile)=='psml' .or. trim(inp%psfile)=='both') then
!
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
   call psmlout(inp,vkb,evkb,mxprj,rr,vpsml,rho,rhomod,rcore,zion,mmax,nrl, &
&               epstot,rc0,fcfact0,fcgrid,rcgrid,log_mesh,irel,nkap,iout)
 end if

 stop

 end subroutine oncvpsp_main

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

 end module oncvpsp_main_m

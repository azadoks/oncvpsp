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
 use build_info_m, only: write_build_info
 use input_text_m, only: read_input_text
 use config_ae_m, only: run_config_ae, aeatom_cached
 use ae_file_m, only: write_ae_file, read_ae_file
 use ps_file_m, only: write_ps_file
#if (defined WITH_TOML)
 use input_toml_m, only: read_input_toml
#endif
 use constants_m, only: dp, mxprj, mxstate, mxcnf, mxl, mxlin, mxkap, nrhomod
 use pieces_m, only: piecewise_t, pw_reach, pw_free
 use modcore_m, only: core_pieces_smooth
 use free_tail_m, only: free_tail_log
 use well_mesh_m, only: well_mesh_t
 use band_ecut_m, only: band_ecut, nband_lev, band_lev, band_width, &
&                       band_vint
 use logmesh_m, only: lm_fh_reach
 use constants_m, only: nlev, nerrdg, nerrlv, well_onset_default
 use input_echo_m, only: write_atom, write_lmax, write_icmod, write_epsh, &
&                        write_mesh, fx, vkb_head, vkb_line
 use relativity_m, only: set_relativity
 implicit none

!Input variables
 integer, intent(in) :: irel,istage

!Local variables
 integer :: ii,ierr,iexc,iexct,ios,iprint,it,icmod,lpopt
 integer :: jj,kk,ll,l1,lloc,lmax,lt,inline
 integer :: mch,mchf,mmax,n1,n2,nc,nlim,nlloc,nlmax,irpsh,nrl
 integer :: nv,ncnf,nvt,upf_mesh
 integer :: iprj,nkap,nnwell,npsh,mfree,mchx
 integer :: ikap,kap,mkap
 integer,allocatable :: npa(:,:)
!
 integer :: dtime(8),na(mxstate),la(mxstate),np(mxlin+1)
 integer :: nacnf(mxstate,mxcnf),lacnf(mxstate,mxcnf),nvcnf(mxcnf)
 integer :: irc(mxlin+1),nodes(mxl+1)
 integer :: nproj(mxlin+1)
 integer :: ncon(mxlin+1),nbas(mxlin+1)
! the ghost tests' margins for the diagnostics block (run_ghosts), per kap
 integer :: ngh(mxl+1,mxkap),ghn(mxstate,mxl+1,mxkap)
 real(dp) :: ghe(2,mxstate,mxl+1,mxkap),ghp(2,mxl+1,mxkap)

 real(dp) :: al,amesh,csc,csc1,depsh,depsht,drl,eeel
 real(dp) :: eeig,eexc
 real(dp) :: emax,epsh1,epsh1t,epsh2,epsh2t,rxpsh
 real(dp) :: et,etest,emin,sls
 real(dp) :: fcfact,fcfact0,rcfact,slack,dvloc0,soscale
 real(dp) :: fcfact_min,fcfact_max,fcfact_step
 real(dp) :: rcfact_min,rcfact_max,rcfact_step
 real(dp) :: fcgrid(3),rcgrid(3)
 real(dp) :: rr1,rcion,rciont,rcmax,rct,rlmax,rpkt
 real(dp) :: amesh_in,rr1_in,rmax_in,rmax,log_mesh(3),c_light_in
 real(dp) :: sf,sw,zz,zion,zval,etot,res
 real(dp) :: cvloc(2,nerrlv),vlchk,wreg
 real(dp) :: eterm(4),rphs(4)
 real(dp) :: fj,cnorm,rcore,eptgt
!
 real(dp) :: cl(mxlin+1),debl(mxlin+1),ea(mxstate,mxkap),ep(mxlin+1,mxkap)
 real(dp) :: fa(mxstate)
 real(dp) :: facnf(mxstate,mxcnf)
 real(dp) :: qcut(mxlin+1),qmsbf(mxlin+1),rc(mxlin+1),rc0(mxlin+1)
 real(dp) :: rpk(mxstate,mxkap)
 real(dp) :: epstot
 real(dp), parameter :: eps=1.0d-8
 real(dp), external :: cinf_taper

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
 logical :: ep_as_given,well_tail_free
!> the cutoff of the pseudo-waves for icmod 6, 7 (budget_cutoff): 0 the
!> largest of the ground-state waves, 1 their weighted mean, 2 the largest
!> of all
 integer :: budget_cut
!> the tail of the well states (0 "well", 1 "free", 2 "confined"), and the
!> n of each l's first well state with "confined" (0 for the default)
 integer :: well_tail,well_n(mxlin+1)
!> the normalization of the well states: 0 that of well_tail, 1 all
!> space, 2 inside rc; and whether they are normalized inside rc
 integer :: well_norm
 logical :: well_norm_rc
!> the weight of the first pseudo-wave of each l with no occupied state in
!> the cutoff of the pseudo-waves (ecut_weight; negative for 2(2l+1))
 real(dp) :: ecut_wt(mxlin+1)
!> the cutoffs at the valence energies (band_ecut_m): the q taken as
!> infinity of each l, the energies, and per energy, l and kappa E_resid at
!> qcut and the cutoffs at the levels band_lev
 real(dp) :: qinf(mxlin+1),eband(3)
 real(dp), allocatable :: bandq(:,:,:),bandc(:,:,:,:)
 integer :: nband,ie
 real(dp) :: elo,ehi
 real(dp) :: well_onset
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

 character(len=2) :: atsym
 character(len=4) :: psfile

!names in the messages: the equation solved and the program
 character(len=19) :: treatment
 character(len=9) :: pname
 character(len=16) :: exe

 integer, parameter :: INPUT_STDIN=1, INPUT_TEXT=2, INPUT_TOML=3
 integer :: input_mode
 integer :: unit
 character(len=1024) :: input_filename
 character(len=1024) :: ae_filename
 character(len=1024) :: ps_filename
 logical :: skip_ld

 input_mode = INPUT_STDIN
 ae_filename=''
 ps_filename=''
 skip_ld=.false.

!the arrays indexed by l have a last index for kappa (see doc/fr_merge_plan.md):
!1 for j=l+1/2 and 2 for j=l-1/2 when fully relativistic, only 1 otherwise
 if(irel==2) then
   nkap=2
   treatment='fully-relativistic'
   pname='oncvpsp_r'
 else
   nkap=1
   if(irel==0) then
     treatment='non-relativistic'
   else
     treatment='scalar-relativistic'
   end if
   pname='oncvpsp'
 end if
 if(irel==0) then
   exe='oncvpspnr'
 else if(irel==1) then
   exe='oncvpsp'
 else
   exe='oncvpspr'
 end if
 if(istage==1) then
   exe=trim(exe)//'_ae'
 else if(istage==2) then
   exe=trim(exe)//'_ps'
 end if
 exe=trim(exe)//'.x'
 allocate(eatcnf(mxstate,nkap,nlev,mxcnf))

 parse_args: block
   integer :: iarg
   character(len=256) :: arg
   do iarg = 1, command_argument_count()
     call get_command_argument(iarg, arg)
     select case(arg)
     case('-h', '--help')
       call write_help()
       stop
     case('-v', '--version')
       write(6,'(3a)') 'ONCVPSP (',trim(treatment),') version 4.0.1'
       call write_build_info(6,trim(exe),trim(treatment))
       stop
     case('-i', '--input')
       if(iarg + 1 > command_argument_count()) then
         write(6,'(/3a)') trim(pname),': ERROR --input requires a file name'
         error stop 1
       end if
       call get_command_argument(iarg + 1, input_filename)
       input_mode = INPUT_TEXT
#if (defined WITH_TOML)
     case('-t', '--toml-input')
       if(iarg + 1 > command_argument_count()) then
         write(6,'(/3a)') trim(pname), &
&              ': ERROR --toml-input requires a file name'
         error stop 1
       end if
       call get_command_argument(iarg + 1, input_filename)
       input_mode = INPUT_TOML
#else
     case('-t', '--toml-input')
       write(6,'(/3a)') trim(pname),': ERROR TOML input is not enabled', &
&            ' in this build'
       error stop 1
#endif
     case('-a', '--ae-file')
!only the stages take an all-electron data file; the whole program ignores it
       if(istage/=0) then
         if(iarg + 1 > command_argument_count()) then
           write(6,'(/3a)') trim(pname), &
&                ': ERROR --ae-file requires a file name'
           error stop 1
         end if
         call get_command_argument(iarg + 1, ae_filename)
       end if
     case('-p', '--ps-file')
!the all-electron stage has no pseudopotential
       if(istage/=1) then
         if(iarg + 1 > command_argument_count()) then
           write(6,'(/3a)') trim(pname), &
&                ': ERROR --ps-file requires a file name'
           error stop 1
         end if
         call get_command_argument(iarg + 1, ps_filename)
       end if
     case('--no-log-derivatives')
       skip_ld=(istage/=1)
     case default
     end select
   end do
 end block parse_args

 if(irel==2) then
   write(6,'(a/a//)') &
&        'ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&        'relativistic version 4.0.1 03/01/2019'
 else
   write(6,'(a/a//)') &
&        'ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&        trim(treatment)//' version 4.0.1 03/01/2019'
 end if

 write(6,'(a/a/a//)') &
&      'While it is not required under the terms of the GNU GPL, it is',&
&      'suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)', &
&      'in any publication utilizing these pseudopotentials.'

 call write_build_info(6,trim(exe),trim(treatment))
 write(6,'(/)')


 if(istage/=0 .and. len_trim(ae_filename) == 0) then
   write(6,'(/4a)') trim(pname),': ERROR ',trim(exe),' requires --ae-file'
   error stop 1
 end if

 ios = 0
 ep_as_given=.false.
 well_onset=well_onset_default
 well_tail=0
 well_n(:)=0
 well_norm=0
 ecut_wt(:)=-1.0d0
 wreg=0.0d0
 select case(input_mode)
 case(INPUT_STDIN)
   call read_input_text(stdin, inline, &
&                       atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
&                       lmax, rc, ep, ncon, nbas, qcut, &
&                       lloc, lpopt, dvloc0, nproj, debl, &
&                       icmod, fcfact, rcfact, slack, budget_cut, &
&                       fcfact_min, fcfact_max, fcfact_step, &
&                       rcfact_min, rcfact_max, rcfact_step, &
&                       epsh1, epsh2, depsh, rxpsh, &
&                       rlmax, drl, upf_mesh, &
&                       ncnf, nvcnf, nacnf, lacnf, facnf, &
&                       amesh_in, rr1_in, rmax_in, c_light_in, ep_as_given, &
&                       well_onset, well_tail, wreg, well_n, well_norm, &
&                       ecut_wt)
 case(INPUT_TEXT)
   open(newunit=unit, file=input_filename, status='old', action='read', &
&       iostat=ios)
   if(ios /= 0) then
     write(6, '(/3a)') trim(pname),': ERROR cannot open input file, file=', &
&       trim(input_filename)
     error stop 1
   end if
   call read_input_text(unit, inline, &
&                       atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
&                       lmax, rc, ep, ncon, nbas, qcut, &
&                       lloc, lpopt, dvloc0, nproj, debl, &
&                       icmod, fcfact, rcfact, slack, budget_cut, &
&                       fcfact_min, fcfact_max, fcfact_step, &
&                       rcfact_min, rcfact_max, rcfact_step, &
&                       epsh1, epsh2, depsh, rxpsh, &
&                       rlmax, drl, upf_mesh, &
&                       ncnf, nvcnf, nacnf, lacnf, facnf, &
&                       amesh_in, rr1_in, rmax_in, c_light_in, ep_as_given, &
&                       well_onset, well_tail, wreg, well_n, well_norm, &
&                       ecut_wt)
   close(unit)
#if (defined WITH_TOML)
 case(INPUT_TOML)
   open(newunit=unit, file=input_filename, status='old', action='read', &
&       iostat=ios)
   if(ios /= 0) then
     write(6, '(/3a)') trim(pname),': ERROR cannot open input file, file=', &
&       trim(input_filename)
     error stop 1
   end if
   call read_input_toml(unit, &
&                       atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
&                       lmax, rc, ep, ncon, nbas, qcut, &
&                       lloc, lpopt, dvloc0, nproj, debl, &
&                       icmod, fcfact, rcfact, slack, budget_cut, &
&                       fcfact_min, fcfact_max, fcfact_step, &
&                       rcfact_min, rcfact_max, rcfact_step, &
&                       epsh1, epsh2, depsh, rxpsh, &
&                       rlmax, drl, upf_mesh, &
&                       ncnf, nvcnf, nacnf, lacnf, facnf, &
&                       amesh_in, rr1_in, rmax_in, c_light_in, ep_as_given, &
&                       well_onset, well_tail, wreg, well_n, well_norm, &
&                       ecut_wt)
   close(unit)
#endif
 case default
   write(6,'(/2a,i0)') trim(pname),': ERROR unknown input mode, input_mode=', &
&        input_mode
   error stop 1
 end select
 well_tail_free=well_tail==1
!the free tail has no finite norm over all space
 if(well_tail==1 .and. well_norm==1) then
   write(6,'(/a)') 'oncvpsp: ERROR well_norm = "space" (1) is not defined '// &
&       'for well_tail = "free": the free tail has no finite norm'
   error stop 1
 end if
 well_norm_rc=well_norm==2 .or. (well_norm==0 .and. well_tail/=0)
 if(well_tail/=2 .and. any(well_n(:)>0)) then
   write(6,'(/a)') 'oncvpsp: ERROR well_n is for well_tail = "confined" only'
   error stop 1
 end if

 nvcnf(1)=nv
 do ii=1,nc+nv
   nacnf(ii,1)=na(ii)
   lacnf(ii,1)=la(ii)
   facnf(ii,1)=fa(ii)
 end do

 do jj=2,ncnf+1
   do ii=1,nc
     nacnf(ii,jj)=na(ii)
     lacnf(ii,jj)=la(ii)
     facnf(ii,jj)=fa(ii)
   end do
 end do

 if(icmod==0) then
   fcfact=0.0d0
 end if

!icmod=4 coarse-search grids: min, max, step
 fcgrid=[fcfact_min,fcfact_max,fcfact_step]
 rcgrid=[rcfact_min,rcfact_max,rcfact_step]

 call check_data(atsym,zz,fcfact,rcfact,slack,fcgrid,rcgrid,epsh1,epsh2,depsh, &
&                rxpsh,rlmax,drl,upf_mesh, &
&                fa,facnf, &
&                rc,ep,qcut,debl,nc,nv,iexc,lmax,lloc,lpopt,icmod, &
&                ncnf,na,la,nvcnf,nacnf,lacnf,ncon,nbas,nproj,psfile, &
&                amesh_in,rr1_in,rmax_in,c_light_in)
 call set_relativity(c_light_in)

 nrl=int((rlmax/drl)-0.5d0)+1

!PWSCF wants an even number of mesh pointe (not made so when fully
!relativistic, see doc/fr_merge_plan.md)
 if(mod(nrl,2)/=0 .and. irel/=2) nrl=nrl+1

!log mesh: the input values, or the built-in ones for those 0 or less
 if(amesh_in>0.0d0) then
   amesh=amesh_in
 else
!Alternative built-in meshes: amesh=1.012d0 (coarser), or 1.003d0 and
!1.0015d0 (finer)
   amesh=1.006d0
 end if

 al=log(amesh)
 if(rr1_in>0.0d0) then
   rr1=rr1_in
 else
   rr1=0.0005d0/zz
   if(irel/=2) rr1=min(rr1,0.0005d0/10)
 end if
 if(rmax_in>0.0d0) then
   rmax=rmax_in
 else
   rmax=45.0d0
 end if
 mmax=log(rmax/rr1)/al
!the log mesh for the echoes of the input, only when the input sets it, so
!that the output of an input without it is unchanged
 log_mesh(:)=0.0d0
 if(amesh_in>0.0d0 .or. rr1_in>0.0d0 .or. rmax_in>0.0d0) then
   log_mesh(:)=[amesh,rr1,rmax]
 end if

!calculate zion for output
 zion=zz
 do ii=1,nc
   zion=zion-fa(ii)
 end do


 allocate(rr(mmax))
 allocate(rho(mmax),rhoc(mmax),rhot(mmax))
 allocate(uu(mmax,2),up(mmax,2),uupsa(mmax,nkap,nv))
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
 allocate(rhoae(mmax,nv),rhops(mmax,nv),rhotae(mmax))
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

 do ii=1,mmax
   rr(ii)=rr1*exp(al*(ii-1))
 end do

!
! full potential atom solution
!
 if(istage==2) then
   call read_ae_file(ae_filename,irel,nkap,atsym,zz,nc,nv,iexc,na,la,fa, &
&                    mmax,rr,ncnf,nvcnf,nacnf,lacnf,facnf, &
&                    it,etot,eterm,res,ea,rhoc,rho,vfull, &
&                    istcnf,ierrcnf,itcnf, &
&                    kerrcnf,nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf)
!aeatom prints the libxc functional description on its first exc_libxc call
   if(iexc<0) then
     block
       real(dp) :: vxct(mmax),exct(mmax)
       call exc_libxc(iexc,al,rho,vxct,exct,rr,mmax)
     end block
   end if
 else
   ierr=0
   call aeatom_cached(irel,nkap,na,la,ea,fa,rpk,nc,nc+nv,it,rhoc,rho, &
&                     rr,vfull,zz,mmax,iexc,etot,eterm,res,ierr)
   if(ierr/=0) then
     write(6,'(/a,i0,a,i0,a,es12.4)') &
&          'oncvpsp: ERROR aeatom failed for the reference configuration'// &
&          ', ierr=',ierr,', it=',it,', res=',res
     error stop 1
   end if
 end if
!
!

! the core radii are used as given (doc/smoothness_plan.md, R-1)
 rcmax=0.0d0
 do l1=1,lmax+1
   rcmax=max(rcmax,rc(l1))
 end do
 do l1=1,lmax+1
   if(rc(l1)==0.0d0) then
     rc(l1)=rcmax
   end if
 end do

 nproj(lloc+1)=0
 rc0(:)=rc(:)
!input value, for the echoes (modcore3 resets fcfact for icmod=4)
 fcfact0=fcfact

! output printing (echos input data, with all-electron eigenvalues added)

 call write_atom(atsym,zz,nc,nv,iexc,psfile)
 write(6,'(a/a)') '#','#   n    l    f        energy (Ha)'
 do ii=1,nc+nv
   if(irel/=2) then
     write(6,'(2i5,f8.2,1pe18.7)') na(ii),la(ii),fa(ii),ea(ii,1)
   else if(la(ii)==0) then
     write(6,'(2i5,f8.2,1p,d18.7)') na(ii),la(ii),fa(ii),ea(ii,1)
     ea(ii,2)=ea(ii,1)
   else
     write(6,'(2i5,f8.2,1p,2d18.7)') na(ii),la(ii),fa(ii),ea(ii,1),ea(ii,2)
   end if
 end do

 call write_lmax(lmax,ep_as_given,well_onset,well_tail,wreg,well_norm)
 write(6,'(a/a)') '#','#   l,   rc,      ep,       ncon, nbas, qcut'
 do l1=1,lmax+1
   write(6,'(i5,2a,2i5,a)') l1-1,fx(rc(l1),10,5),fx(ep(l1,1),10,5), &
&        ncon(l1),nbas(l1),fx(qcut(l1),10,5)
 end do

 write(6,'(a/a/a,a)') '#','# LOCAL POTENTIAL','# lloc, lpopt,  rc(5),', &
&      '   dvloc0'
 write(6,'(2i5,3a)') lloc,lpopt,fx(rc(5),10,5),'   ',fx(dvloc0,10,5)

 write(6,'(a/a/a)') '#','# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', &
&      trim(vkb_head(lmax,well_n,ecut_wt))
 do l1=1,lmax+1
   write(6,'(a)') trim(vkb_line(lmax,l1,nproj,debl,well_n,ecut_wt))
 end do

 call write_icmod(icmod,fcfact,rcfact,slack,budget_cut,fcgrid,rcgrid)

 call write_epsh(epsh1,epsh2,depsh,rxpsh)

 if(upf_mesh==0) then
   write(6,'(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl'
   write(6,'(2a)') fx(rlmax,8,4),fx(drl,8,4)
 else
   write(6,'(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl, upf_mesh'
   write(6,'(2a,i8)') fx(rlmax,8,4),fx(drl,8,4),upf_mesh
 end if

 write(6,'(a/a/a)') '#','# TEST CONFIGURATIONS','# ncnf'
 write(6,'(i5)') ncnf
 write(6,'(a/a)') '# nvcnf','#   n    l    f'
 do jj=2,ncnf+1
   write(6,'(i5)') nvcnf(jj)
   do ii=nc+1,nc+nvcnf(jj)
     write(6,'(2i5,f8.2)') nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
   end do
   write(6,'(a)') '#'
 end do

 call write_mesh(log_mesh)

 write(6,'(//a)') 'Reference configufation results'
 write(6,'(a,i6)') '  iterations',it
 if(it >= 100) then
   write(6,'(2a)') trim(pname), &
&        ': ERROR all-electron reference atom not converged'
   error stop 1
 end if
 write(6,'(a,1p,d18.8)') '  all-electron total energy (Ha)',etot
 write(6,'(a,1p,e10.2)') '  residual of the screening potential (Ha)',res
 write(6,'(a)') '  energy terms (Ha)'
 write(6,'(a,f22.10)') '    total               ',etot
 write(6,'(a,f22.10)') '    kinetic             ',eterm(1)
 write(6,'(a,f22.10)') '    electron-nucleus    ',eterm(2)
 write(6,'(a,f22.10)') '    Hartree             ',eterm(3)
 write(6,'(a,f22.10)') '    exchange-correlation',eterm(4)

 if(istage==1) then
!all-electron halves of the test configurations; their pseudo-atom halves
!and all output about them come from the matching _ps.x
   do jj=1,ncnf+1
     call run_config_ae(jj,nacnf,lacnf,facnf,nc,nvcnf,rr,zz, &
&                       mmax,iexc,ea,irel,nkap,istcnf,ierrcnf,itcnf, &
&                       kerrcnf,nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf)
   end do

   call write_ae_file(ae_filename,irel,nkap,atsym,zz,nc,nv,iexc,na,la,fa, &
&                     mmax,rr,it,etot,eterm,res,ea,rhoc,rho,vfull, &
&                     ncnf,nvcnf,nacnf,lacnf,facnf,istcnf,ierrcnf,itcnf, &
&                     kerrcnf,nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf)

   write(6,'(/2a)') 'All-electron data written to ',trim(ae_filename)
   stop
 end if

!the last log mesh point at or inside each rc, for the bookkeeping of the
!routines that index the mesh (R-2); the radii themselves stay as given
 rcmax=0.0d0
 irc(:)=0
 do l1=1,max(lmax+1,lloc+1)
   irc(l1)=0
   do ii=1,mmax
     if(rr(ii)<=rc(l1)) irc(l1)=ii
   end do
   if(rc(l1)>0.0d0 .and. irc(l1)<2) then
     write(6,'(a,i2,a,es12.4)') &
&         'oncvpsp: ERROR rc inside the mesh start for l=',l1-1,': rc=',rc(l1)
     error stop 1
   end if
   rcmax=max(rcmax,rc(l1))
 end do

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
 write(6,'(/a/a)') 'Begin loop to  construct optimized pseudo wave functions',&
&      'and semi-local pseudopoentials for all angular momenta'

!temporarily set this to 1 so that the pseudo wave function needed for the
!local potential will be generated.  Reset after run_vkb.
 npa(:,:)=0
 nproj(lloc+1)=1
 do l1=1,lmax+1
   vwell(:,:)=0.0d0
   do jj=1,mxprj
     wmesh(jj)%mx=0
   end do
   ll=l1-1
   if(nkap==1 .or. ll==0) then
     mkap=1
   else
     mkap=2
   end if
! loop on J = ll +/- 1/2 when nkap=2; fj is the weight of each J
   do ikap=1,mkap
     if(nkap==1) then
       kap=-(ll+1)
       fj=1.0d0
     else if(ikap==1) then
       kap=-(ll+1)
       fj=dble(ll+1)/dble(2*ll+1)
     else
       kap=  ll
       fj=dble(ll)/dble(2*ll+1)
     end if

     uu(:,:)=0.0d0; qq(:,:)=0.0d0
     iprj=0
     if(ikap==1) inwell(:)=.false.

!get principal quantum number for the highest core state for this l
     npa(1,l1)=l1
     do kk=1,nc
       if(la(kk)==l1-1) npa(1,l1)=na(kk)+1
     end do !kk

!get all-electron bound states for projectors
     if(nv/=0) then
       do kk=nc+1,nc+nv
         if(la(kk)==l1-1) then
           iprj=iprj+1
           et=ea(kk,ikap)
           call boundstate(irel,na(kk),ll,kap,ierr,et, &
&                          rr,zz,vfull,uu,up,mmax,mch)
           if(ierr /= 0) then
             if(irel==2) then
               write(6,'(a,i2,a,i2,a,i2)') &
&                    'oncvpsp_r 409: ldiracfb convergence ERROR, n=',&
&                     na(kk),'  l=',ll, 'kap-',kap
             else
               write(6,'(/a,3i4)') &
&                    'oncvpsp-387: lschfb convergence ERROR n,l,iter=', &
&                    na(kk),la(kk),it
             end if
             error stop 1
           end if

           epa(iprj,l1,ikap)=ea(kk,ikap)
           npa(iprj,l1)=na(kk)
           if(irel==2) then
             call renorm_r(uu,rr,ll,kap,zz,mmax,cnorm)
           else
             cnorm=1.0d0
           end if
           uua(:,iprj,ikap)=uu(:,1)
           upa(:,iprj,ikap)=cnorm*up(:,1)
         end if !la(kk)==l1-1
         if(iprj==nproj(l1)) exit
       end do !kk
     end if !nv/=0

!get all-electron well states for projectors
!if there were no valence states, use ep from input data for 1st well state
!otherwise shift up by input debl

     if(iprj==0) epa(1,l1,ikap)=ep(l1,ikap)
     if(iprj<nproj(l1))then
       do kk=1,nproj(l1)-iprj
         iprj=iprj+1
         if(iprj>1 .and. debl(l1)<=0.0d0) then
           if(irel==2) then
             write(6,'(a,f8.3,a,i2/a/a)') 'oncvpsp: ERROR debl =',debl(l1), &
&                  ' for l=',ll, &
&                  ' not allowed with 2 or more scattering states', &
&                  'program will stop'
           else
             write(6,'(a,f8.3,a,i2/a/a)') 'oncvpsp: ERROR debl =',debl(l1), &
&                  ' for l=',ll, &
&                  ' ERROR not allowed with 2 or more scattering states', &
&                  'program will stop'
           end if
           error stop 1
         end if

! npa is not reset for ikap=2 to retain same node count for consistency with
! keeping same well potential as for ikap=1.
         if(iprj>1) then
           epa(iprj,l1,ikap)=epa(iprj-1,l1,ikap)+debl(l1)
           if(ikap==1) then
             npa(iprj,l1)=npa(iprj-1,l1)+1
           end if
         end if
!with the confined well, the n of the channel's first well state is the
!input's well_n, if set; the next ones follow by one each
         if(well_tail==2 .and. kk==1 .and. ikap==1 .and. well_n(l1)>0) then
           if(well_n(l1)<npa(iprj,l1)) then
             write(6,'(/a,i2,a,i3,a,i3)') 'oncvpsp: ERROR well_n for l=',ll, &
&                 ' is',well_n(l1),', below the first well state''s n,', &
&                 npa(iprj,l1)
             error stop 1
           end if
           npa(iprj,l1)=well_n(l1)
         end if

! the well potential is calculated for the less-bound ikap=1 state and kept
! for ikap=2, which is more deeply bound and should always have a bound
! state for this combination of npa and vwell.  For irel=2 npa takes the n
! of the well state found, which can have more nodes than asked for; else
! npa keeps the input n.
         if(ikap==1 .and. irel==2) then
           call wellstate(npa(iprj,l1),ll,kap,rc(l1),epa(iprj,l1,ikap), &
&                         ep_as_given,well_onset,well_tail==2,.false.,rr, &
&                         vfull, &
&                         vwell(1,iprj),uu,up,zz,mmax,mch,irel, &
&                         inwell(iprj),eptgt,wmesh(iprj))
         else if(ikap==1) then
           nnwell=npa(iprj,l1)
           call wellstate(nnwell,ll,kap,rc(l1),epa(iprj,l1,ikap),ep_as_given, &
&                         well_onset,well_tail==2,well_tail_free,rr,vfull, &
&                         vwell(1,iprj),uu,up,zz,mmax,mch,irel, &
&                         inwell(iprj),eptgt,wmesh(iprj))
         else
           et=epa(iprj,l1,ikap)
           if(wmesh(iprj)%mx>0) then
!the j = l - 1/2 state in the well of j = l + 1/2, on its continued mesh
             wmesh(iprj)%uu(:,:)=0.0d0
             wmesh(iprj)%up(:,:)=0.0d0
             call boundstate(irel,npa(iprj,l1),ll,kap,ierr, &
&                            epa(iprj,l1,ikap),wmesh(iprj)%rr,zz, &
&                            wmesh(iprj)%vv,wmesh(iprj)%uu,wmesh(iprj)%up, &
&                            wmesh(iprj)%mx,mchx)
             uu(:,:)=wmesh(iprj)%uu(1:mmax,:)
             up(:,:)=wmesh(iprj)%up(1:mmax,:)
           else
             call boundstate(irel,npa(iprj,l1),ll,kap,ierr, &
&                            epa(iprj,l1,ikap),rr,zz,vwell(1,iprj),uu,up, &
&                            mmax,mch)
           end if

           if(ierr /= 0) then
             write(6,'(a,i2,a,i2,a,i3,a,f12.6)') &
&                  'oncvpsp_r 474: ldiracfb convergence ERROR, n=',&
&                   npa(iprj,l1),'  l=',ll, 'kap=',kap, ' epa=', &
&                   epa(iprj,l1,ikap)

             error stop 1
           end if

           write(6,'(/a,i3,a,i3,a,i3/a,f8.4)') &
&                   '   Wellstate for l =',ll,'  n =',npa(iprj,l1), &
&                   ' kap=',kap, &
&                   '     eigenvalue = ',et
         end if

         if(irel==2 .and. wmesh(iprj)%mx>0) then
           call renorm_r(wmesh(iprj)%uu,wmesh(iprj)%rr,ll,kap,zz, &
&                        wmesh(iprj)%mx,cnorm)
           uu(:,:)=wmesh(iprj)%uu(1:mmax,:)
         else if(irel==2) then
           call renorm_r(uu,rr,ll,kap,zz,mmax,cnorm)
         else
           cnorm=1.0d0
         end if
         uua(:,iprj,ikap)=uu(:,1)
         upa(:,iprj,ikap)=cnorm*up(:,1)
!a state the well confines is normalized inside rc with well_norm "rc",
!and by default with the confined well, as with the free tail: its norm
!over all space is the well's, and would scale its E_resid with the width
!of the well; else over all space, as the solvers give it (the kernels
!see neither)
         if(well_norm_rc .and. wmesh(iprj)%mx>0) then
           call fpovlp(uua(1,iprj,ikap),uua(1,iprj,ikap),mmax,rc(l1),ll, &
&                      kap,zz,cnorm,rr,irel)
           cnorm=1.0d0/sqrt(cnorm)
           uua(:,iprj,ikap)=cnorm*uua(:,iprj,ikap)
           upa(:,iprj,ikap)=cnorm*upa(:,iprj,ikap)
           wmesh(iprj)%uu(:,:)=cnorm*wmesh(iprj)%uu(:,:)
           wmesh(iprj)%up(:,:)=cnorm*wmesh(iprj)%up(:,:)
         end if
!with the free tail, a well state is the atom's solution at exactly the
!energy asked for (the well's eigenvalue only for j = l - 1/2, D-7c),
!integrated outward past the reach of everything taken at or near rc
!(ae_delta's nodes below, and the interpolants'); normalized inside rc,
!as its norm over all space is the well's; and the windowed free tail
!beyond, so that nothing sees the well but that energy
!(doc/smoothness_plan.md, phase 8)
         if(well_tail_free .and. inwell(iprj)) then
           if(ikap==1) epa(iprj,l1,ikap)=eptgt
           mfree=min(mmax,irc(l1)+pw_reach+16+lm_fh_reach()+2)
           uu(:,:)=0.0d0
           up(:,:)=0.0d0
           call scatterstate(irel,nnwell,ll,kap,ierr,epa(iprj,l1,ikap),rr, &
&                            zz,vfull,uu,up,mmax,mfree)
           if(ierr/=0) then
             write(6,'(a,i2,a,f12.8)') &
&                  'oncvpsp: ERROR scattering state failed for l=',ll, &
&                  ' at energy',epa(iprj,l1,ikap)
             error stop 1
           end if
           uua(:,iprj,ikap)=uu(:,1)
           upa(:,iprj,ikap)=up(:,1)
           call fpovlp(uua(1,iprj,ikap),uua(1,iprj,ikap),mmax,rc(l1),ll, &
&                      kap,zz,cnorm,rr,irel)
           uua(:,iprj,ikap)=uua(:,iprj,ikap)/sqrt(cnorm)
           upa(:,iprj,ikap)=upa(:,iprj,ikap)/sqrt(cnorm)
           call free_tail_log(irel,ll,kap,epa(iprj,l1,ikap),rr,vfull, &
&                             uua(1,iprj,ikap),upa(1,iprj,ikap),mmax, &
&                             rc(l1),qcut(l1),mfree)
         end if
       end do !kk
     end if !iprj<nproj(l1)

     do iprj=1,nproj(l1)

!the relativistic term Delta of the all-electron wave about the taper
!window [0.95 rc, rc] (ae_delta, D-4), for the projectors; the local
!potential of lloc<4 takes Delta/(2u) on the window (S-1)
       ii=max(1,1+int(log(0.95d0*rc(l1)/rr(1))/al)-pw_reach-8)
       jj=min(mmax,irc(l1)+pw_reach+16)
       call ae_delta(irel,ll,kap,epa(iprj,l1,ikap),rr,vfull,uua(1,iprj,ikap), &
&                    upa(1,iprj,ikap),mmax,ii,jj,dlta(1,iprj,l1,ikap))
       if(l1==lloc+1 .and. iprj==1) then
         do kk=ii,jj
           if(uua(kk,1,ikap)==0.0d0) then
             write(6,'(a,i2,a,f10.6)') &
&                  'oncvpsp: ERROR node of the local channel l=',ll, &
&                  ' all-electron wave near rc, at r=',rr(kk)
             error stop 1
           end if
           vrloc(kk,ikap)=0.5d0*dlta(kk,1,l1,ikap)/uua(kk,1,ikap)
         end do
       end if
     end do


!get all-electron overlap matrix
     do jj=1,nproj(l1)
       do ii=1,jj
         call fpovlp(uua(1,ii,ikap),uua(1,jj,ikap),mmax,rc(l1),ll,kap, &
&                    zz,qq(ii,jj),rr,irel)
         qq(jj,ii)=qq(ii,jj)
       end do
     end do
     call run_optimize(epa(1,l1,ikap),ll,kap,irel,mmax,mxprj,rr, &
&                      uua(1,1,ikap),upa(1,1,ikap),qq, &
&                      well_tail_free .and. inwell(:),wmesh,rc(l1), &
&                      irc(l1), &
&                      qcut(l1),qmsbf(l1),ncon(l1),nbas(l1),nproj(l1), &
&                      pswf(1,1,l1,ikap),vp(1,l1,ikap),vkb(1,1,l1,ikap), &
&                      vfull,cvgplt(1,1,1,l1,ikap),ncurve, &
&                      ecurve(1,1,1,l1,ikap), &
&                      pswfx(1,1,l1,ikap),chix(1,1,l1,ikap), &
&                      ektail(1,1,l1,ikap),wreg,optdg(1,1,l1,ikap), &
&                      wshare(1,1,l1,ikap),qinf(l1))

!the weights of its pseudo-waves in the cutoff of the pseudo-waves (the
!weighted mean, and which waves the largest is over): the occupation of
!each valence state (its share fj of J), or for the first wave of a channel
!with none occupied ecut_weight, by default 2(2l+1), as if the channel were
!full; in a solid it holds about the charge that hybridization puts inside
!rc, an input as the atom cannot give it
     jj=0
     do kk=nc+1,nc+nv
       if(la(kk)==ll .and. jj<nproj(l1)) then
         jj=jj+1
         wcurve(jj,l1,ikap)=fj*fa(kk)
       end if
     end do
     if(.not. any(wcurve(:,l1,ikap)>0.0d0)) then
       if(ecut_wt(l1)>=0.0d0) then
         wcurve(1,l1,ikap)=fj*ecut_wt(l1)
       else
         wcurve(1,l1,ikap)=fj*dble(2*(2*ll+1))
       end if
     else if(ecut_wt(l1)>=0.0d0) then
       write(6,'(/a,i2,a)') 'oncvpsp: ERROR ecut_weight for l=',ll, &
&           ' is for a channel with no occupied state; this one has one'
       error stop 1
     end if

   end do !ikap
   wellprj(:,l1)=inwell(:)
 end do !l1

! construct Vanderbilt / Kleinman-Bylander projectors

 write(6,'(/a,a)') 'Construct Vanderbilt / Kleinmman-Bylander projectors'

 call run_vkb(lmax,lloc,lpopt,dvloc0,rc,irc,nproj,rr,mmax,mxprj,pswf,pswfx, &
&             chix,vfull,vp,evkb,vkb,nlim,dlta,vrloc,irel,nkap,vkbpw, &
&             vlocpw,bbeig)

!restore this to its proper value
 nproj(lloc+1)=0

 deallocate(uua,upa)

! accumulate charge and eigenvalues
! pseudo wave functions are calculated with VKB projectors for
! maximum consistency of unscreening
! get all-electron and pseudopotential valence-state by valence-state
! charge densities

! null charge and eigenvalue accumulators
 uupsa(:,:,:)=0.0d0
 eeig=0.0d0
 zval=0.0d0
 rhops(:,:)=0.0d0
 rho(:)=0.0d0
 nodes(:)=0
 rhoae(:,:)=0.0d0
 rhotae(:)=0.0d0
 do kk=1,nv
   ll=la(nc+kk)
   l1=ll+1

   if(nkap==1 .or. ll==0) then
     mkap=1
   else
     mkap=2
   end if
   do ikap=1,mkap
     if(nkap==1) then
       kap=-(ll+1)
       fj=1.0d0
     else if(ikap==1) then
       kap=-(ll+1)
       fj=dble(ll+1)/dble(2*ll+1)
     else
       kap=  ll
       fj=dble(ll)/dble(2*ll+1)
     end if

     et=ea(nc+kk,ikap)

     call boundstate(irel,na(nc+kk),ll,kap,ierr,et, &
&                    rr,zz,vfull,uu,up,mmax,mch)
     if(ierr /= 0) then
       if(irel==2) then
         write(6,'(a,i2,a,i2,a,i2,a,i4)') &
&              'oncvpsp_r 541: ldiracfb convergence ERROR, n=',&
&               na(nc+kk),'  l=',ll, '  kap-',kap,'  ierr=',ierr
       else
         write(6,'(/a,3i4)') &
&              'oncvpsp-461: lschfb convergence ERROR n,l,iter=', &
&              na(nc+kk),la(nc+kk),it
       end if
       error stop 1
     end if

     if(irel==2) then
       rhoae(:,kk)=rhoae(:,kk)+fj*((uu(:,1)/rr(:))**2 + (uu(:,2)/rr(:))**2)
       emax=0.90d0*et
       emin=1.10d0*et
     else
       rhoae(:,kk)=(uu(:,1)/rr(:))**2
       emax=0.75d0*et
       emin=1.25d0*et
     end if

     call lschvkbb(ll+nodes(l1)+1,ll,nproj(l1),ierr,et,emin,emax, &
&                  rr,vp(1,lloc+1,ikap),vlocpw,vkbpw(1:nproj(l1),l1,ikap), &
&                  evkb(1,l1,ikap),uu,up,mmax,mch)
     if(ierr /= 0) then
       if(irel==2) then
         write(6,'(a,i2,a,i2)') 'oncvpsp_r: lschvkbb convergence ERROR, n=',&
&         ll+1,'  l=',ll
       else
         write(6,'(a,3i4)') 'oncvpsp: lschvkbb ERROR',ll+nodes(l1)+1,ll,ierr
         flush(6)
       end if
       error stop 1
     end if

! save valence pseudo wave functions for upfout
     uupsa(:,ikap,kk)=uu(:,1)

     if(irel==2) then
       rhops(:,kk)=rhops(:,kk)+fj*(uu(:,1)/rr(:))**2
     else
       rhops(:,kk)=(uu(:,1)/rr(:))**2
     end if
     eeig=eeig+fa(nc+kk)*fj*et

   end do !ikap
   rhotae(:)=rhotae(:)+fa(nc+kk)*rhoae(:,kk)
   rho(:)=rho(:)+fa(nc+kk)*rhops(:,kk)
   zval=zval+fa(nc+kk)
   nodes(l1)=nodes(l1)+1
 end do !kk

 if(irel==2) then
! decomposition into scalar-relativistic and spin-orbit projectors

   call sr_so_r(lmax,nproj,rr,mmax,mxprj,evkb,vkbpw, &
&               vsr,esr,vso,eso,vsrpw,vsopw)
 end if

 allocate(rhomod(mmax,nrhomod))

 rhomod(:,:)=0.0d0

! construct model core charge based on monotonic polynomial fit
! or Teter function fit

 if(icmod==1) then
   call modcore(icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod,corepw, &
&               fcfact,rcfact,rcmax,mmax,rr,nc,nv,la,zion,iexc)

 else if(icmod==2) then
   call modcore2(icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
&               fcfact,rcfact,mmax,rr,nc,nv,la,zion,iexc)

 else if(icmod==6 .or. icmod==7) then
   ecutw=0.0d0
   if(fcfact<=0.0d0) then
     call modcore_budget_table(ecurve,ncurve,mxprj*4*nkap,wcurve, &
&                              budget_cut==1,budget_cut==2,rcfact,ecutw)
   end if
!the probe densities of the hardness of icmod=7
   allocate(probes(mmax,mxprj*(lmax+1)),lprobe(mxprj*(lmax+1)))
   nprobe=0
   if(icmod==7) then
     call modcore_probes(lmax,lloc,nproj,epa,wellprj,vp,vlocpw,vkbpw,evkb, &
&                        rr,mmax,mxprj, &
&                        nkap,probes,lprobe,nprobe)
   end if
   call modcore_budget(icmod,rho,rhoc,rhomod,corepw,fcfact,ecutw, &
&               budget_cut==1,rcfact, &
&               slack, &
&               minval(rc(1:lmax+1)),mmax,rr,iexc,zion, &
&               probes(:,1:nprobe),lprobe(1:nprobe),nprobe)
   deallocate(probes,lprobe)

 else if(icmod>=3) then
   call modcore3(icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
&               fcfact,rcfact,fcgrid,rcgrid,mmax,rr,nc,nv,la,zion,iexc)

 end if
! the pieces of a model with no joint (the Teter function, or none)
 if(corepw(1)%nseg==0) call core_pieces_smooth(rr,rhomod,mmax,corepw)

! screening potential for pseudocharge

 call vout(1,rho,rhomod(1,1),vo,vxc,zval,eeel,eexc,rr,mmax,iexc)

! switch the screened potentials to the ionic tail beyond rcmax, iterating
! the pseudo valence states with it to a fixed point (after the model core,
! which V_ion uses)

 allocate(vion(mmax))
 call ionic_screening(nkap,nc,nv,la,fa,ea,lmax,lloc,nproj,rr,vp,vlocpw, &
&                     vkbpw,evkb,mmax,mxprj,zz,zion,zval,iexc,rcmax,rhoc, &
&                     rhotae,rhomod(1,1),rho,rhops,uupsa,eeig,vo,vxc, &
&                     eeel,eexc,vion)

! total energy output

 epstot= eeig + eexc - 0.5d0*eeel
 write(6,'(/a,f12.6/)') 'Pseudoatom total energy', epstot

 call run_diag(lmax,npa,epa,wellprj,lloc,irc, &
&               vlocpw,vkbpw,evkb,nproj,rr,vfull,vp,zz,mmax,mxprj,irel,nkap)

!the cutoff of each channel at the valence energies, the top band_width
!of the occupied valence eigenvalues (from the highest down to the lowest
!or band_width below it, a valence band's width, which leaves out the
!semicore states): its pseudo scattering
!state there, normalized to 1 inside rc, with the free tail
!(band_ecut_m), for each channel with no occupied state (Si d): the
!measure of the cutoff its charge in a solid needs, which the convergence
!tables of its well states, at energies the input chose, are not (a
!channel with an occupied state has it in the table of that state)
 allocate(bandq(3,mxlin+1,nkap),bandc(nband_lev,3,mxlin+1,nkap))
 bandq(:,:,:)=0.0d0
 bandc(:,:,:,:)=0.0d0
 nband=0
 if(nv>0) then
   elo=minval(ea(nc+1:nc+nv,1))
   ehi=elo
   do kk=nc+1,nc+nv
     if(fa(kk)>0.0d0) ehi=max(ehi,ea(kk,1))
   end do
   elo=max(elo,ehi-band_width)
   if(ehi-elo>1.0d-6) then
     nband=3
     eband(:)=elo+(ehi-elo)*[0.0d0,0.5d0,1.0d0]
   else
     nband=1
     eband(1)=elo
   end if
   write(6,'(/a/a/a)') &
&       'Cutoffs at the valence energies: the pseudo scattering state of', &
&       'each l with no occupied state, normalized to 1 inside rc', &
&       '                e (Ha)  E_resid(qcut)  Ecut (Ha) at 1e-2 ... 1e-5'
   do l1=1,lmax+1
     ll=l1-1
     if(any(la(nc+1:nc+nv)==ll .and. fa(nc+1:nc+nv)>0.0d0)) cycle
     mkap=nkap
     if(ll==0) mkap=1
     do ikap=1,mkap
       do ie=1,nband
         call band_ecut(ll,nproj(l1),eband(ie),elo-band_vint,rr, &
&                       vp(1,lloc+1,ikap),vlocpw, &
&                       vkbpw(1:nproj(l1),l1,ikap),evkb(1:nproj(l1),l1,ikap), &
&                       mmax,rc(l1),irc(l1),qcut(l1),qinf(l1), &
&                       bandq(ie,l1,ikap),bandc(:,ie,l1,ikap))
         write(6,'(a,i1,a,i1,f9.4,1p,e15.4,0p,4f8.2)') ' band l=',ll, &
&             ' ikap=',ikap,eband(ie),bandq(ie,l1,ikap),bandc(:,ie,l1,ikap)
       end do
     end do
   end do
 end if

 if(irel==2) then
   call run_diag_sr_so_r(lmax,npa,epa,wellprj,lloc,irc, &
&                        vlocpw,vsrpw,esr,vsopw,eso,nproj,rr,vfull,vp,zz, &
&                        mmax,mxprj)

! save full info for psml
   allocate(vpsml(mmax,mxl+2,mxkap))
   do l1=1,max(lmax+1,lloc+1)
     vpsml(:,l1,1) = vp(:,l1,1) - vo(:)
     vpsml(:,l1,2) = vp(:,l1,2) - vo(:)
   end do

! ghost testing

   write(6,'(/a)') 'Ghost test for J = L + 1/2'
 end if

 ngh(:,:)=0
 ghn(:,:,:)=0
 ghe(:,:,:,:)=0.0d0
 ghp(:,:,:)=0.0d0
 call run_ghosts(lmax,la,ea(1,1),nc,nv,lloc,rc,qmsbf, &
&                vlocpw,vkbpw(:,:,1),evkb(1,1,1),nproj,rr,vp,mmax,mxprj, &
&                ngh(:,1),ghn(:,:,1),ghe(:,:,:,1),ghp(:,:,1))

 if(irel==2) then
   write(6,'(/a)') 'Ghost test for J = L - 1/2'
   vkb(:,:,1,2)=vkb(:,:,1,1)
   evkb(:,1,2)=evkb(:,1,1)
   do jj=1,nproj(1)
     vkbpw(jj,1,2)=vkbpw(jj,1,1)
   end do
   call run_ghosts(lmax,la,ea(1,2),nc,nv,lloc,rc,qmsbf, &
&                  vlocpw,vkbpw(:,:,2),evkb(1,1,2),nproj,rr,vp,mmax,mxprj, &
&                  ngh(:,2),ghn(:,:,2),ghe(:,:,:,2),ghp(:,:,2))
   vkb(:,:,1,2)=0.0d0
   evkb(:,1,2)=0.0d0
   do jj=1,nproj(1)
     call pw_free(vkbpw(jj,1,2))
   end do
 end if

! unscreen semi-local potentials, taking the scalar-relativistic average
! when fully relativistic

 do l1=1,max(lmax+1,lloc+1)
   ll=l1-1
   if(irel/=2 .or. ll==0 .or. ll==4) then
     vpuns(:,l1)=vp(:,l1,1)-vo(:)
   else
     vpuns(:,l1)=(dble(ll+1)*vp(:,l1,1)+ dble(ll)*vp(:,l1,2))/dble(2*ll+1) &
&                -vo(:)
   end if
 end do

! -zion/r far out, where the all-electron charge outranges rho: switched
! over the outer half of the mesh (O-5), not cut where rho underflows to 0
 do ii=1,mmax
   sw=cinf_taper(rr(ii),0.5d0*rr(mmax),rr(mmax))
   do l1=1,max(lmax+1,lloc+1)
     vpuns(ii,l1)=sw*vpuns(ii,l1)-(1.0d0-sw)*zion/rr(ii)
   end do
 end do

! how closely the local potential is -zion/r beyond rcmax
 call coulomb_tail(vpuns(1,lloc+1),vion,rr,mmax,zion,rcmax)


! loop over reference plus test atom configurations
 rhot(:)=rho(:)
 do jj=1,ncnf+1

   write(6,'(/a,i2)') 'Test configuration',jj-1

! charge density is initialized to that of reference configuration

   rhot(:)=rho(:)

   if(istage/=2) then
     call run_config_ae(jj,nacnf,lacnf,facnf,nc,nvcnf,rr,zz, &
&                       mmax,iexc,ea,irel,nkap,istcnf,ierrcnf,itcnf, &
&                       kerrcnf,nvtcnf,natcnf,latcnf,fatcnf,eatcnf,etstcnf)
   end if

   call run_config_ps(jj,istcnf,ierrcnf,itcnf,kerrcnf,nvtcnf, &
&                     natcnf,latcnf,fatcnf,eatcnf,etstcnf,nc,rhot,rhomod,rr, &
&                     rcmax,mmax,mxprj,iexc,etot,epstot,nproj,vpuns, &
&                     lloc,vlocpw,vkbpw,evkb,irel,nkap,ntcf(jj), &
&                     itcf(:,:,jj),etcf(:,:,jj),dtcf(:,jj))

 end do !jj

 call run_plot(lmax,npa,epa,wellprj,lloc,irc, &
&               vkb,vlocpw,vkbpw,evkb,nproj,rr,vfull,vp,vpuns,zz,mmax,mxprj, &
&               drl,nrl, &
&               rho,rhoc,rhomod,irel,nkap,cvgplt)



!the pseudopotential for oncvpsp_ld.x
 if(len_trim(ps_filename)/=0) then
   call write_ps_file(trim(ps_filename),irel,nkap,atsym,zz,mmax,rr,lmax, &
&                     lloc,nproj,rc,epsh1,epsh2,depsh,rxpsh,epa,vp,vlocpw, &
&                     evkb,vkbpw)
 end if

!the energies of the log derivatives, from epsh2 down by depsh; none with
!--no-log-derivatives
 if(skip_ld) then
   npsh=0
   rphs(:)=0.0d0
   allocate(pshae(npsh,4,nkap),pshps(npsh,4,nkap))
 else
   npsh=int(((epsh2-epsh1)/depsh)-0.5d0)+1
   allocate(pshae(npsh,4,nkap),pshps(npsh,4,nkap))
   call run_phsft(lmax,lloc,nproj,epa,npsh,epsh2,depsh,rxpsh,vlocpw,vkbpw, &
&                 evkb, &
&                 rr,vfull,vp,zz,mmax,mxprj,rc,irel,nkap,rphs,pshae,pshps)
 end if

 call gnu_script(epa,evkb,lmax,lloc,mxprj,nproj,nkap)

 call vloc_cutoff(rr,mmax,rho,vpuns(1,lloc+1),zion,cvloc,vlchk)

 call diagout(lmax,nproj,mxprj,nkap,cvgplt,npsh,epsh2,depsh,rphs,pshae, &
&             pshps,ngh,ghn,ghe,ghp,(icmod==6 .or. icmod==7) &
&             .and. fcfact<=0.0d0 .and. budget_cut==1,rcfact,ecutw,wcurve, &
&             ektail,icmod,rr,mmax,rcmax,rhoc,corepw,lloc,evkb,bbeig, &
&             ncnf+1,ntcf,itcf,etcf,dtcf,cvloc,vlchk,optdg,wshare,nband,eband, &
&             bandq,bandc)
 deallocate(pshae,pshps)

 if(trim(psfile)=='psp8' .or. trim(psfile)=='both') then
   soscale=1.0d0
   call linout(lmax,lloc,rc,vkbpw,evkb,vsrpw,esr,vsopw,eso,nproj,rr,vpuns, &
&              vlocpw,rho,corepw,rhotae,rhoc,zz,zion,mmax,mxprj,iexc,icmod, &
&              nrl,drl, &
&              atsym,soscale, &
&              na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&              fa,rc0,ep,qcut,ep_as_given,well_onset,well_tail,wreg,debl, &
&              well_n,well_norm,ecut_wt,facnf,dvloc0, &
&              fcfact,fcfact0,rcfact, &
&              slack,budget_cut,fcgrid,rcgrid, &
&              epsh1,epsh2,depsh,rxpsh,rlmax,psfile,log_mesh,upf_mesh,nkap)
 end if

 if(trim(psfile)=='upf' .or. trim(psfile)=='both') then
   call upfout(lmax,lloc,rc,vkbpw,evkb,nproj,rr,vpuns,vlocpw,rho,corepw, &
&              zz,zion,mmax,mxprj,iexc,icmod,nrl,drl,atsym,epstot, &
&              na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&              fa,rc0,ep,qcut,ep_as_given,well_onset,well_tail,wreg,debl, &
&              well_n,well_norm,ecut_wt,facnf,dvloc0, &
&              fcfact0,rcfact, &
&              slack,budget_cut,fcgrid,rcgrid, &
&              epsh1,epsh2,depsh,rxpsh,rlmax,psfile,uupsa,ea,log_mesh, &
&              upf_mesh,nkap)
 end if


 if(trim(psfile)=='psml' .or. trim(psfile)=='both') then
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
   call psmlout(lmax,lloc,rc,vkb,evkb,mxprj,nproj,rr,vpsml,rho,rhomod, &
&               rcore, &
&               zz,zion,mmax,iexc,icmod,nrl,drl,atsym,epstot, &
&               na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&               fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact0,rcfact, &
&               slack,budget_cut,fcgrid,rcgrid, &
&               epsh1,epsh2,depsh,rxpsh,rlmax,psfile,log_mesh,upf_mesh, &
&               irel,nkap,ep_as_given,well_onset,well_tail,wreg,well_n, &
&               well_norm,ecut_wt)
 end if

 stop
 contains

 subroutine write_help()
!! the usage of this executable: what it does and its options
 character(len=16) :: full

 full=exe(1:index(exe,'_')-1)
 if(istage==0) full=exe
 write(6,'(3a)') 'Usage: ',trim(exe),' [options] < INPUT.dat'
 write(6,'(a)') ''
 if(istage==0) then
   write(6,'(3a)') 'Generates a ',trim(treatment),' optimized norm-conserving'
   write(6,'(a)') 'Vanderbilt pseudopotential from the input (D. R. Hamann, PRB'
   write(6,'(a)') '88, 085117 (2013)) and writes its tests and the'
   write(6,'(a)') 'pseudopotential file(s) of the input''s psfile to stdout.'
 else if(istage==1) then
   write(6,'(3a)') 'The all-electron stage of ',trim(full),'.x: solves the'
   write(6,'(3a)') trim(treatment),' all-electron atom and the all-electron'
   write(6,'(a)') 'halves of the test configurations, which depend on no'
   write(6,'(a)') 'pseudopotential input, and writes them to --ae-file for'
   write(6,'(3a)') trim(full),'_ps.x. Its stdout is the all-electron part of'
   write(6,'(3a)') 'the output of ',trim(full),'.x.'
 else
   write(6,'(3a)') 'The pseudopotential stage of ',trim(full),'.x: the same'
   write(6,'(a)') 'input and output, but the all-electron atom is read from'
   write(6,'(3a)') 'the --ae-file of ',trim(full),'_ae.x instead of solved,'
   write(6,'(a)') 'so a loop over pseudopotential inputs solves it once.'
 end if
 write(6,'(a)') ''
 write(6,'(a)') 'Options:'
 write(6,'(a)') '  -h, --help'
 write(6,'(a)') '      Show this help and exit'
 write(6,'(a)') '  -v, --version'
 write(6,'(a)') '      Show the version and the build, and exit'
 write(6,'(a)') '  -i, --input FILE'
 write(6,'(a)') '      Read the input, in the .dat text format, from FILE'
 write(6,'(a)') '      instead of stdin'
 write(6,'(a)') '  -t, --toml-input FILE'
 write(6,'(a)') '      Read the input from the TOML file FILE'
#if !(defined WITH_TOML)
 write(6,'(a)') '      (not available: this build has no TOML)'
#endif
 if(istage==1) then
   write(6,'(a)') '  -a, --ae-file FILE'
   write(6,'(a)') '      Write the all-electron data file FILE (required)'
 else if(istage==2) then
   write(6,'(a)') '  -a, --ae-file FILE'
   write(6,'(3a)') '      Read the all-electron data file FILE, written by ', &
&                  trim(full),'_ae.x'
   write(6,'(a)') '      for the same atom, reference configuration and mesh'
   write(6,'(a)') '      (required); the input''s test configurations must'
   write(6,'(a)') '      be among those the file was written for'
 end if
 if(istage/=1) then
   write(6,'(a)') '  -p, --ps-file FILE'
   write(6,'(a)') '      Also write the pseudopotential data file FILE, from'
   write(6,'(a)') '      which oncvpsp_ld.x computes the log derivatives at'
   write(6,'(a)') '      any radii and energies'
   write(6,'(a)') '  --no-log-derivatives'
   write(6,'(a)') '      Leave out the log derivatives (their plot data and'
   write(6,'(a)') '      the DIAGNOSTICS ld lines); oncvpsp_ld.x gives them'
   write(6,'(a)') '      from --ps-file'
 end if
 write(6,'(a)') ''
 write(6,'(a)') 'Errors are reported on stdout, with a nonzero exit status.'
 return
 end subroutine write_help

 end subroutine oncvpsp_main

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

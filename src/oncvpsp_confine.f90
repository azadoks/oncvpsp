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
! The confined atom, oncvpsp_confine.x (doc/smoothness_plan.md,
! diagnostics item 7): the all-electron atom and the pseudopotential in a
! Wigner-Seitz sphere, at the frozen potentials of the all-electron data
! file of an _ae.x (ae_file_m) and the pseudopotential data file of a run
! with --ps-file (ps_file_m), for any relativistic treatment.
!
! For each l and j and each radius rs of the sphere, the energies of the
! band edges of the Wigner-Seitz rule: where r phi'/phi = 0 at rs (the
! bottom of a band), -l-1 (its centre) and phi = 0 (its top), the roots of
! the phases atan(r phi'/phi) - pi (nodes inside rs) of phsft_node_phase,
! which fall monotonically with the energy, so a scan brackets every root
! whatever its step and Brent's method takes it to round-off. The
! pseudopotential's phase counts instead the poles of its r phi'/phi at rs
! below the energy, the states of phi(rs) = 0 in the sphere: its wave
! function can gain a node inside the projectors (Pb p, I s with a large
! debl), but r phi'/phi falls between its poles for any Hermitian
! potential. They are the eigenvalues of ws_eigenvalues in the basis of
! that condition, each taken to the root of u(rs) beside it. The
! pseudopotential's phases are shifted by the multiple of pi by which
! they differ from the all-electron ones at most energies, so an
! all-electron root and the pseudopotential root of the same phase are
! the same band's edge; n is the all-electron state's nodes inside rs plus
! l+1. For each pseudopotential edge, the cutoff that the state needs in
! the sphere, from the variational energies of ws_basis_m: a cutoff for
! every channel at the energies of the valence bands, with no tail beyond
! rs to invent. The output, on stdout:
!   ws_r l ikap rs(bohr)
!   ws l ikap edge n e_ae e_ps            (Ha; edge 1 bottom, 2 centre,
!                                          3 top)
!   ws_ecut l ikap edge n tol ecut        (Ha: the cutoff at which the
!                                          energy error falls to tol)
!   ws_floor l ikap edge n ecut err       (the error at the largest cutoff)
!   ws_conv l ikap edge n ecut err        (--table: every basis size)
!   ws_missing l ikap edge n e_ae         (no pseudopotential root of that
!                                          phase in the scan)
!   ws_core l ikap edge n e_ae            (the same, the band of a core
!                                          state, which the
!                                          pseudopotential leaves out)
!   ws_paired_off l ikap edge n e_ps      (a pseudopotential edge of
!                                          another phase, nearer the
!                                          all-electron one than its
!                                          partner)
!   ws_extra l ikap edge e_ps             (a pseudopotential edge with no
!                                          all-electron partner: a state
!                                          the all-electron atom has not,
!                                          such as a ghost)
! for l=0..3 (or those of --channels) and each kappa of l, ikap 1 for
! j=l+1/2 and 2 for j=l-1/2.
!
! Then the atom compressed: the self-consistent all-electron atom (aeatom)
! and pseudo-atom (psatom) of the reference configuration in a smooth wall
! 1/2 k x**2 exp(-dwall/x), x = r - r0 beyond r0 (by default the largest
! core radius), for a few strengths k, each from the free atom:
!   wall_r0 r0(bohr)
!   wall k etot_ae etot_ps de_ae de_ps    (Ha, k in Ha/bohr**2; de the
!                                          change from the free atom, k=0)
!   wall_eig k n l ikap e_ae e_ps         (the valence eigenvalues, Ha)
!   wall_failed k atom                    (atom ae or ps: no
!                                          self-consistent solution)
! and last a summary of both:
!   confine_summary max_de ecut4 ecut5 n_missing n_extra n_paired_off
!                   n_wall_failed max_dde ecut4_top ecut5_top
! max_de the largest |e_ps - e_ae| of the edges below 0 (Ha), ecut4 and
! ecut5 the largest cutoffs of those edges at 1e-4 and 1e-5 Ha (Infinity
! if one does not reach it), the counts of the ws_missing, ws_extra,
! ws_paired_off and wall_failed lines, max_dde the largest
! |de_ps - de_ae| of the walls (Ha), and ecut4_top and ecut5_top those of
! the top edges alone, whose errors have no tail from the potential's
! slope at rs (ws_basis_m); NaN and -1 for the parts left out.

 program oncvpsp_confine
 use constants_m, only: dp, pi
 use ae_file_m, only: read_ae_potential
 use ps_file_m, only: ps_file_t, read_ps_file, ps_file_free
 use phsft_m, only: phsft_t, phsft_ae_init, phsft_ps_init, phsft_free, &
&                   phsft_node_phase, phsft_node_phases, phsft_mch
 use ws_basis_m, only: ws_edge_d, ws_convergence, ws_eigenvalues, &
&                       ws_grid_init, ws_grid_t
 use ugrid_m, only: ug_brent
 use relativity_m, only: set_relativity
 use constants_m, only: mxstate, mxprj, mxl, mxlin
 use build_info_m, only: write_build_info
 use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
&                                         ieee_positive_inf
 implicit none

!> relative tolerance of the comparison of the two files' meshes and
!> constants, as in ae_file_m
 real(dp), parameter :: rtol=1.0d-12
!> the most radii of a run
 integer, parameter :: mxq=100
!> the most edges of a phase scan
 integer, parameter :: mxedge=400
!> the energy errors of the cutoffs (Ha)
 real(dp), parameter :: tols(5)=[1.0d-2,1.0d-3,1.0d-4,1.0d-5,1.0d-6]
!> the default radii, in units of the largest core radius
 real(dp), parameter :: rsdef(3)=[1.25d0,1.5d0,2.0d0]
!> the default step of the scan (Ha): it only brackets the roots of
!> monotonic phases, which no step can miss
 real(dp), parameter :: descan=0.1d0
!> the scan of the pseudopotential's phases runs this far (Ha) beyond the
!> energies of the all-electron edges, for their partners
 real(dp), parameter :: emargin=1.0d0
!> the points and spacing (Ha) of the phases' alignment about the
!> reference energy
 integer, parameter :: nanch=21
 real(dp), parameter :: danch=0.05d0
!> the default strengths of the wall (Ha/bohr**2)
 real(dp), parameter :: kdef(3)=[0.01d0,0.03d0,0.1d0]
!> the onset of the wall, x**2 exp(-dwall/x): smooth at r0, so that outwf_m
!> takes it as part of the local potential (bohr)
 real(dp), parameter :: dwall=0.5d0
!> the most strengths of a run
 integer, parameter :: mxk=20

 character(len=1024) :: ae_filename,ps_filename,arg
 character(len=2) :: atsym
 character(len=19) :: treatment
 integer :: iarg,irel,mmax,nq,ii,iq,l1,ll,ikap,kap,nvkb,ios,lc,nl
 integer :: npsh,nae,nps,iedge,ie,je,nb,ierr,itol,nn
 integer :: lsel(4)
 integer :: kae(mxedge),kps(mxedge),edae(mxedge),edps(mxedge)
 real(dp) :: zz,c_light_f,emin,emax,de,ecutmax,rcmax,rs,amesh,ee(3)
 real(dp) :: rq(mxq),eae(mxedge),eps(mxedge),thal,ex
 real(dp), allocatable :: rr(:),vfull(:),shae(:,:),shps(:,:),ecut(:),err(:)
 real(dp), allocatable :: egrid(:)
 integer, allocatable :: ishift(:)
 integer :: ianch(nanch)
 logical :: given_radii,given_energies,table,sphere,given_r0
 integer :: nk,nc,nv,iexc,ik,ierrae,ierrps,kk,itae,itps,nkap
 integer :: na(mxstate),la(mxstate),natp(mxstate),latp(mxstate),nav(4)
 integer :: nprj(mxlin+1)
 real(dp) :: fa(mxstate),ea(mxstate,2),r0,kwall(0:mxk)
 real(dp) :: etae(0:mxk),etps(0:mxk),eae_k(mxstate,2),eps_k(mxstate,2)
 real(dp) :: eaw(mxstate,2),rpk(mxstate,2),eterm(4),res,eatp(mxstate,2)
 real(dp) :: fatp(mxstate,2)
 real(dp), allocatable :: vext(:),rhoc(:),rhoa(:),vi(:),vpuns(:,:),rhops(:)
 real(dp), allocatable :: eigae(:,:,:),eigps(:,:,:)
!> the angular momenta of the output
 logical :: chan(0:3)
 type(ps_file_t) :: pf
 type(phsft_t), target :: sae,sps
!> the quadrature of the channel in each sphere, shared by its edges
 type(ws_grid_t), allocatable :: wsg(:)
!> the poles of each sphere (edir, ndir of the radius at hand)
 type :: poles_t
   real(dp), allocatable :: e(:)
 end type poles_t
 type(poles_t), allocatable :: polq(:)
!> the phases of every radius at the anchors
 real(dp), allocatable :: tanae(:,:),tanps(:,:)
!> the solver, the point of the scan and the phase of edges' root
 type(phsft_t), pointer :: psol=>null()
 real(dp) :: ph_e0,ph_tgt,ph_shift
!> the radius of psol's phase
 integer :: ph_iq
!> the poles of the pseudopotential's r phi'/phi at rs (Ha), and whether
!> the phase of psol counts them rather than its nodes
 integer :: ndir
 real(dp), allocatable :: edir(:)
 logical :: ph_poles
!> the partner of each all-electron edge (0 for none), and whether each
!> pseudopotential edge is one
 integer :: jpair(mxedge)
 logical :: used(mxedge)
!> the summary: the counts of the lines, the largest |e_ps - e_ae| and
!> cutoffs (1e-4, 1e-5 Ha; then those of the top edges alone) of the edges
!> below 0, whether one of them does not reach it, and the largest
!> |de_ps - de_ae| of the walls (-1 for none)
 integer :: nmiss,nextra,noff,nwfail
 real(dp) :: smax_de,scut(4),sdde,ec
 logical :: sinf(4)

 ae_filename=''
 ps_filename=''
 given_radii=.false.
 given_energies=.false.
 table=.false.
 sphere=.true.
 given_r0=.false.
 nk=size(kdef)
 kwall(1:nk)=kdef(:)
 ecutmax=200.0d0
 chan(:)=.true.
 nmiss=-1
 nextra=-1
 noff=-1
 nwfail=-1
 smax_de=-1.0d0
 scut(:)=-1.0d0
 sinf(:)=.false.
 sdde=-1.0d0
 nq=0
 iarg=1
 do while(iarg<=command_argument_count())
   call get_command_argument(iarg,arg)
   select case(arg)
   case('-h','--help')
     call write_help()
     stop
   case('-v','--version')
     call write_build_info(6,'oncvpsp_confine.x','any relativistic treatment')
     stop
   case('-a','--ae-file')
     call next_arg(ae_filename)
   case('-p','--ps-file')
     call next_arg(ps_filename)
   case('-r','--radii')
     call next_arg(arg)
     nq=count_list(arg)
     if(nq>mxq) then
       write(6,'(a,i0,a)') 'oncvpsp_confine: ERROR more than ',mxq,' radii'
       error stop 1
     end if
     read(arg,*,iostat=ios) rq(1:nq)
     if(ios/=0 .or. nq<1 .or. any(rq(1:nq)<=0.0d0)) then
       write(6,'(2a)') 'oncvpsp_confine: ERROR --radii must be positive ', &
&            'numbers separated by commas'
       error stop 1
     end if
     given_radii=.true.
   case('-l','--channels')
     call next_arg(arg)
     nl=min(count_list(arg),size(lsel)+1)
     lsel(:)=-1
     ios=1
     if(nl<=size(lsel)) read(arg,*,iostat=ios) lsel(1:nl)
     if(ios/=0) lsel(1)=-1
     if(any(lsel(1:min(nl,size(lsel)))<0) &
&       .or. any(lsel(1:min(nl,size(lsel)))>3)) then
       write(6,'(2a)') 'oncvpsp_confine: ERROR --channels must be ', &
&            'values of l from 0 to 3 separated by commas'
       error stop 1
     end if
     chan(:)=.false.
     chan(lsel(1:nl))=.true.
   case('-e','--energies')
     call next_arg(arg)
     read(arg,*,iostat=ios) ee(1:3)
     if(ios/=0 .or. count_list(arg)/=3 .or. ee(3)<=0.0d0 &
&       .or. ee(2)<=ee(1)) then
       write(6,'(2a)') 'oncvpsp_confine: ERROR --energies must be ', &
&            'emin,emax,de with emin < emax and de > 0'
       error stop 1
     end if
     given_energies=.true.
   case('-c','--ecut-max')
     call next_arg(arg)
     read(arg,*,iostat=ios) ecutmax
     if(ios/=0 .or. ecutmax<=0.0d0) then
       write(6,'(a)') 'oncvpsp_confine: ERROR --ecut-max must be positive'
       error stop 1
     end if
   case('-t','--table')
     table=.true.
   case('-s','--no-sphere')
     sphere=.false.
   case('-w','--walls')
     call next_arg(arg)
     if(trim(arg)=='none') then
       nk=0
     else
       nk=count_list(arg)
       if(nk>mxk) then
         write(6,'(a,i0,a)') 'oncvpsp_confine: ERROR more than ',mxk, &
&             ' walls'
         error stop 1
       end if
       read(arg,*,iostat=ios) kwall(1:nk)
       if(ios/=0 .or. any(kwall(1:nk)<=0.0d0)) then
         write(6,'(2a)') 'oncvpsp_confine: ERROR --walls must be positive ', &
&              'numbers separated by commas, or none'
         error stop 1
       end if
     end if
   case('-R','--wall-radius')
     call next_arg(arg)
     read(arg,*,iostat=ios) r0
     if(ios/=0 .or. r0<=0.0d0) then
       write(6,'(a)') 'oncvpsp_confine: ERROR --wall-radius must be positive'
       error stop 1
     end if
     given_r0=.true.
   case default
     write(6,'(2a)') 'oncvpsp_confine: ERROR unknown argument ',trim(arg)
     error stop 1
   end select
   iarg=iarg+1
 end do
 if(len_trim(ae_filename)==0 .or. len_trim(ps_filename)==0) then
   write(6,'(a)') 'oncvpsp_confine: ERROR requires --ae-file and --ps-file'
   error stop 1
 end if

 call read_ae_potential(trim(ae_filename),irel,atsym,zz,mmax,rr,c_light_f, &
&                       vfull,nc=nc,nv=nv,iexc=iexc,na=na,la=la,fa=fa,ea=ea)
 call read_ps_file(trim(ps_filename),pf)
 if(pf%irel/=irel .or. pf%atsym/=atsym .or. differ(pf%zz,zz) &
&   .or. pf%mmax/=mmax .or. differ(pf%c_light,c_light_f) &
&   .or. pf%iexc/=iexc) then
   write(6,'(2a)') 'oncvpsp_confine: ERROR the files differ in ', &
&       'relativistic treatment, atom, mesh size, c_light or functional'
   error stop 1
 end if
 if(any(abs(pf%rr(:)-rr(:))>rtol*rr(:))) then
   write(6,'(a)') 'oncvpsp_confine: ERROR the files'' radial meshes differ'
   error stop 1
 end if
 call set_relativity(c_light_f)

 if(irel==0) then
   treatment='non-relativistic'
 else if(irel==1) then
   treatment='scalar-relativistic'
 else
   treatment='fully-relativistic'
 end if
 call write_build_info(6,'oncvpsp_confine.x',trim(treatment))

!the radii: beyond every core radius, where the two atoms' potentials are
!the same and the phases comparable
 lc=pf%lloc+1
 rcmax=max(maxval(pf%rc(1:pf%lmax+1)),pf%rc(lc))
 if(.not. given_radii) then
   nq=size(rsdef)
   rq(1:nq)=rsdef(:)*rcmax
 end if
 do iq=1,nq
   if(rq(iq)<=rcmax) then
     write(6,'(a,f10.5,a,f10.5,a)') 'oncvpsp_confine: ERROR the radius', &
&         rq(iq),' is not beyond the largest core radius',rcmax,' bohr'
     error stop 1
   end if
 end do

!the window of the edges, that of the log derivatives
 if(given_energies) then
   emin=ee(1)
   emax=ee(2)
   de=ee(3)
 else
   emin=pf%epsh1
   emax=pf%epsh2
   de=descan
 end if
 npsh=int(((emax-emin+2.0d0*emargin)/de)-0.5d0)+1
 amesh=rr(2)/rr(1)
 allocate(shae(npsh,nq),shps(npsh,nq),egrid(npsh),ishift(npsh))
 do ii=1,npsh
   egrid(ii)=emax+emargin-de*dble(ii-1)
 end do

 if(sphere) then
   write(6,'(a)') '# ws_r l ikap rs(bohr)'
   write(6,'(a)') '# ws l ikap edge n e_ae(Ha) e_ps(Ha)'
   write(6,'(a)') '# ws_ecut l ikap edge n tol(Ha) ecut(Ha)'
   write(6,'(a)') '# ws_floor l ikap edge n ecut(Ha) err(Ha)'
   if(table) write(6,'(a)') '# ws_conv l ikap edge n ecut(Ha) err(Ha)'
   write(6,'(a)') '# ws_missing l ikap edge n e_ae(Ha)'
   write(6,'(a)') '# ws_core l ikap edge n e_ae(Ha)'
   write(6,'(a)') '# ws_paired_off l ikap edge n e_ps(Ha)'
   write(6,'(a)') '# ws_extra l ikap edge e_ps(Ha)'
   nmiss=0
   nextra=0
   noff=0
   do l1=1,4
     ll=l1-1
     if(.not. chan(ll)) cycle
     do ikap=1,merge(1,pf%nkap,l1==1)
       kap=0
       if(irel==2 .and. ikap==1) kap=-(ll+1)
       if(irel==2 .and. ikap==2) kap=ll
       nvkb=pf%nvkb(l1,ikap)
! the solvers of the channel at every radius, and each sphere's
! quadrature and poles (refined to the roots of u(rs))
       call phsft_ae_init(irel,ll,kap,rr,vfull,zz,mmax,nq,rq(1:nq),sae)
       call phsft_ps_init(ll,nvkb,rr,pf%vloc(:,ikap),pf%vlocpw, &
&                         pf%vkbpw(1:nvkb,l1,ikap),pf%evkb(1:nvkb,l1,ikap), &
&                         mmax,nq,rq(1:nq),sps)
       allocate(wsg(nq),polq(nq),tanae(nanch,nq),tanps(nanch,nq))
       do iq=1,nq
         call ws_grid_init(rq(iq),sqrt(2.0d0*ecutmax),rr,mmax, &
&             pf%vloc(:,ikap),pf%vlocpw,nvkb,pf%vkbpw(1:nvkb,l1,ikap), &
&             wsg(iq))
         call ws_eigenvalues(ll,3,wsg(iq),pf%evkb(1:nvkb,l1,ikap),ndir, &
&             edir,ierr)
         if(ierr/=0) ndir=0
         ph_iq=iq
         call refine_poles(sps)
         polq(iq)%e=edir(1:ndir)
         deallocate(edir)
       end do

! the phases at rs on the grid, the pseudopotential's shifted by the
! multiple of pi by which they differ from the all-electron ones near the
! first projector's reference energy, where the pseudopotential matches
! the all-electron log derivative by construction: the core states of l
! below it, and a ghost below it too.  They agree beyond the core radii
! but where they turn steeply, at bound states, where a small difference
! in their energies is a difference by pi in the phases (a narrow band,
! such as a semicore state's in a large sphere, where the reference
! energy of phsft_pair can be): the most frequent multiple on +-0.5 Ha
! of it.  Without projectors (l = lloc or beyond lmax), the most frequent
! on the grid.  Over the grid of a wide window (epsh1 -12 Ha for Cs
! sp), the multiple below a core band and that above several Ha, where
! the pseudopotential no longer follows the atom, outnumbered the right
! one and paired each band with the next one down
! (all the radii from one solution at each energy)
       do ii=1,npsh
         call scan_phases(egrid(ii),shae(ii,:),shps(ii,:))
       end do
       if(nvkb>0) then
         do ii=1,nanch
           ex=pf%epa(l1,ikap)+danch*dble(ii-(nanch+1)/2)
           call scan_phases(ex,tanae(ii,:),tanps(ii,:))
         end do
       end if

       do iq=1,nq
         rs=rq(iq)
         write(6,'(a,2i4,1p,e24.16)') 'ws_r',ll,ikap,rs
         ph_iq=iq
         ndir=size(polq(iq)%e)
         edir=polq(iq)%e
         do ii=1,npsh
           ishift(ii)=nint((shae(ii,iq)-shps(ii,iq))/pi)
         end do
         if(nvkb>0) then
           do ii=1,nanch
             ianch(ii)=nint((tanae(ii,iq)-tanps(ii,iq))/pi)
           end do
           thal=pi*dble(majority(ianch))
         else
           thal=pi*dble(majority(ishift))
         end if
         shps(:,iq)=shps(:,iq)+thal

         call edges(sae,.false.,shae(:,iq),0.0d0,nae,eae,kae,edae)
         call edges(sps,.true.,shps(:,iq),thal,nps,eps,kps,edps)
         deallocate(edir)

! every all-electron edge, lowest first, and its partner: the nearest
! pseudopotential edge of the same phase, over the whole scan, so that an
! edge in the window whose partner is outside it is not without one
         call sort_edges(nae,eae,kae,edae)
         call sort_edges(nps,eps,kps,edps)
         used(1:nps)=.false.
         do ie=1,nae
           jpair(ie)=0
           do ii=1,nps
             if(kps(ii)/=kae(ie) .or. edps(ii)/=edae(ie)) cycle
             if(jpair(ie)==0) then
               jpair(ie)=ii
             else if(abs(eps(ii)-eae(ie))<abs(eps(jpair(ie))-eae(ie))) then
               jpair(ie)=ii
             end if
           end do
           if(jpair(ie)>0) used(jpair(ie))=.true.
         end do

! those in the window
         do ie=1,nae
           if(eae(ie)<emin .or. eae(ie)>emax) cycle
           iedge=edae(ie)
           nn=n_of(eae(ie))
           je=jpair(ie)
           if(je==0) then
! the band of a core state, which the pseudopotential leaves out
             if(any(na(1:nc)==nn .and. la(1:nc)==ll)) then
               write(6,'(a,4i4,1p,e24.16)') 'ws_core',ll,ikap,iedge,nn, &
&                   eae(ie)
             else
               write(6,'(a,4i4,1p,e24.16)') 'ws_missing',ll,ikap,iedge, &
&                   nn,eae(ie)
               nmiss=nmiss+1
             end if
             cycle
           end if
           write(6,'(a,4i4,1p,2e24.16)') 'ws',ll,ikap,iedge,nn,eae(ie),eps(je)
! an edge of another phase nearer to it than its partner: the phases pair
! a band of one atom with another band of the other
           do ii=1,nps
             if(ii==je .or. edps(ii)/=iedge) cycle
             if(abs(eps(ii)-eae(ie))<abs(eps(je)-eae(ie))) then
               write(6,'(a,4i4,1p,e24.16)') 'ws_paired_off',ll,ikap,iedge, &
&                   nn,eps(ii)
               noff=noff+1
               exit
             end if
           end do
           if(eae(ie)<0.0d0) smax_de=max(smax_de,abs(eps(je)-eae(ie)))

! the cutoff of the pseudopotential's state
           call ws_convergence(ll,iedge,wsg(iq),pf%evkb(1:nvkb,l1,ikap), &
&               eps(je),nb,ecut,err,ierr)
           if(ierr/=0) then
             if(eae(ie)<0.0d0) sinf(1:2)=.true.
             if(eae(ie)<0.0d0 .and. iedge==3) sinf(3:4)=.true.
             cycle
           end if
           do itol=1,size(tols)
             call cutoff_at(tols(itol),ec)
! the summary's cutoffs, at 1e-4 and 1e-5 Ha, of all edges and of the top
             if(eae(ie)<0.0d0 .and. (itol==3 .or. itol==4)) then
               do kk=itol-2,merge(itol,itol-2,iedge==3),2
                 if(ec<0.0d0) sinf(kk)=.true.
                 scut(kk)=max(scut(kk),ec)
               end do
             end if
           end do
           write(6,'(a,4i4,1p,2e24.16)') 'ws_floor',ll,ikap,iedge,nn, &
&               ecut(nb),err(nb)
           if(table) then
             do ii=1,nb
               if(err(ii)==0.0d0) cycle
               write(6,'(a,4i4,1p,2e24.16)') 'ws_conv',ll,ikap,iedge,nn, &
&                   ecut(ii),err(ii)
             end do
           end if
         end do

! the pseudopotential's edges in the window that no all-electron edge
! has as its partner
         do ii=1,nps
           if(used(ii) .or. eps(ii)<emin .or. eps(ii)>emax) cycle
           write(6,'(a,3i4,1p,e24.16)') 'ws_extra',ll,ikap,edps(ii),eps(ii)
           nextra=nextra+1
         end do
       end do
       call phsft_free(sae)
       call phsft_free(sps)
       deallocate(wsg,polq,tanae,tanps)
     end do
   end do
 end if

!the atom in walls: k=0 the free atom
 if(nk>0) then
   if(.not. given_r0) r0=rcmax
   kwall(0)=0.0d0
   nkap=pf%nkap
   allocate(vext(mmax),rhoc(mmax),rhoa(mmax),vi(mmax),rhops(mmax))
   allocate(vpuns(mmax,mxl+2),eigae(mxstate,2,0:nk),eigps(mxstate,2,0:nk))
   eigae(:,:,:)=0.0d0
   eigps(:,:,:)=0.0d0
   vpuns(:,:)=0.0d0
   vpuns(:,pf%lloc+1)=pf%vpuns(:)
   nprj(:)=0
   do l1=1,4
     nprj(l1)=pf%nvkb(l1,1)
   end do
! the valence states as run_config_ps numbers them, n from l+1
   do l1=1,4
     nav(l1)=l1-1
   end do
   do kk=1,nv
     latp(kk)=la(nc+kk)
     l1=latp(kk)+1
     nav(l1)=nav(l1)+1
     natp(kk)=nav(l1)
     fatp(kk,:)=fa(nc+kk)
   end do
   write(6,'(a)') '# wall_r0 r0(bohr)'
   write(6,'(2a)') '# wall k(Ha/bohr**2) etot_ae(Ha) etot_ps(Ha) de_ae(Ha)', &
&        ' de_ps(Ha)'
   write(6,'(a)') '# wall_eig k n l ikap e_ae(Ha) e_ps(Ha)'
   write(6,'(a)') '# wall_failed k atom'
   write(6,'(a,1p,e24.16)') 'wall_r0',r0
   nwfail=0
   ierrae=0
   ierrps=0
   do ik=0,nk
     do ii=1,mmax
       vext(ii)=wall(kwall(ik),rr(ii)-r0)
     end do
! each atom from the wall before (the free atom from the data file's
! potential and the run's pseudo-density): the walls ramp the atom in
     if(ik==0) then
       eaw(:,:)=ea(:,:)
       vi(:)=vfull(:)
     else if(ierrae/=0) then
       eaw(:,:)=eigae(:,:,0)
       vi(:)=vfull(:)
     end if
     rpk(:,:)=0.0d0
     ierrae=0
     call aeatom(irel,nkap,na,la,eaw(:,1:nkap),fa,rpk(:,1:nkap),nc,nc+nv, &
&                itae,rhoc,rhoa,rr,vi,zz,mmax,iexc,etae(ik),eterm,res,ierrae, &
&                vext,.true.)
     eigae(:,:,ik)=eaw(:,:)
     if(ik==0 .or. ierrps/=0) then
       rhops(:)=pf%rhov(:)
       eatp(1:nv,:)=ea(nc+1:nc+nv,:)
     end if
     ierrps=0
     call psatom(natp,latp,eatp(:,1:nkap),fatp,nv,itps,pf%rhomod,rhops, &
&                rr,pf%rcmax,mmax,mxprj,iexc,etps(ik),nprj,vpuns,pf%lloc, &
&                pf%vlocpw,pf%vkbpw(:,:,1:nkap),pf%evkb(:,:,1:nkap),ierrps, &
&                nkap,vext)
     eigps(1:nv,:,ik)=eatp(1:nv,:)
     if(ierrae/=0) write(6,'(a,1p,e24.16,a)') 'wall_failed',kwall(ik),' ae'
     if(ierrps/=0) write(6,'(a,1p,e24.16,a)') 'wall_failed',kwall(ik),' ps'
     if(ik==0 .and. (ierrae/=0 .or. ierrps/=0)) then
       write(6,'(a)') 'oncvpsp_confine: ERROR the free atom failed'
       error stop 1
     end if
     if(ierrae/=0 .or. ierrps/=0) then
       nwfail=nwfail+1
       cycle
     end if
     write(6,'(a,1p,5e24.16)') 'wall',kwall(ik),etae(ik),etps(ik), &
&         etae(ik)-etae(0),etps(ik)-etps(0)
     if(ik>0) sdde=max(sdde,abs(etps(ik)-etps(0)-etae(ik)+etae(0)))
     do kk=1,nv
       ll=la(nc+kk)
       do ikap=1,merge(1,nkap,ll==0)
         write(6,'(a,1p,e24.16,3i4,2e24.16)') 'wall_eig',kwall(ik), &
&             na(nc+kk),ll,ikap,eigae(nc+kk,ikap,ik),eigps(kk,ikap,ik)
       end do
     end do
   end do
   deallocate(vext,rhoc,rhoa,vi,rhops,vpuns,eigae,eigps)
 end if

 write(6,'(3a)') '# confine_summary max_de(Ha) ecut4(Ha) ecut5(Ha)', &
&     ' n_missing n_extra n_paired_off n_wall_failed max_dde(Ha)', &
&     ' ecut4_top(Ha) ecut5_top(Ha)'
 write(6,'(a,1p,3e24.16,4i6,3e24.16)') 'confine_summary', &
&     summ(smax_de,.false.),summ(scut(1),sinf(1)),summ(scut(2),sinf(2)), &
&     nmiss,nextra,noff,nwfail,summ(sdde,.false.),summ(scut(3),sinf(3)), &
&     summ(scut(4),sinf(4))

 deallocate(shae,shps,egrid,ishift)
 call ps_file_free(pf)
 deallocate(rr,vfull)
 stop

 contains

 subroutine edges(sol,poles,sh,shift,ned,eed,ked,ied)
!! the edges of a scan: the energies at which the phase (plus shift)
!! crosses atan(D) + k pi for the D of each edge, each to round-off by
!! Brent's method in the interval of the grid that brackets it; poles for
!! the phase that counts the poles edir rather than the nodes, whose top
!! edges (phi(rs) = 0) are the poles themselves, where the count steps: the
!! k-th at the phase -pi/2 - (k-1) pi just below it
 type(phsft_t), intent(inout), target :: sol
 logical, intent(in) :: poles
 real(dp), intent(in) :: sh(npsh),shift
 integer, intent(out) :: ned,ked(mxedge),ied(mxedge)
 real(dp), intent(out) :: eed(mxedge)
 real(dp) :: thlo,thhi,er
 integer :: jj,kk,jed,ierr2
 psol=>sol
 ph_poles=poles
 ph_shift=shift
 ned=0
 do jj=1,npsh-1
   ph_e0=egrid(jj)
   thlo=min(sh(jj),sh(jj+1))
   thhi=max(sh(jj),sh(jj+1))
   do jed=1,3
     if(poles .and. jed==3) cycle
     do kk=ceiling((thlo-target(jed))/pi),floor((thhi-target(jed))/pi)
       ph_tgt=target(jed)+pi*dble(kk)
! a crossing on the point below is the next interval's
       if(ph_tgt==sh(jj+1) .and. jj<npsh-1) cycle
       call ug_brent(phase,egrid(jj+1),ph_e0,er,ierr2)
       if(ierr2/=0) cycle
! a jump (a pair of nodes inside the projectors), not a root
       if(abs(phase(er))>0.25d0*pi) cycle
       if(ned==mxedge) then
         write(6,'(a)') 'oncvpsp_confine: ERROR too many edges in the scan'
         error stop 1
       end if
       ned=ned+1
       eed(ned)=er
       ked(ned)=kk
       ied(ned)=jed
     end do
   end do
 end do
 if(poles) then
   do kk=1,ndir
     if(edir(kk)<egrid(npsh) .or. edir(kk)>egrid(1)) cycle
     if(ned==mxedge) then
       write(6,'(a)') 'oncvpsp_confine: ERROR too many edges in the scan'
       error stop 1
     end if
     ned=ned+1
     eed(ned)=edir(kk)
     ked(ned)=nint(shift/pi)-(kk-1)
     ied(ned)=3
   end do
 end if
 nullify(psol)
 end subroutine edges

 subroutine scan_phases(ex,thae,thps)
!! the phases of every radius at ex from one solution of each atom: the
!! all-electron one's, and the pseudopotential's less pi per pole of each
!! sphere below ex (phase_of)
 real(dp), intent(in) :: ex
 real(dp), intent(out) :: thae(nq),thps(nq)
 integer :: jq
 call phsft_node_phases(sae,ex,thae)
 call phsft_node_phases(sps,ex,thps,nodes=.false.)
 do jq=1,nq
   thps(jq)=thps(jq)-pi*dble(count(polq(jq)%e<ex))
 end do
 end subroutine scan_phases

 real(dp) function phase(ex)
!! the phase of psol at ex plus ph_shift, less ph_tgt
 real(dp), intent(in) :: ex
 phase=phase_of(ex)+ph_shift-ph_tgt
 end function phase

 real(dp) function phase_of(ex)
!! the phase of psol at ex: less pi per node inside rs, or with ph_poles
!! per pole of r phi'/phi at rs below ex
 real(dp), intent(in) :: ex
 if(ph_poles) then
   phase_of=phsft_node_phase(psol,ph_iq,ex,nodes=.false.) &
&           -pi*dble(count(edir(1:ndir)<ex))
 else
   phase_of=phsft_node_phase(psol,ph_iq,ex)
 end if
 end function phase_of

 subroutine refine_poles(sol)
!! the poles edir in the scan, each from its variational value to the root
!! of u(rs) beside it, in a bracket that grows from 1e-6 Ha on both sides
!! until u changes sign, up to halfway to the poles beside it (or 1 Ha). A
!! variational bound lies above the root of the same Hamiltonian, but
!! the basis's quadrature of the pieces and the radial equation's differ,
!! by 3e-8 Ha at 1.5 Ha in Si f, so the root can be on either side
 type(phsft_t), intent(inout), target :: sol
 real(dp) :: e0,u0,lo,hi,ww,er,elow,ehigh
 integer :: kk,ierr2
 logical :: found
 psol=>sol
 do kk=1,ndir
   e0=edir(kk)
   if(e0<egrid(npsh) .or. e0>egrid(1)) cycle
   elow=e0-1.0d0
   if(kk>1) elow=max(elow,0.5d0*(edir(kk-1)+e0))
   ehigh=e0+1.0d0
   if(kk<ndir) ehigh=min(ehigh,0.5d0*(e0+edir(kk+1)))
   u0=u_rs(e0)
   if(u0==0.0d0) cycle
   ww=1.0d-6*max(1.0d0,abs(e0))
   found=.false.
   do
     lo=max(e0-ww,elow)
     hi=min(e0+ww,ehigh)
     if(u_rs(lo)*u0<0.0d0) then
       hi=e0
       found=.true.
     else if(u_rs(hi)*u0<0.0d0) then
       lo=e0
       found=.true.
     end if
     if(found .or. (lo<=elow .and. hi>=ehigh)) exit
     ww=4.0d0*ww
   end do
   if(.not. found) cycle
   call ug_brent(u_rs,lo,hi,er,ierr2)
   if(ierr2==0) edir(kk)=er
 end do
 nullify(psol)
 end subroutine refine_poles

 real(dp) function u_rs(ex)
!! u of psol at rs and the energy ex
 real(dp), intent(in) :: ex
 real(dp) :: th
 th=phsft_node_phase(psol,ph_iq,ex,nodes=.false.,uq_out=u_rs)
 end function u_rs

 integer function majority(iv)
!! the most frequent value, the smallest of a tie
 integer, intent(in) :: iv(:)
 integer :: jj,nc,nbest
 nbest=0
 majority=iv(1)
 do jj=1,size(iv)
   nc=count(iv(:)==iv(jj))
   if(nc>nbest .or. (nc==nbest .and. iv(jj)<majority)) then
     nbest=nc
     majority=iv(jj)
   end if
 end do
 end function majority

 real(dp) function wall(kk,xx)
!! the wall, 1/2 kk xx**2 exp(-dwall/xx) for xx > 0, 0 inside
 real(dp), intent(in) :: kk,xx
 wall=0.0d0
 if(xx>0.0d0) wall=0.5d0*kk*xx**2*exp(-dwall/xx)
 end function wall

 real(dp) function target(jed)
!! atan(D) of an edge: 0, atan(-l-1), -pi/2
 integer, intent(in) :: jed
 if(jed==3) then
   target=-0.5d0*pi
 else
   target=atan(ws_edge_d(ll,jed))
 end if
 end function target

 subroutine sort_edges(ned,eed,ked,ied)
!! ascending in energy, by insertion
 integer, intent(in) :: ned
 integer, intent(inout) :: ked(mxedge),ied(mxedge)
 real(dp), intent(inout) :: eed(mxedge)
 real(dp) :: et
 integer :: jj,kk,kt,it
 do jj=2,ned
   et=eed(jj)
   kt=ked(jj)
   it=ied(jj)
   kk=jj-1
   do while(kk>=1)
     if(eed(kk)<=et) exit
     eed(kk+1)=eed(kk)
     ked(kk+1)=ked(kk)
     ied(kk+1)=ied(kk)
     kk=kk-1
   end do
   eed(kk+1)=et
   ked(kk+1)=kt
   ied(kk+1)=it
 end do
 end subroutine sort_edges

 integer function n_of(ex)
!! the all-electron state's nodes inside rs, plus l+1
 real(dp), intent(in) :: ex
 real(dp), allocatable :: uu(:,:),up(:,:)
 integer :: mch,jj,nnode,nn2,ierr2
 mch=phsft_mch(rr,mmax,rs)
 allocate(uu(mmax,2),up(mmax,2))
 uu(:,:)=0.0d0
 up(:,:)=0.0d0
 nn2=0
 call scatterstate(irel,nn2,ll,kap,ierr2,ex,rr,zz,vfull,uu,up,mmax,mch)
 nnode=0
 do jj=2,mmax
   if(rr(jj)>rs) exit
   if(uu(jj,1)*uu(jj-1,1)<0.0d0) nnode=nnode+1
 end do
 n_of=nnode+ll+1
 deallocate(uu,up)
 end function n_of

 subroutine cutoff_at(tol,ec)
!! the cutoff ec at which the error falls to tol (-1 if it does not):
!! between the first basis size below it and the one before, linear in
!! ln(err) against the cutoff
 real(dp), intent(in) :: tol
 real(dp), intent(out) :: ec
 integer :: jj
 ec=-1.0d0
 do jj=1,nb
   if(err(jj)==0.0d0 .and. jj<nb) cycle
   if(err(jj)>=tol) cycle
   if(jj==1 .or. err(jj)<=0.0d0) then
     ec=ecut(jj)
   else if(err(jj-1)<=0.0d0) then
     ec=ecut(jj)
   else
     ec=ecut(jj-1)+(ecut(jj)-ecut(jj-1)) &
&       *log(err(jj-1)/tol)/log(err(jj-1)/err(jj))
   end if
   write(6,'(a,4i4,1p,2e24.16)') 'ws_ecut',ll,ikap,iedge,nn,tol,ec
   return
 end do
 end subroutine cutoff_at

 real(dp) function summ(xx,inf)
!! a value of the summary: Infinity if inf, NaN if xx < 0 (none)
 real(dp), intent(in) :: xx
 logical, intent(in) :: inf
 if(inf) then
   summ=ieee_value(1.0d0,ieee_positive_inf)
 else if(xx<0.0d0) then
   summ=ieee_value(1.0d0,ieee_quiet_nan)
 else
   summ=xx
 end if
 end function summ

 subroutine write_help()
!! the usage: what the program does and its options
 write(6,'(a)') 'Usage: oncvpsp_confine.x -a AE_FILE -p PS_FILE [options]'
 write(6,'(a)') ''
 write(6,'(a)') 'A pseudopotential and its all-electron atom in a Wigner-Seitz'
 write(6,'(a)') 'sphere of radius rs, at their frozen potentials: for each l'
 write(6,'(a)') '(and j), the band edges of the Wigner-Seitz rule, where'
 write(6,'(a)') 'r phi''/phi = 0 at rs (bottom), -l-1 (centre) and phi = 0'
 write(6,'(a)') '(top), both atoms'' in the energy window; and the cutoff the'
 write(6,'(a)') 'pseudopotential''s state needs in the sphere, from its'
 write(6,'(a)') 'variational energy in the spherical Bessel functions of the'
 write(6,'(a)') 'same boundary condition. The output, on stdout:'
 write(6,'(a)') '  ws_r l ikap rs(bohr)'
 write(6,'(a)') '  ws l ikap edge n e_ae(Ha) e_ps(Ha)'
 write(6,'(a)') '  ws_ecut l ikap edge n tol(Ha) ecut(Ha)'
 write(6,'(a)') '  ws_floor l ikap edge n ecut(Ha) err(Ha)'
 write(6,'(a)') '  ws_conv l ikap edge n ecut(Ha) err(Ha)       (-t)'
 write(6,'(a)') '  ws_missing l ikap edge n e_ae(Ha)'
 write(6,'(a)') 'edge 1 bottom, 2 centre, 3 top; n the all-electron state''s'
 write(6,'(a)') 'nodes inside rs plus l+1; ikap 1 for j = l+1/2 (and every l'
 write(6,'(a)') 'if not fully relativistic), 2 for j = l-1/2. ws_ecut for the'
 write(6,'(a)') 'errors 1e-2 to 1e-6 Ha that are reached; ws_missing for an'
 write(6,'(a)') 'all-electron edge with no pseudopotential partner.'
 write(6,'(a)') '  ws_core l ikap edge n e_ae(Ha)'
 write(6,'(a)') '  ws_paired_off l ikap edge n e_ps(Ha)'
 write(6,'(a)') '  ws_extra l ikap edge e_ps(Ha)'
 write(6,'(a)') 'ws_core for one that is a core state''s band, ws_paired_off'
 write(6,'(a)') 'for a pseudopotential edge of another phase nearer it than'
 write(6,'(a)') 'its partner, ws_extra for a pseudopotential edge with no'
 write(6,'(a)') 'all-electron partner (a ghost). Then the walls:'
 write(6,'(a)') '  wall_r0 r0(bohr)'
 write(6,'(a)') '  wall k etot_ae etot_ps de_ae de_ps     (Ha, Ha/bohr**2)'
 write(6,'(a)') '  wall_eig k n l ikap e_ae(Ha) e_ps(Ha)'
 write(6,'(a)') '  wall_failed k ae|ps'
 write(6,'(a)') 'and last:'
 write(6,'(a)') '  confine_summary max_de ecut4 ecut5 n_missing n_extra'
 write(6,'(a)') '                  n_paired_off n_wall_failed max_dde'
 write(6,'(a)') '                  ecut4_top ecut5_top'
 write(6,'(a)') 'the largest |e_ps - e_ae| and the largest cutoffs at 1e-4'
 write(6,'(a)') 'and 1e-5 Ha of the edges below 0 (Infinity if one does not'
 write(6,'(a)') 'reach it), the counts of those lines, the largest'
 write(6,'(a)') '|de_ps - de_ae| of the walls, and the cutoffs of the top'
 write(6,'(a)') 'edges alone, free of the tail from the potential''s slope at'
 write(6,'(a)') 'rs; NaN and -1 for a part left out.'
 write(6,'(a)') ''
 write(6,'(a)') 'Options:'
 write(6,'(a)') '  -h, --help'
 write(6,'(a)') '      Show this help and exit'
 write(6,'(a)') '  -v, --version'
 write(6,'(a)') '      Show the build and exit'
 write(6,'(a)') '  -a, --ae-file FILE'
 write(6,'(a)') '      The all-electron data file of an _ae.x (required)'
 write(6,'(a)') '  -p, --ps-file FILE'
 write(6,'(a)') '      The pseudopotential data file of a run of oncvpsp.x,'
 write(6,'(a)') '      oncvpspnr.x, oncvpspr.x or their _ps.x with --ps-file,'
 write(6,'(a)') '      for the same atom, treatment and mesh (required)'
 write(6,'(a)') '  -r, --radii R1,R2,...'
 write(6,'(a)') '      The radii of the sphere (bohr), each beyond every core'
 write(6,'(a)') '      radius. Default 1.25, 1.5 and 2 times the largest'
 write(6,'(a)') '  -l, --channels L1,L2,...'
 write(6,'(a)') '      The angular momenta, from 0 to 3 (both j of each when'
 write(6,'(a)') '      fully relativistic). Default all four'
 write(6,'(a)') '  -e, --energies EMIN,EMAX,DE'
 write(6,'(a)') '      The window of the edges and the step of the scan of'
 write(6,'(a)') '      the phases (Ha). Default the input''s epsh1,epsh2, and'
 write(6,'(a)') '      0.1'
 write(6,'(a)') '  -c, --ecut-max ECUT'
 write(6,'(a)') '      The largest cutoff of the basis (Ha). Default 200'
 write(6,'(a)') '  -t, --table'
 write(6,'(a)') '      Also print the error at every basis size'
 write(6,'(a)') '  -s, --no-sphere'
 write(6,'(a)') '      Leave out the Wigner-Seitz sphere'
 write(6,'(a)') '  -w, --walls K1,K2,... | none'
 write(6,'(a)') '      The strengths k of the wall (Ha/bohr**2), or none.'
 write(6,'(a)') '      Default 0.01, 0.03 and 0.1'
 write(6,'(a)') '  -R, --wall-radius R0'
 write(6,'(a)') '      The onset of the wall (bohr). Default the largest'
 write(6,'(a)') '      core radius'
 write(6,'(a)') ''
 write(6,'(a)') 'Errors are reported on stdout, with a nonzero exit status.'
 end subroutine write_help

 subroutine next_arg(val)
!! the value of the option at iarg, which it passes
 character(len=1024), intent(out) :: val
 if(iarg+1>command_argument_count()) then
   call get_command_argument(iarg,val)
   write(6,'(3a)') 'oncvpsp_confine: ERROR ',trim(val),' requires a value'
   error stop 1
 end if
 iarg=iarg+1
 call get_command_argument(iarg,val)
 end subroutine next_arg

 integer function count_list(str)
!! the number of comma-separated items of str
 character(len=*), intent(in) :: str
 integer :: kk
 count_list=1
 do kk=1,len_trim(str)
   if(str(kk:kk)==',') count_list=count_list+1
 end do
 end function count_list

 elemental logical function differ(aa,bb)
!! .true. if aa and bb differ by more than rtol relative to the larger
 real(dp), intent(in) :: aa,bb
 differ=abs(aa-bb)>rtol*max(1.0d0,abs(aa),abs(bb))
 end function differ

 end program oncvpsp_confine

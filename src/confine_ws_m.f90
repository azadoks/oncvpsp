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
!> The Wigner-Seitz sphere of oncvpsp_confine.x (whose header describes
!> the method and the ws lines): for each channel, the band edges of both
!> atoms at each radius, their pairing, and the cutoffs of the
!> pseudopotential's states (ws_basis_m).
 module confine_ws_m
 use constants_m, only: dp, pi, mxstate
 use kappa_m, only: kap_of
 use ps_file_m, only: ps_file_t
 use phsft_m, only: phsft_t, phsft_ae_init, phsft_ps_init, phsft_free, &
&                   phsft_node_phase, phsft_node_phases, phsft_mch
 use ws_basis_m, only: ws_edge_d, ws_convergence, ws_eigenvalues, &
&                       ws_grid_init, ws_grid_t
 use ugrid_m, only: ug_brent
 use confine_summary_m, only: confine_summary_t
 implicit none
 private
 public :: confine_sphere

!> the most edges of a phase scan
 integer, parameter :: mxedge=400
!> the energy errors of the cutoffs (Ha)
 real(dp), parameter :: tols(5)=[1.0d-2,1.0d-3,1.0d-4,1.0d-5,1.0d-6]
!> the scan of the pseudopotential's phases runs this far (Ha) beyond the
!> energies of the all-electron edges, for their partners
 real(dp), parameter :: emargin=1.0d0
!> the points and spacing (Ha) of the phases' alignment about the
!> reference energy
 integer, parameter :: nanch=21
 real(dp), parameter :: danch=0.05d0

!> the poles of each sphere (edir, ndir of the radius at hand)
 type :: poles_t
   real(dp), allocatable :: e(:)
 end type poles_t

!> the edges of a scan: energy, multiple of pi and edge (1 bottom, 2
!> centre, 3 top)
 type :: edges_t
   integer :: n=0
   real(dp) :: e(mxedge)
   integer :: k(mxedge),ied(mxedge)
 end type edges_t

!> the functions of ug_brent: the solver, the point of the scan and the
!> phase of edges' root, at the radius ph_iq
 type(phsft_t), pointer :: psol=>null()
 real(dp) :: ph_e0,ph_tgt,ph_shift
 integer :: ph_iq
!> the poles of the pseudopotential's r phi'/phi at rs (Ha), and whether
!> the phase of psol counts them rather than its nodes
 integer :: ndir
 real(dp), allocatable :: edir(:)
 logical :: ph_poles

 contains

 subroutine confine_sphere(pf,irel,zz,mmax,rr,vfull,nc,na,la,chan,nq,rq, &
&                          emin,emax,de,ecutmax,table,sm,iout)
!! the ws lines of every channel of chan at the radii rq, of the edges
!! in [emin,emax] from a scan of step de, and the sphere's part of the
!! summary
 type(ps_file_t), intent(in) :: pf
 integer, intent(in) :: irel,mmax,nc,na(mxstate),la(mxstate),nq
 real(dp), intent(in) :: zz,rr(mmax),vfull(mmax),rq(nq),emin,emax,de
 real(dp), intent(in) :: ecutmax
 logical, intent(in) :: chan(0:3),table
 integer, intent(in) :: iout
 type(confine_summary_t), intent(inout) :: sm
 real(dp), allocatable :: egrid(:)
 integer :: npsh,ii,l1,ikap

 npsh=int(((emax-emin+2.0d0*emargin)/de)-0.5d0)+1
 allocate(egrid(npsh))
 do ii=1,npsh
   egrid(ii)=emax+emargin-de*dble(ii-1)
 end do

 write(iout,'(a)') '# ws_r l ikap rs(bohr)'
 write(iout,'(a)') '# ws l ikap edge n e_ae(Ha) e_ps(Ha)'
 write(iout,'(a)') '# ws_ecut l ikap edge n tol(Ha) ecut(Ha)'
 write(iout,'(a)') '# ws_floor l ikap edge n ecut(Ha) err(Ha)'
 if(table) write(iout,'(a)') '# ws_conv l ikap edge n ecut(Ha) err(Ha)'
 write(iout,'(a)') '# ws_missing l ikap edge n e_ae(Ha)'
 write(iout,'(a)') '# ws_core l ikap edge n e_ae(Ha)'
 write(iout,'(a)') '# ws_paired_off l ikap edge n e_ps(Ha)'
 write(iout,'(a)') '# ws_extra l ikap edge e_ps(Ha)'
 sm%nmiss=0
 sm%nextra=0
 sm%noff=0
 do l1=1,4
   if(.not. chan(l1-1)) cycle
   do ikap=1,merge(1,pf%nkap,l1==1)
     call ws_channel(pf,irel,zz,mmax,rr,vfull,nc,na,la,l1-1,ikap,nq,rq, &
&                    egrid,emin,emax,ecutmax,table,sm,iout)
   end do
 end do
 deallocate(egrid)
 end subroutine confine_sphere

 subroutine ws_channel(pf,irel,zz,mmax,rr,vfull,nc,na,la,ll,ikap,nq,rq, &
&                      egrid,emin,emax,ecutmax,table,sm,iout)
!! the ws lines of l and ikap at every radius, from one scan of the
!! phases of both atoms on egrid
 type(ps_file_t), intent(in) :: pf
 integer, intent(in) :: irel,mmax,nc,na(mxstate),la(mxstate),ll,ikap,nq
 real(dp), intent(in) :: zz,rr(mmax),vfull(mmax),rq(nq),egrid(:),emin
 real(dp), intent(in) :: emax,ecutmax
 logical, intent(in) :: table
!> unit of the text output
 integer, intent(in) :: iout
 type(confine_summary_t), intent(inout) :: sm
 type(phsft_t), target :: sae,sps
!> the quadrature of the channel in each sphere, shared by its edges
 type(ws_grid_t), allocatable :: wsg(:)
 type(poles_t), allocatable :: polq(:)
!> the phases of every radius on the grid and at the anchors
 real(dp), allocatable :: shae(:,:),shps(:,:),tanae(:,:),tanps(:,:)
 real(dp) :: ex
 integer :: l1,kap,nvkb,iq,ii,ierr,npsh
 integer :: ifail

 l1=ll+1
 npsh=size(egrid)
 kap=0
 if(irel==2) kap=kap_of(ll,ikap)
 nvkb=pf%nvkb(l1,ikap)
! the solvers of the channel at every radius, and each sphere's
! quadrature and poles (refined to the roots of u(rs))
 call phsft_ae_init(irel,ll,kap,rr,vfull,zz,mmax,nq,rq(1:nq),sae,iout,ifail)
 if(ifail/=0) error stop 1
 call phsft_ps_init(ll,nvkb,rr,pf%vloc(:,ikap),pf%vlocpw, &
&                   pf%vkbpw(1:nvkb,l1,ikap),pf%evkb(1:nvkb,l1,ikap), &
&                   mmax,nq,rq(1:nq),sps,iout,ifail)
 if(ifail/=0) error stop 1
 allocate(wsg(nq),polq(nq),tanae(nanch,nq),tanps(nanch,nq))
 allocate(shae(npsh,nq),shps(npsh,nq))
 do iq=1,nq
   call ws_grid_init(rq(iq),sqrt(2.0d0*ecutmax),rr,mmax, &
&       pf%vloc(:,ikap),pf%vlocpw,nvkb,pf%vkbpw(1:nvkb,l1,ikap), &
&       wsg(iq))
   call ws_eigenvalues(ll,3,wsg(iq),pf%evkb(1:nvkb,l1,ikap),ndir, &
&       edir,ierr)
   if(ierr/=0) ndir=0
   ph_iq=iq
   call refine_poles(sps,egrid)
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
   call scan_phases(sae,sps,polq,egrid(ii),shae(ii,:),shps(ii,:),iout)
 end do
 if(nvkb>0) then
   do ii=1,nanch
     ex=pf%epa(l1,ikap)+danch*dble(ii-(nanch+1)/2)
     call scan_phases(sae,sps,polq,ex,tanae(ii,:),tanps(ii,:),iout)
   end do
 end if

 do iq=1,nq
   call ws_radius(pf,irel,zz,mmax,rr,vfull,nc,na,la,ll,kap,ikap,nvkb, &
&                 rq(iq),iq,sae,sps,polq(iq),wsg(iq),egrid,shae(:,iq), &
&                 shps(:,iq),tanae(:,iq),tanps(:,iq),emin,emax,table,sm,iout)
 end do
 call phsft_free(sae)
 call phsft_free(sps)
 deallocate(wsg,polq,tanae,tanps,shae,shps)
 end subroutine ws_channel

 subroutine ws_radius(pf,irel,zz,mmax,rr,vfull,nc,na,la,ll,kap,ikap,nvkb, &
&                     rs,iq,sae,sps,pol,wsg,egrid,shae,shps,tanae,tanps, &
&                     emin,emax,table,sm,iout)
!! the ws lines of l and ikap at the radius rs (the iq-th of the
!! solvers), from the phases shae and shps of the scan and tanae and
!! tanps of the anchors
 type(ps_file_t), intent(in) :: pf
 integer, intent(in) :: irel,mmax,nc,na(mxstate),la(mxstate),ll,kap,ikap
 integer, intent(in) :: nvkb,iq
 real(dp), intent(in) :: zz,rr(mmax),vfull(mmax),rs,egrid(:),shae(:)
 real(dp), intent(in) :: tanae(nanch),tanps(nanch),emin,emax
 real(dp), intent(inout) :: shps(:)
 type(phsft_t), intent(inout), target :: sae,sps
 type(poles_t), intent(in) :: pol
 type(ws_grid_t), intent(in) :: wsg
 logical, intent(in) :: table
!> unit of the text output
 integer, intent(in) :: iout
 type(confine_summary_t), intent(inout) :: sm
 type(edges_t) :: eae,eps
 real(dp), allocatable :: ecut(:),err(:)
 real(dp) :: thal
 integer :: ianch(nanch),ishift(size(egrid))
!> the partner of each all-electron edge (0 for none), and whether each
!> pseudopotential edge is one
 integer :: jpair(mxedge)
 logical :: used(mxedge)
 integer :: npsh,ii,ie,je,iedge,nn,nb,ierr

 npsh=size(egrid)
 write(iout,'(a,2i4,1p,e24.16)') 'ws_r',ll,ikap,rs
 ph_iq=iq
 ndir=size(pol%e)
 edir=pol%e
 do ii=1,npsh
   ishift(ii)=nint((shae(ii)-shps(ii))/pi)
 end do
 if(nvkb>0) then
   do ii=1,nanch
     ianch(ii)=nint((tanae(ii)-tanps(ii))/pi)
   end do
   thal=pi*dble(majority(ianch))
 else
   thal=pi*dble(majority(ishift))
 end if
 shps(:)=shps(:)+thal

 call edges(sae,.false.,ll,egrid,shae,0.0d0,eae)
 call edges(sps,.true.,ll,egrid,shps,thal,eps)
 deallocate(edir)

! every all-electron edge, lowest first, and its partner: the nearest
! pseudopotential edge of the same phase, over the whole scan, so that an
! edge in the window whose partner is outside it is not without one
 call sort_edges(eae)
 call sort_edges(eps)
 used(1:eps%n)=.false.
 do ie=1,eae%n
   jpair(ie)=0
   do ii=1,eps%n
     if(eps%k(ii)/=eae%k(ie) .or. eps%ied(ii)/=eae%ied(ie)) cycle
     if(jpair(ie)==0) then
       jpair(ie)=ii
     else if(abs(eps%e(ii)-eae%e(ie))<abs(eps%e(jpair(ie))-eae%e(ie))) then
       jpair(ie)=ii
     end if
   end do
   if(jpair(ie)>0) used(jpair(ie))=.true.
 end do

! those in the window
 do ie=1,eae%n
   if(eae%e(ie)<emin .or. eae%e(ie)>emax) cycle
   iedge=eae%ied(ie)
   nn=n_of(irel,ll,kap,eae%e(ie),rr,zz,vfull,mmax,rs,iout)
   je=jpair(ie)
   if(je==0) then
! the band of a core state, which the pseudopotential leaves out
     if(any(na(1:nc)==nn .and. la(1:nc)==ll)) then
       write(iout,'(a,4i4,1p,e24.16)') 'ws_core',ll,ikap,iedge,nn,eae%e(ie)
     else
       write(iout,'(a,4i4,1p,e24.16)') 'ws_missing',ll,ikap,iedge,nn,eae%e(ie)
       sm%nmiss=sm%nmiss+1
     end if
     cycle
   end if
   write(iout,'(a,4i4,1p,2e24.16)') 'ws',ll,ikap,iedge,nn,eae%e(ie),eps%e(je)
! an edge of another phase nearer to it than its partner: the phases pair
! a band of one atom with another band of the other
   do ii=1,eps%n
     if(ii==je .or. eps%ied(ii)/=iedge) cycle
     if(abs(eps%e(ii)-eae%e(ie))<abs(eps%e(je)-eae%e(ie))) then
       write(iout,'(a,4i4,1p,e24.16)') 'ws_paired_off',ll,ikap,iedge,nn, &
&           eps%e(ii)
       sm%noff=sm%noff+1
       exit
     end if
   end do
   if(eae%e(ie)<0.0d0) sm%smax_de=max(sm%smax_de,abs(eps%e(je)-eae%e(ie)))

! the cutoff of the pseudopotential's state
   call ws_convergence(ll,iedge,wsg,pf%evkb(1:nvkb,ll+1,ikap),eps%e(je), &
&       nb,ecut,err,ierr)
   if(ierr/=0) then
     if(eae%e(ie)<0.0d0) sm%sinf(1:2)=.true.
     if(eae%e(ie)<0.0d0 .and. iedge==3) sm%sinf(3:4)=.true.
     cycle
   end if
   call ws_cutoffs(ll,ikap,iedge,nn,eae%e(ie)<0.0d0,nb,ecut,err,table,sm,iout)
 end do

! the pseudopotential's edges in the window that no all-electron edge
! has as its partner
 do ii=1,eps%n
   if(used(ii) .or. eps%e(ii)<emin .or. eps%e(ii)>emax) cycle
   write(iout,'(a,3i4,1p,e24.16)') 'ws_extra',ll,ikap,eps%ied(ii),eps%e(ii)
   sm%nextra=sm%nextra+1
 end do
 end subroutine ws_radius

 subroutine ws_cutoffs(ll,ikap,iedge,nn,bound,nb,ecut,err,table,sm,iout)
!! the ws_ecut, ws_floor (and with table ws_conv) lines of an edge from
!! the errors err of its state at the cutoffs ecut, and the summary's
!! cutoffs, at 1e-4 and 1e-5 Ha, of all edges and of the top, of a bound
!! one
 integer, intent(in) :: ll,ikap,iedge,nn,nb
 logical, intent(in) :: bound,table
 real(dp), intent(in) :: ecut(nb),err(nb)
 integer, intent(in) :: iout
 type(confine_summary_t), intent(inout) :: sm
 real(dp) :: ec
 integer :: itol,kk,ii
 do itol=1,size(tols)
   call cutoff_at(tols(itol),ec)
   if(bound .and. (itol==3 .or. itol==4)) then
     do kk=itol-2,merge(itol,itol-2,iedge==3),2
       if(ec<0.0d0) sm%sinf(kk)=.true.
       sm%scut(kk)=max(sm%scut(kk),ec)
     end do
   end if
 end do
 write(iout,'(a,4i4,1p,2e24.16)') 'ws_floor',ll,ikap,iedge,nn, &
&     ecut(nb),err(nb)
 if(table) then
   do ii=1,nb
     if(err(ii)==0.0d0) cycle
     write(iout,'(a,4i4,1p,2e24.16)') 'ws_conv',ll,ikap,iedge,nn, &
&         ecut(ii),err(ii)
   end do
 end if

 contains

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
   write(iout,'(a,4i4,1p,2e24.16)') 'ws_ecut',ll,ikap,iedge,nn,tol,ec
   return
 end do
 end subroutine cutoff_at

 end subroutine ws_cutoffs

 subroutine edges(sol,poles,ll,egrid,sh,shift,ed)
!! the edges of a scan: the energies at which the phase (plus shift)
!! crosses atan(D) + k pi for the D of each edge, each to round-off by
!! Brent's method in the interval of the grid that brackets it; poles for
!! the phase that counts the poles edir rather than the nodes, whose top
!! edges (phi(rs) = 0) are the poles themselves, where the count steps: the
!! k-th at the phase -pi/2 - (k-1) pi just below it
 type(phsft_t), intent(inout), target :: sol
 logical, intent(in) :: poles
 integer, intent(in) :: ll
 real(dp), intent(in) :: egrid(:),sh(:),shift
 type(edges_t), intent(out) :: ed
 real(dp) :: thlo,thhi,er
 integer :: jj,kk,jed,ierr2,npsh
 npsh=size(egrid)
 psol=>sol
 ph_poles=poles
 ph_shift=shift
 ed%n=0
 do jj=1,npsh-1
   ph_e0=egrid(jj)
   thlo=min(sh(jj),sh(jj+1))
   thhi=max(sh(jj),sh(jj+1))
   do jed=1,3
     if(poles .and. jed==3) cycle
     do kk=ceiling((thlo-target(ll,jed))/pi),floor((thhi-target(ll,jed))/pi)
       ph_tgt=target(ll,jed)+pi*dble(kk)
! a crossing on the point below is the next interval's
       if(ph_tgt==sh(jj+1) .and. jj<npsh-1) cycle
       call ug_brent(phase,egrid(jj+1),ph_e0,er,ierr2)
       if(ierr2/=0) cycle
! a jump (a pair of nodes inside the projectors), not a root
       if(abs(phase(er))>0.25d0*pi) cycle
       call add_edge(ed,er,kk,jed)
     end do
   end do
 end do
 if(poles) then
   do kk=1,ndir
     if(edir(kk)<egrid(npsh) .or. edir(kk)>egrid(1)) cycle
     call add_edge(ed,edir(kk),nint(shift/pi)-(kk-1),3)
   end do
 end if
 nullify(psol)
 end subroutine edges

 subroutine add_edge(ed,ee,kk,jed)
!! an edge more
 type(edges_t), intent(inout) :: ed
 real(dp), intent(in) :: ee
 integer, intent(in) :: kk,jed
 if(ed%n==mxedge) then
   write(6,'(a)') 'oncvpsp_confine: ERROR too many edges in the scan'
   error stop 1
 end if
 ed%n=ed%n+1
 ed%e(ed%n)=ee
 ed%k(ed%n)=kk
 ed%ied(ed%n)=jed
 end subroutine add_edge

 subroutine scan_phases(sae,sps,polq,ex,thae,thps,iout)
!! the phases of every radius at ex from one solution of each atom: the
!! all-electron one's, and the pseudopotential's less pi per pole of each
!! sphere below ex (phase_of)
 type(phsft_t), intent(inout) :: sae,sps
 type(poles_t), intent(in) :: polq(:)
 real(dp), intent(in) :: ex
 integer, intent(in) :: iout
 real(dp), intent(out) :: thae(:),thps(:)
 integer :: jq
 integer :: ifail
 call phsft_node_phases(sae,ex,thae,iout=iout,ifail=ifail)
 if(ifail/=0) error stop 1
 call phsft_node_phases(sps,ex,thps,nodes=.false.,iout=iout,ifail=ifail)
 if(ifail/=0) error stop 1
 do jq=1,size(polq)
   thps(jq)=thps(jq)-pi*dble(count(polq(jq)%e<ex))
 end do
 end subroutine scan_phases

 real(dp) function phase(ex)
!! the phase of psol at ex plus ph_shift, less ph_tgt
 real(dp), intent(in) :: ex
 integer, parameter :: iout=6
 phase=phase_of(ex,iout)+ph_shift-ph_tgt
 end function phase

 real(dp) function phase_of(ex,iout)
!! the phase of psol at ex: less pi per node inside rs, or with ph_poles
!! per pole of r phi'/phi at rs below ex
 real(dp), intent(in) :: ex
 integer, intent(in) :: iout
 integer :: ifail
 if(ph_poles) then
   phase_of=phsft_node_phase(psol,ph_iq,ex,nodes=.false., &
&        iout=iout,ifail=ifail) &
&           -pi*dble(count(edir(1:ndir)<ex))
   if(ifail/=0) error stop 1
 else
   phase_of=phsft_node_phase(psol,ph_iq,ex,iout=iout,ifail=ifail)
   if(ifail/=0) error stop 1
 end if
 end function phase_of

 subroutine refine_poles(sol,egrid)
!! the poles edir in the scan, each from its variational value to the root
!! of u(rs) beside it, in a bracket that grows from 1e-6 Ha on both sides
!! until u changes sign, up to halfway to the poles beside it (or 1 Ha). A
!! variational bound lies above the root of the same Hamiltonian, but
!! the basis's quadrature of the pieces and the radial equation's differ,
!! by 3e-8 Ha at 1.5 Ha in Si f, so the root can be on either side
 type(phsft_t), intent(inout), target :: sol
 real(dp), intent(in) :: egrid(:)
 real(dp) :: e0,u0,lo,hi,ww,er,elow,ehigh
 integer :: kk,ierr2
 logical :: found
 psol=>sol
 do kk=1,ndir
   e0=edir(kk)
   if(e0<egrid(size(egrid)) .or. e0>egrid(1)) cycle
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
 integer, parameter :: iout=6
 real(dp) :: th
 integer :: ifail
 th=phsft_node_phase(psol,ph_iq,ex,nodes=.false.,uq_out=u_rs, &
&      iout=iout,ifail=ifail)
 if(ifail/=0) error stop 1
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

 real(dp) function target(ll,jed)
!! atan(D) of an edge of l: 0, atan(-l-1), -pi/2
 integer, intent(in) :: ll,jed
 if(jed==3) then
   target=-0.5d0*pi
 else
   target=atan(ws_edge_d(ll,jed))
 end if
 end function target

 subroutine sort_edges(ed)
!! ascending in energy, by insertion
 type(edges_t), intent(inout) :: ed
 real(dp) :: et
 integer :: jj,kk,kt,it
 do jj=2,ed%n
   et=ed%e(jj)
   kt=ed%k(jj)
   it=ed%ied(jj)
   kk=jj-1
   do while(kk>=1)
     if(ed%e(kk)<=et) exit
     ed%e(kk+1)=ed%e(kk)
     ed%k(kk+1)=ed%k(kk)
     ed%ied(kk+1)=ed%ied(kk)
     kk=kk-1
   end do
   ed%e(kk+1)=et
   ed%k(kk+1)=kt
   ed%ied(kk+1)=it
 end do
 end subroutine sort_edges

 integer function n_of(irel,ll,kap,ex,rr,zz,vfull,mmax,rs,iout)
!! the all-electron state's nodes inside rs, plus l+1
 use scatterstate_m, only: scatterstate
 integer, intent(in) :: irel,ll,kap,mmax
 real(dp), intent(in) :: ex,rr(mmax),zz,vfull(mmax),rs
 integer, intent(in) :: iout
 real(dp), allocatable :: uu(:,:),up(:,:)
 integer :: mch,jj,nnode,nn2,ierr2
 mch=phsft_mch(rr,mmax,rs)
 allocate(uu(mmax,2),up(mmax,2))
 uu(:,:)=0.0d0
 up(:,:)=0.0d0
 nn2=0
 call scatterstate(irel,nn2,ll,kap,ierr2,ex,rr,zz,vfull,uu,up,mmax,mch,iout)
 nnode=0
 do jj=2,mmax
   if(rr(jj)>rs) exit
   if(uu(jj,1)*uu(jj-1,1)<0.0d0) nnode=nnode+1
 end do
 n_of=nnode+ll+1
 deallocate(uu,up)
 end function n_of

 end module confine_ws_m

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
!> the number of points of the phases' alignment about the reference
!> energy
 integer, parameter :: nanch=21
!> the spacing of the points of the phases' alignment (Ha)
 real(dp), parameter :: danch=0.05d0

!> the poles of each sphere (edir, ndir of the radius at hand)
 type :: poles_t
   real(dp), allocatable :: ee(:)  !< the poles, ascending (Ha)
 end type poles_t

!> the edges of a scan: energy, multiple of pi and edge (1 bottom, 2
!> centre, 3 top)
 type :: edges_t
   integer :: ne=0  !< number of edges
   real(dp) :: ee(mxedge)  !< energy of each edge (Ha)
   integer :: iband(mxedge)  !< multiple of pi of each edge
   integer :: ied(mxedge)  !< edge of each: 1 bottom, 2 centre, 3 top
 end type edges_t

 contains

 subroutine confine_sphere(pf,irel,c_light,zz,mmax,adams,rr,vfull,nc,na,la, &
&                          chan,nq, &
&                          rq,emin,emax,de,ecutmax,table,sm,iout,ifail)
!! the ws lines of every channel of chan at the radii rq, of the edges
!! in [emin,emax] from a scan of step de, and the sphere's part of the
!! summary
!!
!! @param pf  the pseudopotential data file
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param rr  log radial mesh (bohr)
!! @param vfull  all-electron potential (Ha)
!! @param nc  number of core states
!! @param na  principal quantum numbers
!! @param la  angular momenta
!! @param chan  whether to treat each l, 0..3
!! @param nq  number of sphere radii
!! @param rq  sphere radii (bohr)
!! @param emin,emax  energy window of the ws lines (Ha); the scan runs
!!                   emargin beyond it
!! @param de  energy step of the scan (Ha)
!! @param ecutmax  largest cutoff of the basis (Ha)
!! @param table  whether to write the ws_conv lines
!! @param sm  the summary; its sphere part set
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use confine_summary_m, only: confine_summary_t
 use ps_file_m, only: ps_file_t
 use adams_m, only: adams_t

!Input variables
 type(ps_file_t), intent(in) :: pf
 integer, intent(in) :: irel,mmax,nc,na(mxstate),la(mxstate),nq
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: zz,rr(mmax),vfull(mmax),rq(nq),emin,emax,de
 real(dp), intent(in) :: ecutmax
 logical, intent(in) :: chan(0:3),table
 integer, intent(in) :: iout

!Input/output variables
 type(confine_summary_t), intent(inout) :: sm

!Output variables
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: egrid(:)
 integer :: npsh,ii,l1,ikap

 ifail=0
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
   do ikap=1,merge(1,pf%pot%nkap,l1==1)
     call ws_channel(pf,irel,c_light,zz,mmax,adams,rr,vfull,nc,na,la,l1-1, &
&                    ikap,nq, &
&                    rq,egrid,emin,emax,ecutmax,table,sm,iout,ifail)
     if(ifail/=0) return
   end do
 end do
 deallocate(egrid)
 end subroutine confine_sphere

 subroutine ws_channel(pf,irel,c_light,zz,mmax,adams,rr,vfull,nc,na,la,ll, &
&                      ikap,nq, &
&                      rq,egrid,emin,emax,ecutmax,table,sm,iout,ifail)
!! the ws lines of l and ikap at every radius, from one scan of the
!! phases of both atoms on egrid
!!
!! @param pf  the pseudopotential data file
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param rr  log radial mesh (bohr)
!! @param vfull  all-electron potential (Ha)
!! @param nc  number of core states
!! @param na  principal quantum numbers
!! @param la  angular momenta
!! @param ll  angular momentum
!! @param ikap  index of the kappa of l, 1 or 2
!! @param nq  number of sphere radii
!! @param rq  sphere radii (bohr)
!! @param egrid  energies of the scan, descending (Ha)
!! @param emin,emax  energy window of the ws lines (Ha)
!! @param ecutmax  largest cutoff of the basis (Ha)
!! @param table  whether to write the ws_conv lines
!! @param sm  the summary; its sphere part updated
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use confine_summary_m, only: confine_summary_t
 use kappa_m, only: kap_of
 use phsft_m, only: phsft_t, phsft_ae_init, phsft_ps_init, phsft_free
 use ps_file_m, only: ps_file_t
 use ws_basis_m, only: ws_eigenvalues, ws_grid_init, ws_grid_t
 use adams_m, only: adams_t

!Input variables
 type(ps_file_t), intent(in) :: pf
 integer, intent(in) :: irel,mmax,nc,na(mxstate),la(mxstate),ll,ikap,nq
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: zz,rr(mmax),vfull(mmax),rq(nq),egrid(:),emin
 real(dp), intent(in) :: emax,ecutmax
 logical, intent(in) :: table
!> unit of the text output
 integer, intent(in) :: iout

!Input/output variables
 type(confine_summary_t), intent(inout) :: sm

!Output variables
 integer, intent(out) :: ifail

!Local variables
 type(phsft_t) :: sae,sps
!> the quadrature of the channel in each sphere, shared by its edges
 type(ws_grid_t), allocatable :: wsg(:)
 type(poles_t), allocatable :: polq(:)
!> the phases of every radius on the grid and at the anchors
 real(dp), allocatable :: shae(:,:),shps(:,:),tanae(:,:),tanps(:,:)
 real(dp) :: ex
!> the poles of the pseudopotential's r phi'/phi at a radius (Ha)
 real(dp), allocatable :: edir(:)
 integer :: l1,kap,nvkb,iq,ii,ierr,npsh,ndir

 l1=ll+1
 npsh=size(egrid)
 kap=0
 if(irel==2) kap=kap_of(ll,ikap)
 nvkb=pf%nvkb(l1,ikap)
! the solvers of the channel at every radius, and each sphere's
! quadrature and poles (refined to the roots of u(rs))
 call phsft_ae_init(irel,c_light,ll,kap,rr,vfull,zz,mmax,adams,nq,rq(1:nq), &
&                   sae,iout, &
&                   ifail)
 if(ifail/=0) return
 call phsft_ps_init(ll,nvkb,rr,pf%vloc(:,ikap),pf%pot%vlocpw, &
&                   pf%pot%vkbpw(1:nvkb,l1,ikap),pf%pot%evkb(1:nvkb,l1,ikap), &
&                   mmax,adams,nq,rq(1:nq),sps,iout,ifail)
 if(ifail/=0) return
 allocate(wsg(nq),polq(nq),tanae(nanch,nq),tanps(nanch,nq))
 allocate(shae(npsh,nq),shps(npsh,nq))
 do iq=1,nq
   call ws_grid_init(rq(iq),sqrt(2.0d0*ecutmax),rr,mmax, &
&       pf%vloc(:,ikap),pf%pot%vlocpw,nvkb,pf%pot%vkbpw(1:nvkb,l1,ikap), &
&       wsg(iq),ifail)
   if(ifail/=0) return
   call ws_eigenvalues(ll,3,wsg(iq),pf%pot%evkb(1:nvkb,l1,ikap),ndir, &
&       edir,ierr,ifail)
   if(ifail/=0) return
   if(ierr/=0) ndir=0
   call refine_poles(sps,iq,egrid,edir(1:ndir),iout,ifail)
   if(ifail/=0) return
   polq(iq)%ee=edir(1:ndir)
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
   call scan_phases(sae,sps,polq,egrid(ii),shae(ii,:),shps(ii,:),iout, &
&                   ifail)
   if(ifail/=0) return
 end do
 if(nvkb>0) then
   do ii=1,nanch
     ex=pf%epa(l1,ikap)+danch*dble(ii-(nanch+1)/2)
     call scan_phases(sae,sps,polq,ex,tanae(ii,:),tanps(ii,:),iout,ifail)
     if(ifail/=0) return
   end do
 end if

 do iq=1,nq
   call ws_radius(pf,irel,c_light,zz,mmax,adams,rr,vfull,nc,na,la,ll,kap,ikap, &
&                 nvkb, &
&                 rq(iq),iq,sae,sps,polq(iq),wsg(iq),egrid,shae(:,iq), &
&                 shps(:,iq),tanae(:,iq),tanps(:,iq),emin,emax,table,sm,iout, &
&                 ifail)
   if(ifail/=0) return
 end do
 call phsft_free(sae)
 call phsft_free(sps)
 deallocate(wsg,polq,tanae,tanps,shae,shps)
 end subroutine ws_channel

 subroutine ws_radius(pf,irel,c_light,zz,mmax,adams,rr,vfull,nc,na,la,ll,kap, &
&                     ikap, &
&                     nvkb,rs,iq,sae,sps,pol,wsg,egrid,shae,shps,tanae,tanps, &
&                     emin,emax,table,sm,iout,ifail)
!! the ws lines of l and ikap at the radius rs (the iq-th of the
!! solvers), from the phases shae and shps of the scan and tanae and
!! tanps of the anchors
!!
!! @param pf  the pseudopotential data file
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param rr  log radial mesh (bohr)
!! @param vfull  all-electron potential (Ha)
!! @param nc  number of core states
!! @param na  principal quantum numbers
!! @param la  angular momenta
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2; 0 unless irel=2
!! @param ikap  index of the kappa of l, 1 or 2
!! @param nvkb  number of projectors of the channel
!! @param rs  radius of the sphere (bohr)
!! @param iq  index of rs among the radii of the solvers
!! @param sae  the all-electron phase solver
!! @param sps  the pseudopotential's phase solver
!! @param pol  the poles of the sphere (Ha)
!! @param wsg  the quadrature of the channel in the sphere
!! @param egrid  energies of the scan, descending (Ha)
!! @param shae  all-electron phases on egrid
!! @param shps  pseudopotential phases on egrid; the multiple of pi that
!!              aligns them with shae added on return
!! @param tanae,tanps  all-electron and pseudopotential phases at the
!!                     nanch alignment energies
!! @param emin,emax  energy window of the ws lines (Ha)
!! @param table  whether to write the ws_conv lines
!! @param sm  the summary; its sphere part updated
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use confine_summary_m, only: confine_summary_t
 use phsft_m, only: phsft_t
 use ps_file_m, only: ps_file_t
 use ws_basis_m, only: ws_convergence, ws_grid_t
 use adams_m, only: adams_t

!Input variables
 type(ps_file_t), intent(in) :: pf
 integer, intent(in) :: irel,mmax,nc,na(mxstate),la(mxstate),ll,kap,ikap
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 integer, intent(in) :: nvkb,iq
 real(dp), intent(in) :: zz,rr(mmax),vfull(mmax),rs,egrid(:),shae(:)
 real(dp), intent(in) :: tanae(nanch),tanps(nanch),emin,emax
 type(poles_t), intent(in) :: pol
 type(ws_grid_t), intent(in) :: wsg
 logical, intent(in) :: table
!> unit of the text output
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: shps(:)
 type(phsft_t), intent(inout) :: sae,sps
 type(confine_summary_t), intent(inout) :: sm

!Output variables
 integer, intent(out) :: ifail

!Local variables
 type(edges_t) :: eae,eps
 real(dp), allocatable :: ecut(:),err(:)
 real(dp) :: thal
 integer :: ianch(nanch),ishift(size(egrid))
!> the partner of each all-electron edge (0 for none), and whether each
!> pseudopotential edge is one
 integer :: jpair(mxedge)
 logical :: used(mxedge)
 integer :: npsh,ii,ie,je,iedge,nn,nb,ierr

 ifail=0
 npsh=size(egrid)
 write(iout,'(a,2i4,1p,e24.16)') 'ws_r',ll,ikap,rs
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

 call edges(sae,iq,.false.,pol%ee,ll,egrid,shae,0.0d0,eae,iout,ifail)
 if(ifail/=0) return
 call edges(sps,iq,.true.,pol%ee,ll,egrid,shps,thal,eps,iout,ifail)
 if(ifail/=0) return

! every all-electron edge, lowest first, and its partner: the nearest
! pseudopotential edge of the same phase, over the whole scan, so that an
! edge in the window whose partner is outside it is not without one
 call sort_edges(eae)
 call sort_edges(eps)
 used(1:eps%ne)=.false.
 do ie=1,eae%ne
   jpair(ie)=0
   do ii=1,eps%ne
     if(eps%iband(ii)/=eae%iband(ie) .or. eps%ied(ii)/=eae%ied(ie)) cycle
     if(jpair(ie)==0) then
       jpair(ie)=ii
     else if(abs(eps%ee(ii)-eae%ee(ie))<abs(eps%ee(jpair(ie))-eae%ee(ie))) then
       jpair(ie)=ii
     end if
   end do
   if(jpair(ie)>0) used(jpair(ie))=.true.
 end do

! those in the window
 do ie=1,eae%ne
   if(eae%ee(ie)<emin .or. eae%ee(ie)>emax) cycle
   iedge=eae%ied(ie)
   call n_of(irel,c_light,ll,kap,eae%ee(ie),rr,zz,vfull,mmax,adams,rs,nn,ifail)
   if(ifail/=0) return
   je=jpair(ie)
   if(je==0) then
! the band of a core state, which the pseudopotential leaves out
     if(any(na(1:nc)==nn .and. la(1:nc)==ll)) then
       write(iout,'(a,4i4,1p,e24.16)') 'ws_core',ll,ikap,iedge,nn,eae%ee(ie)
     else
       write(iout,'(a,4i4,1p,e24.16)') 'ws_missing',ll,ikap,iedge,nn,eae%ee(ie)
       sm%nmiss=sm%nmiss+1
     end if
     cycle
   end if
   write(iout,'(a,4i4,1p,2e24.16)') 'ws',ll,ikap,iedge,nn,eae%ee(ie),eps%ee(je)
! an edge of another phase nearer to it than its partner: the phases pair
! a band of one atom with another band of the other
   do ii=1,eps%ne
     if(ii==je .or. eps%ied(ii)/=iedge) cycle
     if(abs(eps%ee(ii)-eae%ee(ie))<abs(eps%ee(je)-eae%ee(ie))) then
       write(iout,'(a,4i4,1p,e24.16)') 'ws_paired_off',ll,ikap,iedge,nn, &
&           eps%ee(ii)
       sm%noff=sm%noff+1
       exit
     end if
   end do
   if(eae%ee(ie)<0.0d0) sm%smax_de=max(sm%smax_de,abs(eps%ee(je)-eae%ee(ie)))

! the cutoff of the pseudopotential's state
   call ws_convergence(ll,iedge,wsg,pf%pot%evkb(1:nvkb,ll+1,ikap),eps%ee(je), &
&       nb,ecut,err,ierr,ifail)
   if(ifail/=0) return
   if(ierr/=0) then
     if(eae%ee(ie)<0.0d0) sm%sinf(1:2)=.true.
     if(eae%ee(ie)<0.0d0 .and. iedge==3) sm%sinf(3:4)=.true.
     cycle
   end if
   call ws_cutoffs(ll,ikap,iedge,nn,eae%ee(ie)<0.0d0,nb,ecut,err,table,sm,iout)
 end do

! the pseudopotential's edges in the window that no all-electron edge
! has as its partner
 do ii=1,eps%ne
   if(used(ii) .or. eps%ee(ii)<emin .or. eps%ee(ii)>emax) cycle
   write(iout,'(a,3i4,1p,e24.16)') 'ws_extra',ll,ikap,eps%ied(ii),eps%ee(ii)
   sm%nextra=sm%nextra+1
 end do
 end subroutine ws_radius

 subroutine ws_cutoffs(ll,ikap,iedge,nn,bound,nb,ecut,err,table,sm,iout)
!! the ws_ecut, ws_floor (and with table ws_conv) lines of an edge from
!! the errors err of its state at the cutoffs ecut, and the summary's
!! cutoffs, at 1e-4 and 1e-5 Ha, of all edges and of the top, of a bound
!! one
!!
!! @param ll  angular momentum
!! @param ikap  index of the kappa of l, 1 or 2
!! @param iedge  the edge: 1 bottom, 2 centre, 3 top
!! @param nn  all-electron state's nodes inside rs, plus l+1
!! @param bound  whether the all-electron edge is below 0, so that its
!!               cutoffs enter the summary
!! @param nb  number of basis sizes
!! @param ecut  cutoff of each basis size (Ha)
!! @param err  energy error of the state at each basis size (Ha)
!! @param table  whether to write the ws_conv lines
!! @param sm  the summary; scut and sinf updated for a bound edge
!! @param iout  unit of the text output
 use confine_summary_m, only: confine_summary_t

!Input variables
 integer, intent(in) :: ll,ikap,iedge,nn,nb
 logical, intent(in) :: bound,table
 real(dp), intent(in) :: ecut(nb),err(nb)
 integer, intent(in) :: iout

!Input/output variables
 type(confine_summary_t), intent(inout) :: sm

!Local variables
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

 subroutine edges(sol,iq,poles,edir,ll,egrid,sh,shift,ed,iout,ifail)
!! the edges of a scan at the iq-th radius of sol: the energies at which
!! the phase (plus shift) crosses atan(D) + k pi for the D of each edge,
!! each to round-off by Brent's method in the interval of the grid that
!! brackets it; poles for the phase that counts the poles edir rather than
!! the nodes, whose top edges (phi(rs) = 0) are the poles themselves, where
!! the count steps: the k-th at the phase -pi/2 - (k-1) pi just below it
!!
!! @param sol  the phase solver
!! @param iq  index of the radius among those of sol
!! @param poles  whether the phase counts the poles edir, not the nodes
!! @param edir  poles of r phi'/phi at the radius, ascending (Ha)
!! @param ll  angular momentum
!! @param egrid  energies of the scan, descending (Ha)
!! @param sh  phases on egrid, shift included
!! @param shift  multiple of pi added to the phase of sol
!! @param ed  the edges (output)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use phsft_m, only: phsft_t
 use root_m, only: root_brent

!Input/output variables
 type(phsft_t), intent(inout) :: sol

!Input variables
 integer, intent(in) :: iq,ll
 logical, intent(in) :: poles
 real(dp), intent(in) :: edir(:),egrid(:),sh(:),shift
!> unit of the text output
 integer, intent(in) :: iout

!Output variables
 type(edges_t), intent(out) :: ed
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: thlo,thhi,er,e0,tgt
 integer :: jj,kk,jed,ierr2,jfail,npsh
 ifail=0
 npsh=size(egrid)
 ed%ne=0
 do jj=1,npsh-1
   e0=egrid(jj)
   thlo=min(sh(jj),sh(jj+1))
   thhi=max(sh(jj),sh(jj+1))
   do jed=1,3
     if(poles .and. jed==3) cycle
     do kk=ceiling((thlo-target(ll,jed))/pi),floor((thhi-target(ll,jed))/pi)
       tgt=target(ll,jed)+pi*dble(kk)
! a crossing on the point below is the next interval's
       if(tgt==sh(jj+1) .and. jj<npsh-1) cycle
       call root_brent(phase,egrid(jj+1),e0,er,ierr2,jfail)
       if(jfail/=0) ifail=jfail
       if(ifail/=0) return
       if(ierr2/=0) cycle
! a jump (a pair of nodes inside the projectors), not a root
       if(abs(phase(er))>0.25d0*pi) cycle
       if(ifail/=0) return
       call add_edge(ed,er,kk,jed,iout,ifail)
       if(ifail/=0) return
     end do
   end do
 end do
 if(poles) then
   do kk=1,size(edir)
     if(edir(kk)<egrid(npsh) .or. edir(kk)>egrid(1)) cycle
     call add_edge(ed,edir(kk),nint(shift/pi)-(kk-1),3,iout,ifail)
     if(ifail/=0) return
   end do
 end if

 contains

 real(dp) function phase(ex)
!! the phase of sol at ex plus shift, less tgt; 0 after a failure (ifail
!! of the host set), which stops root_brent at once
 real(dp), intent(in) :: ex
 real(dp) :: th
 phase=0.0d0
 if(ifail/=0) return
 call phase_at(sol,iq,poles,edir,ex,th,iout,ifail)
 if(ifail/=0) return
 phase=th+shift-tgt
 end function phase

 end subroutine edges

 subroutine add_edge(ed,ee,kk,jed,iout,ifail)
!! an edge more
!!
!! @param ed  the edges; ee, kk and jed added
!! @param ee  energy of the edge (Ha)
!! @param kk  multiple of pi of the edge
!! @param jed  the edge: 1 bottom, 2 centre, 3 top
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_input

!Input/output variables
 type(edges_t), intent(inout) :: ed

!Input variables
 real(dp), intent(in) :: ee
 integer, intent(in) :: kk,jed
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ifail

 ifail=0
 if(ed%ne==mxedge) then
   write(iout,'(/a,i0)') 'add_edge: ERROR too many edges in the scan, ', &
&        mxedge
   ifail=fail_input
   return
 end if
 ed%ne=ed%ne+1
 ed%ee(ed%ne)=ee
 ed%iband(ed%ne)=kk
 ed%ied(ed%ne)=jed
 end subroutine add_edge

 subroutine scan_phases(sae,sps,polq,ex,thae,thps,iout,ifail)
!! the phases of every radius at ex from one solution of each atom: the
!! all-electron one's, and the pseudopotential's less pi per pole of each
!! sphere below ex (phase_of)
!!
!! @param sae  the all-electron phase solver
!! @param sps  the pseudopotential's phase solver
!! @param polq  the poles of each sphere (Ha)
!! @param ex  energy (Ha)
!! @param thae  all-electron phase at each radius (output)
!! @param thps  pseudopotential phase at each radius (output)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use phsft_m, only: phsft_t, phsft_node_phases

!Input/output variables
 type(phsft_t), intent(inout) :: sae,sps

!Input variables
 type(poles_t), intent(in) :: polq(:)
 real(dp), intent(in) :: ex
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: thae(:),thps(:)
 integer, intent(out) :: ifail

!Local variables
 integer :: jq
 call phsft_node_phases(sae,ex,thae,iout=iout,ifail=ifail)
 if(ifail/=0) return
 call phsft_node_phases(sps,ex,thps,nodes=.false.,iout=iout,ifail=ifail)
 if(ifail/=0) return
 do jq=1,size(polq)
   thps(jq)=thps(jq)-pi*dble(count(polq(jq)%ee<ex))
 end do
 end subroutine scan_phases

 subroutine phase_at(sol,iq,poles,edir,ex,th,iout,ifail)
!! the phase of sol at ex at its iq-th radius: less pi per node inside rs,
!! or with poles per pole of r phi'/phi at rs (edir) below ex
!!
!! @param sol  the phase solver
!! @param iq  index of the radius among those of sol
!! @param poles  whether to count the poles edir, not the nodes
!! @param edir  poles of r phi'/phi at the radius (Ha)
!! @param ex  energy (Ha)
!! @param th  the phase
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use phsft_m, only: phsft_t, phsft_node_phase

!Input/output variables
 type(phsft_t), intent(inout) :: sol

!Input variables
 integer, intent(in) :: iq
 logical, intent(in) :: poles
 real(dp), intent(in) :: edir(:),ex
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: th
 integer, intent(out) :: ifail

 if(poles) then
   call phsft_node_phase(sol,iq,ex,th,nodes=.false.,iout=iout,ifail=ifail)
   th=th-pi*dble(count(edir<ex))
 else
   call phsft_node_phase(sol,iq,ex,th,iout=iout,ifail=ifail)
 end if
 end subroutine phase_at

 subroutine refine_poles(sol,iq,egrid,edir,iout,ifail)
!! the poles edir of the iq-th radius of sol in the scan, each from its
!! variational value to the root of u(rs) beside it, in a bracket that
!! grows from 1e-6 Ha on both sides until u changes sign, up to halfway to
!! the poles beside it (or 1 Ha). A variational bound lies above the root
!! of the same Hamiltonian, but the basis's quadrature of the pieces and
!! the radial equation's differ, by 3e-8 Ha at 1.5 Ha in Si f, so the
!! root can be on either side
!!
!! @param sol  the pseudopotential's phase solver
!! @param iq  index of the radius among those of sol
!! @param egrid  energies of the scan, descending (Ha)
!! @param edir  poles of r phi'/phi at the radius, ascending (Ha); those
!!              in the scan refined on return
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use phsft_m, only: phsft_t, phsft_node_phase
 use root_m, only: root_brent

!Input/output variables
 type(phsft_t), intent(inout) :: sol
 real(dp), intent(inout) :: edir(:)

!Input variables
 integer, intent(in) :: iq
 real(dp), intent(in) :: egrid(:)
!> unit of the text output
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: e0,u0,lo,hi,ww,er,elow,ehigh
 integer :: kk,ierr2,jfail,ndir
 logical :: found
 ifail=0
 ndir=size(edir)
 do kk=1,ndir
   e0=edir(kk)
   if(e0<egrid(size(egrid)) .or. e0>egrid(1)) cycle
   elow=e0-1.0d0
   if(kk>1) elow=max(elow,0.5d0*(edir(kk-1)+e0))
   ehigh=e0+1.0d0
   if(kk<ndir) ehigh=min(ehigh,0.5d0*(e0+edir(kk+1)))
   u0=u_rs(e0)
   if(ifail/=0) return
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
   if(ifail/=0) return
   if(.not. found) cycle
   call root_brent(u_rs,lo,hi,er,ierr2,jfail)
   if(jfail/=0) ifail=jfail
   if(ifail/=0) return
   if(ierr2==0) edir(kk)=er
 end do

 contains

 real(dp) function u_rs(ex)
!! u of sol at rs and the energy ex; 0 after a failure (ifail of the host
!! set), which stops root_brent at once
 real(dp), intent(in) :: ex
 real(dp) :: th
 u_rs=0.0d0
 if(ifail/=0) return
 call phsft_node_phase(sol,iq,ex,th,nodes=.false.,uq_out=u_rs, &
&      iout=iout,ifail=ifail)
 if(ifail/=0) u_rs=0.0d0
 end function u_rs

 end subroutine refine_poles

 integer function majority(iv)
!! the most frequent value, the smallest of a tie
!!
!! @param iv  the values

!Input variables
 integer, intent(in) :: iv(:)

!Local variables
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
!!
!! @param ll  angular momentum
!! @param jed  the edge: 1 bottom, 2 centre, 3 top
 use ws_basis_m, only: ws_edge_d

!Input variables
 integer, intent(in) :: ll,jed

 if(jed==3) then
   target=-0.5d0*pi
 else
   target=atan(ws_edge_d(ll,jed))
 end if
 end function target

 subroutine sort_edges(ed)
!! ascending in energy, by insertion
!!
!! @param ed  the edges, sorted on return

!Input/output variables
 type(edges_t), intent(inout) :: ed

!Local variables
 real(dp) :: et
 integer :: jj,kk,kt,it
 do jj=2,ed%ne
   et=ed%ee(jj)
   kt=ed%iband(jj)
   it=ed%ied(jj)
   kk=jj-1
   do while(kk>=1)
     if(ed%ee(kk)<=et) exit
     ed%ee(kk+1)=ed%ee(kk)
     ed%iband(kk+1)=ed%iband(kk)
     ed%ied(kk+1)=ed%ied(kk)
     kk=kk-1
   end do
   ed%ee(kk+1)=et
   ed%iband(kk+1)=kt
   ed%ied(kk+1)=it
 end do
 end subroutine sort_edges

 subroutine n_of(irel,c_light,ll,kap,ex,rr,zz,vfull,mmax,adams,rs,nn,ifail)
!! the all-electron state's nodes inside rs, plus l+1
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2; 0 unless irel=2
!! @param ex  energy of the state (Ha)
!! @param rr  log radial mesh (bohr)
!! @param zz  atomic number
!! @param vfull  all-electron potential (Ha)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param rs  radius of the sphere (bohr)
!! @param nn  the nodes inside rs, plus l+1
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use scatterstate_m, only: scatterstate
 use phsft_m, only: phsft_mch
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: irel,ll,kap,mmax
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: ex,rr(mmax),zz,vfull(mmax),rs

!Output variables
 integer, intent(out) :: nn,ifail

!Local variables
 real(dp), allocatable :: uu(:,:),up(:,:)
 integer :: mch,jj,nnode,nn2
 mch=phsft_mch(rr,mmax,rs)
 allocate(uu(mmax,2),up(mmax,2))
 uu(:,:)=0.0d0
 up(:,:)=0.0d0
 nn2=0
 nn=0
 call scatterstate(irel,c_light,nn2,ll,kap,ex,rr,zz,vfull,uu,up,mmax,adams, &
&                  mch,ifail)
 if(ifail/=0) return
 nnode=0
 do jj=2,mmax
   if(rr(jj)>rs) exit
   if(uu(jj,1)*uu(jj-1,1)<0.0d0) nnode=nnode+1
 end do
 nn=nnode+ll+1
 deallocate(uu,up)
 end subroutine n_of

 end module confine_ws_m

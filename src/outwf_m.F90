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
!> Outward integration of the radial equation
!>   u'' = [l(l+1)/r**2 + 2(V - E)] u + 2 beta_j
!> for the homogeneous solution and one inhomogeneous solution per source
!> beta_j (the projectors of vkboutwf), where V and the beta_j are functions
!> of pieces (pieces_m) with kinks at the same few radii (doc/
!> smoothness_plan.md, R-10). A multistep method whose stencil spans a
!> kink has an error that depends on where in the cell the kink falls, so
!> no integration on the mesh nodes is smooth in the radii. Here each
!> segment between kinks is integrated on its own grid, uniform in ln(r)
!> with the mesh step, that starts exactly at the segment's lower bound
!> and runs og_ext points past its upper one (the pieces hold their smooth
!> extensions there): the first segment on the nodes themselves from the
!> power series at the first node, each later one from (u, u') at its
!> lower bound interpolated from the grid before, every one primed by
!> Runge-Kutta sub-steps and continued by the implicit Adams-Moulton steps
!> of the other integrators (adams_m), with V and the sources interpolated
!> onto the grid from their pieces. The overlaps
!> with the sources are taken per segment on its grid by the 8th-order
!> cumulative rule (ugrid_m) to the segment's bound, and the solution is
!> interpolated back onto the mesh nodes. Every step then has the same
!> form wherever the kinks sit between the nodes, and the result is a
!> smooth function of their radii, to the method's truncation.

 module outwf_m
 use constants_m, only: dp
 use pieces_m, only: piecewise_t, pw_maxseg, pw_segment, pw_node_values
 use adams_m, only: adams_t
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: outwf_t, og_ext, outwf_setup, outwf_integrate, outwf_overlap
 public :: outwf_nodes, outwf_count, outwf_index
 public :: outwf_index_at
 public :: outwf_solve, outwf_free

!> the most segments: the breakpoints of the sources and of the potential
!> are the same radii (0.95 rc, rc, 0.95 rcloc, rcloc), pw_maxseg at most
 integer, parameter :: og_maxseg=pw_maxseg+1
!> grid points past a segment's end: the reach of the interpolant that
!> restarts the next segment there and of the cumulative rule's
!> interpolation at the bound (both lm_fh_reach=27), and the cell rule's
!> stencil
 integer, parameter :: og_ext=32
!> Runge-Kutta sub-steps of each primed step
 integer, parameter :: nsub=4

 type og_seg_t
!> grid points (0..npt-1)
   integer :: npt=0
!> the last grid point integrated (-1 for none)
   integer :: jmax=-1
!> the piece of the potential
   integer :: ipv=0
!> the piece of the sources (0 beyond their support)
   integer :: ipb=0
!> mesh position (units of the step, 0 at rr(1)) of grid point 0
   real(dp) :: t0=0.0d0
!> the position of the segment's end on the grid
   real(dp) :: tend=0.0d0
!> radius at the grid points (bohr)
   real(dp), allocatable :: rg(:)
!> al**2 (l(l+1) + 2 V r**2) at the grid points
   real(dp), allocatable :: cv(:)
!> 2 al**2 r**2 at the grid points, the coefficient of -E
   real(dp), allocatable :: r2(:)
!> 2 al**2 beta r**2 at the grid points, one column per source
   real(dp), allocatable :: src(:,:)
!> beta at the grid points, one column per source
   real(dp), allocatable :: beta(:,:)
!> cv at the sub-step points of the primed steps (from the first grid
!> point)
   real(dp), allocatable :: cvs(:)
!> r2 at the sub-step points of the primed steps
   real(dp), allocatable :: r2s(:)
!> src at the sub-step points of the primed steps
   real(dp), allocatable :: srcs(:,:)
!> the solutions at the grid points, homogeneous first, then one per source
   real(dp), allocatable :: uu(:,:)
!> their first derivatives in the position
   real(dp), allocatable :: up(:,:)
!> their second derivatives in the position
   real(dp), allocatable :: upp(:,:)
 end type og_seg_t

 type outwf_t
   integer :: ll=0  !< angular momentum
   integer :: nsrc=0  !< number of sources
   integer :: nseg=0  !< number of segments
   integer :: mmax=0  !< size of the log radial mesh
!> order and coefficients of the Adams-Moulton steps; a segment's first
!> order-1 steps are primed by Runge-Kutta
   type(adams_t) :: adams
!> the mesh step in ln(r)
   real(dp) :: al=0.0d0
!> the end of the sources' support (bohr)
   real(dp) :: rsrc=0.0d0
!> the position (as a mesh index) of the last breakpoint
   real(dp) :: tlast=0.0d0
!> the integration limit of the last outwf_integrate (mesh position)
   real(dp) :: plim=0.0d0
!> the merged bounds, bound(0)=0, bound(nseg)=rr(mmax)
   real(dp) :: bound(0:og_maxseg)=0.0d0
   real(dp), allocatable :: rr(:)  !< log radial mesh (bohr)
!> the weights of the first segment's overlaps (lm_integral_weights), nw1
!> of them
   real(dp), allocatable :: w1(:)
!> the weights times the first segment's sources, wb1(i,j) = w1(i)
!> beta_j(r_i) (outwf_overlaps with BLAS)
   real(dp), allocatable :: wb1(:,:)
!> the upper limit of the weights w1 (bohr), -1 for none
   real(dp) :: rb1=-1.0d0
   integer :: nw1=0  !< number of the weights w1
   type(og_seg_t) :: seg(og_maxseg)  !< the segments, nseg of them
 end type outwf_t

 contains

 subroutine outwf_setup(ll,rr,mmax,adams,vloc,vlocpw,nsrc,srcpw,ow,iout, &
&                       ifail)

!! the segments and their grids, with the potential and the sources on
!! them: everything that does not depend on the energy
!!
!! @param ll  angular momentum
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param vloc  the potential on the nodes; its difference from the node
!!              values of vlocpw is a smooth function added to every piece
!!              (the screening of the pseudo-atom's iteration)
!! @param vlocpw  the screened local potential as pieces, its last bound
!!                the mesh end
!! @param nsrc  number of sources
!! @param srcpw  the sources' pieces, all with the same bounds
!! @param ow  the setup
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use logmesh_m, only: lm_fh_weights, lm_fh_reach, lm_step, lm_position

!Input variables
 integer, intent(in) :: ll,mmax,nsrc
 real(dp), intent(in) :: rr(mmax),vloc(mmax)
 type(piecewise_t), intent(in) :: vlocpw,srcpw(nsrc)
 type(adams_t), intent(in) :: adams
 integer, intent(in) :: iout

!Output variables
 type(outwf_t), intent(out), target :: ow
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: vs(:),vn(:),fpiece(:)
 real(dp) :: al,als,sls,brk(2*pw_maxseg+2),tt,rm,plast
 integer :: kk,nbrk,nseg,npt,nsol,nprime
 type(og_seg_t), pointer :: sg

 ifail=0
 allocate(vs(mmax),vn(mmax),fpiece(mmax))
 ow%ll=ll
 ow%mmax=mmax
 ow%adams=adams
 nprime=adams%order-1
 ow%nsrc=nsrc
 allocate(ow%rr(mmax))
 ow%rr(:)=rr(:)
 al=lm_step(rr,mmax)
 als=al**2
 sls=dble(ll*(ll+1))
 ow%al=al
 nsol=nsrc+1
 plast=dble(mmax-1)

 call pw_node_values(rr,mmax,vlocpw,vn)
 vs(:)=vloc(:)-vn(:)

! the breakpoints: the bounds of the potential and of the sources inside
! the mesh, in order, each once, and the segments between them
 call breakpoints(ifail)
 if(ifail/=0) return

 do kk=1,nseg
   sg=>ow%seg(kk)
   if(kk<nseg) then
     rm=0.5d0*(ow%bound(kk-1)+ow%bound(kk))
   else
     rm=max(ow%bound(kk-1),rr(1))*exp(0.5d0*al)
   end if
   sg%ipv=pw_segment(vlocpw,rm)
   if(sg%ipv==0) then
     write(iout,'(/a)') &
&          'outwf_setup: ERROR the potential''s pieces end inside the mesh'
     ifail=1
     return
   end if
   sg%ipb=0
   if(nsrc>0) sg%ipb=pw_segment(srcpw(1),rm)
   fpiece(:)=vlocpw%vals(:,sg%ipv)+vs(:)

   if(kk==1) then
     call first_grid(ifail)
     if(ifail/=0) return
   else
     call later_grid(kk,ifail)
     if(ifail/=0) return
   end if
   allocate(sg%uu(0:npt-1,nsol),sg%up(0:npt-1,nsol),sg%upp(0:npt-1,nsol))
   sg%uu(:,:)=0.0d0
   sg%up(:,:)=0.0d0
   sg%upp(:,:)=0.0d0
   sg%jmax=-1
 end do
 deallocate(vs,vn,fpiece)
 return

 contains

 subroutine breakpoints(ifail)

!! the bounds of the potential and of the sources inside the mesh, in
!! order, each once (brk, nbrk), and the segments between them

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: ii,jj,kk

 ifail=0
 nbrk=0
 do kk=1,vlocpw%nseg
   call add_break(vlocpw%bound(kk))
 end do
 ow%rsrc=0.0d0
 if(nsrc>0) then
   do jj=2,nsrc
     if(srcpw(jj)%nseg/=srcpw(1)%nseg .or. &
&       any(srcpw(jj)%bound(0:srcpw(1)%nseg) &
&           /=srcpw(1)%bound(0:srcpw(1)%nseg))) then
       write(error_unit,'(/a)') &
&           'outwf_setup: ERROR the sources'' segments differ'
       error stop 1
     end if
   end do
   do kk=1,srcpw(1)%nseg
     call add_break(srcpw(1)%bound(kk))
   end do
   ow%rsrc=srcpw(1)%bound(srcpw(1)%nseg)
 end if
 do ii=1,nbrk-1
   do jj=ii+1,nbrk
     if(brk(jj)<brk(ii)) then
       tt=brk(ii)
       brk(ii)=brk(jj)
       brk(jj)=tt
     end if
   end do
 end do
 if(nbrk>og_maxseg-1) then
   write(iout,'(/a,i0)') 'outwf_setup: ERROR too many breakpoints ',nbrk
   ifail=1
   return
 end if
 nseg=nbrk+1
 ow%nseg=nseg
 ow%bound(:)=0.0d0
 ow%bound(1:nbrk)=brk(1:nbrk)
 ow%bound(nseg)=rr(mmax)
 ow%tlast=0.0d0
 if(nbrk>0) ow%tlast=1.0d0+lm_position(rr,mmax,brk(nbrk))
 return
 end subroutine breakpoints

 subroutine first_grid(ifail)

!! the first segment's grid: the nodes up to the first bound and og_ext
!! beyond (all of them for one segment)

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: jj,i1

 ifail=0
 sg%t0=0.0d0
 if(nseg>1) then
   sg%tend=lm_position(rr,mmax,ow%bound(1))
   i1=floor(sg%tend)+1
   if(i1<4) then
     write(iout,'(/a,es12.4)') &
&              'outwf_setup: ERROR breakpoint below the 4th node, r=', &
&              ow%bound(1)
     ifail=1
     return
   end if
   npt=min(mmax,ceiling(sg%tend)+og_ext+1)
   if(npt-1<ceiling(sg%tend)+og_ext) then
     write(iout,'(/a,es12.4)') &
&              'outwf_setup: ERROR breakpoint too close to the mesh end, r=', &
&              ow%bound(1)
     ifail=1
     return
   end if
 else
   sg%tend=plast
   npt=mmax
 end if
 sg%npt=npt
 allocate(sg%rg(0:npt-1),sg%cv(0:npt-1),sg%r2(0:npt-1))
 allocate(sg%src(0:npt-1,nsrc),sg%beta(0:npt-1,nsrc))
 sg%rg(0:npt-1)=rr(1:npt)
 sg%cv(0:npt-1)=als*sls+2.0d0*als*fpiece(1:npt)*rr(1:npt)**2
 sg%r2(0:npt-1)=2.0d0*als*rr(1:npt)**2
 sg%beta(:,:)=0.0d0
 sg%src(:,:)=0.0d0
 do jj=1,nsrc
   if(sg%ipb==0) exit
   sg%beta(0:npt-1,jj)=srcpw(jj)%vals(1:npt,sg%ipb)
   sg%src(0:npt-1,jj)=2.0d0*als*sg%beta(0:npt-1,jj)*rr(1:npt)**2
 end do
 call substeps(0.0d0)
 return
 end subroutine first_grid

 subroutine later_grid(kk,ifail)

!! the grid of segment kk>1: from the bound, the mesh step, to og_ext past
!! the next bound or to the mesh end

!Input variables
 integer, intent(in) :: kk

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: jj

 ifail=0
 sg%t0=lm_position(rr,mmax,ow%bound(kk-1))
 if(kk<nseg) then
   sg%tend=lm_position(rr,mmax,ow%bound(kk))-sg%t0
   npt=ceiling(sg%tend)+og_ext+1
 else
   sg%tend=plast-sg%t0
   npt=floor(sg%tend)+1
 end if
 npt=min(npt,floor(plast-sg%t0)+1)
 if(kk<nseg .and. npt-1<ceiling(sg%tend)+og_ext) then
   write(iout,'(/a,es12.4)') &
&            'outwf_setup: ERROR breakpoint too close to the mesh end, r=', &
&            ow%bound(kk)
   ifail=1
   return
 end if
 if(sg%t0+nprime>plast) then
   write(iout,'(/a,es12.4)') &
&            'outwf_setup: ERROR breakpoint too close to the mesh end, r=', &
&            ow%bound(kk-1)
   ifail=1
   return
 end if
 sg%npt=npt
 allocate(sg%rg(0:npt-1),sg%cv(0:npt-1),sg%r2(0:npt-1))
 allocate(sg%src(0:npt-1,nsrc),sg%beta(0:npt-1,nsrc))
 do jj=0,npt-1
   sg%rg(jj)=rr(1)*exp(al*(sg%t0+dble(jj)))
 end do
 sg%rg(0)=ow%bound(kk-1)
 sg%r2(0:npt-1)=2.0d0*als*sg%rg(0:npt-1)**2
 call fh_offset(fpiece,mmax,sg%t0,0,npt-1,sg%cv(0:npt-1))
 sg%cv(0:npt-1)=als*sls+sg%cv(0:npt-1)*sg%r2(0:npt-1)
 sg%beta(:,:)=0.0d0
 sg%src(:,:)=0.0d0
 do jj=1,nsrc
   if(sg%ipb==0) exit
   call fh_offset(srcpw(jj)%vals(:,sg%ipb),mmax,sg%t0,0,npt-1, &
&                     sg%beta(0:npt-1,jj))
   sg%src(0:npt-1,jj)=sg%beta(0:npt-1,jj)*sg%r2(0:npt-1)
 end do
 call substeps(sg%t0)
 return
 end subroutine later_grid

 subroutine substeps(tq)

!! the potential and the sources at the sub-step points of the primed
!! steps of the segment sg, from the mesh position tq

 real(dp), intent(in) :: tq
 real(dp) :: wts(2*lm_fh_reach()+2)
 integer :: is,js,nsp,ilo,ihi,nw

 nsp=2*nsub*nprime
 allocate(sg%cvs(0:nsp),sg%r2s(0:nsp),sg%srcs(0:nsp,nsrc))
 sg%srcs(:,:)=0.0d0
 do is=0,nsp
   tt=tq+dble(is)/dble(2*nsub)
   sg%r2s(is)=2.0d0*als*(rr(1)*exp(al*tt))**2
! the interpolant's weights at tt, shared by the potential and the sources
   call lm_fh_weights(mmax,tt,ilo,ihi,wts)
   nw=ihi-ilo+1
   sg%cvs(is)=als*sls+sum(wts(1:nw)*fpiece(ilo:ihi))*sg%r2s(is)
   do js=1,nsrc
     if(sg%ipb==0) exit
     sg%srcs(is,js)=sum(wts(1:nw)*srcpw(js)%vals(ilo:ihi,sg%ipb)) &
&                  *sg%r2s(is)
   end do
 end do
 end subroutine substeps

 subroutine add_break(bb)

!! a bound inside the mesh, if not there already

 real(dp), intent(in) :: bb
 integer :: ib

 if(bb<=0.0d0 .or. bb>=rr(mmax)) return
 do ib=1,nbrk
   if(brk(ib)==bb) return
 end do
 nbrk=nbrk+1
 brk(nbrk)=bb
 end subroutine add_break

 end subroutine outwf_setup

 subroutine outwf_free(ow)

!! releases the grids
!!
!! @param ow  the setup; its arrays deallocated

!Input/output variables
 type(outwf_t), intent(inout), target :: ow

!Local variables
 integer :: kk

 do kk=1,og_maxseg
   if(allocated(ow%seg(kk)%rg)) deallocate(ow%seg(kk)%rg,ow%seg(kk)%cv, &
&     ow%seg(kk)%r2,ow%seg(kk)%src,ow%seg(kk)%beta)
   if(allocated(ow%seg(kk)%cvs)) deallocate(ow%seg(kk)%cvs,ow%seg(kk)%r2s, &
&     ow%seg(kk)%srcs)
   if(allocated(ow%seg(kk)%uu)) deallocate(ow%seg(kk)%uu,ow%seg(kk)%up, &
&     ow%seg(kk)%upp)
   ow%seg(kk)%npt=0
   ow%seg(kk)%jmax=-1
 end do
 if(allocated(ow%rr)) deallocate(ow%rr)
 if(allocated(ow%w1)) deallocate(ow%w1)
 if(allocated(ow%wb1)) deallocate(ow%wb1)
 ow%rb1=-1.0d0
 ow%nw1=0
 ow%nseg=0
 return
 end subroutine outwf_free

 subroutine outwf_integrate(ow,ee,mlast)

!! the homogeneous solution (1) and the inhomogeneous one of each source
!! (1+j) at the energy ee, to the mesh point mlast (and og_ext points past it
!! on its segment's grid)
!!
!! @param ow  the setup
!! @param ee  energy (Ha)
!! @param mlast  index of the radius to which the solutions are wanted

 use logmesh_m, only: lm_integral_weights
 use ugrid_m, only: ug_interp

!Input variables
 real(dp), intent(in) :: ee
 integer, intent(in) :: mlast

!Input/output variables
 type(outwf_t), intent(inout), target :: ow

!Local variables
 real(dp) :: al,als,pm,rb
 integer :: kk,nsol,jmax,ll,nprime
 type(og_seg_t), pointer :: sg

 al=ow%al
 als=al**2
 ll=ow%ll
 nsol=ow%nsrc+1
 nprime=ow%adams%order-1
 pm=dble(mlast-1)
 ow%plim=pm
! the first segment's overlaps end at its bound, or at the limit, so their
! weights change only with the limit, if at all
 rb=min(ow%bound(1),ow%rr(1)*exp(al*pm))
 if(rb/=ow%rb1) then
   call lm_integral_weights(ow%rr,ow%mmax,2*ow%ll+2,rb,ow%nw1,ow%w1)
   ow%rb1=rb
#if defined WITH_BLAS_DOT
   if(allocated(ow%wb1)) deallocate(ow%wb1)
   sg=>ow%seg(1)
   if(ow%nsrc>0 .and. sg%ipb/=0) then
     m0=min(ow%nw1,sg%npt)
     allocate(ow%wb1(m0,ow%nsrc))
     do jj=1,ow%nsrc
       ow%wb1(:,jj)=ow%w1(1:m0)*sg%beta(0:m0-1,jj)
     end do
   end if
#endif
 end if

 do kk=1,ow%nseg
   sg=>ow%seg(kk)
   sg%jmax=-1
   if(sg%t0>pm) cycle
   if(kk<ow%nseg .and. ow%seg(kk+1)%t0<=pm) then
     jmax=sg%npt-1
   else
     jmax=min(sg%npt-1,ceiling(pm-sg%t0)+og_ext)
   end if
   jmax=max(jmax,nprime)
   sg%jmax=jmax

   if(kk==1) then
     call series_start()
   else
     call restart(kk)
   end if

! prime the first steps by classical Runge-Kutta in sub-steps
   call prime()

! the implicit Adams-Moulton steps, as the other integrators
   call march()
 end do
 return

 contains

 subroutine series_start()

!! the first segment's start: the power series at the first node, as
!! lschps and lschkb start: r**(l+1) for the homogeneous solution, the
!! particular solution of a source going as r**(l+1) for the others

!Local variables
 real(dp) :: akb
 integer :: jj

 sg%uu(0,1)=sg%rg(0)**(ll+1)
 sg%up(0,1)=al*dble(ll+1)*sg%uu(0,1)
 do jj=1,ow%nsrc
   akb=2.0d0*sg%beta(0,jj)/sg%rg(0)**(ll+1)/dble(6+4*ll)
   sg%uu(0,jj+1)=akb*sg%rg(0)**(ll+3)
   sg%up(0,jj+1)=al*dble(ll+3)*sg%uu(0,jj+1)
 end do
 return
 end subroutine series_start

 subroutine restart(kk)

!! segment kk>1 restarts from the grid before, at its lower bound

!Input variables
 integer, intent(in) :: kk

!Local variables
 real(dp) :: y1,y2
 integer :: isol
 type(og_seg_t), pointer :: prev

 prev=>ow%seg(kk-1)
 do isol=1,nsol
   call ug_interp(prev%uu(0:prev%jmax,isol),prev%jmax+1,prev%tend,y1)
   call ug_interp(prev%up(0:prev%jmax,isol),prev%jmax+1,prev%tend,y2)
   sg%uu(0,isol)=y1
   sg%up(0,isol)=y2
 end do
 return
 end subroutine restart

 subroutine prime()

!! the first nprime steps of every solution by classical Runge-Kutta in
!! nsub sub-steps each, and u'' at their points

!Local variables
 real(dp) :: cf,hh,y1,y2,k1(2),k2(2),k3(2),k4(2)
 integer :: ii,isol,istep,m0,mh,m1

 hh=1.0d0/dble(nsub)
 do isol=1,nsol
   y1=sg%uu(0,isol)
   y2=sg%up(0,isol)
   do istep=1,nprime*nsub
     m0=2*(istep-1)
     mh=2*istep-1
     m1=2*istep
     k1=rhs(m0,isol,y1,y2)
     k2=rhs(mh,isol,y1+0.5d0*hh*k1(1),y2+0.5d0*hh*k1(2))
     k3=rhs(mh,isol,y1+0.5d0*hh*k2(1),y2+0.5d0*hh*k2(2))
     k4=rhs(m1,isol,y1+hh*k3(1),y2+hh*k3(2))
     y1=y1+hh*(k1(1)+2.0d0*k2(1)+2.0d0*k3(1)+k4(1))/6.0d0
     y2=y2+hh*(k1(2)+2.0d0*k2(2)+2.0d0*k3(2)+k4(2))/6.0d0
     if(mod(istep,nsub)==0) then
       sg%uu(istep/nsub,isol)=y1
       sg%up(istep/nsub,isol)=y2
     end if
   end do
   do ii=0,nprime
     cf=sg%cv(ii)-ee*sg%r2(ii)
     sg%upp(ii,isol)=al*sg%up(ii,isol)+cf*sg%uu(ii,isol)
     if(isol>1) sg%upp(ii,isol)=sg%upp(ii,isol)+sg%src(ii,isol-1)
   end do
 end do
 return
 end subroutine prime

 subroutine march()

!! the implicit Adams-Moulton steps from nprime to jmax

!Local variables
 real(dp) :: cf,c0,b1,b2
 integer :: ii,isol

 c0=ow%adams%cm(0)
 do ii=nprime,jmax-1
   cf=sg%cv(ii+1)-ee*sg%r2(ii+1)
   do isol=1,nsol
     b1=sg%uu(ii,isol)+am_out(ow%adams,sg%up(0,isol),ii+1)
     b2=sg%up(ii,isol)+am_out(ow%adams,sg%upp(0,isol),ii+1)
     if(isol>1) b2=b2+c0*sg%src(ii+1,isol-1)
     call am_solve(c0,0.0d0,1.0d0,cf,al,b1,b2,sg%uu(ii+1,isol), &
&                  sg%up(ii+1,isol))
     sg%upp(ii+1,isol)=al*sg%up(ii+1,isol)+cf*sg%uu(ii+1,isol)
     if(isol>1) sg%upp(ii+1,isol)=sg%upp(ii+1,isol)+sg%src(ii+1,isol-1)
   end do
 end do
 return
 end subroutine march

 function rhs(mm,isol,u1,u2)

!! the right-hand side (u', u'') at the sub-step point mm for the solution
!! isol

 integer, intent(in) :: mm,isol
 real(dp), intent(in) :: u1,u2
 real(dp) :: rhs(2)

 rhs(1)=u2
 rhs(2)=al*u2+(sg%cvs(mm)-ee*sg%r2s(mm))*u1
 if(isol>1) rhs(2)=rhs(2)+sg%srcs(mm,isol-1)
 end function rhs

 end subroutine outwf_integrate

 subroutine outwf_overlap(ow,isol,jsrc,ss,skip1)

!! int_0^rsrc phi_isol beta_jsrc dr, per segment on its grid, cut at the
!! integration limit if that is inside the sources' support
!!
!! @param ow  the setup, integrated
!! @param isol  the solution, 1 the homogeneous one, 1+j that of source j
!! @param jsrc  the source
!! @param ss  the integral
!! @param skip1  if present and true, the first segment's part is left out
!!               (outwf_overlaps with BLAS)

 use ugrid_m, only: ug_interp, ug_cumint_down

!Input variables
 type(outwf_t), intent(in), target :: ow
 integer, intent(in) :: isol,jsrc
 logical, intent(in), optional :: skip1

!Output variables
 real(dp), intent(out) :: ss

!Local variables
 real(dp), allocatable :: ff(:),cum(:)
 real(dp) :: sk,cq,tq,s4(4),c4(4),y4(4),t4(4)
 integer :: kk,nn,jj
 type(og_seg_t), pointer :: sg

 ss=0.0d0
 do kk=1,ow%nseg
   sg=>ow%seg(kk)
   if(sg%ipb==0 .or. sg%jmax<0) cycle
   if(kk==1 .and. present(skip1)) then
     if(skip1) cycle
   end if
   if(kk==1) then
! lm_integral to rb1 by its weights (outwf_integrate), the sum compensated
! in four lanes (independent chains of additions, so four times faster);
! the fast math of Release folds the compensation away, which is kept
! (5-7 % of the run time; the integrals of logmesh_m and ugrid_m keep it)
     nn=min(ow%nw1,sg%jmax+1)
     s4(:)=0.0d0
     c4(:)=0.0d0
     do jj=1,nn-3,4
       y4(:)=ow%w1(jj:jj+3)*sg%uu(jj-1:jj+2,isol)*sg%beta(jj-1:jj+2,jsrc) &
&            - c4(:)
       t4(:)=s4(:) + y4(:)
       c4(:)=(t4(:) - s4(:)) - y4(:)
       s4(:)=t4(:)
     end do
     do jj=4*(nn/4)+1,nn
       y4(1)=ow%w1(jj)*sg%uu(jj-1,isol)*sg%beta(jj-1,jsrc) - c4(1)
       t4(1)=s4(1) + y4(1)
       c4(1)=(t4(1) - s4(1)) - y4(1)
       s4(1)=t4(1)
     end do
     sk=(s4(1) + s4(2)) + (s4(3) + s4(4))
   else
     nn=sg%jmax+1
     allocate(ff(nn),cum(nn))
     do jj=0,sg%jmax
       ff(jj+1)=ow%al*sg%rg(jj)*sg%uu(jj,isol)*sg%beta(jj,jsrc)
     end do
     call ug_cumint_down(ff,nn,1.0d0,cum)
     tq=min(sg%tend,ow%plim-sg%t0)
     call ug_interp(cum,nn,tq,cq)
     sk=cum(1)-cq
     deallocate(ff,cum)
   end if
   ss=ss+sk
 end do
 return
 end subroutine outwf_overlap

 subroutine outwf_overlaps(ow,ovl)

!! every overlap of outwf_overlap at once. Built with WITH_BLAS_DOT (the
!! CMake option ONCVPSP_BLAS_DOT), the first segment's, most of the work, by
!! one dgemm of the solutions with the weighted sources: faster, but rounded
!! as the linked BLAS sums, so the results depend on it at round-off
!!
!! @param ow  the setup, integrated
!! @param ovl  ovl(isol,jsrc), int phi_isol beta_jsrc dr

 use lapack_m, only: dgemm

!Input variables
 type(outwf_t), intent(in), target :: ow

!Output variables
 real(dp), intent(out) :: ovl(ow%nsrc+1,ow%nsrc)

!Local variables
 integer :: isol,jsrc
#if defined WITH_BLAS_DOT
 real(dp) :: ss
 integer :: nn
 type(og_seg_t), pointer :: sg

 ovl(:,:)=0.0d0
 sg=>ow%seg(1)
 if(ow%nsrc>0 .and. sg%ipb/=0 .and. sg%jmax>=0) then
   nn=min(ow%nw1,sg%jmax+1)
   call dgemm('T','N',ow%nsrc+1,ow%nsrc,nn,1.0d0,sg%uu,size(sg%uu,1), &
&             ow%wb1,size(ow%wb1,1),0.0d0,ovl,ow%nsrc+1)
 end if
 do jsrc=1,ow%nsrc
   do isol=1,ow%nsrc+1
     call outwf_overlap(ow,isol,jsrc,ss,.true.)
     ovl(isol,jsrc)=ovl(isol,jsrc)+ss
   end do
 end do
#else
 do jsrc=1,ow%nsrc
   do isol=1,ow%nsrc+1
     call outwf_overlap(ow,isol,jsrc,ovl(isol,jsrc))
   end do
 end do
#endif
 return
 end subroutine outwf_overlaps

 subroutine outwf_combine(ow,coef)

!! the homogeneous solution becomes phi_1 + sum_j coef(j) phi_(1+j)
!!
!! @param ow  the setup, integrated
!! @param coef  the coefficients of the inhomogeneous solutions

!Input/output variables
 type(outwf_t), intent(inout), target :: ow

!Input variables
 real(dp), intent(in) :: coef(ow%nsrc)

!Local variables
 integer :: kk,jj
 type(og_seg_t), pointer :: sg

 do kk=1,ow%nseg
   sg=>ow%seg(kk)
   if(sg%jmax<0) cycle
   do jj=1,ow%nsrc
     sg%uu(0:sg%jmax,1)=sg%uu(0:sg%jmax,1)+coef(jj)*sg%uu(0:sg%jmax,jj+1)
     sg%up(0:sg%jmax,1)=sg%up(0:sg%jmax,1)+coef(jj)*sg%up(0:sg%jmax,jj+1)
   end do
 end do
 return
 end subroutine outwf_combine

 subroutine outwf_nodes(ow,isol,mlast,uu,up)

!! a solution on the mesh nodes 1..mlast (0 beyond), each node from the grid
!! of its segment
!!
!! @param ow  the setup, integrated
!! @param isol  the solution
!! @param mlast  the last node wanted
!! @param uu  the solution
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr), its
!!            derivative in the mesh index, al r du/dr

!Input variables
 type(outwf_t), intent(in), target :: ow
 integer, intent(in) :: isol,mlast

!Output variables
 real(dp), intent(out) :: uu(ow%mmax),up(ow%mmax)

!Local variables
 integer :: ii,kk,i1,i2
 type(og_seg_t), pointer :: sg

 uu(:)=0.0d0
 up(:)=0.0d0
 i2=0
 do kk=1,ow%nseg
   sg=>ow%seg(kk)
   if(sg%jmax<0) exit
! the nodes of this segment: bound(kk-1) < r <= bound(kk)
   i1=i2+1
   if(kk<ow%nseg) then
     i2=floor(log(ow%bound(kk)/ow%rr(1))/ow%al)+1
     do while(i2<ow%mmax .and. ow%rr(i2+1)<=ow%bound(kk))
       i2=i2+1
     end do
     do while(i2>=1 .and. ow%rr(i2)>ow%bound(kk))
       i2=i2-1
     end do
   else
     i2=ow%mmax
   end if
   i2=min(i2,mlast)
   if(kk==1) then
     do ii=i1,min(i2,sg%jmax+1)
       uu(ii)=sg%uu(ii-1,isol)
       up(ii)=sg%up(ii-1,isol)
     end do
   else
     call fh_offset(sg%uu(0:sg%jmax,isol),sg%jmax+1,dble(i1-1)-sg%t0,0, &
&                   i2-i1,uu(i1:i2),sg%up(0:sg%jmax,isol),up(i1:i2))
   end if
   if(i2>=mlast) exit
 end do
 return
 end subroutine outwf_nodes

 integer function outwf_count(ow,isol,rn,mnode)

!! the sign changes of a solution between consecutive points of the grids
!! (the nodes of the first segment, then each grid over its own segment)
!! at radii above rn, up to the mesh point mnode
!!
!! @param ow  the setup, integrated
!! @param isol  the solution
!! @param rn  radius below which sign changes are not counted
!! @param mnode  the last mesh point counted

!Input variables
 type(outwf_t), intent(in), target :: ow
 integer, intent(in) :: isol,mnode
 real(dp), intent(in) :: rn

!Local variables
 real(dp) :: uprev,pm,tt
 integer :: kk,jj,jlast
 logical :: have
 type(og_seg_t), pointer :: sg

 outwf_count=0
 pm=dble(mnode-1)
 have=.false.
 uprev=0.0d0
 do kk=1,ow%nseg
   sg=>ow%seg(kk)
   if(sg%jmax<0) exit
   if(kk<ow%nseg) then
     tt=min(sg%tend,pm-sg%t0)
   else
     tt=pm-sg%t0
   end if
   jlast=min(sg%jmax,floor(tt))
   do jj=0,jlast
     if(have .and. sg%rg(jj)>rn .and. uprev*sg%uu(jj,isol)<0.0d0) then
       outwf_count=outwf_count+1
     end if
     uprev=sg%uu(jj,isol)
     have=.true.
   end do
   if(sg%t0+tt>=pm) exit
 end do
 return
 end function outwf_count

 subroutine outwf_index(ow,evkb,ovl,mnode,neig,iout,ifail)

!! the number of eigenvalues below the energy of the last outwf_integrate
!! of the radial equation with the projectors of coefficients evkb on
!! (0, R], u(R) = 0, R the last grid point of outwf_count at or below the
!! mesh point mnode, which must lie beyond the sources' support. A node
!! count of the solution with the projectors is not this number: the
!! projectors' terms start as r**(l+1) as well, so where their sum and the
!! homogeneous solution cancel at small r a node enters from r = 0 without
!! an eigenvalue passing, and moves out with the energy. The index changes
!! only where an eigenvalue passes. With A = H_loc - E on (0, R] and the
!! sources B (H = H_loc + B D B^T, D = diag(evkb)), the inertia of the
!! Schur complements (Haynsworth) gives
!!   n_-(H - E) = n_-(A) + n_-(-D^-1 - B^T A^-1 B) - n_-(-D^-1)
!! where n_-(A) is the node count of the homogeneous solution phi_1 on
!! (0, R] (no node enters at 0 with a local potential), and
!! A^-1 beta_j = -(phi_(1+j) - phi_(1+j)(R)/phi_1(R) phi_1), so that
!! (-D^-1 - B^T A^-1 B)_ij = -delta_ij/evkb(i) + <beta_i|phi_(1+j)>
!!   - phi_(1+j)(R)/phi_1(R) <beta_i|phi_1>
!! (symmetric to the integration's truncation, and symmetrized). Sources
!! with evkb = 0 drop out. Called before outwf_combine
!!
!! @param ow  the setup, integrated
!! @param evkb  coefficients of the projectors
!! @param ovl  the overlaps of outwf_overlaps
!! @param mnode  the last mesh point counted
!! @param neig  the number of eigenvalues below the energy
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use lapack_m, only: dsyev

!Input variables
 type(outwf_t), intent(in), target :: ow
 real(dp), intent(in) :: evkb(ow%nsrc),ovl(ow%nsrc+1,ow%nsrc)
 integer, intent(in) :: mnode
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: neig
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: mm(:,:),ev(:),fr(:),work(:)
 real(dp) :: uprev,pm,tt,rlast
 integer, allocatable :: isrc(:)
 integer :: kk,jj,jlast,kl,jl,nloc,nact,ii,info
 logical :: have
 type(og_seg_t), pointer :: sg

 ifail=0
! the node count of phi_1 to R, as outwf_count, and the point R
 nloc=0
 pm=dble(mnode-1)
 have=.false.
 uprev=0.0d0
 kl=1
 jl=0
 do kk=1,ow%nseg
   sg=>ow%seg(kk)
   if(sg%jmax<0) exit
   if(kk<ow%nseg) then
     tt=min(sg%tend,pm-sg%t0)
   else
     tt=pm-sg%t0
   end if
   jlast=min(sg%jmax,floor(tt))
   do jj=0,jlast
     if(have .and. uprev*sg%uu(jj,1)<0.0d0) nloc=nloc+1
     uprev=sg%uu(jj,1)
     have=.true.
     kl=kk
     jl=jj
   end do
   if(sg%t0+tt>=pm) exit
 end do
! a zero of phi_1 at R: R at the point before, with no node between
 if(ow%seg(kl)%uu(jl,1)==0.0d0 .and. jl>0) jl=jl-1
 rlast=ow%seg(kl)%rg(jl)
 if(rlast<ow%rsrc) then
   write(iout,'(/a,2es12.4)') 'outwf_index: ERROR R inside the sources,'// &
&       ' R, rsrc=',rlast,ow%rsrc
   ifail=1
   return
 end if

 nact=count(evkb/=0.0d0)
 neig=nloc
 if(nact==0) return
 allocate(mm(nact,nact),ev(nact),fr(ow%nsrc+1),isrc(nact),work(3*nact+8))
 isrc=pack([(jj,jj=1,ow%nsrc)],evkb/=0.0d0)
 fr(:)=ow%seg(kl)%uu(jl,:)
 do jj=1,nact
   do ii=1,nact
     mm(ii,jj)=ovl(isrc(jj)+1,isrc(ii)) &
&             -fr(isrc(jj)+1)/fr(1)*ovl(1,isrc(ii))
   end do
 end do
 mm=0.5d0*(mm+transpose(mm))
 do jj=1,nact
   mm(jj,jj)=mm(jj,jj)-1.0d0/evkb(isrc(jj))
 end do
 call dsyev('N','U',nact,mm,nact,ev,work,size(work),info)
 if(info/=0) then
   write(iout,'(/a,i4)') 'outwf_index: dsyev ERROR, stopping info =',info
   ifail=1
   return
 end if
 neig=nloc+count(ev<0.0d0)-count(evkb(isrc)>0.0d0)
 deallocate(mm,ev,fr,isrc,work)
 return
 end subroutine outwf_index

 subroutine outwf_index_at(ll,rr,mmax,adams,vloc,vlocpw,nsrc,srcpw,evkb,ee, &
&                          mnode,neig,iout,ifail)

!! outwf_index of the solution at the energy ee, on its own setup
!!
!! @param ll  angular momentum
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param vloc  the local potential on the mesh
!! @param vlocpw  the screened local potential as pieces
!! @param nsrc  number of sources
!! @param srcpw  the sources' pieces
!! @param evkb  coefficients of the projectors
!! @param ee  energy (Ha)
!! @param mnode  the last mesh point counted
!! @param neig  the number of eigenvalues below ee
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

!Input variables
 integer, intent(in) :: ll,mmax,nsrc,mnode
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: rr(mmax),vloc(mmax),evkb(nsrc),ee
 type(piecewise_t), intent(in) :: vlocpw,srcpw(nsrc)
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: neig
 integer, intent(out) :: ifail

!Local variables
 type(outwf_t) :: ow
 real(dp), allocatable :: uu(:),up(:)

 ifail=0
 allocate(uu(mmax),up(mmax))
 call outwf_setup(ll,rr,mmax,adams,vloc,vlocpw,nsrc,srcpw,ow,iout,ifail)
 if(ifail/=0) return
 call outwf_solve(ow,ee,evkb,mnode,mnode,uu,up,neig, &
&      iout=iout,ifail=ifail)
 if(ifail/=0) return
 call outwf_free(ow)
 deallocate(uu,up)
 return
 end subroutine outwf_index_at

 subroutine outwf_solve(ow,ep,evkb,mlast,mnode,uu,up,neig,iout,ifail)

!! the outward solution with the sources as Vanderbilt-Kleinman-Bylander
!! projectors of coefficients evkb: the homogeneous solution plus the
!! inhomogeneous ones with the coefficients that make the projections
!! self-consistent (vkboutwf), on the nodes to mlast, and optionally its
!! eigenvalue index. Its node count is no measure of the energy: a node can
!! enter at r = 0 without an eigenvalue passing (outwf_index)
!!
!! @param ow  the setup
!! @param ep  energy (Ha)
!! @param evkb  coefficients of the projectors
!! @param mlast  index of the radius to which the solution is wanted
!! @param mnode  the point R of the index, <= mlast
!! @param uu  the solution
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr), its
!!            derivative in the mesh index
!! @param neig  optional, the number of eigenvalues below ep on
!!              (0, r(mnode)] (outwf_index)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use lapack_m, only: dgesv

!Input/output variables
 type(outwf_t), intent(inout), target :: ow

!Input variables
 real(dp), intent(in) :: ep,evkb(ow%nsrc)
 integer, intent(in) :: mlast,mnode
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: uu(ow%mmax),up(ow%mmax)
 integer, intent(out), optional :: neig
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: gg(:,:),gg0(:),ovl(:,:)
 integer, allocatable :: ipiv(:)
 integer :: ii,jj,info,nsrc

 ifail=0
 nsrc=ow%nsrc
 call outwf_integrate(ow,ep,mlast)
 if(nsrc>0) then
   allocate(gg(nsrc,nsrc),gg0(nsrc),ipiv(nsrc),ovl(nsrc+1,nsrc))
   call outwf_overlaps(ow,ovl)
   if(present(neig)) then
     call outwf_index(ow,evkb,ovl,mnode,neig,iout,ifail)
     if(ifail/=0) return
   end if
   do jj=1,nsrc
     gg0(jj)=evkb(jj)*ovl(1,jj)
     do ii=1,nsrc
       gg(jj,ii)=-evkb(jj)*ovl(ii+1,jj)
     end do
     gg(jj,jj)=1.0d0+gg(jj,jj)
   end do
   call dgesv(nsrc,1,gg,nsrc,ipiv,gg0,nsrc,info)
   if(info/=0) then
     write(iout,'(/a,i4)') 'outwf_solve: dgesv ERROR, stopping info =',info
     ifail=1
     return
   end if
   call outwf_combine(ow,gg0)
   deallocate(gg,gg0,ipiv,ovl)
 end if
 if(nsrc==0 .and. present(neig)) then
   allocate(ovl(1,0))
   call outwf_index(ow,evkb,ovl,mnode,neig,iout,ifail)
   if(ifail/=0) return
   deallocate(ovl)
 end if
 call outwf_nodes(ow,1,mlast,uu,up)
 return
 end subroutine outwf_solve

 subroutine fh_offset(ff,nn,tt0,j1,j2,out,ff2,out2)

!! the interpolant (lm_fh_value) of ff, and of ff2 if present, at the
!! positions tt0+j, j=j1..j2: away from the ends of the node range all share
!! the weights of the offset tt0 (lm_fh_weights), applied as one dot product
!! each; near them each position has its own, shared by ff and ff2
!!
!! @param ff  values on the nodes
!! @param nn  number of nodes
!! @param tt0  the offset, 0<=tt0
!! @param j1,j2  the positions wanted
!! @param out  the values, out(1:j2-j1+1)
!! @param ff2,out2  a second function on the same nodes, and its values

 use logmesh_m, only: lm_fh_weights, lm_fh_reach

!Input variables
 integer, intent(in) :: nn,j1,j2
 real(dp), intent(in) :: ff(nn),tt0
 real(dp), intent(in), optional :: ff2(nn)

!Output variables
 real(dp), intent(out) :: out(j2-j1+1)
 real(dp), intent(out), optional :: out2(j2-j1+1)

!Local variables
 real(dp) :: wts(2*lm_fh_reach()+2),we(2*lm_fh_reach()+2),frac,tt
 integer :: jj,ilo,ihi,nreach,jref,i0,nw,ie,je,ioff

 nreach=lm_fh_reach()
 frac=tt0-floor(tt0)
 jref=nreach+1
 nw=0
 if(nn>=2*nreach+3) then
   call lm_fh_weights(nn,frac+dble(jref),ilo,ihi,wts)
   nw=ihi-ilo+1
 end if
 do jj=j1,j2
   tt=tt0+dble(jj)
   i0=floor(tt)
   if(nw>0 .and. i0-nreach>=0 .and. i0+nreach+1<=nn-1) then
     ioff=ilo-jref+i0
     out(jj-j1+1)=dot(wts,ff(ioff:ioff+nw-1),nw)
     if(present(ff2)) out2(jj-j1+1)=dot(wts,ff2(ioff:ioff+nw-1),nw)
   else
     call lm_fh_weights(nn,tt,ie,je,we)
     out(jj-j1+1)=dot(we,ff(ie:je),je-ie+1)
     if(present(ff2)) out2(jj-j1+1)=dot(we,ff2(ie:je),je-ie+1)
   end if
 end do
 return

 contains

 real(dp) function dot(aa,bb,nq)

!! sum_i aa(i) bb(i), i=1..nq, in order

 integer, intent(in) :: nq
 real(dp), intent(in) :: aa(nq),bb(nq)
 integer :: iq

 dot=0.0d0
 do iq=1,nq
   dot=dot+aa(iq)*bb(iq)
 end do
 end function dot

 end subroutine fh_offset

 include 'adams.inc'

 end module outwf_m

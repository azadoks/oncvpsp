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
!> The plane-wave cutoff of a pseudopotential state in a Wigner-Seitz
!> sphere (oncvpsp_confine.x, doc/smoothness_plan.md, diagnostics item 7).
!> The state of l in the sphere of radius rs with r phi'/phi = D at rs (D =
!> 0 at the bottom of a band, -l-1 at its centre and infinite at its top,
!> the Wigner-Seitz rule) is expanded in the functions phi_n = r j_l(q_n
!> r), whose q_n are those of the same condition: they are the eigenstates
!> of the free radial Hamiltonian in the sphere, with eigenvalues q_n**2/2,
!> orthogonal, so the energy of the pseudo-Hamiltonian in the first N of
!> them is a variational bound (Rayleigh-Ritz) that falls to the state's
!> energy as N grows, and q_N**2/2 is the cutoff of the N functions. The
!> boundary needs no tail beyond the sphere: the basis and the state meet
!> it alike. The error has a q**-5 tail from the potential's slope at rs
!> at the bottom and centre (integrating <phi_n|V psi> by parts leaves
!> V'(rs) psi(rs) phi_n(rs)), none at the top (psi(rs) = 0); in Si it is
!> below the tails of the pseudopotential's own kinks, as slow at the top.
!> The integrals are by Gauss-Legendre panels between the kinks of the
!> potential's and the projectors' pieces, each panel shorter than a
!> radian of the largest q, so they are exact to round-off.
!> All the eigenvalues in the full basis (ws_eigenvalues) at the top's
!> condition, phi(rs) = 0, count the poles of the pseudopotential's
!> r phi'/phi at rs below an energy, the count that makes its phase
!> continuous where its nodes do not (oncvpsp_confine.x).

 module ws_basis_m
 use constants_m, only: dp, pi
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: ws_edge_d, ws_convergence, ws_eigenvalues, ws_grid_init
 public :: ws_grid_t
 public :: ws_no_basis

!> the quadrature of a channel in the sphere (ws_grid_init): its points,
!> weights, and the local potential and the projectors on them
 type :: ws_grid_t
   real(dp) :: rs=0.0d0  !< radius of the sphere (bohr)
   real(dp) :: qmax=0.0d0  !< the largest q of the basis (bohr**-1)
   integer :: ng=0  !< number of quadrature points
   integer :: nvkb=0  !< number of projectors
   real(dp), allocatable :: rg(:)  !< the quadrature points (bohr)
   real(dp), allocatable :: wg(:)  !< their weights
   real(dp), allocatable :: vg(:)  !< the local potential at them (Ha)
!> the projectors at them, vkbg(nvkb,ng)
   real(dp), allocatable :: vkbg(:,:)
 end type ws_grid_t

!> ierr of ws_convergence, ws_eigenvalues: no basis function below qmax
 integer, parameter :: ws_no_basis=1
!> Gauss-Legendre nodes per panel
 integer, parameter :: ngl=16
!> the step of the search for the roots of the basis's condition (in q rs)
 real(dp), parameter :: dxroot=0.05d0

 contains

 real(dp) function ws_edge_d(ll,iedge)

!! the logarithmic derivative r phi'/phi at the sphere of each edge: 1
!! bottom (0), 2 centre (-l-1), 3 top (phi = 0: huge)
!!
!! @param ll  angular momentum
!! @param iedge  1, 2 or 3

!Input variables
 integer, intent(in) :: ll,iedge

 select case(iedge)
 case(1)
   ws_edge_d=0.0d0
 case(2)
   ws_edge_d=-dble(ll+1)
 case default
   ws_edge_d=huge(1.0d0)
 end select
 return
 end function ws_edge_d

 subroutine ws_roots(ll,iedge,rs,qmax,nq,qq,ifail)

!! the q of the basis functions r j_l(q r) of an edge's condition at rs,
!! from 0 to qmax: x j_l'(x) = 0 (bottom), j_(l-1)(x) = 0 (centre, cos x
!! for l=0) or j_l(x) = 0 (top), x = q rs; with x = 0 for the bottom of l=0
!! (a constant). Each is a pole-free function of x, so every sign change on
!! the grid brackets one root, which Brent's method takes to round-off
!!
!! @param ll  angular momentum, 0..3
!! @param iedge  1 bottom, 2 centre, 3 top
!! @param rs  radius of the sphere (bohr)
!! @param qmax  the largest q (bohr**-1)
!! @param nq  number of roots (output)
!! @param qq  the roots (output)
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use sbf_m, only: sbf_trig
 use root_m, only: root_brent
 use failure_m, only: fail_internal, failed, msg_ws_roots_brent

!Input variables
 integer, intent(in) :: ll,iedge
 real(dp), intent(in) :: rs,qmax

!Output variables
 integer, intent(out) :: nq
 real(dp), allocatable, intent(out) :: qq(:)
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: xmax,xa,xb,fa,fb,xr
 real(dp), allocatable :: work(:)
 integer :: ii,nx,ierr

 ifail=0
 xmax=qmax*rs
 nx=int(xmax/dxroot)+1
 allocate(work(nx+1))
 nq=0
 if(iedge==1 .and. ll==0) then
   nq=1
   work(1)=0.0d0
 end if
 xa=0.5d0*dxroot
 fa=gg(xa)
 do ii=1,nx
   xb=xa+dxroot
   fb=gg(xb)
   if(sign(1.0d0,fa)/=sign(1.0d0,fb)) then
     call root_brent(gg,xa,xb,xr,ierr,ifail)
     if(ifail/=0) return
     if(ierr/=0) then
       write(error_unit,'(/a,i0)') &
&           'ws_roots: ERROR root_brent failed, ierr=',ierr
       ifail=failed(fail_internal,msg_ws_roots_brent,ierr)
       return
     end if
     if(xr<=xmax) then
       nq=nq+1
       work(nq)=xr
     end if
   end if
   xa=xb
   fa=fb
 end do
 qq=work(1:nq)/rs
 deallocate(work)
 return

 contains

 real(dp) function gg(xx)
!! the edge's condition, pole-free
 real(dp), intent(in) :: xx
 real(dp) :: ss,cc
 ss=sin(xx)
 cc=cos(xx)
 select case(iedge)
 case(1)
! x j_l' = x j_(l-1) - (l+1) j_l, and -x j_1 for l=0
   if(ll==0) then
     gg=sbf_trig(1,xx,ss,cc)
   else
     gg=xx*sbf_trig(ll-1,xx,ss,cc)-dble(ll+1)*sbf_trig(ll,xx,ss,cc)
   end if
 case(2)
   if(ll==0) then
     gg=cc
   else
     gg=sbf_trig(ll-1,xx,ss,cc)
   end if
 case default
   gg=sbf_trig(ll,xx,ss,cc)
 end select
 end function gg

 end subroutine ws_roots

 subroutine ws_grid_init(rs,qmax,rr,mmax,vloc,vlocpw,nvkb,vkbpw,grid,ifail)

!! the quadrature of a channel in the sphere and the pseudopotential on it,
!! shared by the Hamiltonians of every edge's condition at rs: Gauss-
!! Legendre panels between the kinks of the pieces inside the sphere, none
!! longer than a radian of the largest q (or 0.1 bohr)
!!
!! @param rs  radius of the sphere (bohr), beyond the projectors and the
!!            local potential's core radius
!! @param qmax  the largest q of the basis (bohr**-1)
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param vloc  local potential on the mesh
!! @param vlocpw  the screened local potential as pieces
!! @param nvkb  number of projectors
!! @param vkbpw  the projectors as pieces
!! @param grid  the quadrature and the potentials on it (output)
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use logmesh_m, only: lm_interp
 use pieces_m, only: piecewise_t, pw_value, pw_node_values
 use quad_m, only: quad_gauss_legendre

!Input variables
 integer, intent(in) :: mmax,nvkb
 real(dp), intent(in) :: rs,qmax,rr(mmax),vloc(mmax)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(nvkb)

!Output variables
 type(ws_grid_t), intent(out) :: grid
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: dv(:),bounds(:)
 real(dp) :: xg(ngl),wxg(ngl),dvq,ra,rb,hmax
 integer :: ng,nbound,ii,jj,kk,np,ip

 grid%rs=rs
 grid%qmax=qmax
 grid%nvkb=nvkb

 allocate(bounds(2+vlocpw%nseg+sum(vkbpw(:)%nseg)))
 ifail=0
 nbound=1
 bounds(1)=rr(1)
 call add_bounds(vlocpw)
 do jj=1,nvkb
   call add_bounds(vkbpw(jj))
 end do
 nbound=nbound+1
 bounds(nbound)=rs
 call sort(bounds(1:nbound))
 hmax=min(0.1d0,1.0d0/max(qmax,1.0d0))
 ng=0
 do ii=1,nbound-1
   if(bounds(ii+1)<=bounds(ii)) cycle
   ng=ng+ngl*max(1,ceiling((bounds(ii+1)-bounds(ii))/hmax))
 end do
 grid%ng=ng
 allocate(grid%rg(ng),grid%wg(ng),grid%vg(ng),grid%vkbg(nvkb,ng))
 ng=0
 do ii=1,nbound-1
   if(bounds(ii+1)<=bounds(ii)) cycle
   np=max(1,ceiling((bounds(ii+1)-bounds(ii))/hmax))
   do ip=1,np
     ra=bounds(ii)+(bounds(ii+1)-bounds(ii))*dble(ip-1)/dble(np)
     rb=bounds(ii)+(bounds(ii+1)-bounds(ii))*dble(ip)/dble(np)
     call quad_gauss_legendre(ra,rb,ngl,xg,wxg,ifail)
     if(ifail/=0) return
     grid%rg(ng+1:ng+ngl)=xg(:)
     grid%wg(ng+1:ng+ngl)=wxg(:)
     ng=ng+ngl
   end do
 end do

!the local potential: its pieces, plus the smooth difference of the mesh
!values from theirs (as outwf_m)
 allocate(dv(mmax))
 call pw_node_values(rr,mmax,vlocpw,dv)
 dv(:)=vloc(:)-dv(:)
 do kk=1,ng
   call lm_interp(rr,dv,mmax,grid%rg(kk),dvq)
   grid%vg(kk)=pw_value(rr,mmax,vlocpw,grid%rg(kk))+dvq
 end do

!the projectors
 do kk=1,ng
   do ii=1,nvkb
     grid%vkbg(ii,kk)=pw_value(rr,mmax,vkbpw(ii),grid%rg(kk))
   end do
 end do
 deallocate(dv,bounds)
 return

 contains

 subroutine add_bounds(pw)
!! the inner bounds of a piecewise function inside the sphere
 type(piecewise_t), intent(in) :: pw
 integer :: kb
 do kb=1,pw%nseg
   if(pw%bound(kb)>rr(1) .and. pw%bound(kb)<rs) then
     nbound=nbound+1
     bounds(nbound)=pw%bound(kb)
   end if
 end do
 end subroutine add_bounds

 end subroutine ws_grid_init

 subroutine ws_convergence(ll,iedge,grid,evkb,eex,nb,ecut,err,ierr,ifail)

!! the energy error of the pseudopotential state of an edge in the first N
!! basis functions, N = 1..nb, against its energy eex from the radial
!! equation: the eigenvalue of the same index in each, the index that of
!! the eigenvalue nearest eex in all nb (0 for the N below it)
!!
!! @param ll  angular momentum, 0..3
!! @param iedge  1 bottom, 2 centre, 3 top
!! @param grid  the channel's quadrature in the sphere (ws_grid_init)
!! @param evkb  coefficients of the projectors
!! @param eex  the state's energy (Ha)
!! @param nb  number of basis functions (output)
!! @param ecut  q_N**2/2 of each N (Ha, output)
!! @param err  the energy error of each N (Ha, output), 0 where N is below
!!             the state's index
!! @param ierr  0, or ws_no_basis if no basis function is below qmax
!!              (output)
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use lapack_m, only: dsyev
 use failure_m, only: fail_numerical, failed, msg_ws_convergence_dsyev

!Input variables
 integer, intent(in) :: ll,iedge
 type(ws_grid_t), intent(in) :: grid
 real(dp), intent(in) :: evkb(grid%nvkb),eex

!Output variables
 integer, intent(out) :: nb,ierr,ifail
 real(dp), allocatable, intent(out) :: ecut(:),err(:)

!Local variables
 real(dp), allocatable :: hmat(:,:),hh(:,:),ev(:),work(:)
 integer :: np,idx,lwork,info

 call ws_hamiltonian(ll,iedge,grid,evkb,nb,ecut,hmat,ierr,ifail)
 if(ifail/=0) return
 allocate(err(nb))
 err(:)=0.0d0
 if(ierr/=0) return

!the first N of them, N = nb first for the index
 allocate(hh(nb,nb),ev(nb))
 lwork=max(1,3*nb)
 allocate(work(lwork))
 idx=0
 do np=nb,1,-1
   hh(1:np,1:np)=hmat(1:np,1:np)
   call dsyev('N','U',np,hh,nb,ev,work,lwork,info)
   if(info/=0) then
     write(error_unit,'(/a,i0)') 'ws_convergence: ERROR dsyev info ',info
     ifail=failed(fail_numerical,msg_ws_convergence_dsyev)
     return
   end if
   if(np==nb) idx=minloc(abs(ev(1:nb)-eex),1)
   if(np>=idx) err(np)=ev(idx)-eex
 end do
 deallocate(hmat,hh,ev,work)
 return
 end subroutine ws_convergence

 subroutine ws_eigenvalues(ll,iedge,grid,evkb,nb,ev,ierr,ifail)

!! the eigenvalues of the pseudo-Hamiltonian in all the basis functions of
!! an edge's condition below qmax, lowest first: variational bounds on the
!! states of that condition in the sphere, each above its state's energy
!!
!! @param ll  angular momentum, 0..3
!! @param iedge  1 bottom, 2 centre, 3 top
!! @param grid  the channel's quadrature in the sphere (ws_grid_init)
!! @param evkb  coefficients of the projectors
!! @param nb  number of basis functions and eigenvalues (output)
!! @param ev  the eigenvalues (Ha, output)
!! @param ierr  0, or ws_no_basis if no basis function is below qmax
!!              (output)
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use lapack_m, only: dsyev
 use failure_m, only: fail_numerical, failed, msg_ws_eigenvalues_dsyev

!Input variables
 integer, intent(in) :: ll,iedge
 type(ws_grid_t), intent(in) :: grid
 real(dp), intent(in) :: evkb(grid%nvkb)

!Output variables
 integer, intent(out) :: nb,ierr,ifail
 real(dp), allocatable, intent(out) :: ev(:)

!Local variables
 real(dp), allocatable :: ecut(:),hmat(:,:),work(:)
 integer :: lwork,info

 call ws_hamiltonian(ll,iedge,grid,evkb,nb,ecut,hmat,ierr,ifail)
 if(ifail/=0) return
 allocate(ev(nb))
 if(ierr/=0) return
 lwork=max(1,3*nb)
 allocate(work(lwork))
 call dsyev('N','U',nb,hmat,nb,ev,work,lwork,info)
 if(info/=0) then
   write(error_unit,'(/a,i0)') 'ws_eigenvalues: ERROR dsyev info ',info
   ifail=failed(fail_numerical,msg_ws_eigenvalues_dsyev)
   return
 end if
 deallocate(ecut,hmat,work)
 return
 end subroutine ws_eigenvalues

 subroutine ws_hamiltonian(ll,iedge,grid,evkb,nb,ecut,hmat,ierr,ifail)

!! the pseudo-Hamiltonian in the nb basis functions of an edge's condition
!! below qmax, and their cutoffs q_n**2/2 (the arguments as
!! ws_convergence's)
!!
!! @param ll  angular momentum, 0..3
!! @param iedge  1 bottom, 2 centre, 3 top
!! @param grid  the channel's quadrature in the sphere (ws_grid_init)
!! @param evkb  coefficients of the projectors (Ha)
!! @param nb  number of basis functions (output)
!! @param ecut  q_n**2/2 of each basis function (Ha, output)
!! @param hmat  the pseudo-Hamiltonian in the basis (Ha, output)
!! @param ierr  0, or ws_no_basis if no basis function is below qmax
!!              (output)
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use sbf_m, only: sbf_trig

!Input variables
 integer, intent(in) :: ll,iedge
 type(ws_grid_t), intent(in) :: grid
 real(dp), intent(in) :: evkb(grid%nvkb)

!Output variables
 integer, intent(out) :: nb,ierr,ifail
 real(dp), allocatable, intent(out) :: ecut(:),hmat(:,:)

!Local variables
 real(dp), allocatable :: qq(:),phi(:,:),proj(:,:),vphi(:,:)
 real(dp) :: ss,cc,xx
 integer :: ii,jj,kk

 ierr=0
 call ws_roots(ll,iedge,grid%rs,grid%qmax,nb,qq,ifail)
 if(ifail/=0) return
 if(nb==0) then
   ierr=ws_no_basis
   allocate(ecut(0),hmat(0,0))
   return
 end if
 allocate(ecut(nb))
 ecut(:)=0.5d0*qq(:)**2

 associate(ng=>grid%ng,rg=>grid%rg,wg=>grid%wg,vg=>grid%vg)
!the normalized basis functions and their projections
   allocate(phi(ng,nb),proj(grid%nvkb,nb))
   do jj=1,nb
     do kk=1,ng
       xx=qq(jj)*rg(kk)
       ss=sin(xx)
       cc=cos(xx)
       phi(kk,jj)=rg(kk)*sbf_trig(ll,xx,ss,cc)
     end do
     phi(:,jj)=phi(:,jj)/sqrt(sum(wg(:)*phi(:,jj)**2))
   end do
   proj(:,:)=0.0d0
   do ii=1,grid%nvkb
     do kk=1,ng
       proj(ii,:)=proj(ii,:)+wg(kk)*grid%vkbg(ii,kk)*phi(kk,:)
     end do
   end do

!the Hamiltonian in all nb
   allocate(hmat(nb,nb),vphi(ng,nb))
   do ii=1,nb
     vphi(:,ii)=wg(:)*vg(:)*phi(:,ii)
   end do
   hmat(:,:)=matmul(transpose(vphi),phi)
   do jj=1,nb
     do ii=1,jj
       hmat(ii,jj)=hmat(ii,jj)+sum(evkb(:)*proj(:,ii)*proj(:,jj))
       hmat(jj,ii)=hmat(ii,jj)
     end do
     hmat(jj,jj)=hmat(jj,jj)+ecut(jj)
   end do
 end associate
 deallocate(qq,phi,proj,vphi)
 return
 end subroutine ws_hamiltonian

 subroutine sort(aa)
!! ascending, by insertion (a few values)
!!
!! @param aa  the values, sorted on return

!Input/output variables
 real(dp), intent(inout) :: aa(:)

!Local variables
 real(dp) :: tt
 integer :: ii,jj
 do ii=2,size(aa)
   tt=aa(ii)
   jj=ii-1
   do while(jj>=1)
     if(aa(jj)<=tt) exit
     aa(jj+1)=aa(jj)
     jj=jj-1
   end do
   aa(jj+1)=tt
 end do
 end subroutine sort

 end module ws_basis_m

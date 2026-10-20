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
!> Values, derivatives, integrals and roots of a function given on the log
!> radial mesh, at any radius, not only at mesh points: the operations that
!> let a radius (a core radius, a matching or crossover radius) vary
!> continuously instead of being rounded to the mesh. Every result is a
!> smooth (C-infinity, to round-off) function of the radii it is evaluated
!> at, so a quantity built from them is as smooth in its inputs as the
!> mathematics allows. The method is in doc/smoothness_plan.md: the
!> Floater-Hormann rational interpolant of blending degree lm_blend in
!> x = ln(r), which is uniform on the mesh, of the values (lm_interp), of
!> derivative fields made with fixed stencils (lm_derivs) and of the
!> cumulative integral (lm_cumint, lm_integral).
!>
!> The interpolant sums over the whole mesh. Its cardinal functions fall
!> off only as 1/|x-x_i|, so cutting the sum off near the radius would make
!> a jump of that order for rough data; the whole sum costs mmax operations
!> and is C-infinity. Its round-off is about eps*sum_i |f_i|/(1+|x-x_i|/al).
!> The mesh is rr(i)=rr(1)*amesh**(i-1), as everywhere in the programs,
!> with mmax>=11.

 module logmesh_m
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: qp => real128
 implicit none
 private
 public :: lm_blend, lm_interp, lm_derivs, lm_cumint, lm_integral, lm_root

!> blending degree d of the Floater-Hormann interpolant: it reproduces
!> polynomials of degree d in ln(r), and its error is O(al**(d+1))
 integer, parameter :: lm_blend=6

!> half width of the derivative stencils: 11 points, of 8th order or better
!> for the derivatives 1..4
 integer, parameter :: nhalf=5

!> points of the cell-integration stencils: degree 7, 8th order
 integer, parameter :: ncell=8

 contains

 subroutine lm_interp(rr,ff,mmax,rq,fq)

!! value at the radius rq of the function given on the log mesh, from its
!! Floater-Hormann interpolant in ln(r)
!!
!! @param rr  log radial mesh
!! @param ff  function on the mesh
!! @param mmax  size of the log radial mesh
!! @param rq  radius, rr(1)<=rq<=rr(mmax) (bohr)
!! @param fq  value at rq

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),ff(mmax),rq

!Output variables
 real(dp), intent(out) :: fq

!Local variables
 real(dp) :: tt

 tt=mesh_position(rr,mmax,rq)
 fq=fh_value(ff,mmax,tt,node_index(rr,mmax,rq,tt))
 return
 end subroutine lm_interp

 subroutine lm_derivs(rr,ff,mmax,rq,nder,dd)

!! value and radial derivatives up to order nder at the radius rq of a smooth
!! function given on the log mesh: the derivative fields in ln(r) on the
!! mesh, made with fixed centred stencils, converted to d/dr and interpolated
!! to rq, so that the result is smooth in rq
!!
!! @param rr  log radial mesh
!! @param ff  function on the mesh
!! @param mmax  size of the log radial mesh
!! @param rq  radius, away from the ends of the mesh by the stencils' reach
!!            (bohr)
!! @param nder  highest derivative wanted, 0<=nder<=4
!! @param dd  d^k ff/dr^k at rq in dd(k+1), k=0..nder

!Input variables
 integer, intent(in) :: mmax,nder
 real(dp), intent(in) :: rr(mmax),ff(mmax),rq

!Output variables
 real(dp), intent(out) :: dd(nder+1)

!Local variables
 real(dp), allocatable :: field(:)
 real(dp) :: cc(0:2*nhalf,0:4,0:2*nhalf),dx(0:4),hh,tt
 integer :: ii,jj,kk,inode,ishift
! Stirling numbers of the first kind: r**k d^k/dr^k = sum_j s(k,j) d^j/dx^j
 real(dp), parameter :: stirling(4,4)=reshape([ &
&   1.0d0, -1.0d0,  2.0d0, -6.0d0, &
&   0.0d0,  1.0d0, -3.0d0, 11.0d0, &
&   0.0d0,  0.0d0,  1.0d0, -6.0d0, &
&   0.0d0,  0.0d0,  0.0d0,  1.0d0],[4,4])

 hh=mesh_step(rr,mmax)
 tt=mesh_position(rr,mmax,rq)
 inode=node_index(rr,mmax,rq,tt)

! stencil weights for each position of the node in the 11-point window
 do ii=0,2*nhalf
   call stencil_weights(ii,cc(:,:,ii))
 end do

! d^k f/dx^k on the whole mesh, interpolated to rq; the window is centred
! on the node, and shifted by ishift to fit at the ends of the mesh
 dx(0)=fh_value(ff,mmax,tt,inode)
 allocate(field(mmax))
 do kk=1,nder
   do ii=1,mmax
     ishift=max(0,nhalf+1-ii) - max(0,ii+nhalf-mmax)
     field(ii)=0.0d0
     do jj=0,2*nhalf
       field(ii)=field(ii) + cc(jj,kk,nhalf-ishift)*ff(ii-nhalf+ishift+jj)
     end do
     field(ii)=field(ii)/hh**kk
   end do
   dx(kk)=fh_value(field,mmax,tt,inode)
 end do
 deallocate(field)

! to d^k f/dr^k
 dd(1)=dx(0)
 do kk=1,nder
   dd(kk+1)=0.0d0
   do jj=1,kk
     dd(kk+1)=dd(kk+1) + stirling(kk,jj)*dx(jj)
   end do
   dd(kk+1)=dd(kk+1)/rq**kk
 end do
 return
 end subroutine lm_derivs

 subroutine lm_cumint(rr,ff,mmax,mm,cum)

!! cumulative integral of a function given on the log mesh, from 0 to each
!! mesh point, the part below rr(1) from its power law r**mm
!!
!! @param rr  log radial mesh
!! @param ff  integrand on the mesh
!! @param mmax  size of the log radial mesh
!! @param mm  power of r with which ff goes to 0 at the origin
!! @param cum  int_0^rr(i) ff dr in cum(i)

!Input variables
 integer, intent(in) :: mmax,mm
 real(dp), intent(in) :: rr(mmax),ff(mmax)

!Output variables
 real(dp), intent(out) :: cum(mmax)

!Local variables
 real(dp) :: ww(ncell,0:ncell-2),hh,cell,sum_hi,sum_lo,yy,tmp
 integer :: ii,jj,istart,ishift

 hh=mesh_step(rr,mmax)
 do ishift=0,ncell-2
   call cell_weights(ishift,ww(:,ishift))
 end do

! each cell [x_i,x_i+1] integrates the degree-7 interpolant of ff*r (dr=r dx)
! on the 8 points i-3..i+4, shifted to fit at the ends; the running sum is
! compensated (Kahan), so the round-off does not grow with mmax
 sum_hi=ff(1)*rr(1)/dble(mm+1)
 sum_lo=0.0d0
 cum(1)=sum_hi
 do ii=1,mmax-1
   istart=min(max(ii-(ncell/2-1),1),mmax-ncell+1)
   ishift=ii-istart
   cell=0.0d0
   do jj=1,ncell
     cell=cell + ww(jj,ishift)*ff(istart+jj-1)*rr(istart+jj-1)
   end do
   yy=hh*cell - sum_lo
   tmp=sum_hi + yy
   sum_lo=(tmp - sum_hi) - yy
   sum_hi=tmp
   cum(ii+1)=sum_hi
 end do
 return
 end subroutine lm_cumint

 subroutine lm_integral(rr,ff,mmax,mm,ra,rb,ss)

!! integral of a function given on the log mesh between any two radii: the
!! interpolated cumulative integral at rb less that at ra, so smooth in both
!!
!! @param rr  log radial mesh
!! @param ff  integrand on the mesh
!! @param mmax  size of the log radial mesh
!! @param mm  power of r with which ff goes to 0 at the origin
!! @param ra,rb  limits, 0<=ra, rb<=rr(mmax) (bohr)
!! @param ss  int_ra^rb ff dr

!Input variables
 integer, intent(in) :: mmax,mm
 real(dp), intent(in) :: rr(mmax),ff(mmax),ra,rb

!Output variables
 real(dp), intent(out) :: ss

!Local variables
 real(dp), allocatable :: cum(:)

 allocate(cum(mmax))
 call lm_cumint(rr,ff,mmax,mm,cum)
 ss=cum_at(rb) - cum_at(ra)
 deallocate(cum)
 return

 contains

 real(dp) function cum_at(rq)

!! the cumulative integral at rq: the power law below rr(1), else the
!! interpolant of cum
!!
!! @param rq  radius (bohr)

!Input variables
 real(dp), intent(in) :: rq

!Local variables
 real(dp) :: tt

 if(rq<=0.0d0) then
   cum_at=0.0d0
 else if(rq<rr(1)) then
   cum_at=cum(1)*(rq/rr(1))**(mm+1)
 else
   tt=mesh_position(rr,mmax,rq)
   cum_at=fh_value(cum,mmax,tt,node_index(rr,mmax,rq,tt))
 end if
 return
 end function cum_at

 end subroutine lm_integral

 subroutine lm_root(rr,ff,mmax,ra,rb,rq,ierr)

!! root, between ra and rb, of the interpolant (lm_interp) of a function
!! given on the log mesh, converged to round-off by Brent's method in ln(r)
!!
!! @param rr  log radial mesh
!! @param ff  function on the mesh
!! @param mmax  size of the log radial mesh
!! @param ra,rb  bracket (bohr); the interpolant must change sign between
!!               them
!! @param rq  the root (bohr)
!! @param ierr  0, or 1 if the interpolant has the same sign at ra and rb

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),ff(mmax),ra,rb

!Output variables
 real(dp), intent(out) :: rq
 integer, intent(out) :: ierr

!Local variables
 real(dp) :: xa,xb,xc,fa,fb,fc,de,dlast,tol1,xm,pp,qq,rs,ss
 integer :: iter

 xa=log(ra)
 xb=log(rb)
 fa=value_at(xa)
 fb=value_at(xb)
 ierr=0
 rq=ra
 if(fa==0.0d0) return
 rq=rb
 if(fb==0.0d0) return
 if(sign(1.0d0,fa)==sign(1.0d0,fb)) then
   ierr=1
   return
 end if

! Brent (zeroin) with no tolerance of its own: it stops when the bracket
! is a few rounding units of ln(r) wide
 xc=xa
 fc=fa
 de=xb-xa
 dlast=de
 do iter=1,200
   if(sign(1.0d0,fb)==sign(1.0d0,fc)) then
     xc=xa
     fc=fa
     de=xb-xa
     dlast=de
   end if
   if(abs(fc)<abs(fb)) then
     xa=xb
     xb=xc
     xc=xa
     fa=fb
     fb=fc
     fc=fa
   end if
   tol1=epsilon(1.0d0)*(2.0d0*abs(xb)+1.0d-3)
   xm=0.5d0*(xc-xb)
   if(abs(xm)<=tol1 .or. fb==0.0d0) exit
   if(abs(dlast)>=tol1 .and. abs(fa)>abs(fb)) then
! inverse quadratic interpolation, or the secant
     ss=fb/fa
     if(xa==xc) then
       pp=2.0d0*xm*ss
       qq=1.0d0-ss
     else
       qq=fa/fc
       rs=fb/fc
       pp=ss*(2.0d0*xm*qq*(qq-rs)-(xb-xa)*(rs-1.0d0))
       qq=(qq-1.0d0)*(rs-1.0d0)*(ss-1.0d0)
     end if
     if(pp>0.0d0) qq=-qq
     pp=abs(pp)
     if(2.0d0*pp<min(3.0d0*xm*qq-abs(tol1*qq),abs(dlast*qq))) then
       dlast=de
       de=pp/qq
     else
       de=xm
       dlast=de
     end if
   else
! bisection
     de=xm
     dlast=de
   end if
   xa=xb
   fa=fb
   if(abs(de)>tol1) then
     xb=xb+de
   else
     xb=xb+sign(tol1,xm)
   end if
   fb=value_at(xb)
 end do
 rq=exp(xb)
 return

 contains

 real(dp) function value_at(xx)

!! the interpolant at x=ln(r)
!!
!! @param xx  ln(r)

!Input variables
 real(dp), intent(in) :: xx

!Local variables
 real(dp) :: rx,tt

 rx=exp(xx)
 tt=mesh_position(rr,mmax,rx)
 value_at=fh_value(ff,mmax,tt,node_index(rr,mmax,rx,tt))
 return
 end function value_at

 end subroutine lm_root

 real(dp) function mesh_step(rr,mmax)

!! the step al of the log mesh in ln(r), from its ends
!!
!! @param rr  log radial mesh
!! @param mmax  size of the log radial mesh

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax)

 mesh_step=log(rr(mmax)/rr(1))/dble(mmax-1)
 return
 end function mesh_step

 real(dp) function mesh_position(rr,mmax,rq)

!! position of rq on the mesh in units of the step, 0 at rr(1):
!! (ln(rq)-ln(rr(1)))/al
!!
!! @param rr  log radial mesh
!! @param mmax  size of the log radial mesh
!! @param rq  radius (bohr)

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),rq

 mesh_position=log(rq/rr(1))/mesh_step(rr,mmax)
 return
 end function mesh_position

 integer function node_index(rr,mmax,rq,tt)

!! the index of the mesh point that rq is on exactly, or 0
!!
!! @param rr  log radial mesh
!! @param mmax  size of the log radial mesh
!! @param rq  radius (bohr)
!! @param tt  its mesh position (mesh_position)

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),rq,tt

!Local variables
 integer :: ii

 node_index=0
 ii=nint(tt)+1
 if(ii>=1 .and. ii<=mmax) then
   if(rr(ii)==rq) node_index=ii
 end if
 return
 end function node_index

 real(dp) function fh_value(ff,nn,tt,inode)

!! Floater-Hormann interpolant of blending degree lm_blend through the
!! values ff on equispaced nodes, at the position tt, in barycentric form
!! (Floater and Hormann, Numer. Math. 107, 315 (2007)); node i has the weight
!! (-1)**i sum_k binomial(d,i-k) over the windows k..k+d of the nodes that
!! hold it, which is +/-2**d away from the ends
!!
!! @param ff  values on the nodes
!! @param nn  number of nodes
!! @param tt  position in units of the node spacing, 0 at the first node
!! @param inode  index of the node that tt is on, or 0

!Input variables
 integer, intent(in) :: nn,inode
 real(dp), intent(in) :: ff(nn),tt

!Local variables
 real(dp) :: ww,dt,num,den
 integer :: ii,kk,id

 if(inode>0) then
   fh_value=ff(inode)
   return
 end if
 id=min(lm_blend,nn-1)
 num=0.0d0
 den=0.0d0
 do ii=1,nn
   dt=tt-dble(ii-1)
   if(dt==0.0d0) then
     fh_value=ff(ii)
     return
   end if
   if(ii>id .and. ii<=nn-id) then
     ww=dble(2**id)
   else
     ww=0.0d0
     do kk=max(1,ii-id),min(ii,nn-id)
       ww=ww+binomial(id,ii-kk)
     end do
   end if
   if(mod(ii,2)==0) ww=-ww
   num=num+ww*ff(ii)/dt
   den=den+ww/dt
 end do
 fh_value=num/den
 return
 end function fh_value

 real(dp) function binomial(nn,kk)

!! the binomial coefficient nn over kk
!!
!! @param nn,kk  its arguments, 0<=kk<=nn

!Input variables
 integer, intent(in) :: nn,kk

!Local variables
 integer :: ii

 binomial=1.0d0
 do ii=1,kk
   binomial=binomial*dble(nn-kk+ii)/dble(ii)
 end do
 return
 end function binomial

 subroutine stencil_weights(iz,cc)

!! weights of the 11-point stencils for d^k/dx^k, k=0..4, at the node iz of
!! the window of nodes 0..2*nhalf (unit spacing), by Fornberg's algorithm
!! (Math. Comp. 51, 699 (1988)) in quadruple precision
!!
!! @param iz  the node the derivatives are at, 0..2*nhalf (nhalf: centred)
!! @param cc  weight of node j for d^k in cc(j,k)

!Input variables
 integer, intent(in) :: iz

!Output variables
 real(dp), intent(out) :: cc(0:2*nhalf,0:4)

!Local variables
 real(qp) :: wq(0:2*nhalf,0:4),xq(0:2*nhalf),zq,c1,c2,c3,c4,c5
 integer :: ii,jj,kk,mn

 do ii=0,2*nhalf
   xq(ii)=real(ii,qp)
 end do
 zq=real(iz,qp)
 wq(:,:)=0.0_qp
 c1=1.0_qp
 c4=xq(0)-zq
 wq(0,0)=1.0_qp
 do ii=1,2*nhalf
   mn=min(ii,4)
   c2=1.0_qp
   c5=c4
   c4=xq(ii)-zq
   do jj=0,ii-1
     c3=xq(ii)-xq(jj)
     c2=c2*c3
     if(jj==ii-1) then
       do kk=mn,1,-1
         wq(ii,kk)=c1*(real(kk,qp)*wq(ii-1,kk-1)-c5*wq(ii-1,kk))/c2
       end do
       wq(ii,0)=-c1*c5*wq(ii-1,0)/c2
     end if
     do kk=mn,1,-1
       wq(jj,kk)=(c4*wq(jj,kk)-real(kk,qp)*wq(jj,kk-1))/c3
     end do
     wq(jj,0)=c4*wq(jj,0)/c3
   end do
   c1=c2
 end do
 cc(:,:)=real(wq,dp)
 return
 end subroutine stencil_weights

 subroutine cell_weights(ishift,ww)

!! weights of the integral over the cell [0,1] of the degree-7 interpolant
!! through the nodes -ishift..ncell-1-ishift (unit spacing): the integrals
!! of the Lagrange basis polynomials, in quadruple precision
!!
!! @param ishift  how far the window starts left of the cell
!! @param ww  the weights

!Input variables
 integer, intent(in) :: ishift

!Output variables
 real(dp), intent(out) :: ww(ncell)

!Local variables
 real(qp) :: poly(0:ncell-1),xq(ncell),den,acc
 integer :: ii,jj,kk,ndeg

 do ii=1,ncell
   xq(ii)=real(ii-1-ishift,qp)
 end do
 do ii=1,ncell
! coefficients of prod_{j/=i} (t-x_j), and its value at x_i
   poly(:)=0.0_qp
   poly(0)=1.0_qp
   ndeg=0
   den=1.0_qp
   do jj=1,ncell
     if(jj==ii) cycle
     ndeg=ndeg+1
     do kk=ndeg,1,-1
       poly(kk)=poly(kk-1)-xq(jj)*poly(kk)
     end do
     poly(0)=-xq(jj)*poly(0)
     den=den*(xq(ii)-xq(jj))
   end do
   acc=0.0_qp
   do kk=0,ncell-1
     acc=acc+poly(kk)/real(kk+1,qp)
   end do
   ww(ii)=real(acc/den,dp)
 end do
 return
 end subroutine cell_weights

 end module logmesh_m

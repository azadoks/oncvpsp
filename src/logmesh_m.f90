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
!> The interpolant is made local without a switch. Floater-Hormann blends
!> the degree-d interpolants p_k on the windows of nodes k..k+d with
!> weights lambda_k(x) ~ |x-x_k|**-(d+1); the far p_k grow as |x-x_k|**d,
!> so its cardinal functions fall off only as 1/|x-x_i|, and anything rough
!> anywhere in the array (a kink, unset values past a solver's last point)
!> leaks everywhere at that rate. Here each window's weight is also
!> multiplied by a C-infinity bump omega of the distance from x to the
!> window's centre, 1 within nwin1 steps and 0 beyond nwin2: the result is
!> still C-infinity in x, interpolates, reproduces polynomials of degree d,
!> and has no poles (omega does not increase away from x, so the
!> alternating-sum argument of Floater and Hormann holds), and a node more
!> than about nwin2+d steps away has no influence at all. The mesh is
!> rr(i)=rr(1)*amesh**(i-1), as everywhere in the programs, with mmax>=11.

 module logmesh_m
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: lm_blend, lm_interp, lm_derivs, lm_cumint, lm_integral, lm_root
 public :: lm_integral_weights, lm_head, lm_origin
 public :: lm_ncell, lm_cell_weights, lm_fh_value, lm_fh_weights, lm_fh_reach

!> blending degree d of the Floater-Hormann interpolant: it reproduces
!> polynomials of degree d in ln(r), and its error is O(al**(d+1))
 integer, parameter :: lm_blend=6

!> half width of the derivative stencils: 11 points, of 8th order or better
!> for the derivatives 1..4
 integer, parameter :: nhalf=5

!> points of the cell-integration stencils: degree 7, 8th order
 integer, parameter :: lm_ncell=8
 integer, parameter :: ncell=lm_ncell

!> nodes of the interpolant in r at the origin (lm_origin), a factor of
!> about 2 apart
 integer, parameter :: norig=10

!> the window bump of the interpolant: 1 for a window centre within nwin1
!> steps of x, falling smoothly to 0 at nwin2 steps
 real(dp), parameter :: nwin1=10.0d0, nwin2=20.0d0

!> the weights of the 11-point stencils for d^k/dx^k, k=0..4, at each node
!> of the window 0..2*nhalf (unit spacing), stencil_table(:,k,iz), by
!> Fornberg's algorithm (Math. Comp. 51, 699 (1988)); and of the integral
!> over the cell [0,1] of the degree-7 interpolant through the nodes
!> -ishift..ncell-1-ishift, cell_table(:,ishift), the integrals of its
!> Lagrange basis polynomials: rational numbers, rounded once to double
!> precision by scripts/logmesh_weights.py
 include 'logmesh_weights.inc'

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
 real(dp) :: dx(0:4),hh,tt
 integer :: ii,jj,kk,inode,ishift,ilo,ihi,nreach
! Stirling numbers of the first kind: r**k d^k/dr^k = sum_j s(k,j) d^j/dx^j
 real(dp), parameter :: stirling(4,4)=reshape([ &
&   1.0d0, -1.0d0,  2.0d0, -6.0d0, &
&   0.0d0,  1.0d0, -3.0d0, 11.0d0, &
&   0.0d0,  0.0d0,  1.0d0, -6.0d0, &
&   0.0d0,  0.0d0,  0.0d0,  1.0d0],[4,4])

 hh=mesh_step(rr,mmax)
 tt=mesh_position(rr,mmax,rq)
 inode=node_index(rr,mmax,rq,tt)

! d^k f/dx^k on the nodes the interpolant reaches (fh_value), interpolated
! to rq; the stencil is centred on the node, and shifted by ishift to fit
! at the ends of the mesh
 dx(0)=fh_value(ff,mmax,tt,inode)
 nreach=ceiling(nwin2)+lm_blend+1
 ilo=max(1,floor(tt)+1-nreach)
 ihi=min(mmax,floor(tt)+1+nreach)
 allocate(field(mmax))
 field(:)=0.0d0
 do kk=1,nder
   do ii=ilo,ihi
     ishift=max(0,nhalf+1-ii) - max(0,ii+nhalf-mmax)
     field(ii)=0.0d0
     do jj=0,2*nhalf
       field(ii)=field(ii) &
&                + stencil_table(jj,kk,nhalf-ishift)*ff(ii-nhalf+ishift+jj)
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

 call cumint_to(rr,ff,mmax,mm,mmax,cum)
 return
 end subroutine lm_cumint

 subroutine cumint_to(rr,ff,mmax,mm,nlim,cum)

!! lm_cumint to the mesh point nlim only: each value is the one of the whole
!! mesh, since the running sum goes up from the origin and every cell keeps
!! the stencil it has on the whole mesh
!!
!! @param rr  log radial mesh
!! @param ff  integrand on the mesh
!! @param mmax  size of the log radial mesh
!! @param mm  power of r with which ff goes to 0 at the origin
!! @param nlim  last mesh point wanted, 1..mmax
!! @param cum  int_0^rr(i) ff dr in cum(i), i=1..nlim

!Input variables
 integer, intent(in) :: mmax,mm,nlim
 real(dp), intent(in) :: rr(mmax),ff(mmax)

!Output variables
 real(dp), intent(out) :: cum(nlim)

!Local variables
 real(dp) :: hh,cell,sum_hi,sum_lo,yy
! volatile: the fast math of Release (-fassociative-math) would fold the
! compensation (tmp-sum_hi)-yy to 0
 real(dp), volatile :: tmp
 integer :: ii,jj,istart,ishift

 hh=mesh_step(rr,mmax)

! each cell [x_i,x_i+1] integrates the degree-7 interpolant of ff*r (dr=r dx)
! on the 8 points i-3..i+4, shifted to fit at the ends; the running sum is
! compensated (Kahan), so the round-off does not grow with mmax
 sum_hi=ff(1)*rr(1)/dble(mm+1)
 sum_lo=0.0d0
 cum(1)=sum_hi
 do ii=1,nlim-1
   istart=min(max(ii-(ncell/2-1),1),mmax-ncell+1)
   ishift=ii-istart
   cell=0.0d0
   do jj=1,ncell
     cell=cell + cell_table(jj,ishift)*ff(istart+jj-1)*rr(istart+jj-1)
   end do
   yy=hh*cell - sum_lo
   tmp=sum_hi + yy
   sum_lo=(tmp - sum_hi) - yy
   sum_hi=tmp
   cum(ii+1)=sum_hi
 end do
 return
 end subroutine cumint_to

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
 integer :: nlim

! the cumulative integral only as far as the interpolant at max(ra,rb) reads
! it (its nodes, and the windows they lie in), so the same as on the whole
! mesh
 nlim=mmax
 if(max(ra,rb)>rr(1)) then
   nlim=min(mmax,floor(mesh_position(rr,mmax,max(ra,rb))) + 2 &
&       + lm_fh_reach() + lm_blend)
 end if
 allocate(cum(nlim))
 call cumint_to(rr,ff,mmax,mm,nlim,cum)
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
   cum_at=fh_value(cum,nlim,tt,node_index(rr,mmax,rq,tt))
 end if
 return
 end function cum_at

 end subroutine lm_integral

 subroutine lm_integral_weights(rr,mmax,mm,rb,nw,wts)

!! the weights of lm_integral from 0 to rb: the integral of any ff is
!! sum_{i=1}^{nw} wts(i) ff(i), the linear functional lm_integral applies
!! (the two agree to round-off), so that many integrals to one radius cost
!! a product each (outwf_m)
!!
!! @param rr  log radial mesh
!! @param mmax  size of the log radial mesh
!! @param mm  power of r with which the integrands go to 0 at the origin
!! @param rb  upper limit, 0<=rb<=rr(mmax) (bohr)
!! @param nw  number of weights, the last mesh point they reach
!! @param wts  the weights, wts(1:nw)

!Input variables
 integer, intent(in) :: mmax,mm
 real(dp), intent(in) :: rr(mmax),rb

!Output variables
 integer, intent(out) :: nw
 real(dp), allocatable, intent(out) :: wts(:)

!Local variables
 real(dp) :: fw(2*lm_fh_reach()+2),hh,tt,tail,total
 integer :: nlim,ilo,ihi,ii,jj,istart,ishift,inode

 if(rb<rr(1)) then
   nw=1
   allocate(wts(1))
   wts(1)=0.0d0
   if(rb>0.0d0) wts(1)=rr(1)/dble(mm+1)*(rb/rr(1))**(mm+1)
   return
 end if

! the interpolant's weights on the cumulative integral at rb, cut as in
! lm_integral
 tt=mesh_position(rr,mmax,rb)
 nlim=min(mmax,floor(tt) + 2 + lm_fh_reach() + lm_blend)
 inode=node_index(rr,mmax,rb,tt)
 if(inode>0) then
   ilo=inode
   ihi=inode
   fw(1)=1.0d0
 else
   call lm_fh_weights(nlim,tt,ilo,ihi,fw)
 end if

! cum(k) is the start rr(1) ff(1)/(mm+1) plus the cells below k, so cell c
! enters with the weights at k > c, the tail of fw beyond c
 hh=mesh_step(rr,mmax)
 nw=min(mmax,ihi+ncell/2)
 allocate(wts(nw))
 wts(:)=0.0d0
 total=sum(fw(1:ihi-ilo+1))
 wts(1)=rr(1)/dble(mm+1)*total
 do ii=1,ihi-1
   tail=total
   if(ii>=ilo) tail=sum(fw(ii+2-ilo:ihi-ilo+1))
   istart=min(max(ii-(ncell/2-1),1),mmax-ncell+1)
   ishift=ii-istart
   do jj=1,ncell
     wts(istart+jj-1)=wts(istart+jj-1) &
&                     + hh*tail*cell_table(jj,ishift)*rr(istart+jj-1)
   end do
 end do
 return
 end subroutine lm_integral_weights

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
!!              (the iteration limit is an error)

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
 if(iter>200) then
   write(error_unit,'(a,es12.4,a,es12.4)') &
&       'lm_root: ERROR not converged in 200 steps, r=',exp(xb),', f=',fb
   error stop 1
 end if
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

 subroutine lm_origin(rr,ff,mmax,f0,ierr)

!! value at r = 0 of a function smooth in r, given on the log mesh: the
!! Floater-Hormann interpolant in r (not ln(r), in which the origin is at
!! infinity) of blending degree lm_blend through norig nodes a factor of
!! about 2 apart from rr(1), at r = 0. There each window's weight is
!! about 2**-(lm_blend+1) of that of the window before, so the window
!! nearest the origin dominates and the denominator cannot vanish; the
!! error is that of a polynomial through nodes within rr(1)*2**lm_blend
!! of the origin, below round-off for the functions of the programs
!!
!! @param rr  log radial mesh
!! @param ff  function on the mesh
!! @param mmax  size of the log radial mesh
!! @param f0  value at r = 0
!! @param ierr  0, or 1 if the mesh is too short for the nodes

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),ff(mmax)

!Output variables
 real(dp), intent(out) :: f0
 integer, intent(out) :: ierr

!Local variables
 real(dp) :: xn(norig),ww,pp,num,den
 integer :: mstep,ii,jj,kk

 f0=0.0d0
 ierr=0
 mstep=max(1,nint(log(2.0d0)/mesh_step(rr,mmax)))
 if(1+mstep*(norig-1)>mmax) then
   ierr=1
   return
 end if
! the nodes in units of rr(1)
 do ii=1,norig
   xn(ii)=rr(1+mstep*(ii-1))/rr(1)
 end do
! barycentric form, node i with the weight sum_k (-1)**k / prod_{j/=i}
! (x_i-x_j) over the windows k..k+lm_blend that hold it
 num=0.0d0
 den=0.0d0
 do ii=1,norig
   ww=0.0d0
   do kk=max(1,ii-lm_blend),min(ii,norig-lm_blend)
     pp=1.0d0
     do jj=kk,kk+lm_blend
       if(jj/=ii) pp=pp/(xn(ii)-xn(jj))
     end do
     if(mod(kk,2)==0) pp=-pp
     ww=ww+pp
   end do
   ww=-ww/xn(ii)
   num=num+ww*ff(1+mstep*(ii-1))
   den=den+ww
 end do
 f0=num/den
 return
 end subroutine lm_origin

 real(dp) function lm_head(f1,rho1,rho2,al,pw)

!! the integral over x = ln(r) below rr(1) of an integrand f = rho*g that
!! goes as a power of r there, g as r**pw: f1/p, with p = pw + the power
!! of the density through its first two nodes, ln(rho2/rho1)/al. With
!! the density taken as constant (p = pw) the head is wrong where the
!! density is relativistic: an s state goes as r**gamma, gamma =
!! sqrt(1 - (Z/c)**2), and the density as r**(2 gamma - 2)
!! (doc/numerics_audit.md, C-3). The density is used rather than f
!! itself: a GGA potential is not smooth enough at the first nodes to
!! give a power. Where rho1 or rho2 is not positive, p = pw
!!
!! @param f1  integrand at rr(1)
!! @param rho1  density at rr(1)
!! @param rho2  density at rr(2)
!! @param al  step of the log mesh
!! @param pw  power of g

!Input variables
 real(dp), intent(in) :: f1,rho1,rho2,al,pw

 if(rho1>0.0d0 .and. rho2>0.0d0) then
   lm_head=f1/(pw+log(rho2/rho1)/al)
 else
   lm_head=f1/pw
 end if
 return
 end function lm_head

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

!! the windowed Floater-Hormann interpolant (see the module description) of
!! blending degree lm_blend through the values ff on equispaced nodes, at
!! the position tt, in barycentric form (Floater and Hormann, Numer. Math.
!! 107, 315 (2007)): node i has the weight
!! (-1)**i sum_k omega(|tt-c_k|) binomial(d,i-k) over the windows k..k+d
!! that hold it, c_k their centres; omega=1 gives Floater-Hormann
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
 real(dp) :: wk(2*(ceiling(nwin2)+lm_blend+1)+lm_blend+1),bin(0:lm_blend)
 integer :: ii,kk,id,ilo,ihi,nreach,klo,khi

 if(inode>0) then
   fh_value=ff(inode)
   return
 end if
 id=min(lm_blend,nn-1)
! the nodes of the windows whose centres are within nwin2 of tt
 nreach=ceiling(nwin2)+id+1
 ilo=max(1,floor(tt)+1-nreach)
 ihi=min(nn,floor(tt)+1+nreach)
! the window weights of the windows those nodes lie in, and the binomials,
! once (the node loop below would evaluate each window id+1 times)
 klo=max(1,ilo-id)
 khi=min(ihi,nn-id)
 do kk=klo,khi
   wk(kk-klo+1)=window(abs(tt-dble(kk-1)-0.5d0*dble(id)))
 end do
 do kk=0,id
   bin(kk)=binomial(id,kk)
 end do
 num=0.0d0
 den=0.0d0
 do ii=ilo,ihi
   dt=tt-dble(ii-1)
   if(dt==0.0d0) then
     fh_value=ff(ii)
     return
   end if
   ww=0.0d0
   do kk=max(1,ii-id),min(ii,nn-id)
     ww=ww+wk(kk-klo+1)*bin(ii-kk)
   end do
   if(ww==0.0d0) cycle
   if(mod(ii,2)==0) ww=-ww
   ww=ww/dt
   num=num+ww*ff(ii)
   den=den+ww
 end do
 fh_value=num/den
 return
 end function fh_value

 real(dp) function window(ss)

!! the C-infinity bump of the window weights: 1 for ss<=nwin1, 0 for
!! ss>=nwin2, and the smooth step exp(-1/z)/(exp(-1/z)+exp(-1/(1-z)))
!! between, z=(nwin2-ss)/(nwin2-nwin1); that is 1/(1+exp(t)), t=1/z-
!! 1/(1-z), one exponential, taken as exp(-t)/(1+exp(-t)) for t>0 so
!! that it cannot overflow
!!
!! @param ss  distance of a window's centre from the position (steps)

!Input variables
 real(dp), intent(in) :: ss

!Local variables
 real(dp) :: zz,tt,ee

 if(ss<=nwin1) then
   window=1.0d0
 else if(ss>=nwin2) then
   window=0.0d0
 else
   zz=(nwin2-ss)/(nwin2-nwin1)
   tt=1.0d0/zz-1.0d0/(1.0d0-zz)
   if(tt>0.0d0) then
     ee=exp(-tt)
     window=ee/(1.0d0+ee)
   else
     window=1.0d0/(1.0d0+exp(tt))
   end if
 end if
 return
 end function window

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

 subroutine lm_cell_weights(ishift,ww)

!! the weights of the cell rule (cell_table), for the uniform grids of
!! ugrid_m

!Input variables
 integer, intent(in) :: ishift

!Output variables
 real(dp), intent(out) :: ww(ncell)

 ww(:)=cell_table(:,ishift)
 return
 end subroutine lm_cell_weights

 real(dp) function lm_fh_value(ff,nn,tt)

!! the interpolant (fh_value) of values on equispaced nodes at the position
!! tt, in units of the spacing, 0 at the first node, for the uniform grids
!! of ugrid_m
!!
!! @param ff  values on the nodes
!! @param nn  number of nodes
!! @param tt  position, 0<=tt<=nn-1

!Input variables
 integer, intent(in) :: nn
 real(dp), intent(in) :: ff(nn),tt

!Local variables
 integer :: inode

 inode=0
 if(tt==dble(nint(tt))) inode=nint(tt)+1
 if(inode<1 .or. inode>nn) inode=0
 lm_fh_value=fh_value(ff,nn,tt,inode)
 return
 end function lm_fh_value

 pure integer function lm_fh_reach()

!! the nodes on each side of a position that its interpolant (fh_value)
!! can use: those of the windows whose centres are within nwin2 of it

 lm_fh_reach=ceiling(nwin2)+lm_blend+1
 return
 end function lm_fh_reach

 subroutine lm_fh_weights(nn,tt,ilo,ihi,wts)

!! the weights of the interpolant (fh_value) at the position tt: its value
!! is sum_{i=ilo}^{ihi} wts(i-ilo+1) ff(i). Away from the ends of the
!! node range they depend on tt only through its fractional part, so one
!! set serves every position with that offset (outwf_m)
!!
!! @param nn  number of nodes
!! @param tt  position in units of the node spacing, 0 at the first node
!! @param ilo,ihi  the range of nodes with weights
!! @param wts  the weights, wts(1:ihi-ilo+1), of size 2*lm_fh_reach()+2
!!             at least

!Input variables
 integer, intent(in) :: nn
 real(dp), intent(in) :: tt

!Output variables
 integer, intent(out) :: ilo,ihi
 real(dp), intent(out) :: wts(2*lm_fh_reach()+2)

!Local variables
 real(dp) :: ww,dt,den
 real(dp) :: wk(2*(ceiling(nwin2)+lm_blend+1)+lm_blend+1),bin(0:lm_blend)
 integer :: ii,kk,id,nreach,klo,khi,inode

 inode=0
 if(tt==dble(nint(tt))) inode=nint(tt)+1
 if(inode>=1 .and. inode<=nn) then
   ilo=inode
   ihi=inode
   wts(1)=1.0d0
   return
 end if
 id=min(lm_blend,nn-1)
 nreach=ceiling(nwin2)+id+1
 ilo=max(1,floor(tt)+1-nreach)
 ihi=min(nn,floor(tt)+1+nreach)
 klo=max(1,ilo-id)
 khi=min(ihi,nn-id)
 do kk=klo,khi
   wk(kk-klo+1)=window(abs(tt-dble(kk-1)-0.5d0*dble(id)))
 end do
 do kk=0,id
   bin(kk)=binomial(id,kk)
 end do
 den=0.0d0
 do ii=ilo,ihi
   dt=tt-dble(ii-1)
   ww=0.0d0
   do kk=max(1,ii-id),min(ii,nn-id)
     ww=ww+wk(kk-klo+1)*bin(ii-kk)
   end do
   if(mod(ii,2)==0) ww=-ww
   wts(ii-ilo+1)=ww/dt
   den=den+wts(ii-ilo+1)
 end do
 do ii=ilo,ihi
   wts(ii-ilo+1)=wts(ii-ilo+1)/den
 end do
 return
 end subroutine lm_fh_weights

 end module logmesh_m

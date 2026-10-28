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
!> Integrals, values and roots of a function given on a uniform grid, the
!> counterpart of logmesh_m for the linear meshes in r and in q of the
!> residual kinetic energy (doc/smoothness_plan.md, Q-1 to Q-4): the same
!> 8th-order cell rule (each cell integrates the degree-7 interpolant on
!> the 8 nodes about it) and the same windowed Floater-Hormann interpolant,
!> here in the grid's index. Positions are in units of the spacing, 0 at
!> the first node. A grid whose spacing is a continuous function of the
!> inputs, with a node count that is not, gives results whose only
!> discontinuity is the rule's own error, below round-off at 8th order on
!> the grids used.

 module ugrid_m
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: ug_weights, ug_cumint_down, ug_interp, ug_root, ug_cross_down
 public :: ug_brent, ug_fmin, ug_fn

!> a function of one variable, for ug_brent
 abstract interface
   real(dp) function ug_fn(xx)
   import :: dp
   implicit none
   real(dp), intent(in) :: xx
   end function ug_fn
 end interface

 contains

 subroutine ug_weights(nn,hh,ww)

!! weights of the integral from the first node to the last: the sum over
!! the cells of the cell rule's weights, hh inside and corrected over the
!! first and last lm_ncell-1 nodes
!!
!! @param nn  number of nodes, >=lm_ncell
!! @param hh  spacing
!! @param ww  the weights

 use logmesh_m, only: lm_ncell, lm_cell_weights

!Input variables
 integer, intent(in) :: nn
 real(dp), intent(in) :: hh

!Output variables
 real(dp), intent(out) :: ww(nn)

!Local variables
 real(dp) :: wc(lm_ncell,0:lm_ncell-2)
 integer :: ii,jj,istart,ishift

 if(nn<lm_ncell) then
   write(error_unit,'(/a,i0)') &
&       'ug_weights: ERROR fewer nodes than the rule, nn=',nn
   error stop 1
 end if
 do ishift=0,lm_ncell-2
   call lm_cell_weights(ishift,wc(:,ishift))
 end do
 ww(:)=0.0d0
 do ii=1,nn-1
   istart=min(max(ii-(lm_ncell/2-1),1),nn-lm_ncell+1)
   ishift=ii-istart
   do jj=1,lm_ncell
     ww(istart+jj-1)=ww(istart+jj-1) + hh*wc(jj,ishift)
   end do
 end do
 return
 end subroutine ug_weights

 subroutine ug_cumint_down(ff,nn,hh,cum)

!! cumulative integral from each node to the last, accumulated from the
!! last node (so that a small value far from it is not a difference of
!! two large ones), compensated (Kahan)
!!
!! @param ff  integrand on the nodes
!! @param nn  number of nodes, >=lm_ncell
!! @param hh  spacing
!! @param cum  int_{x_i}^{x_nn} ff dx in cum(i)

 use logmesh_m, only: lm_ncell, lm_cell_weights

!Input variables
 integer, intent(in) :: nn
 real(dp), intent(in) :: ff(nn),hh

!Output variables
 real(dp), intent(out) :: cum(nn)

!Local variables
 real(dp) :: wc(lm_ncell,0:lm_ncell-2),cell,sum_hi,sum_lo,yy
! volatile: the fast math of Release (-fassociative-math) would fold the
! compensation (tmp-sum_hi)-yy to 0
 real(dp), volatile :: tmp
 integer :: ii,jj,istart,ishift

 if(nn<lm_ncell) then
   write(error_unit,'(/a,i0)') &
&       'ug_cumint_down: ERROR fewer nodes than the rule, nn=', &
&        nn
   error stop 1
 end if
 do ishift=0,lm_ncell-2
   call lm_cell_weights(ishift,wc(:,ishift))
 end do
 sum_hi=0.0d0
 sum_lo=0.0d0
 cum(nn)=0.0d0
 do ii=nn-1,1,-1
   istart=min(max(ii-(lm_ncell/2-1),1),nn-lm_ncell+1)
   ishift=ii-istart
   cell=0.0d0
   do jj=1,lm_ncell
     cell=cell + wc(jj,ishift)*ff(istart+jj-1)
   end do
   yy=hh*cell - sum_lo
   tmp=sum_hi + yy
   sum_lo=(tmp - sum_hi) - yy
   sum_hi=tmp
   cum(ii)=sum_hi
 end do
 return
 end subroutine ug_cumint_down

 subroutine ug_interp(ff,nn,tt,fq)

!! the interpolant of ff at the position tt
!!
!! @param ff  values on the nodes
!! @param nn  number of nodes
!! @param tt  position in units of the spacing, 0 at the first node,
!!            0<=tt<=nn-1
!! @param fq  value at tt

 use logmesh_m, only: lm_fh_value

!Input variables
 integer, intent(in) :: nn
 real(dp), intent(in) :: ff(nn),tt

!Output variables
 real(dp), intent(out) :: fq

 if(tt<0.0d0 .or. tt>dble(nn-1)) then
   write(error_unit,'(/a,es12.4,a,i0)') &
&       'ug_interp: ERROR position off the grid, tt=', &
&        tt,' nn=',nn
   error stop 1
 end if
 fq=lm_fh_value(ff,nn,tt)
 return
 end subroutine ug_interp

 subroutine ug_root(ff,nn,ta,tb,tq,ierr)

!! root, between the positions ta and tb, of the interpolant of ff, by
!! Brent's method with no tolerance of its own (ug_brent)
!!
!! @param ff  values on the nodes
!! @param nn  number of nodes
!! @param ta,tb  bracket (positions); the interpolant must change sign
!!               between them
!! @param tq  the root (position)
!! @param ierr  0, or 1 if the interpolant has the same sign at ta and tb

 use logmesh_m, only: lm_fh_value

!Input variables
 integer, intent(in) :: nn
 real(dp), intent(in) :: ff(nn),ta,tb

!Output variables
 real(dp), intent(out) :: tq
 integer, intent(out) :: ierr

 call ug_brent(interpolant,ta,tb,tq,ierr)
 return

 contains

 real(dp) function interpolant(tt)
!! the interpolant of ff at the position tt
 real(dp), intent(in) :: tt
 interpolant=lm_fh_value(ff,nn,tt)
 end function interpolant

 end subroutine ug_root

 subroutine ug_brent(fun,ta,tb,tq,ierr,fq,fo,fends,converged)

!! root of fun between ta and tb by Brent's method with no tolerance of
!! its own: it stops when the bracket is a few rounding units wide; after
!! 200 steps it reports through converged, or without it stops in error
!!
!! @param fun  the function
!! @param ta,tb  bracket; fun must change sign between them
!! @param tq  the root
!! @param ierr  0, or 1 if fun has the same sign at ta and tb
!! @param fq  fun at tq (optional output)
!! @param fo  fun at the other end of the last bracket (optional output)
!! @param fends  fun at ta and tb (optional output)
!! @param converged  whether the bracket closed before the step limit
!!                   (optional output)

!Input variables
 procedure(ug_fn) :: fun
 real(dp), intent(in) :: ta,tb

!Output variables
 real(dp), intent(out) :: tq
 integer, intent(out) :: ierr
 real(dp), intent(out), optional :: fq,fo,fends(2)
 logical, intent(out), optional :: converged

!Local variables
 real(dp) :: xa,xb,xc,fa,fb,fc,de,dlast,tol1,xm,pp,qq,rs,ss
 integer :: iter

 xa=ta
 xb=tb
 fa=fun(xa)
 fb=fun(xb)
 if(present(fends)) fends(:)=[fa,fb]
 if(present(converged)) converged=.true.
 ierr=0
 tq=ta
 if(present(fq)) fq=fa
 if(present(fo)) fo=fa
 if(fa==0.0d0) return
 tq=tb
 if(present(fq)) fq=fb
 if(present(fo)) fo=fb
 if(fb==0.0d0) return
 if(sign(1.0d0,fa)==sign(1.0d0,fb)) then
   ierr=1
   return
 end if

 xc=xa
 fc=fa
 de=xb-xa
 dlast=de
 if(present(converged)) converged=.false.
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
   if(abs(xm)<=tol1 .or. fb==0.0d0) then
     if(present(converged)) converged=.true.
     exit
   end if
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
   fb=fun(xb)
 end do
 if(iter>200 .and. .not. present(converged)) then
   write(error_unit,'(/a,es12.4,a,es12.4)') &
&       'ug_brent: ERROR not converged in 200 steps, t=',xb,', f=',fb
   error stop 1
 end if
 tq=xb
 if(present(fq)) fq=fb
 if(present(fo)) fo=fc
 return
 end subroutine ug_brent

 real(dp) function ug_fmin(fun,ax,bx)

!! the least value of fun on [ax, bx], by Brent's method (golden section
!! and parabolic steps) to a relative step of 1e-8, and at the ends; the
!! iteration limit is an error
!!
!! @param fun  the function
!! @param ax,bx  the interval

!Input variables
 procedure(ug_fn) :: fun
 real(dp), intent(in) :: ax,bx

!Local variables
 real(dp), parameter :: cgold=0.3819660112501051d0, rtol=1.0d-8
 real(dp) :: aa,bb,vv,ww,xx,ee,dd,fvv,fww,fxx,xm,tol1,tol2,rr1,qq,pp,uu,fu
 integer :: iter
 aa=min(ax,bx)
 bb=max(ax,bx)
 xx=aa + cgold*(bb-aa)
 vv=xx
 ww=xx
 ee=0.0d0
 dd=0.0d0
 fxx=fun(xx)
 fvv=fxx
 fww=fxx
 do iter=1,100
   xm=0.5d0*(aa+bb)
   tol1=rtol*abs(xx) + 1.0d-12
   tol2=2.0d0*tol1
   if(abs(xx-xm)<=tol2-0.5d0*(bb-aa)) exit
   if(iter==100) then
     write(error_unit,'(/a,es12.4,a,es12.4)') &
&          'ug_fmin: ERROR not converged in 100 steps, x=',xx, &
&          ', bracket width=',bb-aa
     error stop 1
   end if
   if(abs(ee)>tol1) then
     rr1=(xx-ww)*(fxx-fvv)
     qq=(xx-vv)*(fxx-fww)
     pp=(xx-vv)*qq - (xx-ww)*rr1
     qq=2.0d0*(qq-rr1)
     if(qq>0.0d0) pp=-pp
     qq=abs(qq)
     if(abs(pp)>=abs(0.5d0*qq*ee) .or. pp<=qq*(aa-xx) &
&       .or. pp>=qq*(bb-xx)) then
       if(xx>=xm) then
         ee=aa-xx
       else
         ee=bb-xx
       end if
       dd=cgold*ee
     else
       ee=dd
       dd=pp/qq
       uu=xx+dd
       if(uu-aa<tol2 .or. bb-uu<tol2) dd=sign(tol1,xm-xx)
     end if
   else
     if(xx>=xm) then
       ee=aa-xx
     else
       ee=bb-xx
     end if
     dd=cgold*ee
   end if
   if(abs(dd)>=tol1) then
     uu=xx+dd
   else
     uu=xx+sign(tol1,dd)
   end if
   fu=fun(uu)
   if(fu<=fxx) then
     if(uu>=xx) then
       aa=xx
     else
       bb=xx
     end if
     vv=ww
     fvv=fww
     ww=xx
     fww=fxx
     xx=uu
     fxx=fu
   else
     if(uu<xx) then
       aa=uu
     else
       bb=uu
     end if
     if(fu<=fww .or. ww==xx) then
       vv=ww
       fvv=fww
       ww=uu
       fww=fu
     else if(fu<=fvv .or. vv==xx .or. vv==ww) then
       vv=uu
       fvv=fu
     end if
   end if
 end do
 ug_fmin=min(fxx,fun(ax),fun(bx))
 return
 end function ug_fmin

 subroutine ug_cross_down(ff,nn,hh,x1,level,xq,found)

!! the largest x at which a positive, decreasing curve given on the grid
!! falls to level: the last cell (from the top) in which ff crosses level
!! from above, then the root in it of the interpolant of ln ff (ln ff is
!! the smoother function; the convergence curves of run_optimize are near
!! exponential). The curve is taken up to its first value that is not
!! positive. found is false if it never crosses level
!!
!! @param ff  values on the nodes
!! @param nn  number of nodes
!! @param hh  spacing
!! @param x1  the first node's coordinate
!! @param level  the level, >0
!! @param xq  the crossing (coordinate)
!! @param found  whether there is one

!Input variables
 integer, intent(in) :: nn
 real(dp), intent(in) :: ff(nn),hh,x1,level

!Output variables
 real(dp), intent(out) :: xq
 logical, intent(out) :: found

!Local variables
 real(dp) :: lf(nn),tq
 integer :: ii,np,ierr

 found=.false.
 xq=0.0d0
 np=nn
 do ii=1,nn
   if(.not. (ff(ii)>0.0d0)) then
     np=ii-1
     exit
   end if
 end do
 if(np<2) return
 lf(1:np)=log(ff(1:np))-log(level)
 do ii=np,2,-1
   if(lf(ii-1)>0.0d0 .and. lf(ii)<=0.0d0) then
     call ug_root(lf(1:np),np,dble(ii-2),dble(ii-1),tq,ierr)
     if(ierr/=0) then
       write(error_unit,'(/a)') &
&           'ug_cross_down: ERROR the bracket lost its sign change'
       error stop 1
     end if
     xq=x1+hh*tq
     found=.true.
     return
   end if
 end do
 return
 end subroutine ug_cross_down

 end module ugrid_m

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
!> One-dimensional searches: the root of a function by Brent's method, its
!> least value on an interval by Brent's minimizer, and the stopping width
!> every root search and polish in the programs uses (root_tol). A leaf
!> module, so that logmesh_m and ugrid_m, and the searches outside them,
!> share one Brent.

 module root_m
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: root_fn, root_brent, root_fmin, root_tol

!> a function of one variable, for root_brent and root_fmin
 abstract interface
   real(dp) function root_fn(xx)
   import :: dp
   implicit none
   real(dp), intent(in) :: xx
   end function root_fn
 end interface

 contains

 pure real(dp) function root_tol(xx,xabs)

!! the width at which a root search or polish stops: four rounding units
!! of |xx| plus xabs, the absolute part for an xx that can be near 0
!!
!! @param xx  the root, or the scale of a function's values
!! @param xabs  the absolute part (0 where xx is bounded away from 0)

!Input variables
 real(dp), intent(in) :: xx,xabs

 root_tol=4.0d0*epsilon(1.0d0)*(abs(xx)+xabs)
 return
 end function root_tol

 subroutine root_brent(fun,ta,tb,tq,ierr,ifail,fq,fo,fends,converged)

!! root of fun between ta and tb by Brent's method with no tolerance of
!! its own: it stops when the bracket is a few rounding units wide
!! (root_tol, with an absolute part of 5e-4 for a root near 0); after
!! 200 steps it reports through converged, or without it fails
!!
!! @param fun  the function
!! @param ta,tb  bracket; fun must change sign between them
!! @param tq  the root
!! @param ierr  0, or 1 if fun has the same sign at ta and tb
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
!!               (fail_convergence after 200 steps without converged)
!! @param fq  fun at tq (optional output)
!! @param fo  fun at the other end of the last bracket (optional output)
!! @param fends  fun at ta and tb (optional output)
!! @param converged  whether the bracket closed before the step limit
!!                   (optional output)

 use failure_m, only: fail_convergence

!Input variables
 procedure(root_fn) :: fun
 real(dp), intent(in) :: ta,tb

!Output variables
 real(dp), intent(out) :: tq
 integer, intent(out) :: ierr,ifail
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
 ifail=0
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
   tol1=0.5d0*root_tol(xb,5.0d-4)
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
&       'root_brent: ERROR not converged in 200 steps, x=',xb,', f=',fb
   ifail=fail_convergence
   return
 end if
 tq=xb
 if(present(fq)) fq=fb
 if(present(fo)) fo=fc
 return
 end subroutine root_brent

 subroutine root_fmin(fun,ax,bx,fmin,ifail)

!! the least value of fun on [ax, bx], by Brent's method (golden section
!! and parabolic steps) to a relative step of 1e-8, and at the ends; the
!! iteration limit is a failure
!!
!! @param fun  the function
!! @param ax,bx  the interval
!! @param fmin  the least value
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
!!               (fail_convergence after 100 steps)

 use failure_m, only: fail_convergence

!Input variables
 procedure(root_fn) :: fun
 real(dp), intent(in) :: ax,bx

!Output variables
 real(dp), intent(out) :: fmin
 integer, intent(out) :: ifail

!Local variables
 real(dp), parameter :: cgold=0.3819660112501051d0, rtol=1.0d-8
 real(dp) :: aa,bb,vv,ww,xx,ee,dd,fvv,fww,fxx,xm,tol1,tol2,rr1,qq,pp,uu,fu
 integer :: iter
 ifail=0
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
&          'root_fmin: ERROR not converged in 100 steps, x=',xx, &
&          ', bracket width=',bb-aa
     ifail=fail_convergence
     return
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
 fmin=min(fxx,fun(ax),fun(bx))
 return
 end subroutine root_fmin

 end module root_m

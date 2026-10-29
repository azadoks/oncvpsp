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
!> The matching of the outward and inward integrations of a bound-state
!> solver at a continuous point (doc/smoothness_plan.md, S-8). Matched at
!> a mesh point, the turning-point index, the eigenvalue of the discrete
!> equations jumps each time that index moves (their outward and inward
!> truncation errors differ, so their Wronskian is not quite constant
!> across the mesh): about 1e-13 Ha, which a well state carries into its
!> well width and everything built on it. Here the matching point is the
!> turning point itself, a continuous mesh index xm (the root of the
!> interpolant of the turning-point coefficient), and the two solutions
!> are compared and joined over the points within match_k of it, with
!> C-infinity weights of the distance to xm:
!>   the energy correction is the bump-weighted average over those points
!>   of the pointwise perturbative correction, and the inward solution's
!>   scale the ratio of the bump-weighted sums of the two;
!>   the returned function is the outward one below xm-match_k, the scaled
!>   inward one above xm+match_k, and the smooth-step blend of the two in
!>   between (they differ there by the truncation error, about 1e-13).
!> Every result is then a smooth function of the energy and the potential.
!> The solutions must be valid on mch-1-match_k..mch+match_k, mch the
!> mesh point just beyond the turning point: the outward integration is
!> carried to mch+match_k, the inward one to mch-1-match_k.

 module match_m
 implicit none
 private
 public :: match_k, match_point, match_weights, match_smax, match_blend

!> half width, in mesh steps, of the matching window
 integer, parameter :: match_k=4

 contains

 subroutine match_point(rr,cf,mmax,mch,xm,iout,ifail)

!! continuous mesh index xm in [mch-1,mch] of the root of the interpolant
!! of cf, which changes sign from <=0 to >0 between rr(mch-1) and rr(mch)
!!
!! @param rr  log radial mesh (bohr)
!! @param cf  the coefficient whose sign change is the turning point
!! @param mmax  size of the log radial mesh
!! @param mch  the mesh point just beyond the sign change
!! @param xm  the continuous index of the root (1 at rr(1))
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_input, fail_not_found, failed, &
&                    msg_match_point_sign, msg_match_point_window
 use constants_m, only: dp
 use logmesh_m, only: lm_root, lm_step

!Input variables
 integer, intent(in) :: mmax,mch
 real(dp), intent(in) :: rr(mmax),cf(mmax)
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: xm
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: rq,al
 integer :: ierr

 ifail=0
 if(mch-1-match_k<2 .or. mch+match_k>mmax) then
   write(iout,'(/a,i0,a,i0)') 'match_point: ERROR matching window at mch=', &
&        mch,' leaves the mesh, mmax=',mmax
   ifail=failed(fail_input,msg_match_point_window)
   return
 end if
 call lm_root(rr,cf,mmax,rr(mch-1),rr(mch),rq,ierr,ifail)
 if(ifail/=0) return
 if(ierr/=0) then
   write(iout,'(/a,i0)') 'match_point: ERROR no sign change below mch=',mch
   ifail=failed(fail_not_found,msg_match_point_sign)
   return
 end if
 al=lm_step(rr,mmax)
 xm=dble(mch-1) + log(rq/rr(mch-1))/al
 return
 end subroutine match_point

 subroutine match_weights(xm,lo,hi,ww,ss)

!! the weights over the mesh points lo..hi about the continuous index xm:
!! ww a C-infinity bump of (n-xm)/match_k summing to 1, ss a C-infinity
!! step of it, 0 below xm-match_k and 1 above xm+match_k
!!
!! @param xm  continuous mesh index of the matching point
!! @param lo,hi  the points, covering xm-match_k..xm+match_k
!! @param ww  bump weights
!! @param ss  step weights

 use constants_m, only: dp

!Input variables
 integer, intent(in) :: lo,hi
 real(dp), intent(in) :: xm

!Output variables
 real(dp), intent(out) :: ww(lo:hi),ss(lo:hi)

!Local variables
 real(dp) :: tt
 integer :: nn

 do nn=lo,hi
   tt=(dble(nn)-xm)/dble(match_k)
   if(tt<=-1.0d0) then
     ww(nn)=0.0d0
     ss(nn)=0.0d0
   else if(tt>=1.0d0) then
     ww(nn)=0.0d0
     ss(nn)=1.0d0
   else
     ww(nn)=exp(-1.0d0/(1.0d0-tt**2))
     ss(nn)=edge(1.0d0+tt)/(edge(1.0d0+tt)+edge(1.0d0-tt))
   end if
 end do
 ww(:)=ww(:)/sum(ww)
 return
 end subroutine match_weights

 subroutine match_blend(lo,hi,mch,ss,fo,ff)

!! the outward solution fo and the inward one ff blended across the
!! window lo..hi by the step ss of match_weights: ff = (1-ss) fo + ss ff
!!
!! @param lo,hi  the points of the window
!! @param mch  the matching point, at which fo and ss are indexed 0
!! @param ss  step weights
!! @param fo  outward solution on the window
!! @param ff  inward solution, the blend on the window on return

 use constants_m, only: dp

!Input variables
 integer, intent(in) :: lo,hi,mch
 real(dp), intent(in) :: ss(lo-mch:hi-mch),fo(lo-mch:hi-mch)

!Input/output variables
 real(dp), intent(inout) :: ff(hi)

!Local variables
 integer :: ii

 do ii=lo,hi
   ff(ii)=(1.0d0-ss(ii-mch))*fo(ii-mch)+ss(ii-mch)*ff(ii)
 end do
 return
 end subroutine match_blend

 real(dp) function match_smax(aa,bb)

!! a C-infinity maximum of two mesh indices, within half a step of the
!! larger: (a + b + sqrt((a-b)**2 + 1))/2
!!
!! @param aa,bb  the two mesh indices, as reals

 use constants_m, only: dp

!Input variables
 real(dp), intent(in) :: aa,bb

 match_smax=0.5d0*(aa+bb+sqrt((aa-bb)**2+1.0d0))
 end function match_smax

 real(dp) function edge(xx)

!! exp(-1/x) for x>0, 0 otherwise
!!
!! @param xx  the argument x

 use constants_m, only: dp

!Input variables
 real(dp), intent(in) :: xx

 if(xx>0.0d0) then
   edge=exp(-1.0d0/xx)
 else
   edge=0.0d0
 end if
 end function edge

 end module match_m

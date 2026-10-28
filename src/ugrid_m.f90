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
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: ug_weights, ug_cumint_down, ug_interp, ug_root, ug_cross_down

 contains

 subroutine ug_weights(nn,hh,ww)

!! weights of the integral from the first node to the last: the sum over
!! the cells of the cell rule's weights, hh inside and corrected over the
!! first and last lm_ncell-1 nodes
!!
!! @param nn  number of nodes, >=lm_ncell
!! @param hh  spacing
!! @param ww  the weights

 use constants_m, only: dp
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

 use constants_m, only: dp
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

 use constants_m, only: dp
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
!! Brent's method with no tolerance of its own (root_brent)
!!
!! @param ff  values on the nodes
!! @param nn  number of nodes
!! @param ta,tb  bracket (positions); the interpolant must change sign
!!               between them
!! @param tq  the root (position)
!! @param ierr  0, or 1 if the interpolant has the same sign at ta and tb

 use constants_m, only: dp
 use logmesh_m, only: lm_fh_value
 use root_m, only: root_brent

!Input variables
 integer, intent(in) :: nn
 real(dp), intent(in) :: ff(nn),ta,tb

!Output variables
 real(dp), intent(out) :: tq
 integer, intent(out) :: ierr

 call root_brent(interpolant,ta,tb,tq,ierr)
 return

 contains

 real(dp) function interpolant(tt)
!! the interpolant of ff at the position tt
 real(dp), intent(in) :: tt
 interpolant=lm_fh_value(ff,nn,tt)
 end function interpolant

 end subroutine ug_root

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

 use constants_m, only: dp
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

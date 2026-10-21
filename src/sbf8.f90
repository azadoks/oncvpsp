!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
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
! calculates spherical Bessel functions using recursion algorithm

 subroutine sbf8(nm,xx,sb_out)

!! @param nm  maximum angular momentum wanted plus 1
!! @param xx  argument of spherical bessel function
!! @param sb_out  output of sbf_l(xx) for l=0,...,nm-1

 use constants_m, only: dp
 use sbf_m, only: sbf_trig
 implicit none

!Input variables
 integer, intent(in) :: nm
 real(dp), intent(in) :: xx

!Input/output variables
 real(dp), intent(inout) :: sb_out(nm)

!Local variables
 integer :: nlim,nn
 real(dp) :: cc,fn,sn,ss,xi,xn,xs
 real(dp),allocatable :: sb(:)


 if(xx<= 1.0d-36) then
!  zero argument section
   sb_out(:)=0.0d0
   sb_out(1)=1.0d0
 else if(xx<1.0d-3) then
!  small argument section
   xn=1.0d0
   xs=0.5d0*xx**2
   do nn=1,nm
     sb_out(nn)=xn*(1.0d0 - xs*(1.0d0 - xs/(4*nn+6))/(2*nn+1))
     xn=xx*xn/(2*nn+1)
   end do
!  trigonometric section (small l values), every order from sin and cos
!  (sbf_trig: the series where the closed forms cancel)
 else if(nm<=4) then
   ss=sin(xx)
   cc=cos(xx)
   do nn=1,nm
     sb_out(nn)=sbf_trig(nn-1,xx,ss,cc)
   end do
 else
!  recursion method (accurate for large l, slow for large arguments),
!  started at order nm+1.36 xx+15 at every xx (the start's error is then
!  below round-off; at small xx a start at nm+1 had cost the top orders
!  1e-9 and jumped where it moved), normalized by
!  sum_l (2l+1) j_l**2 = 1. Its round-off grows with the xx oscillatory
!  orders it passes, about 0.6 eps xx
   nlim=nm+int(1.36d0*xx)+15
   allocate(sb(nlim+1))
   xi=1.0d0/xx
   sb(nlim+1)=0.0d0
   sb(nlim)=1.0d-18
   do nn=nlim-1,1,-1
     sb(nn)=dble(2*nn+1)*xi*sb(nn+1) - sb(nn+2)
     if(abs(sb(nn))>1.0d150) then
       sb(nn:nlim+1)=sb(nn:nlim+1)*1.0d-150
     end if
   end do
   sn=0.0d0
   do nn=1,nlim
     sn=sn + dble(2*nn-1)*sb(nn)*sb(nn)
   end do
   fn=1.0d0/sqrt(sn)
   sb_out(:)=fn*sb(1:nm)
   deallocate(sb)
 end if
 return
 end subroutine sbf8

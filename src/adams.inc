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
! The Adams-Moulton step of order k (adams_m): the sums over the given
! points, and the exact solution of the implicit step for a linear
! equation of two components
!
! Included as internal functions (after "contains") by the radial
! integrators that use them, so the compiler can inline them

 function am_out(adams,yy,jj)

!! the given part of the outward step from jj to jj+1,
!! sum_{i=1}^{k-1} cm(i) yy(jj+1-i)
!!
!! @param adams  the order and coefficients
!! @param yy  derivative of the integrated function in the mesh index
!! @param jj  mesh index of the current point

 use constants_m, only: dp
 use adams_m, only: adams_t
 implicit none

!Input variables
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: yy(*)
 integer, intent(in) :: jj

!Output variables
 real(dp) :: am_out

!Local variables
 integer :: ii

 am_out=0.0d0
 do ii=1,adams%order-1
   am_out=am_out+adams%cm(ii)*yy(jj+1-ii)
 end do
 return
 end function am_out

 function am_in(adams,yy,jj)

!! the given part of the inward step from jj to jj-1,
!! -sum_{i=1}^{k-1} cm(i) yy(jj-1+i)
!!
!! @param adams  the order and coefficients
!! @param yy  derivative of the integrated function in the mesh index
!! @param jj  mesh index of the current point

 use constants_m, only: dp
 use adams_m, only: adams_t
 implicit none

!Input variables
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: yy(*)
 integer, intent(in) :: jj

!Output variables
 real(dp) :: am_in

!Local variables
 integer :: ii

 am_in=0.0d0
 do ii=1,adams%order-1
   am_in=am_in-adams%cm(ii)*yy(jj-1+ii)
 end do
 return
 end function am_in

 subroutine am_solve(c0,a11,a12,a21,a22,b1,b2,y1,y2)

!! the implicit step: y = b + c0 (A y + s) solved for y, the right-hand
!! side b holding the current point, the given part of the step and
!! c0 s. That is (1 - c0 A) y = b; c0 is cm(0) outward, -cm(0) inward
!!
!! @param c0  coefficient of the derivative at the new point
!! @param a11  A(1,1) at the new point
!! @param a12  A(1,2)
!! @param a21  A(2,1)
!! @param a22  A(2,2)
!! @param b1  first component of the right-hand side
!! @param b2  second component
!! @param y1  first component at the new point
!! @param y2  second component at the new point

 use constants_m, only: dp
 implicit none

!Input variables
 real(dp), intent(in) :: c0,a11,a12,a21,a22,b1,b2

!Output variables
 real(dp), intent(out) :: y1,y2

!Local variables
 real(dp) :: d11,d12,d21,d22,det

 d11=1.0d0-c0*a11
 d12=-c0*a12
 d21=-c0*a21
 d22=1.0d0-c0*a22
 det=d11*d22-d12*d21
 y1=(d22*b1-d12*b2)/det
 y2=(d11*b2-d21*b1)/det
 return
 end subroutine am_solve

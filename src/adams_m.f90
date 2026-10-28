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
!> The Adams-Moulton formulas of the radial integrators, of order 4 to 8
!> (adams_order of the input, 4 by default). The step of order k from
!> point n to n+1 (or n-1 inward) is
!>   y(n+1) = y(n) + cm(0) y'(n+1) + sum_{i=1}^{k-1} cm(i) y'(n+1-i),
!> y' the derivative in the mesh index. The equations integrated are
!> linear, y' = A y + s for two components, so the step is solved exactly
!> by a 2x2 linear solve (am_solve in adams.inc) rather than iterated
!> from an Adams-Bashforth predictor. The integrators start from k given
!> points.

 module adams_m
 use constants_m, only: dp
 implicit none
 private
 public :: adams_t, adams_init, adams_min, adams_max

!> the lowest order
 integer, parameter :: adams_min=4
!> the highest order
 integer, parameter :: adams_max=8

!> the coefficients of the order. The default is order 4, so that a
!> declared adams_t needs no adams_init for it
 type adams_t
!> the order k: the step uses the derivatives at k points
   integer :: order=4
!> the coefficients cm(0:k-1), 0 beyond
   real(dp) :: cm(0:adams_max-1)= &
&    [9.0d0, 19.0d0, -5.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]/24.0d0
 end type adams_t

 contains

 function adams_init(order) result(adams)

!! the coefficients of the given order, adams_min to adams_max (the caller
!! checks it): the integrals over [0,1] of the Lagrange polynomials on the
!! nodes 1, 0, -1, ..., 2-k (Abramowitz and Stegun, p. 896), as exact
!! rationals
!!
!! @param order  the order k

!Input variables
 integer, intent(in) :: order

!Output variables
 type(adams_t) :: adams

!Local variables
 real(dp), parameter :: num(0:adams_max-1,adams_min:adams_max)=reshape([ &
&   9.0d0, 19.0d0, -5.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
&   251.0d0, 646.0d0, -264.0d0, 106.0d0, -19.0d0, 0.0d0, 0.0d0, 0.0d0, &
&   475.0d0, 1427.0d0, -798.0d0, 482.0d0, -173.0d0, 27.0d0, 0.0d0, &
&   0.0d0, &
&   19087.0d0, 65112.0d0, -46461.0d0, 37504.0d0, -20211.0d0, 6312.0d0, &
&   -863.0d0, 0.0d0, &
&   36799.0d0, 139849.0d0, -121797.0d0, 123133.0d0, -88547.0d0, &
&   41499.0d0, -11351.0d0, 1375.0d0], [adams_max,adams_max-adams_min+1])
 real(dp), parameter :: den(adams_min:adams_max)= &
&   [24.0d0, 720.0d0, 1440.0d0, 60480.0d0, 120960.0d0]

 adams%order=order
 adams%cm(:)=num(:,order)/den(order)

 return
 end function adams_init

 end module adams_m

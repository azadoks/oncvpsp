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
! Adams extrapolation and interpolation formulas for
! outward and inward integration, Abramowitz and
! Stegun, p. 896
!
! Included as internal functions (after "contains") by the radial
! integrators that use them, so the compiler can inline them

 function aeo(yy,jj)

!Adams outward extrapolation (predictor) to point jj+1:
!increment of the integrated function per mesh step

!yy  derivative of the integrated function with respect to the mesh index
!jj  mesh index of the current point

 use constants_m, only: dp
 implicit none

!Input variables
 real(dp), intent(in) :: yy(*)
 integer, intent(in) :: jj

!Output variables
 real(dp) :: aeo

 aeo=(1.0d0/24.0d0)*(55.0d0*yy(jj)-59.0d0*yy(jj-1) &
& +37.0d0*yy(jj-2)-9.0d0*yy(jj-3))
 return
 end function aeo

 function aio(yy,jj)

!Adams outward interpolation (corrector) to point jj+1:
!increment of the integrated function per mesh step

!yy  derivative of the integrated function with respect to the mesh index
!jj  mesh index of the current point

 use constants_m, only: dp
 implicit none

!Input variables
 real(dp), intent(in) :: yy(*)
 integer, intent(in) :: jj

!Output variables
 real(dp) :: aio

 aio=(1.0d0/24.0d0)*(9.0d0*yy(jj+1)+19.0d0*yy(jj) &
& -5.0d0*yy(jj-1)+yy(jj-2))
 return
 end function aio

 function aei(yy,jj)

!Adams inward extrapolation (predictor) to point jj-1:
!increment of the integrated function per mesh step

!yy  derivative of the integrated function with respect to the mesh index
!jj  mesh index of the current point

 use constants_m, only: dp
 implicit none

!Input variables
 real(dp), intent(in) :: yy(*)
 integer, intent(in) :: jj

!Output variables
 real(dp) :: aei

 aei=-(1.0d0/24.0d0)*(55.0d0*yy(jj)-59.0d0*yy(jj+1) &
& +37.0d0*yy(jj+2)-9.0d0*yy(jj+3))
 return
 end function aei

 function aii(yy,jj)

!Adams inward interpolation (corrector) to point jj-1:
!increment of the integrated function per mesh step

!yy  derivative of the integrated function with respect to the mesh index
!jj  mesh index of the current point

 use constants_m, only: dp
 implicit none

!Input variables
 real(dp), intent(in) :: yy(*)
 integer, intent(in) :: jj

!Output variables
 real(dp) :: aii

 aii=-(1.0d0/24.0d0)*(9.0d0*yy(jj-1)+19.0d0*yy(jj) &
& -5.0d0*yy(jj+1)+yy(jj+2))
 return
 end function aii

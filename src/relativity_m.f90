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
!> The speed of light of the relativistic equations and of the relativistic
!> exchange (iexc=7): its default, the value of the code before it could be
!> set, and the value of a run from its input (c_light_of). The run passes
!> it to the solvers and the exchange as an argument (c_light), so that two
!> runs at once, in threads, can have different values.

 module relativity_m
 use constants_m, only: dp
 implicit none
 private
 public :: c_light_default, c_light_of

!> default speed of light in atomic units (Hartree), 1/alpha.  The 4.0.1
!> value is kept, rather than CODATA 137.035999..., so results are unchanged
 real(dp), parameter :: c_light_default=137.036d0

 contains

 pure real(dp) function c_light_of(c_light_in)

!! the speed of light of a run from the input's
!!
!! @param c_light_in  speed of light of the input, 0 or less for
!!                    c_light_default

!Input variables
 real(dp), intent(in) :: c_light_in

 if(c_light_in>0.0d0) then
   c_light_of=c_light_in
 else
   c_light_of=c_light_default
 end if

 return
 end function c_light_of

 end module relativity_m

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
!> exchange (iexc=6), set once from the input (set_relativity) and fixed
!> afterwards.  The default is the value of the code before it could be set.

 module relativity_m
 use constants_m, only: dp
 implicit none
 private
 public :: c_light_default, c_light, set_relativity

!> default speed of light in atomic units (Hartree), 1/alpha.  The 4.0.1
!> value is kept, rather than CODATA 137.035999..., so results are unchanged
 real(dp), parameter :: c_light_default=137.036_dp

!> speed of light of the scalar- and fully-relativistic equations and of
!> the relativistic exchange
 real(dp), protected :: c_light=c_light_default

 contains

 subroutine set_relativity(c_light_in)

!! sets the speed of light from the input
!!
!! @param c_light_in  speed of light, 0 or less for c_light_default

!Input variables
 real(dp), intent(in) :: c_light_in

 if(c_light_in>0.0_dp) then
   c_light=c_light_in
 else
   c_light=c_light_default
 end if

 return
 end subroutine set_relativity

 end module relativity_m

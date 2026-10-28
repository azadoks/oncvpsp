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

 module cinf_taper_m

 implicit none
 private
 public :: cinf_taper

 contains

! The C-infinity taper from 1 to 0 over a window, used wherever a function
! must be cut off without a kink: the ionic tail, the Teter function of the
! model cores, the unscreened local potential where there is no charge,
! the sphere of the error of icmod 6/7.

 real(dp) function cinf_taper(rad,ra,rb)

!! C-infinity switch from 1 (r <= ra) to 0 (r >= rb): f(1-x)/(f(1-x)+f(x)),
!! f(t) = exp(-1/t), x = (r-ra)/(rb-ra); its transform decays faster than
!! any power of q, so it adds no plane-wave tail of its own
!!
!! @param rad  radius
!! @param ra  start of the window
!! @param rb  end of the window
!! @return the switch at rad

 use constants_m, only: dp

!Input variables
 real(dp), intent(in) :: rad,ra,rb

!Local variables
 real(dp) :: xx,f0,f1

 if(rad<=ra) then
   cinf_taper=1.0d0
 else if(rad>=rb) then
   cinf_taper=0.0d0
 else
   xx=(rad-ra)/(rb-ra)
   f0=exp(-1.0d0/(1.0d0-xx))
   f1=exp(-1.0d0/xx)
   cinf_taper=f0/(f0+f1)
 end if
 return
 end function cinf_taper

 end module cinf_taper_m

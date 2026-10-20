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
 subroutine sbf_overlap(ll,q1,q2,rc,ss)

!! overlap of two spherical Bessel functions of one l inside a radius, in
!! closed form (Lommel's integrals), so exact and smooth in q1, q2 and rc:
!! int_0^rc r**2 j_l(q1 r) j_l(q2 r) dr, with the q1=q2 limit (and q1 near
!! q2) taken without cancellation (see doc/smoothness_plan.md)
!!
!! STUB: returns huge() until implemented (tests/unit/test_logmesh.f90)
!!
!! @param ll  angular momentum l
!! @param q1,q2  wave vectors (bohr^-1), >0
!! @param rc  radius (bohr), >0
!! @param ss  the overlap

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: ll
 real(dp), intent(in) :: q1,q2,rc

!Output variables
 real(dp), intent(out) :: ss

 ss=huge(1.0_dp)
 return
 end subroutine sbf_overlap

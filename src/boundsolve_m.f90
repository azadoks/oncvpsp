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
!> The energy iteration the bound-state solvers share (lschfb, ldiracfb,
!> lschvkbb, lschvkbbe, lschpsbar): the convergence test on the energy
!> correction, and the steps of the energy inside its bracket, by the
!> correction or, where the node count is wrong, by bisection.
 module boundsolve_m
 use constants_m, only: dp
 implicit none
 private
 public :: bs_converged, bs_correct, bs_bisect

 contains

 logical function bs_converged(de,deprev,ee,eps,polish)

!! whether the energy ee has converged: its correction de is below
!! eps*max(|ee|,0.2) and either at round-off (8 ulps) or no longer halving
!! from the one before, deprev
!!
!! @param de  the energy correction
!! @param deprev  the correction before
!! @param ee  the energy
!! @param eps  the relative convergence threshold
!! @param polish  set .true. where de is below eps but not converged
!!                (optional)

!Input variables
 real(dp), intent(in) :: de,deprev,ee,eps

!Input/output variables
 logical, intent(inout), optional :: polish

 bs_converged=.false.
 if(abs(de)<max(abs(ee),0.2d0)*eps) then
   if(abs(de)<=8.0d0*epsilon(1.0d0)*max(abs(ee),0.2d0) .or. &
&     abs(de)>=0.5d0*abs(deprev)) then
     bs_converged=.true.
     return
   end if
   if(present(polish)) polish=.true.
 end if
 return
 end function bs_converged

 subroutine bs_correct(de,ee,emin,emax,half)

!! the energy moved by its correction de, after the bracket's end on its
!! side moves to ee (half: halfway to ee); bisection of the bracket if the
!! step leaves it
!!
!! @param de  the energy correction
!! @param ee  the energy
!! @param emin,emax  the bracket
!! @param half  whether the end moves only halfway to ee

!Input variables
 real(dp), intent(in) :: de
 logical, intent(in) :: half

!Input/output variables
 real(dp), intent(inout) :: ee,emin,emax

 if(de>0.0d0) then
   if(half) then
     emin=0.5d0*(emin+ee)
   else
     emin=ee
   end if
 else
   if(half) then
     emax=0.5d0*(emax+ee)
   else
     emax=ee
   end if
 end if
 ee=ee+de
 if(ee>emax .or. ee<emin) ee=0.5d0*(emax+emin)
 return
 end subroutine bs_correct

 subroutine bs_bisect(ee,emin,emax,too_low,half)

!! the bracket's lower end (too_low: too few nodes) or upper end moved to
!! ee (half: halfway to ee), and ee to the middle of the bracket
!!
!! @param ee  the energy
!! @param emin,emax  the bracket
!! @param too_low  whether ee is below the state
!! @param half  whether the end moves only halfway to ee

!Input variables
 logical, intent(in) :: too_low,half

!Input/output variables
 real(dp), intent(inout) :: ee,emin,emax

 if(too_low) then
   if(half) then
     emin=0.5d0*(emin+ee)
   else
     emin=ee
   end if
 else
   if(half) then
     emax=0.5d0*(emax+ee)
   else
     emax=ee
   end if
 end if
 ee=0.5d0*(emin+emax)
 return
 end subroutine bs_bisect

 end module boundsolve_m

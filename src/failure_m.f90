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
!> The kinds of failure that a routine's ifail reports. No library routine
!> stops: each failure writes its message and returns one of these, and
!> every caller returns the ifail it got, unchanged, up to the program,
!> which stops with exit status 1 for any of them (oncvpsp_main). A caller
!> from another language reads the kind from ifail; the message says which
!> failure of that kind it was.
!>
!> ierr, where a routine has one, is a status of its own that its caller
!> interprets (no root in a bracket, too few nodes, ...), documented with
!> the routine; it is not a failure until the caller makes it one.

 module failure_m
 implicit none
 private
 public :: fail_none, fail_input, fail_convergence, fail_not_found, &
&          fail_construction, fail_numerical, fail_io, fail_internal, &
&          fail_name

!> no failure
 integer, parameter :: fail_none=0
!> the input is invalid, inconsistent or asks for something not built in
!> (check_data, the readers, a functional without Libxc)
 integer, parameter :: fail_input=1
!> an iteration reached its limit (a root search, a self-consistent loop,
!> a bound state, a Newton polish)
 integer, parameter :: fail_convergence=2
!> something the construction needs does not exist for these parameters:
!> a bound state, a root in a bracket, a crossing
 integer, parameter :: fail_not_found=3
!> the pseudopotential cannot be constructed as asked: a node in a pseudo
!> wave function, a ghost, a norm that cannot be kept, a parameter out of
!> its range for this atom
 integer, parameter :: fail_construction=4
!> a numerical operation failed: a singular or non-positive matrix, a
!> LAPACK routine's info, a value that is not finite
 integer, parameter :: fail_numerical=5
!> a file cannot be opened, read or written, or its contents do not match
 integer, parameter :: fail_io=6
!> arguments that the calling code should never pass: a bug, not the input
 integer, parameter :: fail_internal=7

 contains

 pure function fail_name(kfail) result(name)

!! the name of a kind of failure, for messages and bindings
!!
!! @param kfail  the kind of failure

!Input variables
 integer, intent(in) :: kfail

!Output variables
 character(len=12) :: name

 select case(kfail)
 case(fail_none)
   name='none'
 case(fail_input)
   name='input'
 case(fail_convergence)
   name='convergence'
 case(fail_not_found)
   name='not_found'
 case(fail_construction)
   name='construction'
 case(fail_numerical)
   name='numerical'
 case(fail_io)
   name='io'
 case(fail_internal)
   name='internal'
 case default
   name='unknown'
 end select
 return
 end function fail_name

 end module failure_m

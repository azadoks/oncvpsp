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
!> ierr, where a routine has one, is a status that its caller interprets
!> (no root in a bracket, a box state, ...), documented with the routine;
!> it is not a failure until the caller makes it one. The statuses that
!> pass between modules, those of the radial solvers and the
!> self-consistent atoms, are named here (ierr_*).

 module failure_m
 implicit none
 private
 public :: fail_none, fail_input, fail_convergence, fail_not_found, &
&          fail_construction, fail_numerical, fail_io, fail_internal, &
&          fail_name
 public :: ierr_box_state, ierr_iterations, ierr_not_converged, &
&          ierr_ps_scf, ierr_ghost
 public :: ierr_lxc_unknown, ierr_lxc_no_vxc, ierr_lxc_family

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

! The statuses (ierr) of the radial solvers and of the self-consistent
! atoms, which their callers interpret; the values are those of 4.0.1,
! which some messages print. 0 is success.

!> no bound state: the box state stands in for it (bs_solve, lschfb,
!> ldiracfb); the self-consistent atoms accept it during their iteration
 integer, parameter :: ierr_box_state=10
!> the energy iteration of a bound-state solver did not converge (bs_solve,
!> lschfb, ldiracfb, lschvkbbe, lschpsbar)
 integer, parameter :: ierr_iterations=60
!> lschvkbb's energy iteration, or sratom's and relatom's self-consistency,
!> did not converge
 integer, parameter :: ierr_not_converged=100
!> psatom's self-consistency did not converge
 integer, parameter :: ierr_ps_scf=101
!> lschvkbb: every energy tried had more states of l below it than the
!> node count expects, so a state it does not expect (a bound ghost of the
!> projectors) lies below emin
 integer, parameter :: ierr_ghost=200

! The statuses of lxc_init and exc_libxc_check, which check_data reports

!> Libxc does not know the functional
 integer, parameter :: ierr_lxc_unknown=1
!> the functional gives no potential
 integer, parameter :: ierr_lxc_no_vxc=2
!> the functional is neither an LDA nor a GGA (meta-GGA, hybrid, ...)
 integer, parameter :: ierr_lxc_family=3

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

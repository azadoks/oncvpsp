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
!> The assertions of the Fortran unit tests. A check writes what it
!> compared to stderr and counts a failure in nfail; check_done ends the
!> test, with exit status 1 after any failure.
 module unit_check_m
 use constants_m, only: dp, qp
 use, intrinsic :: iso_fortran_env, only: stderr => error_unit
 implicit none
 private
 public :: check, check_done, nfail

!> number of failed checks, counted by check and by the tests' own checks
 integer :: nfail = 0

!> check(label, err, tol) with err and tol real(dp) or real(qp): a failure
!> when err > tol or err is NaN; check(label, ok): a failure when not ok
 interface check
   module procedure check_dp, check_qp, check_ok
 end interface check

 contains

 subroutine check_dp(label, err, tol)
!! writes err, and counts a failure when err is not <= tol (NaN fails)
!!
!! @param label  what was checked
!! @param err  its error
!! @param tol  the largest error that passes

!Input variables
 character(len=*), intent(in) :: label
 real(dp), intent(in) :: err, tol

 write(stderr, '(2a,es12.4)') label, ': error ', err
 if(.not. (err <= tol)) then
   write(stderr, '(3a,es12.4,a,es12.4)') 'FAILED: ', label, ': error ', &
&    err, ' > ', tol
   nfail = nfail + 1
 end if
 return
 end subroutine check_dp

 subroutine check_qp(label, err, tol)
!! writes err, and counts a failure when err is not <= tol (NaN fails)
!!
!! @param label  what was checked
!! @param err  its error, in quadruple precision
!! @param tol  the largest error that passes, in quadruple precision

!Input variables
 character(len=*), intent(in) :: label
 real(qp), intent(in) :: err, tol

 write(stderr, '(2a,es12.4)') label, ': error ', real(err, dp)
 if(.not. (err <= tol)) then
   write(stderr, '(3a,es12.4,a,es12.4)') 'FAILED: ', label, ': error ', &
&    real(err, dp), ' > ', real(tol, dp)
   nfail = nfail + 1
 end if
 return
 end subroutine check_qp

 subroutine check_ok(label, ok)
!! counts a failure, and writes it, when ok is false
!!
!! @param label  what was checked
!! @param ok  the result of the check

!Input variables
 character(len=*), intent(in) :: label
 logical, intent(in) :: ok

 if(.not. ok) then
   write(stderr, '(2a)') 'FAILED: ', label
   nfail = nfail + 1
 end if
 return
 end subroutine check_ok

 subroutine check_done(name)
!! ends the checks of test name: stops with exit status 1 if any failed
!!
!! @param name  the test, for the last line

!Input variables
 character(len=*), intent(in) :: name

 if(nfail > 0) then
   write(stderr, '(2a,i0,a)') name, ': ', nfail, ' check(s) failed'
   error stop 1
 end if
 write(stderr, '(2a)') name, ': all checks passed'
 return
 end subroutine check_done

 end module unit_check_m

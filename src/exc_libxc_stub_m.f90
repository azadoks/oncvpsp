!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University, and M. Verstratte, University of Liege
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

 module exc_libxc_m

 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: exc_libxc, exc_libxc_info, have_libxc

!> whether this build has libxc, for check_data
 logical, parameter :: have_libxc=.false.

 contains

 subroutine exc_libxc(iexc,al,rho,vxc,exc,rr,mmax)
!!  dummy substitute for interface to libxc library to allow successful
!!  build without this library; check_data has rejected the libxc
!!  functionals (have_libxc), so this is not reached
!!
!! @param iexc  exchange-correlation function to be used: -(1000*id1+id2) for
!!              the libxc functionals id1 and id2
!! @param al  log of the mesh ratio amesh of the log radial mesh
!! @param rho  charge density
!! @param vxc  exchange-correlation potential
!! @param exc  exchange-correlation energy density
!! @param rr  log radial mesh
!! @param mmax  dimension of log radial grid

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: mmax, iexc
 real(dp), intent(in) :: al
 real(dp), intent(in) :: rho(mmax),rr(mmax)

!Output variables
 real(dp), intent(out) :: vxc(mmax),exc(mmax)

 write(error_unit,'(a,i8,a)') 'exc_libxc_stub: ERROR iexc = ',iexc, &
&     ' requires libxc'
 write(error_unit,'(a)') &
&     'The present oncvpsp executable was built without libxc.'
 write(error_unit,'(a)') 'Program will stop.'

 error stop 1
 end subroutine exc_libxc

 subroutine exc_libxc_info(iexc,iout)
!!  dummy substitute, not reached either (exc_libxc)
!!
!! @param iexc  exchange-correlation function to be used
!! @param iout  unit of the text output

 implicit none

!Input variables
 integer, intent(in) :: iexc,iout

 write(iout,'(a,i8,a)') 'exc_libxc_stub: ERROR iexc = ',iexc, &
&     ' requires libxc'
 error stop 1
 end subroutine exc_libxc_info

 end module exc_libxc_m

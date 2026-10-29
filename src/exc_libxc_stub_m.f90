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
 public :: exc_libxc, exc_libxc_fxc, exc_libxc_check, exc_libxc_info
 public :: have_libxc

!> whether this build has libxc, for check_data
 logical, parameter :: have_libxc=.false.

 contains

 subroutine exc_libxc(iexc,al,rho,vxc,exc,rr,mmax,ifail)
!!  dummy substitute for interface to libxc library to allow successful
!!  build without this library; check_data has rejected the libxc
!!  functionals (have_libxc), so this is not reached
!!
!! @param iexc  exchange-correlation functional: -(1000*id1+id2) for the libxc
!!              functionals id1 and id2
!! @param al  log of the mesh ratio amesh of the log radial mesh
!! @param rho  charge density
!! @param vxc  exchange-correlation potential (Ha)
!! @param exc  exchange-correlation energy per electron (Ha)
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_input, failed, msg_libxc_stub_exc
 use constants_m, only: dp

!Input variables
 integer, intent(in) :: mmax, iexc
 real(dp), intent(in) :: al
 real(dp), intent(in) :: rho(mmax),rr(mmax)

!Output variables
 real(dp), intent(out) :: vxc(mmax),exc(mmax)
 integer, intent(out) :: ifail

 write(error_unit,'(a,i8,a)') 'exc_libxc_stub: ERROR iexc = ',iexc, &
&     ' requires libxc'
 write(error_unit,'(a)') &
&     'The present oncvpsp executable was built without libxc.'
 write(error_unit,'(a)') 'Program will stop.'

 vxc(:)=0.0d0
 exc(:)=0.0d0
 ifail=failed(fail_input,msg_libxc_stub_exc)
 end subroutine exc_libxc

 subroutine exc_libxc_check(iexc,ierr)
!!  dummy substitute, not reached (exc_libxc)
!!
!! @param iexc  exchange-correlation functional
!! @param ierr  0 on success; the non-zero values differ by routine:
!!              the reasons above

 use failure_m, only: ierr_lxc_unknown

!Input variables
 integer, intent(in) :: iexc

!Output variables
 integer, intent(out) :: ierr

 ierr=ierr_lxc_unknown
 write(error_unit,'(/a,i0)') 'exc_libxc_check: ERROR without libxc, iexc = ', &
&     iexc
 end subroutine exc_libxc_check

 subroutine exc_libxc_fxc(iexc,al,rho,drho,dvxc,rr,mmax,have,ifail)
!!  dummy substitute, not reached either (exc_libxc)
!!
!! @param iexc  exchange-correlation functional
!! @param al  log of the mesh ratio amesh of the log radial mesh
!! @param rho  charge density
!! @param drho  change of the charge density
!! @param dvxc  change of the exchange-correlation potential (Ha)
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param have  whether dvxc was evaluated
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_input, failed, msg_libxc_stub_fxc
 use constants_m, only: dp

!Input variables
 integer, intent(in) :: mmax, iexc
 real(dp), intent(in) :: al
 real(dp), intent(in) :: rho(mmax),drho(mmax),rr(mmax)

!Output variables
 real(dp), intent(out) :: dvxc(mmax)
 logical, intent(out) :: have
 integer, intent(out) :: ifail

 dvxc(:)=0.0d0
 have=.false.
 write(error_unit,'(/a,i0)') 'exc_libxc_fxc: ERROR without libxc, iexc = ', &
&     iexc
 ifail=failed(fail_input,msg_libxc_stub_fxc)
 end subroutine exc_libxc_fxc

 subroutine exc_libxc_info(iexc,iout,ifail)
!!  dummy substitute, not reached either (exc_libxc)
!!
!! @param iexc  exchange-correlation functional
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_input, failed, msg_libxc_stub_info

!Input variables
 integer, intent(in) :: iexc,iout

!Output variables
 integer, intent(out) :: ifail

 write(iout,'(a,i8,a)') 'exc_libxc_stub: ERROR iexc = ',iexc, &
&     ' requires libxc'
 ifail=failed(fail_input,msg_libxc_stub_info)
 end subroutine exc_libxc_info

 end module exc_libxc_m

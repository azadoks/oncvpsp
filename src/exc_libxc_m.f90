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

 implicit none
 private
 public :: exc_libxc, exc_libxc_info, have_libxc

!> whether this build has libxc, for check_data
 logical, parameter :: have_libxc=.true.

 contains

 subroutine exc_libxc(iexc,al,rho,vxc,exc,rr,mmax)
!!  interface to libxc library for potential and energy density calculation
!!
!! @param iexc  exchange-correlation function to be used: -(1000*id1+id2) for
!!              the libxc functionals id1 and id2
!! @param al  log of the mesh ratio amesh of the log radial mesh
!! @param rho  charge density
!! @param vxc  exchange-correlation potential
!! @param exc  exchange-correlation energy density
!! @param rr  log radial mesh
!! @param mmax  dimension of log radial grid

 use functionals_m, only: xc_functl_t, xc_functl_init_functl, &
&                         xc_functl_get_vxc, &
&                         xc_functl_end, XC_DERIV_NUMERICAL
 use constants_m, only: dp, pi
 use derivs_m, only: derivs

!Input variables
 integer, intent(in) :: mmax, iexc
 real(dp), intent(in) :: al
 real(dp), intent(in) :: rho(mmax),rr(mmax)

!Output variables
 real(dp), intent(out) :: vxc(mmax),exc(mmax)

!Local variables
 integer :: ii, unit
 integer ::  id(2), ifunc
 type(xc_functl_t) :: functls(2)
 integer :: nspin, irel, deriv_method
 real(dp),allocatable :: ip(:)

 real(dp), allocatable :: vv(:), ee(:)
 real(dp), allocatable :: dpr(:), dppr(:), dlap(:)
 real(dp), allocatable :: tau_dummy(:), vtau_dummy(:)
 real(dp), allocatable :: rho_loc(:)

 allocate(vv(mmax), ee(mmax))
 allocate(dpr(mmax), dppr(mmax), dlap(mmax))
 allocate(tau_dummy(mmax), vtau_dummy(mmax))
 allocate(rho_loc(mmax))

 rho_loc = rho(1:mmax) /4.0d0 / pi

 call derivs(mmax, rho_loc, al, rr, dpr, dppr, dlap)

!
! determine xc to be used from input code (negative by convention)

 id(1) = abs(iexc)/1000
 id(2) = abs(iexc) - id(1)*1000
 nspin = 1
 irel = 0 !XC_NON_RELATIVISTIC
!numerical derivatives (XC_DERIV_ANALYTICAL is the alternative)
 deriv_method = XC_DERIV_NUMERICAL
 call xc_functl_init_functl(functls(1),id(1),3,0.0d0,nspin,deriv_method)
 call xc_functl_init_functl(functls(2),id(2),3,0.0d0,nspin,deriv_method)

! no use of vtau kin en density optional argument
! ip is an ionization potential parameter for certain xc functionals, not
! supported here. Set to 0.5 Ha
 allocate(ip(nspin))
 ip = 0.5d0
 tau_dummy = 0.0d0
 vtau_dummy = 0.0d0
 vxc(1:mmax) = 0.0d0
 exc(1:mmax) = 0.0d0
 do ifunc = 1, 2
   vv = 0.0d0
   ee = 0.0d0
   call xc_functl_get_vxc(functls(ifunc), mmax, al, rr, rho_loc, dpr, dlap, &
&         tau_dummy, ip, vv, ee, vtau_dummy)
   vxc(1:mmax) = vxc(1:mmax) + vv(1:mmax)
   exc(1:mmax) = exc(1:mmax) + ee(1:mmax)
 end do

! convert to the convention of D. R. Hamann's built-in routines (excwig
! etc.).  Why this is needed is not understood: in principle libxc
! outputs the energy per particle, yet dividing by rho once more here
! gives an LDA exc identical to that of the built-in routines.
 do ii=1,mmax
   exc(ii) = exc(ii) / max(rho_loc(ii),1.0d-20)
 end do

! finish
 call xc_functl_end(functls(1))
 call xc_functl_end(functls(2))
 deallocate(ip)
 deallocate(vv, ee)
 deallocate(dpr, dppr, dlap)
 deallocate(tau_dummy, vtau_dummy)
 deallocate(rho_loc)
 return

 end subroutine exc_libxc

 subroutine exc_libxc_info(iexc,iout)
!!  the description of the libxc functionals of iexc, written once per run
!!  (oncvpsp_run)
!!
!! @param iexc  exchange-correlation function to be used: -(1000*id1+id2) for
!!              the libxc functionals id1 and id2
!! @param iout  unit of the text output

 use functionals_m, only: xc_functl_t, xc_functl_init_functl, &
&                         xc_functl_write_info, xc_functl_end, &
&                         XC_DERIV_NUMERICAL

!Input variables
 integer, intent(in) :: iexc,iout

!Local variables
 type(xc_functl_t) :: functls(2)
 integer :: id(2)

 id(1) = abs(iexc)/1000
 id(2) = abs(iexc) - id(1)*1000
 call xc_functl_init_functl(functls(1),id(1),3,0.0d0,1,XC_DERIV_NUMERICAL)
 call xc_functl_init_functl(functls(2),id(2),3,0.0d0,1,XC_DERIV_NUMERICAL)
 call xc_functl_write_info(functls(1), iout)
 call xc_functl_write_info(functls(2), iout)
 write(iout,'(a)') ''
 call xc_functl_end(functls(1))
 call xc_functl_end(functls(2))
 return
 end subroutine exc_libxc_info

 end module exc_libxc_m

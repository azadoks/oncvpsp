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
!> The exchange-correlation energy per electron and potential of a density
!> for the functional iexc: the built-in ones and libxc's (vout,
!> modcore_budget), and the change of the potential for a change of the
!> density (der2exc).
 module xc_dispatch_m
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: exc_dispatch, fxc_dispatch

 contains

 subroutine exc_dispatch(iexc,c_light,al,rho,vxc,exc,rr,nn,who)

!! the xc potential and energy per electron of rho on the log mesh rr
!!
!! @param iexc  exchange-correlation functional (1-4, 6, 7 built in, <0
!!              libxc)
!! @param c_light  speed of light (Hartree atomic units)
!! @param al  log of the mesh ratio
!! @param rho  charge density on rr
!! @param vxc  exchange-correlation potential (Ha) (output)
!! @param exc  exchange-correlation energy per electron (Ha), output
!! @param rr  log radial mesh (bohr)
!! @param nn  its size
!! @param who  the caller, for the error message

 use exc_libxc_m, only: exc_libxc
 use exc_off_pbe_m, only: excggc
 use exc_relx_m, only: exc_relx
 use exchdl_m, only: exchdl
 use excpzca_m, only: excpzca
 use excvwn_m, only: excvwn
 use excwig_m, only: excwig
 use constants_m, only: dp

!Input variables
 integer, intent(in) :: iexc,nn
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: al,rho(nn),rr(nn)
 character(len=*), intent(in) :: who

!Input/output variables (exc_relx adds to the local ones)
 real(dp), intent(inout) :: vxc(nn),exc(nn)

 if(iexc == 1) then
   call excwig(rho,vxc,exc,nn)
 else if(iexc == 2) then
   call exchdl(rho,vxc,exc,nn)
 else if(iexc == 3) then
   call excpzca(rho,vxc,exc,nn)
 else if(iexc == 4) then
   call excggc(rho,vxc,exc,rr,nn)
 else if(iexc == 6 .or. iexc == 7) then
   call excvwn(rho,vxc,exc,nn)
   if(iexc == 7) call exc_relx(rho,vxc,exc,nn,c_light)
 else if(iexc < 0) then
   call exc_libxc(iexc,al,rho,vxc,exc,rr,nn)
 else
   write(error_unit,'(/a,i4)') who//': ERROR bad input iexc =',iexc
   error stop 1
 end if
 return
 end subroutine exc_dispatch

 subroutine fxc_dispatch(iexc,c_light,al,rho,drho,dvxc,rr,mmax,have)

!! the change of exc_dispatch's potential at rho for the change drho of
!! the density, from the functional's analytic second derivatives: for
!! the built-in LDAs (iexc 1-3, 6, 7) and for libxc's LDAs and GGAs that
!! have them. have is false, and dvxc zero, for the others (the built-in
!! PBE, iexc 4), which der2exc differentiates numerically
!!
!! @param iexc  exchange-correlation functional (1-4, 6, 7 built in, <0
!!              libxc)
!! @param c_light  speed of light (Hartree atomic units)
!! @param al  log of the mesh ratio
!! @param rho  charge density on rr
!! @param drho  change of the charge density
!! @param dvxc  change of the exchange-correlation potential (Ha)
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param have  whether dvxc was evaluated

 use exc_libxc_m, only: exc_libxc_fxc
 use exc_relx_m, only: fxc_relx
 use exchdl_m, only: fxchdl
 use excpzca_m, only: fxcpzca
 use excvwn_m, only: fxcvwn
 use excwig_m, only: fxcwig
 use constants_m, only: dp

!Input variables
 integer, intent(in) :: iexc,mmax
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: al,rho(mmax),drho(mmax),rr(mmax)

!Output variables
 real(dp), intent(out) :: dvxc(mmax)
 logical, intent(out) :: have

!Local variables
 real(dp), allocatable :: fxc(:)

 if(iexc<0) then
   call exc_libxc_fxc(iexc,al,rho,drho,dvxc,rr,mmax,have)
   return
 end if

 allocate(fxc(mmax))
 have=.true.
 if(iexc == 1) then
   call fxcwig(rho,fxc,mmax)
 else if(iexc == 2) then
   call fxchdl(rho,fxc,mmax)
 else if(iexc == 3) then
   call fxcpzca(rho,fxc,mmax)
 else if(iexc == 6 .or. iexc == 7) then
   call fxcvwn(rho,fxc,mmax)
   if(iexc == 7) call fxc_relx(rho,fxc,mmax,c_light)
 else
   fxc(:)=0.0d0
   have=.false.
 end if
 dvxc(:)=fxc(:)*drho(:)
 deallocate(fxc)
 return
 end subroutine fxc_dispatch

 end module xc_dispatch_m

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
!> modcore_budget).
 module xc_dispatch_m
 use constants_m, only: dp
 implicit none
 private
 public :: exc_dispatch

 contains

 subroutine exc_dispatch(iexc,al,rho,vxc,exc,rr,nn,who)

!! the xc potential and energy per electron of rho on the log mesh rr
!!
!! @param iexc  exchange-correlation functional (1-4, 6, 7 built in, <0
!!              libxc)
!! @param al  log of the mesh's step
!! @param rho  charge density on rr
!! @param vxc  exchange-correlation potential (output)
!! @param exc  exchange-correlation energy per electron (output)
!! @param rr  log radial mesh
!! @param nn  its size
!! @param who  the caller, for the error message

!Input variables
 integer, intent(in) :: iexc,nn
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
   if(iexc == 7) call exc_relx(rho,vxc,exc,nn)
 else if(iexc < 0) then
   call exc_libxc(iexc,al,rho,vxc,exc,rr,nn)
 else
   write(6,'(/a,i4)') who//': ERROR bad input iexc =',iexc
   error stop 1
 end if
 return
 end subroutine exc_dispatch

 end module xc_dispatch_m

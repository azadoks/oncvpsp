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


!> The pseudopotential files of a run: psp8, UPF and PSML
 module psp_files_m

 implicit none
 private
 public :: write_psp_files, write_psp_file

 contains

 subroutine write_psp_files(inp,psp,iout,ifail)
!! writes the files of the pseudopotential that the input's psfile asks
!! for ("psp8", "upf", "psml", "both" for all three, "none"), in the
!! output on iout
!!
!! @param inp  the input values (input_m), as oncvpsp_run completed them
!! @param psp  the results of oncvpsp_run
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use input_m, only: input_t
 use results_m, only: results_t

!Input variables
 type(input_t), intent(in) :: inp
 type(results_t), intent(in) :: psp
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ifail

!Local variables
 character(len=4), parameter :: kinds(3)=['psp8','upf ','psml']
 integer :: ii

 ifail=0
 do ii=1,3
   if(trim(inp%psfile)==trim(kinds(ii)) .or. trim(inp%psfile)=='both') then
     call write_psp_file(kinds(ii),inp,psp,iout,ifail)
     if(ifail/=0) return
   end if
 end do

 return
 end subroutine write_psp_files

 subroutine write_psp_file(kind,inp,psp,iunit,ifail)
!! writes the pseudopotential as one file format on iunit, its input echo
!! that of inp (psfile as given)
!!
!! @param kind  the format: "psp8", "upf" or "psml"
!! @param inp  the input values (input_m), as oncvpsp_run completed them
!! @param psp  the results of oncvpsp_run
!! @param iunit  unit of the file
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written
 use constants_m, only: dp
 use input_m, only: input_t
 use results_m, only: results_t
 use linout_m, only: linout
 use upfout_m, only: upfout
 use psmlout_m, only: psmlout
 use failure_m, only: fail_internal, failed, msg_psp_file_kind
 use, intrinsic :: iso_fortran_env, only: error_unit

!Input variables
 character(len=*), intent(in) :: kind
 type(input_t), intent(in) :: inp
 type(results_t), intent(in) :: psp
 integer, intent(in) :: iunit

!Output variables
 integer, intent(out) :: ifail

!Local variables
!the spin-orbit strength of linout: 1, not in use yet
 real(dp), parameter :: soscale=1.0d0

 ifail=0

 select case(trim(kind))
 case('psp8')
   call linout(inp,psp%pot%c_light,psp%pot%vkbpw,psp%pot%evkb,psp%pot%vsrpw, &
&              psp%pot%esr,psp%pot%vsopw,psp%pot%eso,psp%pot%rr,psp%pot%vpuns, &
&              psp%pot%vlocpw,psp%pot%rho,psp%pot%corepw,psp%rhotae,psp%rhoc, &
&              psp%pot%zion,psp%pot%mmax,psp%nrl,soscale,psp%echo%rc_in, &
&              psp%echo%fcfact_in,psp%echo%fcgrid,psp%echo%rcgrid, &
&              psp%echo%log_mesh,psp%pot%nkap,iunit,ifail)
 case('upf')
   call upfout(inp,psp%pot%c_light,psp%pot%vkbpw,psp%pot%evkb,psp%pot%rr, &
&              psp%pot%vpuns,psp%pot%vlocpw,psp%pot%rho,psp%pot%corepw, &
&              psp%pot%zion,psp%pot%mmax,psp%nrl,psp%epstot, &
&              psp%echo%rc_in,psp%echo%fcfact_in,psp%echo%fcgrid, &
&              psp%echo%rcgrid, &
&              psp%uupsa,psp%ea,psp%echo%log_mesh,psp%pot%nkap,iunit,ifail)
 case('psml')
   call psmlout(inp,psp%pot%vkb,psp%pot%evkb,psp%pot%rr,psp%vpsml, &
&               psp%pot%rho,psp%pot%rhomod,psp%rcore,psp%pot%zion, &
&               psp%pot%mmax,psp%nrl,psp%epstot,psp%echo%rc_in, &
&               psp%echo%fcfact_in, &
&               psp%echo%fcgrid,psp%echo%rcgrid,psp%echo%log_mesh, &
&               psp%pot%irel,psp%pot%c_light,psp%pot%nkap,iunit,ifail)
 case default
   write(error_unit,'(/2a)') 'write_psp_file: ERROR no file format ',kind
   ifail=failed(fail_internal,msg_psp_file_kind)
 end select

 return
 end subroutine write_psp_file

 end module psp_files_m

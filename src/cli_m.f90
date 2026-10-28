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
!> The command lines of oncvpsp_ld.x and oncvpsp_confine.x: an option's
!> value, comma-separated lists, and the options both share (--radii,
!> --channels, --energies), each an error on stdout with a nonzero exit
!> status that names the program, prog.
 module cli_m
 implicit none
 private
 public :: cli_next_arg, cli_count_list, cli_radii, cli_channels, &
&          cli_energies

 contains

 subroutine cli_next_arg(prog,iarg,val)
!! the value of the option at iarg, which it passes
!!
!! @param prog  the program name, for the error message
!! @param iarg  index of the option among the command arguments; on
!!              return that of its value
!! @param val  the value, the argument after the option

!Input variables
 character(len=*), intent(in) :: prog

!Input/output variables
 integer, intent(inout) :: iarg

!Output variables
 character(len=1024), intent(out) :: val
 if(iarg+1>command_argument_count()) then
   call get_command_argument(iarg,val)
   write(6,'(4a)') prog,': ERROR ',trim(val),' requires a value'
   error stop 1
 end if
 iarg=iarg+1
 call get_command_argument(iarg,val)
 end subroutine cli_next_arg

 integer function cli_count_list(str)
!! the number of comma-separated items of str
!!
!! @param str  the comma-separated list

!Input variables
 character(len=*), intent(in) :: str

!Local variables
 integer :: kk
 cli_count_list=1
 do kk=1,len_trim(str)
   if(str(kk:kk)==',') cli_count_list=cli_count_list+1
 end do
 end function cli_count_list

 subroutine cli_radii(prog,arg,mxq,nq,rq)
!! --radii: nq (at most mxq) positive radii rq
!!
!! @param prog  the program name, for the error messages
!! @param arg  the value of --radii, radii separated by commas
!! @param mxq  the largest number of radii
!! @param nq  the number of radii read
!! @param rq  the radii (bohr), in rq(1:nq)
 use constants_m, only: dp

!Input variables
 character(len=*), intent(in) :: prog,arg
 integer, intent(in) :: mxq

!Output variables
 integer, intent(out) :: nq
 real(dp), intent(out) :: rq(mxq)

!Local variables
 integer :: ios
 nq=cli_count_list(arg)
 if(nq>mxq) then
   write(6,'(2a,i0,a)') prog,': ERROR more than ',mxq,' radii'
   error stop 1
 end if
 read(arg,*,iostat=ios) rq(1:nq)
 if(ios/=0 .or. nq<1 .or. any(rq(1:nq)<=0.0d0)) then
   write(6,'(3a)') prog,': ERROR --radii must be positive ', &
&       'numbers separated by commas'
   error stop 1
 end if
 end subroutine cli_radii

 subroutine cli_channels(prog,arg,chan)
!! --channels: the values of l from 0 to 3 of the output
!!
!! @param prog  the program name, for the error message
!! @param arg  the value of --channels, values of l separated by commas
!! @param chan  .true. for each l listed

!Input variables
 character(len=*), intent(in) :: prog,arg

!Output variables
 logical, intent(out) :: chan(0:3)

!Local variables
 integer :: lsel(4),nl,ios
 nl=min(cli_count_list(arg),size(lsel)+1)
 lsel(:)=-1
 ios=1
 if(nl<=size(lsel)) read(arg,*,iostat=ios) lsel(1:nl)
 if(ios/=0) lsel(1)=-1
 if(any(lsel(1:min(nl,size(lsel)))<0) &
&   .or. any(lsel(1:min(nl,size(lsel)))>3)) then
   write(6,'(3a)') prog,': ERROR --channels must be values ', &
&       'of l from 0 to 3 separated by commas'
   error stop 1
 end if
 chan(:)=.false.
 chan(lsel(1:nl))=.true.
 end subroutine cli_channels

 subroutine cli_energies(prog,arg,ee)
!! --energies: emin, emax and the step de
!!
!! @param prog  the program name, for the error message
!! @param arg  the value of --energies, emin,emax,de
!! @param ee  emin, emax and de (Ha)
 use constants_m, only: dp

!Input variables
 character(len=*), intent(in) :: prog,arg

!Output variables
 real(dp), intent(out) :: ee(3)

!Local variables
 integer :: ios
 read(arg,*,iostat=ios) ee(1:3)
 if(ios/=0 .or. cli_count_list(arg)/=3 .or. ee(3)<=0.0d0 &
&   .or. ee(2)<=ee(1)) then
   write(6,'(3a)') prog,': ERROR --energies must be ', &
&       'emin,emax,de with emin < emax and de > 0'
   error stop 1
 end if
 end subroutine cli_energies

 end module cli_m

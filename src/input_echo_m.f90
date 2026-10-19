!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
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
! sections of the input echo whose form depends on the input
! optional inputs are written only when they are set, so inputs without
! them echo as they always have. The *_lines routines return the lines,
! for the PSML echo; write_icmod and write_epsh write them to unit 6.

 module input_echo_m

 use constants_m, only: dp
 implicit none
 private
 public :: icmod_lines, epsh_lines, write_icmod, write_epsh

!length and number of the lines returned by the *_lines routines
 integer, parameter, public :: echo_len=100, echo_max=5

 contains

! MODEL CORE CHARGE section
! the icmod=4 search grid is written only when it differs from the default
 subroutine icmod_lines(icmod,fcfact,rcfact,fcgrid,rcgrid,text,ntext)

!icmod  model core charge option
!fcfact  prefactor for model amplitude
!rcfact  prefactor for model scale
!fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!text  the lines of the section
!ntext  number of lines

!Input variables
 integer, intent(in) :: icmod
 real(dp), intent(in) :: fcfact,rcfact
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)

!Output variables
 character(len=echo_len), intent(out) :: text(echo_max)
 integer, intent(out) :: ntext

!Local variables
!defaults, as set by the input readers
 real(dp), parameter :: fcgrid_def(3)=[1.5d0,6.0d0,0.5d0]
 real(dp), parameter :: rcgrid_def(3)=[1.0d0,1.9d0,0.1d0]

 text(1)='#'
 text(2)='# MODEL CORE CHARGE'
 if(icmod==4 .and. (any(fcgrid/=fcgrid_def) .or. any(rcgrid/=rcgrid_def))) then
   text(3)='# icmod, fcfact, rcfact,'
   text(4)='#   fcfact_min, fcfact_max, fcfact_step, rcfact_min, rcfact_max, rcfact_step'
   write(text(5),'(i5,8f10.5)') icmod,fcfact,rcfact,fcgrid(:),rcgrid(:)
   ntext=5
 else
   text(3)='# icmod, fcfact, rcfact'
!icmod=6 rcfact is an error tolerance
   if(icmod==6) then
     write(text(4),'(i5,f10.5,es12.4)') icmod,fcfact,rcfact
   else
     write(text(4),'(i5,2f10.5)') icmod,fcfact,rcfact
   end if
   ntext=4
 end if

 return
 end subroutine icmod_lines

! LOG DERIVATIVE ANALYSIS section
! rxpsh is written only when the input sets it (the readers set -1 otherwise)
 subroutine epsh_lines(epsh1,epsh2,depsh,rxpsh,text,ntext)

!epsh1, epsh2  energy range of the log derivative analysis
!depsh  energy step
!rxpsh  radius of the analysis, negative for the default
!text  the lines of the section
!ntext  number of lines

!Input variables
 real(dp), intent(in) :: epsh1,epsh2,depsh,rxpsh

!Output variables
 character(len=echo_len), intent(out) :: text(echo_max)
 integer, intent(out) :: ntext

 text(1)='#'
 text(2)='# LOG DERIVATIVE ANALYSIS'
 if(rxpsh>=0.0d0) then
   text(3)='# epsh1, epsh2, depsh, rxpsh'
   write(text(4),'(3f8.2,f10.5)') epsh1,epsh2,depsh,rxpsh
 else
   text(3)='# epsh1, epsh2, depsh'
   write(text(4),'(3f8.2)') epsh1,epsh2,depsh
 end if
 ntext=4

 return
 end subroutine epsh_lines

! writes the MODEL CORE CHARGE section to unit 6
 subroutine write_icmod(icmod,fcfact,rcfact,fcgrid,rcgrid)

!arguments as in icmod_lines

!Input variables
 integer, intent(in) :: icmod
 real(dp), intent(in) :: fcfact,rcfact
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)

!Local variables
 character(len=echo_len) :: text(echo_max)
 integer :: ii,ntext

 call icmod_lines(icmod,fcfact,rcfact,fcgrid,rcgrid,text,ntext)
 write(6,'(a)') (trim(text(ii)),ii=1,ntext)

 return
 end subroutine write_icmod

! writes the LOG DERIVATIVE ANALYSIS section to unit 6
 subroutine write_epsh(epsh1,epsh2,depsh,rxpsh)

!arguments as in epsh_lines

!Input variables
 real(dp), intent(in) :: epsh1,epsh2,depsh,rxpsh

!Local variables
 character(len=echo_len) :: text(echo_max)
 integer :: ii,ntext

 call epsh_lines(epsh1,epsh2,depsh,rxpsh,text,ntext)
 write(6,'(a)') (trim(text(ii)),ii=1,ntext)

 return
 end subroutine write_epsh

 end module input_echo_m

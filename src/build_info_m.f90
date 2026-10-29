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
! Git commit, state of the source relative to it, and build time of the
! executables, generated at build time by cmake/BuildInfo.cmake into
! build_info.inc
!
 module build_info_m
 implicit none
 private
 public :: write_build_info, build_info_text

 include 'build_info.inc'

 contains

 subroutine write_build_info(iout,program,treatment)
!! writes the one-line build description
!!
!! @param iout  unit of the text output
!! @param program  name of the executable
!! @param treatment  relativistic treatment (non-, scalar- or
!!                   fully-relativistic)

!Input variables
 integer, intent(in) :: iout
 character(len=*), intent(in) :: program,treatment

 write(iout,'(10a)') 'build: ',program,' (',treatment,'), commit ', &
&      build_commit,', source ',build_source,', compiled ',build_time
 return

 end subroutine write_build_info

 function build_info_text() result(text)
!! the build's commit, the state of the source and the build time, a line
!! each: 'commit ', 'source ', 'compiled ' and the value

!Output variables
 character(len=:), allocatable :: text

 text='commit '//build_commit//new_line('a')//'source '//build_source// &
&     new_line('a')//'compiled '//build_time//new_line('a')
 return

 end function build_info_text

 end module build_info_m

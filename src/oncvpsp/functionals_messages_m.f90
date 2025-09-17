!! Taken from Octopus (2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch)
!! Adapted to oncvpsp by A. Castaneda M. (2019)
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
module xc_messages
    implicit none
    character(1200) :: message(2)
contains

subroutine messages_info(lmess, iunit)
    integer, intent(in) :: lmess, iunit
    integer :: i
    do i=1,lmess
        write(iunit,'("LXC ",A)') trim(message(i))
    end do
end subroutine messages_info

subroutine messages_fatal(lmess)
    integer, intent(in) :: lmess
    integer :: i
    do i=1,lmess
        write(*,'("Fatal ",A)') trim(message(i))
    end do
    stop
end subroutine messages_fatal

subroutine messages_warning(lmess)
    integer, intent(in) :: lmess
    integer :: i
    do i=1,lmess
        write(*,'("Warning ",A)') trim(message(i))
    end do
end subroutine messages_warning

subroutine messages_input_error(str1,str2)
    character(*), intent(in) :: str1,str2
    write(*,'("Input error ",A,A)') trim(str1),trim(str2)
end subroutine messages_input_error

subroutine messages_experimental(str)
    character(*), intent(in) :: str
    write(*,'("Experimental xc ",A)') trim(str)
end subroutine messages_experimental

subroutine messages_not_implemented(str)
    character(*), intent(in) :: str
    write(*,'("Not implemented ",A)') trim(str)
end subroutine messages_not_implemented

end module xc_messages

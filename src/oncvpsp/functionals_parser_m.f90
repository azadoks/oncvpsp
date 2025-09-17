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
module parser
    implicit none
    interface parse_variable
        module procedure real_parser
        module procedure logical_parser
    end interface parse_variable
contains

subroutine real_parser(varval, var)
    real(8), intent(in) :: varval
    real(8), intent(out) :: var
    var = varval
end subroutine real_parser

subroutine logical_parser(varval, var)
    logical, intent(in) :: varval
    logical, intent(out) :: var
    var = varval
end subroutine logical_parser

end module parser

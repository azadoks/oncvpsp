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
! writes the MODEL CORE CHARGE section of the input echo
! the icmod=4 search grid is written only when it differs from the default,
! so inputs without an explicit grid echo as they always have

 subroutine write_icmod(icmod,fcfact,rcfact,fcgrid,rcgrid)

!icmod  model core charge option
!fcfact  prefactor for model amplitude
!rcfact  prefactor for model scale
!fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 integer, intent(in) :: icmod
 real(dp), intent(in) :: fcfact,rcfact
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)

!Local variables
!defaults, as set by the input readers
 real(dp), parameter :: fcgrid_def(3)=[1.5d0,6.0d0,0.5d0]
 real(dp), parameter :: rcgrid_def(3)=[1.0d0,1.9d0,0.1d0]

 if(icmod==4 .and. (any(fcgrid/=fcgrid_def) .or. any(rcgrid/=rcgrid_def))) then
   write(6,'(a/a/a/a)') '#','# MODEL CORE CHARGE', &
&        '# icmod, fcfact, rcfact,', &
&        '#   fcfact_min, fcfact_max, fcfact_step, rcfact_min, rcfact_max, rcfact_step'
   write(6,'(i5,8f10.5)') icmod,fcfact,rcfact,fcgrid(:),rcgrid(:)
 else
   write(6,'(a/a/a)') '#','# MODEL CORE CHARGE', &
&        '# icmod, fcfact, rcfact'
!icmod=6 rcfact is an error tolerance
   if(icmod==6) then
     write(6,'(i5,f10.5,es12.4)') icmod,fcfact,rcfact
   else
     write(6,'(i5,2f10.5)') icmod,fcfact,rcfact
   end if
 end if

 return
 end subroutine write_icmod

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
!> The inward Adams integration of the second-order radial equation
!> u'' = al u' + cf u in the mesh index, shared by the bound-state solvers
!> of the pseudopotential (lschvkbb, lschpsbar)

 module inward2_m

 use constants_m, only: dp
 implicit none
 private
 public :: inward2

 contains

 subroutine inward2(al,cf,mmax,ifrom,ito,uu,up,upp)

!! inward Adams integration from ifrom down to ito, two corrector passes;
!! uu, up and upp are given on ifrom..ifrom+3. upp at ito is that of the
!! corrector before the last
!!
!! @param al  log of the mesh ratio
!! @param cf  coefficient of uu in upp
!! @param mmax  size of the arrays
!! @param ifrom  first point, given
!! @param ito  last point integrated
!! @param uu  radial wave function (*rr)
!! @param up  d(uu)/dn
!! @param upp  d(up)/dn

!Input variables
 integer, intent(in) :: mmax,ifrom,ito
 real(dp), intent(in) :: al,cf(mmax)

!Input/output variables
 real(dp), intent(inout) :: uu(mmax),up(mmax),upp(mmax)

!Local variables
 integer :: ii,it

 do ii=ifrom,ito+1,-1
   uu(ii-1)=uu(ii)+aei(up,ii)
   up(ii-1)=up(ii)+aei(upp,ii)
   do it=1,2
     upp(ii-1)=al*up(ii-1)+cf(ii-1)*uu(ii-1)
     up(ii-1)=up(ii)+aii(upp,ii)
     uu(ii-1)=uu(ii)+aii(up,ii)
   end do
 end do

 return
 end subroutine inward2

 include 'adams.inc'

 end module inward2_m

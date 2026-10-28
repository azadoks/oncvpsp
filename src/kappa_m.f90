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
!> The two spin-orbit components of a channel l when nkap=2 (fully
!> relativistic): how many there are, their kappa, and the share of an
!> occupation each takes, 2j+1 of 2(2l+1).
 module kappa_m
 implicit none
 private
 public :: kap_count, kap_of, j_share

 contains

 integer function kap_count(nkap,ll)
!! the components of l: nkap, and one for l=0
 integer, intent(in) :: nkap,ll
 kap_count=nkap
 if(ll==0) kap_count=1
 end function kap_count

 integer function kap_of(ll,ikap)
!! kappa of component ikap of l: -(l+1) for j=l+1/2 (ikap 1), l for
!! j=l-1/2 (ikap 2)
 integer, intent(in) :: ll,ikap
 if(ikap==1) then
   kap_of=-(ll+1)
 else
   kap_of=ll
 end if
 end function kap_of

 real(dp) function j_share(ll,ikap,ff)
!! the share of the occupation ff of l in component ikap
 use constants_m, only: dp

 integer, intent(in) :: ll,ikap
 real(dp), intent(in) :: ff
 if(ikap==1) then
   j_share=dble(ll+1)*ff/dble(2*ll+1)
 else
   j_share=dble(ll)*ff/dble(2*ll+1)
 end if
 end function j_share

 end module kappa_m

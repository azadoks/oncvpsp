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
 subroutine lschvkbs(ll,ivkb,ee,rr,vloc,vlocpw,vkbpw,evkb,uu,up,mmax,mch)

!! Normalized scattering state for pseudopotential, fully non-local unless
!! ivkb=0, which should only happen if ll=lloc
!!
!! @param ll  angular-momentum quantum number
!! @param ivkb  = 0, 1 or 2 VKB proectors to be used
!! @param ee  scattering state energy
!! @param rr  log radial mesh
!! @param vloc  local pseudopotential on the mesh
!! @param vlocpw  its pieces
!! @param vkbpw  VKB projectors as functions of pieces
!! @param evkb  coefficients of VKB projectors
!! @param uu  output radial wave function (*rr)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of log grid
!! @param mch  index of radius to which uu is computed

 use constants_m, only: dp
 use pieces_m, only: piecewise_t
 use quad_m, only: log_norm
 implicit none

!Input variables
 integer, intent(in) :: mmax,mch
 real(dp), intent(in) :: rr(mmax),vloc(mmax),evkb(ivkb)
 real(dp), intent(in) :: ee
 integer, intent(in) :: ivkb,ll
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(ivkb)

!Output variables
 real(dp), intent(out) :: uu(mmax),up(mmax)


!Local variables
 real(dp) :: amesh,al
 real(dp) :: cn,sn
 integer :: ii


 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

! null arrays to remove leftover garbage

 uu(:)=0.0d0
 up(:)=0.0d0

 call vkboutwf(ll,ivkb,ee,vkbpw,evkb,rr,vloc,vlocpw,uu,up,mmax,mch)

!perform normalization sum

 sn=log_norm(rr,uu(1:mch)**2,mch,ll,dble(ll+1))

!normalize u

 cn=1.0d0/sqrt(sn)

 do ii=1,mch
   up(ii)=cn*up(ii)
   uu(ii)=cn*uu(ii)
 end do
 do ii=mch+1,mmax
   uu(ii)=0.0d0
 end do

 return
 end  subroutine lschvkbs

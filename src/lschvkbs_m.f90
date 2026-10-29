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

 module lschvkbs_m

 implicit none
 private
 public :: lschvkbs

 contains

 subroutine lschvkbs(ll,ivkb,ee,rr,vloc,vlocpw,vkbpw,evkb,uu,up,mmax,adams, &
&      mlast, &
&      iout,ifail)

!! Normalized scattering state for pseudopotential, fully non-local unless
!! ivkb=0, which should only happen if ll=lloc
!!
!! @param ll  angular momentum
!! @param ivkb  = 0, 1 or 2 VKB proectors to be used
!! @param ee  scattering state energy
!! @param rr  log radial mesh (bohr)
!! @param vloc  local pseudopotential on the mesh
!! @param vlocpw  the screened local potential as pieces
!! @param vkbpw  the projectors as pieces
!! @param evkb  coefficients of the projectors
!! @param uu  output radial wave function (*rr)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param mlast  index of radius to which uu is computed
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: dp
 use pieces_m, only: piecewise_t
 use quad_m, only: quad_log_norm
 use vkboutwf_m, only: vkboutwf
 use logmesh_m, only: lm_step
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: mmax,mlast
 type(adams_t), intent(in) :: adams
 integer, intent(in) :: ivkb,ll
 real(dp), intent(in) :: rr(mmax),vloc(mmax),evkb(ivkb)
 real(dp), intent(in) :: ee
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(ivkb)
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: uu(mmax),up(mmax)
 integer, intent(out) :: ifail


!Local variables
 real(dp) :: amesh,al
 real(dp) :: cn,sn
 integer :: ii


 ifail=0
 al=lm_step(rr,mmax)
 amesh = exp(al)

! null arrays to remove leftover garbage

 uu(:)=0.0d0
 up(:)=0.0d0

 call vkboutwf(ll,ivkb,ee,vkbpw,evkb,rr,vloc,vlocpw,uu,up,mmax,adams,mlast, &
&              iout,ifail)
 if(ifail/=0) return

!perform normalization sum

 sn=quad_log_norm(rr,uu(1:mlast)**2,mlast,ll,dble(ll+1))

!normalize u

 cn=1.0d0/sqrt(sn)

 do ii=1,mlast
   up(ii)=cn*up(ii)
   uu(ii)=cn*uu(ii)
 end do
 do ii=mlast+1,mmax
   uu(ii)=0.0d0
 end do

 return
 end  subroutine lschvkbs

 end module lschvkbs_m

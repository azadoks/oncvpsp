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
 subroutine lschkb(nn,ll,ierr,ee,vkb,rr,vv,uu,up,mmax,mch)

!! outward integration of the inhomogeneous radial Schroedinger equation
!! on a logarithmic mesh with local potential and one proector term
!!
!! @param nn  principal quantum number (not used)
!! @param ll  angular-momentum quantum number
!! @param ierr  non-zero return if error
!! @param ee  bound-state energy, input guess and output calculated value
!! @param vkb  Vanderbilt-Kleinman-bylander projector
!! @param rr  log radial mesh
!! @param vv  local pseudopotential
!! @param uu  output radial wave function (*rr)
!! @param up  d(uu)/dr
!! @param mmax  size of log grid
!! @param mch  matching mesh point for inward-outward integrations

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: mmax,mch
 real(dp), intent(in) :: rr(mmax),vv(mmax),vkb(mmax)
 integer, intent(in) :: nn,ll
 real(dp), intent(in) :: ee

!Output variables
 real(dp), intent(out) :: uu(mmax),up(mmax)
 integer, intent(out) :: ierr


!Local variables
 real(dp) :: amesh,al
 real(dp) :: als, cn
 real(dp) :: sls, uout, upout
 real(dp) :: akb,ckb
 integer :: ii, it

 real(dp), allocatable :: upp(:),cf(:)

 allocate(upp(mmax),cf(mmax))


 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

 ierr = 0

 sls=ll*(ll+1)

! null arrays to remove leftover garbage

 uu(:)=0.0d0
 up(:)=0.0d0
 upp(:)=0.0d0

 als=al**2

! coefficient array for uu in differential eq.
 do ii=1,mmax
   cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
 end do

! start wavefunction with series based on projector

 ckb = 2.0d0*vkb(1)/rr(1)**(ll+1)
 akb = ckb/(6.0d0+4.0d0*ll)
 do ii=1,4
   uu(ii)=akb*rr(ii)**(ll+3)
   up(ii)= al*(ll+3)*uu(ii)
   upp(ii)= als*(ll+3)**2*uu(ii)
 end do

! outward integration using predictor once, corrector
! twice

 do ii=4,mch-1
   uu(ii+1)=uu(ii)+aeo(up,ii)
   up(ii+1)=up(ii)+aeo(upp,ii)
   do it=1,2
     upp(ii+1)=al*up(ii+1)+cf(ii+1)*uu(ii+1) &
&            + 2.0d0*als*vkb(ii+1)*rr(ii+1)**2
     up(ii+1)=up(ii)+aio(upp,ii)
     uu(ii+1)=uu(ii)+aio(up,ii)
   end do
 end do

 uout=uu(mch)
 upout=up(mch)

 deallocate(upp,cf)

 return
 contains

 include 'adams.inc'

 end  subroutine lschkb

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
 subroutine fphsft(ll,kap,epsh2,depsh,pshf,rr,vv,zz,mmax,mch,npsh,irel)

!! computes full potential scattering log derivatives
!! returns atan(rr(mch) * du/dr / u) which is sort-of like a phase shift
!! and easier to compare in plots than the log derivatives themselves
!! Pauli-type scalar-relativistic or non-relativistic calculation, or the
!! Dirac equation with the log derivatives of the large component only
!!
!! @param ll  angular momentum
!! @param kap  =l, -(l+1) for j=l -/+ 1/2 (irel=2 only)
!! @param epsh2  upper limit of energy scan
!! @param depsh  increment of scan
!! @param pshf  log derivatives "angles", as above
!! @param rr  radial log grid
!! @param vv  all-electron potential
!! @param zz  atomic number
!! @param mmax  dimension of rr, etc.
!! @param mch  index of radius for log der test
!! @param npsh  number of energy points in scan
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic

 use constants_m, only: dp, pi
 implicit none

!Input variables
 integer, intent(in) :: ll,kap,mmax,npsh,mch,irel
 real(dp), intent(in) :: rr(mmax),vv(mmax)
 real(dp), intent(in) :: depsh,epsh2,zz

!Output variables
 real(dp), intent(out) :: pshf(npsh)

!Local variables
 real(dp), parameter :: pi2=2.0d0*pi
 real(dp) :: al,epsh,phi,phip,pshoff
 integer :: ii,ierr,nn

 real(dp), allocatable :: uu(:,:),up(:,:)

 allocate(uu(mmax,2),up(mmax,2))

 al = 0.01d0 * log(rr(101) / rr(1))

 pshoff = 0.0d0

 do ii = 1,npsh
   epsh = epsh2-(ii-1)*depsh

   call scatstate(irel,nn,ll,kap,ierr,epsh,rr,zz,vv,uu,up,mmax,mch)

   phi = uu(mch,1)/rr(mch)
   phip = (up(mch,1)-al*uu(mch,1))/(al*rr(mch)**2)
   pshf(ii) = atan2(rr(mch)*phip,phi) + pshoff
   if(ii>1) then
     if(pshf(ii)<pshf(ii-1)) then
       pshoff = pshoff+pi2
       pshf(ii) = pshf(ii)+pi2
     end if
   end if
 end do

 deallocate(uu,up)
 return
 end subroutine fphsft

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
! calculates all-electron wave function overlap intergrals

 subroutine fpovlp(gg,hh,nn,ll,zz,ss,rr,irel)

!! product of all-elecctron scalar-relativistic wave functions gg*hh
!! goes like rr**(2*gamma) as rr -> 0
!! integral on usual log mesh from rr=0 to rr(nn)
!!
!! @param gg  wave function one
!! @param hh  wave function two
!! @param nn  outer integration limit
!! @param ll  angular momentum
!! @param zz  nuclear charge
!! @param ss  overlap output
!! @param rr  log radial mesh
!! @param irel  0 non-relativistic, 1 scalar-, 2 fully-relativistic (the
!!              scalar-relativistic overlap is used for both)

 use constants_m, only: dp
 use relativity_m, only: c_light
 implicit none

!Input variables
 integer, intent(in) :: nn,ll
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: gg(nn),hh(nn),rr(nn)
 integer, intent(in) :: irel

!Output variables
 real(dp), intent(out) :: ss

!Local variables
 real(dp) :: r0,amesh,al,fss,gamma
 integer :: ii

 al = 0.01d0 * log(rr(101)/rr(1))
 amesh = exp(al)

 if(irel/=0) then
   fss=(1.0d0/c_light)**2
 else
   fss=1.0d-12
 end if

 if(ll==0) gamma=sqrt(1.0d0-fss*zz**2)
 if(ll>0) gamma=(ll*sqrt(ll**2-fss*zz**2) + &
& (ll+1)*sqrt((ll+1)**2-fss*zz**2))/(2*ll+1)

 r0=rr(1)/sqrt(amesh)
 ss=r0**(2.0d0*gamma+1.0d0)/(2.d0*gamma+1.0d0)
 ss=(ss*gg(1)*hh(1))/rr(1)**(2.d0*gamma)

 do ii = 1, nn - 3
   ss =  ss + al*gg(ii)*hh(ii)*rr(ii)
 end do

 ss=ss + al*(23.d0*rr(nn-2)*gg(nn-2)*hh(nn-2) &
&        + 28.d0*rr(nn-1)*gg(nn-1)*hh(nn-1) &
&        +  9.d0*rr(nn  )*gg(nn  )*hh(nn  ))/24.d0


 return
 end subroutine fpovlp

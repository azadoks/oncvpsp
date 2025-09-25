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

!> calculates all-electron wave function overlap intergrals
!> product of all-electron scalar-relativistic wave functions gg*hh
!> goes like rr**(2*gamma) as rr -> 0
!> integral on usual log mesh from rr=0 to rr(nn)
subroutine fpovlp(gg,hh,nn,ll,zz,ss,rr,srel)
   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables
   !> nuclear charge
   real(dp), intent(in) :: zz
   !> wave function one
   real(dp), intent(in) :: gg(nn)
   !> wave function two
   real(dp), intent(in) :: hh(nn)
   !> log radial mesh
   real(dp), intent(in) :: rr(nn)
   !> outer integration limit
   integer, intent(in) :: nn
   !> angular momentum
   integer, intent(in) :: ll
   !> .true. for scalar-relativistic, .false. for non-relativistic
   logical, intent(in) :: srel

!Output variable
   !> overlap output
   real(dp), intent(out) :: ss

!Local variables
   real(dp) :: r0,amesh,al,fss,gamma
   integer :: ii

   al = 0.01d0 * dlog(rr(101)/rr(1))
   amesh = exp(al)

   if(srel) then
      fss=(1.0d0/137.036d0)**2
   else
      fss=1.0d-12
   end if

   if(ll==0) gamma=dsqrt(1.0d0-fss*zz**2)
   if(ll>0) gamma=(ll*dsqrt(ll**2-fss*zz**2) + &
   & (ll+1)*dsqrt((ll+1)**2-fss*zz**2))/(2*ll+1)

   r0=rr(1)/dsqrt(amesh)
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

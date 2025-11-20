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
module renormalize_m
   implicit none
   private
   public :: renorm_r
contains
!> renormalize Dirac wave function so that large component is normalized
subroutine renorm_r(uu, rr, ll, kap, zz, mmax, cnorm)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> mmax  size of radial grid
   integer, intent(in) :: mmax
   !> rr  log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> zz  atomic number
   real(dp), intent(in) :: zz
   !> ll  angular momentum
   integer, intent(in) :: ll
   !> kap  Dirac kappa
   integer, intent(in) :: kap


   !Output variable
   !> cnorm  renormalization coefficient
   real(dp), intent(out) :: cnorm

   !> Input/Output variable
   !> uu  Dirac wave function
   real(dp), intent(in out) :: uu(mmax, 2)

   !Local variables
   real(dp) :: sn
   real(dp) :: al, amesh
   real(dp) :: r0, gam, cc, cci
   integer :: ii, nin

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = exp(al)

   cc = 137.036d0
   cci = 1.0d0/cc

   gam = sqrt(kap**2 - (zz*cci)**2)

   do ii = mmax, 1, -1
      if (dabs(uu(ii, 1)) > 0.0d0) then
         nin = ii
         exit
      end if
   end do

   r0 = rr(1)/dsqrt(amesh)
   sn = r0**(2.0d0*gam + 1.0d0)/(2.d0*gam + 1.0d0)
   sn = sn*(uu(1, 1)**2)/rr(1)**(2.d0*gam)

   do ii = 1, nin - 3
      sn = sn + al*rr(ii)*uu(ii, 1)**2
   end do

   sn = sn + al*(23.0d0*rr(nin - 2)*uu(nin - 2, 1)**2 &
   &          + 28.0d0*rr(nin - 1)*uu(nin - 1, 1)**2 &
   &          + 9.0d0*rr(nin)*uu(nin, 1)**2)/24.0d0

   cnorm = sqrt(1.0d0/sn)

   uu(:, :) = cnorm*uu(:, :)

   return
end subroutine renorm_r
end module renormalize_m

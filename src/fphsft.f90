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
!> computes full potential scattering log derivatives
!> returns atan(rr(mch) * du/dr / u) which is sort-of like a phase shift
!> and easier to compare in plots than the log derivatives themselves
!> Pauli-type scalar-relativistic calculation
subroutine fphsft(ll, epsh2, depsh, pshf, rr, vv, zz, mmax, mch, npsh, srel)
   implicit none
   integer, parameter :: dp = kind(1.0d0)
   real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
   real(dp), parameter :: pi2 = 2.0d0*pi

   !Input variables
   !> ll  angular momentum
   integer, intent(in) :: ll
   !> mmax  dimension of rr, etc.
   integer, intent(in) :: mmax
   !> npsh  number of energy points in scan
   integer, intent(in) :: npsh
   !> mch  index of radius for log der test
   integer, intent(in) :: mch
   !> rr  radial log grid
   real(dp), intent(in) :: rr(mmax)
   !> vv  all-electron potential
   real(dp), intent(in) :: vv(mmax)
   !> depsh  increment of scan
   real(dp), intent(in) :: depsh
   !> epsh2  upper limit of energy scan
   real(dp), intent(in) :: epsh2
   !> zz  atomic number
   real(dp), intent(in) :: zz
   !> srel .true. for scalar-relativistic, .false. for non-relativistic
   logical, intent(in) :: srel

   !Output variables
   !> pshf  log derivatives "angles", as above
   real(dp), intent(out) :: pshf(npsh)

   !Local variables
   real(dp) :: al, epsh, phi, phip, pshoff
   integer :: ii, ierr, nn

   real(dp), allocatable :: uu(:), up(:)

   allocate (uu(mmax), up(mmax))

   al = 0.01d0*dlog(rr(101)/rr(1))

   pshoff = 0.0d0

   do ii = 1, npsh
      epsh = epsh2 - (ii - 1)*depsh

      call lschfs(nn, ll, ierr, epsh, rr, vv, uu, up, zz, mmax, mch, srel)

      phi = uu(mch)/rr(mch)
      phip = (up(mch) - al*uu(mch))/(al*rr(mch)**2)
      pshf(ii) = atan2(rr(mch)*phip, phi) + pshoff
      if (ii > 1) then
         if (pshf(ii) < pshf(ii - 1)) then
            pshoff = pshoff + pi2
            pshf(ii) = pshf(ii) + pi2
         end if
      end if
   end do

   deallocate (uu, up)
   return
end subroutine fphsft

!
! Copyright (c) 1989-2026 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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
! interpolates various arrays onto linear radial mesh

subroutine linear_interpolate(nc, nv, lmax, lloc, la, mxprj, nproj, icmod, mmax, &
   nrl, drl, rr, vkb, vpuns, rhotae, rhoc, rho, rhomod, uuaea, uupsa, rl, dr, vkbl, &
   vpunsl, rhotael, rhocl, rhol, rhomodl, uuaeal, uupsal)
   implicit none
   integer, parameter :: dp=kind(1.0d0)

   ! Input variables
   !> Number of core states
   integer, intent(in) :: nc
   !> Number of valence states
   integer, intent(in) :: nv
   !> Maximum angular momentum
   integer, intent(in) :: lmax
   !> l for local potential
   integer, intent(in) :: lloc
   !> Angular momentum for each state
   integer, intent(in) :: la(nc + nv)
   !> Maximum number of projectors
   integer, intent(in) :: mxprj
   !> Number of vkb projectors for each l
   integer, intent(in) :: nproj(lmax+1)
   !> Model core charge type
   integer, intent(in) :: icmod
   !> Size of log radial grid
   integer, intent(in) :: mmax
   !> Size of linear radial grid
   integer, intent(in) :: nrl
   ! Linear radial grid spacing
   real(dp), intent(in) :: drl
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> VKB projectors
   real(dp), intent(in) :: vkb(mmax, mxprj, 4)
   !> Unscreened semi-local pseudopotentials
   real(dp), intent(in) :: vpuns(mmax, 5)
   !> Total all-electron charge
   real(dp), intent(in) :: rhotae(mmax)
   !> All-electron core charge
   real(dp), intent(in) :: rhoc(mmax)
   !> Valence pseudocharge
   real(dp), intent(in) :: rho(mmax)
   !> Model core charge
   real(dp), intent(in) :: rhomod(mmax, 5)
   !> All-electron wavefunctions
   real(dp), intent(in) :: uuaea(mmax, nc + nv)
   !> Pseudo wavefunctions
   real(dp), intent(in) :: uupsa(mmax, nv)

   ! Output variables
   !> Linear radial grid
   real(dp), intent(out) :: rl(nrl)
   !> Linear radial grid derivative
   real(dp), intent(out) :: dr(nrl)
   !> VKB projectors on linear grid
   real(dp), intent(out) :: vkbl(nrl, mxprj, 4)
   !> Unscreened semi-local pseudopotentials on linear grid
   real(dp), intent(out) :: vpunsl(nrl, 5)
   !> All-electron charge on linear grid
   real(dp), intent(out) :: rhotael(nrl)
   !> All-electron core charge on linear grid
   real(dp), intent(out) :: rhocl(nrl)
   !> Valence pseudocharge on linear grid
   real(dp), intent(out) :: rhol(nrl)
   !> Model core charge on linear grid
   real(dp), intent(out) :: rhomodl(nrl, 5)
   !> All-electron wavefunctions on linear grid
   real(dp), intent(out) :: uuaeal(nrl, nc + nv)
   !> Pseudo wavefunctions on linear grid
   real(dp), intent(out) :: uupsal(nrl, nv)

   ! Local variables
   integer :: ii, l1
   integer :: iproj

   rl(:) = 0.0_dp
   vkbl(:, :, :) = 0.0_dp
   vpunsl(:, :) = 0.0_dp
   rhotael(:) = 0.0_dp
   rhocl(:) = 0.0_dp
   rhol(:) = 0.0_dp
   rhomodl(:, :) = 0.0_dp
   uuaeal(:, :) = 0.0_dp
   uupsal(:, :) = 0.0_dp

   do ii = 1, nrl
      rl(ii) = drl * real(ii - 1, dp)
      dr(ii) = drl
   end do

   do l1 = 1, max(lmax + 1, lloc + 1)
      ! Interpolate the unscreened semi-local pseudopotentials
      write(6, *) 'Interpolating vpuns for l=', l1 - 1
      call dpnint(rr, vpuns(1, l1), mmax, rl, vpunsl(1, l1), nrl)
      ! If no model core charge, override dpnint extrapolation to zero
      ! for the unscreened semi-local pseudopotentials and local potential
      if (icmod == 0) then
         vpunsl(1, l1) = vpuns(1, l1)
      end if
      ! Interpolate the VKB projectors
      if (l1 /= lloc + 1) then
         do iproj = 1, nproj(l1)
            write(6, *) 'Interpolating vkb for l=', l1 - 1, ' iproj=', iproj
            call dpnint(rr, vkb(1, iproj, l1), mmax, rl, vkbl(1, iproj, l1), nrl)
         end do
      end if
   end do

   ! Interpolate charge densities
   write(6, *) 'Interpolating rhotae, rhoc, rho'
   call dpnint(rr, rhotae, mmax, rl, rhotael, nrl)
   call dpnint(rr, rhoc, mmax, rl, rhocl, nrl)
   call dpnint(rr, rho, mmax, rl, rhol, nrl)

   ! Interpolate the (ii - 1)th derivatives of the model core charge
   write(6, *) 'Interpolating rhomod'
   do ii = 1, 5
      call dpnint(rr, rhomod(1, ii), mmax, rl, rhomodl(1, ii), nrl)
   end do

   ! Interpolate the all-electron wavefunctions
   write(6, *) 'Interpolating uuaea'
   do ii = 1, nc + nv
      call dpnint(rr, uuaea(1, ii), mmax, rl, uuaeal(1, ii), nrl)
   end do

   ! Interpolate the pseudo wavefunctions
   write(6, *) 'Interpolating uupsa'
   do ii = 1, nv
      call dpnint(rr, uupsa(1, ii), mmax, rl, uupsal(1, ii), nrl)
   end do

   return
end subroutine linear_interpolate

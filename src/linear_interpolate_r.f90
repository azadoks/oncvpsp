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

subroutine linear_interpolate_r(nc, nv, lmax, lloc, la, mxprj, nproj, icmod, mmax, &
   nrl, drl, rr, vkb, vpuns, rho, rhomod, uuaea, uupsa, rl, vkbl, vpunsl, rhol,  &
   rhomodl, uuaeal, uupsal)
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
   real(dp), intent(in) :: vkb(mmax, mxprj, 4, 2)
   !> Unscreened semi-local pseudopotentials
   real(dp), intent(in) :: vpuns(mmax, 5)
   !> Valence pseudocharge
   real(dp), intent(in) :: rho(mmax)
   !> Model core charge
   real(dp), intent(in) :: rhomod(mmax, 5)
   !> All-electron wavefunctions
   real(dp), intent(in) :: uuaea(mmax, 2, nc + nv)
   !> Pseudo wavefunctions
   real(dp), intent(in) :: uupsa(mmax, 2, nv)

   ! Output variables
   !> Linear radial grid
   real(dp), intent(out) :: rl(nrl)
   !> VKB projectors on linear grid
   real(dp), intent(out) :: vkbl(nrl, mxprj, 4, 2)
   !> Unscreened semi-local pseudopotentials on linear grid
   real(dp), intent(out) :: vpunsl(nrl, 5)
   !> Valence pseudocharge on linear grid
   real(dp), intent(out) :: rhol(nrl)
   !> Model core charge on linear grid
   real(dp), intent(out) :: rhomodl(nrl, 5)
   !> All-electron wavefunctions on linear grid
   real(dp), intent(out) :: uuaeal(nrl, 2, nc + nv)
   !> Pseudo wavefunctions on linear grid
   real(dp), intent(out) :: uupsal(nrl, 2, nv)

   ! Local variables
   integer :: ii, l1, ll
   integer :: iproj, ikap, mkap

   rl(:) = 0.0_dp
   vkbl(:, :, :, :) = 0.0_dp
   vpunsl(:, :) = 0.0_dp
   rhol(:) = 0.0_dp
   rhomodl(:, :) = 0.0_dp
   uuaeal(:, :, :) = 0.0_dp
   uupsal(:, :, :) = 0.0_dp

   do ii = 1, nrl
      rl(ii) = drl * real(ii - 1, dp)
   end do

   do l1 = 1, max(lmax + 1, lloc + 1)
      ! Interpolate the unscreened semi-local pseudopotentials
      call dpnint(rr, vpuns(1, l1), mmax, rl, vpunsl(1, l1), nrl)
      ! If no model core charge, override dpnint extrapolation to zero
      ! for the unscreened semi-local pseudopotentials and local potential
      if (icmod == 0) then
         vpunsl(1, l1) = vpuns(1, l1)
      end if
   end do

   ! Interpolate the VKB projectors
   do l1 = 1, lmax + 1
      ! Skip local potential angular momentum channel
      if (l1 == lloc + 1) cycle
      ! J = ll + 1/2 for ll = 0
      ! J = ll ± 1/2 for ll > 0
      ll = l1 - 1
      if (ll == 0) then
         mkap = 1
      else
         mkap = 2
      end if
      ! Loop on J = ll ± 1/2
      do ikap = 1, mkap
         do iproj = 1, nproj(l1)
            call dpnint(rr, vkb(1, iproj, l1, ikap), mmax, rl, vkbl(1, iproj, l1, ikap), nrl)
         end do  ! jj
      end do  ! ikap
   end do  ! l1

   ! Interpolate the valence pseudocharge
   call dpnint(rr, rho, mmax, rl, rhol, nrl)

   ! Interpolate the (ii - 1)th derivatives of the model core charge
   do ii = 1, 5
      call dpnint(rr, rhomod(1, ii), mmax, rl, rhomodl(1, ii), nrl)
   end do

   ! Interpolate the all-electron wavefunctions
   do ii = 1, nc + nv
      l1 = la(ii)
      if (l1 == 0) then
         mkap = 1
      else
         mkap = 2
      end if
      do ikap = 1, mkap
         call dpnint(rr, uuaea(1, ikap, ii), mmax, rl, uuaeal(1, ikap, ii), nrl)
      end do
   end do

   ! Interpolate the pseudo wavefunctions
   do ii = 1, nv
      l1 = la(ii + nc)
      if (l1 == 0) then
         mkap = 1
      else
         mkap = 2
      end if
      do ikap = 1, mkap
         call dpnint(rr, uupsa(1, ikap, ii), mmax, rl, uupsal(1, ikap, ii), nrl)
      end do
   end do

   return
end subroutine linear_interpolate_r

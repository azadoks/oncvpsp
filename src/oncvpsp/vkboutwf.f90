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
subroutine vkboutwf(ll, nvkb, ep, vkb, evkb, rr, vloc, uu, up, node, mmax, mch)

   ! computes Vanderbilt / Kleinman-Bylander outward-integrated wave functions

   !ll  angular momentum
   !nvkb  switch for 1 or 2 projedtors
   !ep  energy at which wave function is to be calculated
   !vkb  Vanderbilt-Kleinman-Bylander projectors for this l
   !evkb  projector coefficients
   !rr  log radial mesh
   !vloc  local pseudopotential
   !uu  wave function
   !up  1st derivative of uu
   !node  count of number of nodes from 0 to rr(mch)
   !mmax  dimension of log mesh
   !mch  index of radius to which wave function is to be integrated
   use precision_m, only: dp
   use lsch_m, only: lschps, lschkb
   implicit none

   !Input variables
   integer, intent(in) :: nvkb
   integer, intent(in) :: ll
   integer, intent(in) :: mmax
   integer, intent(in) :: mch
   real(dp), intent(in) :: rr(mmax)
   real(dp), intent(in) :: vloc(mmax)
   real(dp), intent(in) :: vkb(mmax, nvkb)
   real(dp), intent(in) :: evkb(nvkb)
   !> Energy at which wave function is to be calculated
   real(dp), intent(in) :: ep

   !Output variables
   real(dp), intent(out) :: uu(mmax)
   real(dp), intent(out) :: up(mmax)
   integer, intent(out) :: node

   !Local variables
   real(dp) :: rc, rn, ep_tmp
   real(dp), allocatable :: phi(:, :), phip(:, :)
   real(dp), allocatable :: gg0(:), gg(:, :)
   integer, allocatable :: ipiv(:)

   integer :: ii, jj, ierr, info

   uu(:) = 0.0_dp
   up(:) = 0.0_dp
   ep_tmp = ep

   ! homogeneous solution
   call lschps(ll, ierr, ep_tmp, rr, vloc, uu, up, mmax, mch)

   rc = 0.0_dp
   if (nvkb /= 0) then
      ! find cutoff radius for projectors
      do ii = mmax, 1, -1
         if (abs(vkb(ii, 1)) > 0.0_dp) then
            rc = rr(ii)
            exit
         end if
      end do

      allocate (phi(mmax, nvkb), phip(mmax, nvkb))
      allocate (gg(nvkb, nvkb), gg0(nvkb))
      allocate (ipiv(nvkb))

      phi(:, :) = 0.0_dp
      phip(:, :) = 0.0_dp
      gg(:, :) = 0.0_dp
      gg0(:) = 0.0_dp

      ! inhomogeneous solutions
      do ii = 1, nvkb
         call lschkb(ll, ierr, ep_tmp, vkb(1, ii), rr, vloc, phi(1, ii), phip(1, ii), mmax, mch)
      end do

      ! projector matrix elements and coefficient matrix
      do jj = 1, nvkb
         call vpinteg(uu, vkb(1, jj), mch, 2 * ll + 2, gg0(jj), rr)
         gg0(jj) = evkb(jj) * gg0(jj)
         do ii = 1, nvkb
            call vpinteg(phi(1, ii), vkb(1, jj), mch, 2 * ll + 2, gg(jj, ii), rr)
            gg(jj, ii) = -evkb(jj) * gg(jj, ii)
         end do
         gg(jj, jj) = 1.0_dp + gg(jj, jj)
      end do

      ! solve linear equations for coefficients

      !    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      call dgesv(nvkb, 1, gg, nvkb, ipiv, gg0, nvkb, info)
      if (info /= 0) then
         write (6, '(/a,i4)') 'vkbout: dgesv ERROR, stopping info =', info
         stop
      end if

      ! output wave functions
      do jj = 1, nvkb
         uu(:) = uu(:) + gg0(jj) * phi(:, jj)
         up(:) = up(:) + gg0(jj) * phip(:, jj)
      end do

      deallocate (phi, phip)
      deallocate (gg, gg0)
      deallocate (ipiv)
   end if

   ! lower cutoff for node counting to avoid small-r noise
   rn = max(0.1_dp * rc, 0.05_dp)

   node = 0
   do ii = 6, mch
      ! note historic evolution!
      !  if(rr(ii)>0.5d0 .and. uu(ii-1)*uu(ii)<0.0_dp) then
      !  if(rr(ii)>0.25d0 .and. uu(ii-1)*uu(ii)<0.0_dp) then
      !  if(rr(ii)>0.10d0 .and. uu(ii-1)*uu(ii)<0.0_dp) then
      !  if(rr(ii)>0.07d0 .and. uu(ii-1)*uu(ii)<0.0_dp) then
      if (rr(ii) > rn .and. uu(ii - 1) * uu(ii) < 0.0_dp) then
         node = node + 1
      end if
   end do

   return
end subroutine vkboutwf

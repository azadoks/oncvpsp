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
!> computes log derivatives, actually atan(r * ((d psi(r)/dr)/psi(r)))
!> at rr(irphs) comparing all-electron with Vanderbilt-Kleinman-Bylander
!> results for 1 and 2 projectors, or the semi-local pseudpotential
!> when that is the local potential for some l
!> the computed quantity is reminiscent of a scattering phase shift, but isn't
subroutine run_phsft(lmax, lloc, nproj, epa, epsh1, epsh2, depsh, vkb, evkb, &
&                     rr, vfull, vp, zz, mmax, mxprj, irc, srel)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> lmax  maximum angular momentum
   integer, intent(in) :: lmax
   !> lloc  l for local potential
   integer, intent(in) :: lloc
   !> mmax  size of radial grid
   integer, intent(in) :: mmax
   !> mxprj dimension of number of projectors
   integer, intent(in) :: mxprj
   !> nproj  number ov V / KB projectors for  each l
   integer, intent(in) :: nproj(6)
   integer, intent(in) :: irc(6)
   !> epsh1  low energy limit for "phase shift" calculation
   real(dp), intent(in) :: epsh1
   !> epsh2  high energy limit for "phase shift" calculation
   real(dp), intent(in) :: epsh2
   !> depsh  energy increment
   real(dp), intent(in) :: depsh
   !> zz  atomic number
   real(dp), intent(in) :: zz
   !> rr  log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> vp  semi-local pseudopotentials (vp(:,5) is local potential if linear comb.)
   real(dp), intent(in) :: vp(mmax, 5)
   !> ep  bound-state or scattering state reference energies for vkb potentials
   real(dp), intent(in) :: epa(mxprj, 6)
   !> vfull  all-electron potential
   real(dp), intent(in) :: vfull(mmax)
   !> vkb  VKB projectors
   real(dp), intent(in) :: vkb(mmax, mxprj, 4)
   !> evkb  coefficients of VKB projectors
   real(dp), intent(in) :: evkb(mxprj, 4)
   !> srel .true. for scalar-relativistic, .false. for non-relativistic
   logical, intent(in) :: srel

   !Output variables - printing only

   !Local variables
   !> irphs  index of rr beyond which all vp==vlocal
   integer :: irphs
   integer :: ii, ll, l1, npsh
   real(dp) :: epsh

   real(dp), allocatable :: pshf(:), pshp(:)

   npsh = int(((epsh2 - epsh1)/depsh) - 0.5d0) + 1

   allocate (pshf(npsh), pshp(npsh))

   ! loop for phase shift calculation -- full, then local or Kleinman-
   ! Bylander / Vanderbilt

   do l1 = 1, 4

      ll = l1 - 1
      if (ll <= lmax) then
         irphs = irc(l1) + 2
      else
         irphs = irc(lloc + 1)
      end if

      call fphsft(ll, epsh2, depsh, pshf, rr, vfull, zz, mmax, irphs, npsh, srel)
      if (ll .eq. lloc) then
         call vkbphsft(ll, 0, epsh2, depsh, epa(1, l1), pshf, pshp, &
         &                   rr, vp(1, lloc + 1), vkb(1, 1, l1), evkb(1, l1), &
         &                   mmax, irphs, npsh)
      else
         call vkbphsft(ll, nproj(l1), epsh2, depsh, epa(1, l1), pshf, pshp, &
         &                   rr, vp(1, lloc + 1), vkb(1, 1, l1), evkb(1, l1), &
         &                   mmax, irphs, npsh)
      end if

      write (6, '(/a,i2)') 'log derivativve data for plotting, l=', ll
      write (6, '(a,f6.2)') 'atan(r * ((d psi(r)/dr)/psi(r))), r=', rr(irphs)
      write (6, '(a/)') 'l, energy, all-electron, pseudopotential'
      do ii = 1, npsh
         epsh = epsh2 - depsh*dfloat(ii - 1)
         write (6, '(a, i6, 3f12.6)') '! ', ll, epsh, pshf(ii), pshp(ii)
      end do
   end do
   deallocate (pshf, pshp)
   return
end subroutine run_phsft

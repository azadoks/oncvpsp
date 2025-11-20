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
module test_log_der_m
   use schroedinger_m, only: lschfs, vkboutwf
   use dirac_m, only: ldiracfs
   implicit none
   private
   public :: run_phsft
   public :: run_phsft_r
contains
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
   integer :: lmax
   !> lloc  l for local potential
   integer :: lloc
   !> mmax  size of radial grid
   integer :: mmax
   !> mxprj dimension of number of projectors
   integer :: mxprj
   !> nproj  number ov V / KB projectors for  each l
   integer :: nproj(6)
   integer :: irc(6)
   !> epsh1  low energy limit for "phase shift" calculation
   real(dp) :: epsh1
   !> epsh2  high energy limit for "phase shift" calculation
   real(dp) :: epsh2
   !> depsh  energy increment
   real(dp) :: depsh
   !> zz  atomic number
   real(dp) :: zz
   !> rr  log radial grid
   real(dp) :: rr(mmax)
   !> vp  semi-local pseudopotentials (vp(:,5) is local potential if linear comb.)
   real(dp) :: vp(mmax, 5)
   !> ep  bound-state or scattering state reference energies for vkb potentials
   real(dp) :: epa(mxprj, 6)
   !> vfull  all-electron potential
   real(dp) :: vfull(mmax)
   !> vkb  VKB projectors
   real(dp) :: vkb(mmax, mxprj, 4)
   !> evkb  coefficients of VKB projectors
   real(dp) :: evkb(mxprj, 4)
   !> srel .true. for scalar-relativistic, .false. for non-relativistic
   logical :: srel

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
   integer :: ll
   !> mmax  dimension of rr, etc.
   integer :: mmax
   !> npsh  number of energy points in scan
   integer :: npsh
   !> mch  index of radius for log der test
   integer :: mch
   !> rr  radial log grid
   real(dp) :: rr(mmax)
   !> vv  all-electron potential
   real(dp) :: vv(mmax)
   !> depsh  increment of scan
   real(dp) :: depsh
   !> epsh2  upper limit of energy scan
   real(dp) :: epsh2
   !> zz  atomic number
   real(dp) :: zz
   !> srel .true. for scalar-relativistic, .false. for non-relativistic
   logical :: srel

   !Output variables
   !> pshf  log derivatives "angles", as above
   real(dp) :: pshf(npsh)

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
!> computes Vanderfilt / Kleinman-Bylander scattering log derivatives
!> (or semi-local if ivkb==0)
!> returns atan(rr(mch) * du/dr / u) which is sort-of like a phase shift
!> and easier to compare in plots than the log derivatives themselves
!> Pauli-type scalar-relativistic calculation
subroutine vkbphsft(ll, ivkb, epsh2, depsh, ep, pshf, pshp, &
& rr, vloc, vkb, evkb, mmax, mch, npsh)
   implicit none
   integer, parameter :: dp = kind(1.0d0)
   real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
   real(dp), parameter :: pi2 = 2.0d0*pi
   real(dp), parameter :: eps = 1.0d-8

   !Input variables
   !> ll  angular momentum
   integer :: ll
   !> ivkb  number of projectors
   integer :: ivkb
   !> mmax  dimension of rr, etc.
   integer :: mmax
   !> npsh  number of energy points in scan
   integer :: npsh
   !> mch  index of radius for log der test
   integer :: mch
   real(dp) :: rr(mmax)
   !> vloc  local part of psp
   real(dp) :: vloc(mmax)
   !> vkb  VKB projectorsx
   real(dp) :: vkb(mmax, *)
   real(dp) :: evkb(*)
   !> pshf  all-electron log derivatives "angles", as above (input)
   real(dp) :: pshf(npsh)
   !> depsh  increment of scan
   real(dp) :: depsh
   !> epsh2  upper limit of energy scan
   real(dp) :: epsh2
   !> ep  reference energy for psp creation (bound or scattering)
   real(dp) :: ep

   !Output variables
   !> pshp  pseudopotential log derivatives "angles", as above (output)
   real(dp) :: pshp(npsh)

   !Local variables
   real(dp) :: al, dnpi, epsh, phi, phip, pshoff
   real(dp) :: dmin, dmax, dtst, emin, emax, et, psmin, psmax, pst, shift, shift2
   integer :: ii, jj, ierr, node, n2
   logical :: jump

   real(dp), allocatable :: uu(:), up(:)

   allocate (uu(mmax), up(mmax))

   al = 0.01d0*dlog(rr(101)/rr(1))

   pshoff = 0.0d0

   do ii = 1, npsh
      epsh = epsh2 - (ii - 1)*depsh

      call vkboutwf(ll, ivkb, epsh, vkb, evkb, rr, vloc, uu, up, node, mmax, mch)

      phi = uu(mch)/rr(mch)
      phip = (up(mch) - al*uu(mch))/(al*rr(mch)**2)
      pshp(ii) = atan2(rr(mch)*phip, phi) + pshoff

      ! Shift for continuity and to avoid false jumps suggesting ghost states

      ! -2 pi jump is harmless cycling from pi to -pi
      if (ii > 1) then
         if (pshp(ii) < pshp(ii - 1) - 4.0d0) then
            pshoff = pshoff + pi2
            pshp(ii) = pshp(ii) + pi2
            ! +/- pi jump can be actual bound semi-core state, ghost state, or spurious
            ! jump from sudden change of sign of both uu and up with no real change
            ! of shape.
         else if (abs(pshp(ii) - pshp(ii - 1)) > 2.0d0) then
            shift = sign(pi, pshp(ii) - pshp(ii - 1))

            ! interval-halving search to determine if this is a discontinuous pi jump
            ! or varies continuously on some scale, indicating a real or ghost
            ! bound state / resonance
            jump = .true.
            emin = epsh
            emax = epsh + depsh
            psmin = pshp(ii)
            psmax = pshp(ii - 1)
            !       if((psmax-psmin)>pi) psmin=psmin+pi2
            if ((psmax - psmin) > 2.9d0) psmin = psmin + pi2

            do jj = 1, 25

               if (.not. jump) cycle
               et = 0.5d0*(emin + emax)
               shift2 = 0.0d0

               call vkboutwf(ll, ivkb, et, vkb, evkb, rr, vloc, uu, up, node, mmax, mch)
               phi = uu(mch)/rr(mch)
               phip = (up(mch) - al*uu(mch))/(al*rr(mch)**2)
               pst = atan2(rr(mch)*phip, phi) + pshoff
               if ((psmax - pst) > 2.9d0) then
                  pst = pst + pi2
                  shift2 = pi2
               end if

               dmin = abs(pst - psmin)
               dmax = abs(psmax - pst)

               if (dmin > dmax) then
                  emax = et
                  psmax = pst
               else
                  emin = et
                  psmin = pst
               end if

               ! this is the test to see if an interval of change less than pi had been
               ! reached
               if (dmax < 0.5d0 .and. dmin < 0.5d0) jump = .false.

            end do

            ! if this is a spurious abrupt jump, restore +/- pi
            if (jump) then
               pshoff = pshoff - shift
               pshp(ii) = pshp(ii) - shift
               ! if this is a real continuous transition, check for 2 pi issue and fix
            else if (pshp(ii) < pshp(ii - 1) - 2.5d0) then
               pshoff = pshoff + pi2
               pshp(ii) = pshp(ii) + pi2
            end if

         end if
      end if

      ! calculate shift to align with all-electron results
      if (abs(ep - epsh) - 0.5d0*depsh < eps) then
         dnpi = pi*(nint((pshf(ii) - pshp(ii))/pi))
      end if

   end do

   pshp(:) = pshp(:) + dnpi

   deallocate (uu, up)
   return
end subroutine vkbphsft
!> computes log derivatives, actually atan(r * ((d psi(r)/dr)/psi(r)))
!> at rr(irphs) comparing all-electron with Vanderbilt-Kleinman-Bylander
!> results for 1 and 2 projectors, or the semi-local pseudpotential
!> when that is the local potential for some l
!> the computed quantity is reminiscent of a scattering phase shift, but isn't
!> This version is for fully-relativistic pseudopotentials
subroutine run_phsft_r(lmax, lloc, nproj, ep, epsh1, epsh2, depsh, vkb, evkb, &
&                     rr, vfull, vp, zz, mmax, mxprj, irc)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> lmax  maximum angular momentum
   integer :: lmax
   !> lloc  l for local potential
   integer :: lloc
   !> mmax  size of radial grid
   integer :: mmax
   !> mxprj  dimension of number of projectors
   integer :: mxprj
   !> nproj  number ov VKB projectors for  each l
   integer :: nproj(6)
   !> irc  core radii
   integer :: irc(6)
   !> epsh1  low energy limit for "phase shift" calculation
   real(dp) :: epsh1
   !> epsh2  high energy limit for "phase shift" calculation
   real(dp) :: epsh2
   !> depsh  energy increment
   real(dp) :: depsh
   !> zz  atomic number
   real(dp) :: zz
   !> rr  log radial grid
   real(dp) :: rr(mmax)
   !> vp  semi-local pseudopotentials (vp(:,5) is local potential if linear comb.)
   real(dp) :: vp(mmax, 5, 2)
   !> ep  bound-state or scattering state reference energies for vkb potentials
   real(dp) :: ep(6, 2)
   !> vfull  all-electron potential
   real(dp) :: vfull(mmax)
   !> vkb  VKB projectors
   real(dp) :: vkb(mmax, mxprj, 4, 2)
   !> evkb  coefficients of VKB projectors
   real(dp) :: evkb(mxprj, 4, 2)

   !Output variables - printing only

   !Local variables
   !> irphs  index of rr beyond which all vp==vlocal
   integer :: irphs
   integer :: ii, ll, l1, npsh
   integer :: ikap, kap, mkap
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

      if (l1 == 1) then
         mkap = 1
      else
         mkap = 2
      end if
      ! loop on J = ll +/- 1/2
      do ikap = 1, mkap
         if (ikap == 1) kap = -(ll + 1)
         if (ikap == 2) kap = ll
         call fphsft_r(ll, kap, epsh2, depsh, pshf, rr, vfull, zz, mmax, irphs, npsh)
         if (ll .eq. lloc) then
            call vkbphsft(ll, 0, epsh2, depsh, ep(l1, ikap), pshf, pshp, &
            &                   rr, vp(1, lloc + 1, ikap), vkb(1, 1, l1, ikap), evkb(1, l1, ikap), &
            &                   mmax, irphs, npsh)
         else
            call vkbphsft(ll, nproj(l1), epsh2, depsh, ep(l1, ikap), pshf, pshp, &
            &                   rr, vp(1, lloc + 1, ikap), vkb(1, 1, l1, ikap), evkb(1, l1, ikap), &
            &                   mmax, irphs, npsh)
         end if

         write (6, '(/a,i2)') 'log derivativve data for plotting, l=', ll
         write (6, '(a,f6.2)') 'atan(r * ((d psi(r)/dr)/psi(r))), r=', rr(irphs)
         write (6, '(a/)') 'l, energy, all-electron, pseudopotential'
         do ii = 1, npsh
            epsh = epsh2 - depsh*dfloat(ii - 1)
            if (ikap == 1) then
               write (6, '(a, i6, 3f12.6)') '! ', -ll, epsh, pshf(ii), pshp(ii)
            else
               write (6, '(a, i6, 3f12.6)') '! ', ll, epsh, pshf(ii), pshp(ii)
            end if
         end do
      end do !ikap
   end do !l1
   deallocate (pshf, pshp)
   return
end subroutine run_phsft_r
!> computes full potential scattering log derivatives
!> returns atan(rr(mch) * du/dr / u) which is sort-of like a phase shift
!> and easier to compare in plots than the log derivatives themselves
!> Dirac equation, log derivatives based on large component only
subroutine fphsft_r(ll, kap, epsh2, depsh, pshf, rr, vv, zz, mmax, mch, npsh)
   implicit none
   integer, parameter :: dp = kind(1.0d0)
   real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
   real(dp), parameter :: pi2 = 2.0d0*pi

   !Input variables
   !> ll  angular momentum
   integer :: ll
   !> kap =l, -(l+1) for j=l -/+ 1/2
   integer :: kap
   !> mmax  dimension of rr, etc.
   integer :: mmax
   !> npsh  number of energy points in scan
   integer :: npsh
   !> mch  index of radius for log der test
   integer :: mch
   !> rr  radial log grid
   real(dp) :: rr(mmax)
   !> vv  all-electron potential
   real(dp) :: vv(mmax)
   !> depsh  increment of scan
   real(dp) :: depsh
   !> epsh2  upper limit of energy scan
   real(dp) :: epsh2
   !> zz  atomic number
   real(dp) :: zz

   !Output variables
   !> pshf  log derivatives "angles", as above
   real(dp) :: pshf(npsh)

   !Local variables
   real(dp) :: al, epsh, phi, phip, pshoff
   integer :: ii, ierr, nnae

   real(dp), allocatable :: uu(:, :), up(:, :)

   allocate (uu(mmax, 2), up(mmax, 2))

   al = 0.01d0*dlog(rr(101)/rr(1))

   pshoff = 0.0d0

   do ii = 1, npsh
      epsh = epsh2 - (ii - 1)*depsh

      call ldiracfs(nnae, ll, kap, ierr, epsh, rr, zz, vv, uu, up, mmax, mch)

      phi = uu(mch, 1)/rr(mch)
      phip = (up(mch, 1) - al*uu(mch, 1))/(al*rr(mch)**2)
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
end subroutine fphsft_r
end module test_log_der_m

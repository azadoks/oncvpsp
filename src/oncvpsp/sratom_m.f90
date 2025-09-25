!
! Copyright (c) 1989-201r by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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

module sratom_m
   use precision_m, only: dp
   use lsch_m, only: lschfb
   use atom_config_m, only: atom_config_t
   use mesh_m, only: log_mesh_t
   use ae_atom_m, only: ae_atom_t
   implicit none
   private
   public :: sratom
contains
!> self-consistent scalar-relativistic all-electron atom
!> calculation using log mesh (non-relativistic when srel=.false.)
subroutine sratom(na, la, ea, fa, rpk, nc, ncv, it, rhoc, rho, &
                  rr, vi, zz, mmax, iexc, etot, ierr, srel)
   !na  principal quantum number array, dimension ncv
   !la  angular-momenta
   !ea  eigenvalues (output)
   !fa  occupancies
   !rpk  radius of outermost peak of wave function
   !nc  number of core states
   !ncv  number of core+valence states
   !it  number of iterations (output)
   !rr  log radial mesh
   !vi  all-electron potential (output)
   !zz  atomic number
   !mmax  size of log grid
   !iexc  exchange-correlation function to be used
   !etot  all-electron total energy (output)
   !ierr  error flag
   !srel  .true. for scalar-relativistic, .false. for non-relativistic

   !Input variables
   integer, intent(in) :: mmax, iexc, nc, ncv
   integer, intent(in) :: na(ncv), la(ncv)
   real(dp), intent(in) :: zz
   real(dp), intent(in) :: fa(ncv), rr(mmax)
   logical, intent(in) :: srel

   !Output variables
   integer, intent(out) :: it, ierr
   real(dp), intent(out) :: etot
   real(dp), intent(out) :: ea(ncv), rpk(ncv)
   real(dp), intent(out) :: rho(mmax), rhoc(mmax), vi(mmax)

   !Local function
   real(dp) :: tfapot

   !Local variables
   integer :: nin, mch
   real(dp) :: amesh, al
   real(dp) :: dr, eeel, eexc, et, rl, rl1, sd, sf, sn, eeig
   real(dp) :: thl, vn, zion
   integer :: ii, jj
   logical :: convg

   real(dp), allocatable :: u(:), up(:)
   real(dp), allocatable :: vo(:), vi1(:), vo1(:), vxc(:)

   ! blend parameter for Anderson iterative potential mixing
   real(dp), parameter :: bl = 0.5d0

   allocate (u(mmax), up(mmax))
   allocate (vo(mmax), vi1(mmax), vo1(mmax), vxc(mmax))

   ! why all this is necessary is unclear, but it seems to be
   u(:) = 0.d0; up(:) = 0.d0; vo(:) = 0.d0; vi1(:) = 0.d0; vo1(:) = 0.d0; vxc(:) = 0.d0
   dr = 0.d0; eeel = 0.d0; eexc = 0.d0; et = 0.d0; rl = 0.d0; rl1 = 0.d0
   sd = 0.d0; sf = 0.d0; sn = 0.d0; eeig = 0.d0; thl = 0.d0; vn = 0.d0; zion = 0.d0
   nin = 0; mch = 0

   al = 0.01d0 * dlog(rr(101) / rr(1))
   amesh = dexp(al)

   do ii = 1, mmax
      vi(ii) = tfapot(rr(ii), zz)
   end do

   ! starting approximation for energies
   sf = 0.0d0
   do ii = 1, ncv
      sf = sf + fa(ii)
      zion = zz + 1.0d0 - sf
      ea(ii) = -0.5d0 * (zion / na(ii))**2
      if (ea(ii) > vi(mmax)) ea(ii) = 2.0d0 * vi(mmax)
   end do

   ! big self  self-consietency loop

   do it = 1, 100
      convg = .true.

      rhoc(:) = 0.0d0
      rho(:) = 0.0d0

      ! solve for bound states in turn
      eeig = 0.0d0
      do ii = 1, ncv

         ! skip unoccupied states
         if (fa(ii) == 0.0d0) then
            ea(ii) = 0.0d0
            cycle
         end if
         et = ea(ii)
         ierr = 0
         call lschfb(na(ii), la(ii), ierr, et, rr, vi, u, up, zz, mmax, mch, srel)
         if (ierr /= 0) then
            write (6, '(/a,3i4)') 'sratom123: lschfb convergence ERROR n,l,iter=', &
            &       na(ii), la(ii), it
            stop
         end if

         ! overall convergence criterion based on eps within lschfb
         if (ea(ii) /= et) convg = .false.
         ea(ii) = et

         ! accumulate charge and eigenvalues
         eeig = eeig + fa(ii) * ea(ii)
         rho(:) = rho(:) + fa(ii) * (u(:) / rr(:))**2
         if (ii <= nc) then
            rhoc(:) = rhoc(:) + fa(ii) * (u(:) / rr(:))**2
         end if

         ! find outermost peak of wavefunction
         do jj = mch - 1, 1, -1
            if (up(jj) * up(jj + 1) < 0.0d0) then
               rpk(ii) = rr(jj)
               exit
            end if
         end do

      end do

      if (ierr /= 0) then
         exit
      end if

      ! output potential
      call vout(0, rho, rhoc, vo, vxc, sf - zz, eeel, eexc, rr, mmax, iexc)

      etot = eeig + eexc - 0.5d0 * eeel

      ! generate next iteration using d. g. anderson''s
      ! method
      thl = 0.0d0
      if (it > 1) then
         sn = 0.0d0
         sd = 0.0d0
         do ii = 1, mmax
            rl = vo(ii) - vi(ii)
            rl1 = vo1(ii) - vi1(ii)
            dr = rl - rl1
            sn = sn + rl * dr * rr(ii)**2
            sd = sd + dr * dr * rr(ii)**2
         end do
         thl = sn / sd
      end if

      do ii = 1, mmax
         vn = (1.0d0 - bl) * ((1.0d0 - thl) * vi(ii) + thl * vi1(ii)) &
         &   + bl * ((1.0d0 - thl) * vo(ii) + thl * vo1(ii))
         vi1(ii) = vi(ii)
         vo1(ii) = vo(ii)
         vi(ii) = vn
      end do

      if (convg) exit

      if (it == 100 .and. .not. convg) then
         write (6, '(/a)') 'sratom: WARNING failed to converge'
      end if

   end do  !it

   if (.not. convg .and. ierr == 0) then
      ierr = 100
   end if

   ! total energy output

   ! output potential for e-e interactions

   call vout(0, rho, rhoc, vo, vxc, sf, eeel, eexc, rr, mmax, iexc)

   etot = eeig + eexc - 0.5d0 * eeel

   deallocate (u, up)
   deallocate (vo, vi1, vo1, vxc)
   return

end subroutine sratom

subroutine sratom_structured(config, mesh, iexc, atom, ierr, srel)
   ! Constants
   !> Parameter for Anderson iterative potential mixing
   real(dp), parameter :: MIX_BETA = 0.5_dp
   !> Maximum number of self-consistent iterations
   integer, parameter :: ITER_MAX = 100
   !> Tolerance for eigenvalue convergence (Ha)
   real(dp), parameter :: EIG_TOL = 1.0e-22_dp

   ! Input variables
   type(atom_config_t), intent(in) :: config
   type(log_mesh_t), intent(in) :: mesh
   integer, intent(in) :: iexc
   logical, optional, intent(in) :: srel

   ! Output variables
   type(ae_atom_t), intent(out) :: atom
   integer, intent(out) :: ierr

   ! Local variables
   !> Total number of states
   integer :: n_states
   !> Number of occupied states
   integer :: n_occupied_states
   !> Current input and output potentials
   real(dp), allocatable :: vin(:), vout(:)
   !> Previous input and output potentials
   real(dp), allocatable :: vin_old(:), vout_old(:)
   !> Exchange-correlation potential
   real(dp), allocatable :: vxc(:)
   !> Occupation cumulative sum
   real(dp) :: occ_cumsum
   !> Convergence flag
   logical :: convgerged
   !> Previous eigenvalue
   real(dp) :: eig_old
   !> Matching points for inward-outward integration
   integer, allocatable :: in_out_matches(:)

   n_states = size(config%ns)
   n_occupied_states = count(config%occs > tiny(1.0_dp))
   if (n_occupied_states == 0) then
      error stop 'sratom_structured: No occupied states in configuration'
   end if

   atom%etot = 0.0_dp
   atom%ehart = 0.0_dp
   atom%exc = 0.0_dp
   allocate(atom%ns(n_states))
   atom%ns = config%ns(1:n_states)
   allocate(atom%ells(n_states))
   atom%ells = config%ells(1:n_states)
   allocate(atom%occs(n_states))
   atom%occs = config%occs(1:n_states)
   allocate(atom%eigs(n_states))
   atom%eigs = 0.0_dp
   allocate(atom%peaks(n_states))
   atom%peaks = 0.0_dp
   allocate(atom%chg_den_tot(mesh%len))
   atom%chg_den_tot = 0.0_dp
   allocate(atom%pot(mesh%len))
   atom%pot = 0.0_dp
   allocate(atom%wfs(mesh%len, n_states))
   atom%wfs = 0.0_dp
   allocate(atom%wfs_deriv(mesh%len, n_states))
   atom%wfs_deriv = 0.0_dp
   allocate(atom%chg_dens(mesh%len, n_states))
   atom%chg_dens = 0.0_dp

   allocate(vin(mesh%len))
   vin(:) = 0.0_dp
   allocate(vout(mesh%len))
   vout(:) = 0.0_dp
   allocate(vin_old(mesh%len))
   vin_old(:) = 0.0_dp
   allocate(vout_old(mesh%len))
   vout_old(:) = 0.0_dp
   allocate(vxc(mesh%len))
   vxc(:) = 0.0_dp
   allocate(in_out_matches(n_states))
   in_out_matches(:) = 0

   ierr = 0

   ! Use Thomas-Fermi potential as initial guess
   do i = 1, mesh%len
      vin(i) = tfapot(mesh%r(i), config%z, dp)
   end do

   ! Initial eigenvalue estimates
   occ_cumsum = 0.0_dp
   do i = 1, n_states
      if (atom%occs(i) < tiny(1.0_dp)) then
         atom%eigs(i) = 0.0_dp
         cycle
      end if
      occ_cumsum = occ_cumsum + atom%occs(i)
      atom%eigs(i) = -0.5_dp * ((config%z + 1.0_dp - occ_cumsum) / atom%ns(i))**2
      if (atom%eigs(i) > vin(mesh%len)) then
         atom%eigs(i) = 2.0d0 * vin(mesh%len)
      end if
   end do

   ! Solve self-consistent equations
   do iter = 1, ITER_MAX
      convgerged = .true.
      ! Solve for bound states in turn
      do i = 1, n_states
         if (atom%occs(i) < tiny(1.0_dp)) then
            cycle
         end if
         ! Copy the current eigenvalue for the convergence check below
         eig_old = atom%eigs(i)
         ! Compute the bound state eigenvalue, wavefunction, and wavefunction derivative
         call lschfb(atom%ns(i), atom%ells(i), ierr, atom%eigs(i), mesh%r, vin, &
                     atom%wfs(:, i), atom%wfs_deriv(:, i), config%z, mesh%len, in_out_matches(i), &
                     merge(srel, .false., present(srel)))
         ! Errors in lschfb are unrecoverable here
         if (ierr /= 0) then
            error stop 'sratom_structured: lschfb failed'
         end if
         atom%chg_den_tot(:) = atom%chg_den_tot(:) + atom%chg_dens(:, i)
         if (abs(atom%eigs(i) - eig_old) > EIG_TOL) then
            convgerged = .false.
         end if
      end do  ! n_states
      if (converged) exit
      call compute_hartree(atom%chg_den_tot, atom%pot_hartree, atom%ehartree, atom%z - sum(atom%occs), mesh%r, mesh%len)
      call compute_xc(atom%chg_den_tot, atom%pot_xc, atom%exc, iexc, mesh%r, mesh%len)
      vout(:) = atom%pot_hartree(:) + atom%pot_xc(:)
      call mix_anderson(vin, vout, vin_old, vout_old, MIX_BETA, mesh%r, mesh%len)
   end do  ! iter

   if (.not. converged) then
      ierr = 1
      write (6, '(/a)') 'sratom_structured: WARNING failed to converge'
   end if

   ! Fill in atom contents with final results
   atom%pot_tot(:) = vin(:)
   atom%etot = sum(atom%eigs * atom%occs) + atom%exc - 0.5_dp * atom%ehartree
   do i = 1, n_states
      atom%chg_dens(:, i) = atom%occs(i) * (atom%wfs(:, i) / mesh%r(:))**2
      atom%peaks(i) = find_outermost_peak(mesh%r, atom%wfs(:, i), mesh%len, in_out_matches(i))
   end do  ! n_states

   return
end subroutine sratom_structured

function find_outermost_peak(r, wf, len, match_point) result(r_peak)
   ! Input variables
   real(dp), intent(in) :: r(len)
   real(dp), intent(in) :: wf(len)
   integer, intent(in) :: len
   integer, intent(in) :: match_point
   ! Output variable
   real(dp) :: r_peak

   integer :: i

   r_peak = 0.0_dp
   do i = match_point - 1, 1, -1
      if (wf(i) * wf(i + 1) < 0.0_dp) then
         r_peak = r(i)
         exit
      end if
   end do

   return
end function find_outermost_peak

subroutine mix_anderson(vin, vout, vin_old, vout_old, beta, r, len, first_iter)
   ! Input variables
   real(dp), intent(in) :: beta
   real(dp), intent(in) :: r(len)
   integer, intent(in) :: len
   logical, intent(in) :: first_iter

   ! Input/Output variables
   real(dp), intent(in out) :: vin(len)
   real(dp), intent(in out) :: vout(len)
   real(dp), intent(in out) :: vin_old(len)
   real(dp), intent(in out) :: vout_old(len)

   ! Local variables
   integer :: i
   real(dp) :: thl
   real(dp) :: sn
   real(dp) :: sd
   real(dp) :: rl
   real(dp) :: rl1
   real(dp) :: dr
   real(dp) :: vtemp

   thl = 0.0_dp
   if (.not. first_iter) then
      sn = 0.0_dp
      sd = 0.0_dp
      do i = 1, len
         rl = vout(i) - vin(i)
         rl1 = vout_old(i) - vin_old(i)
         dr = rl - rl1
         sn = sn + rl * dr * r(i)**2
         sd = sd + dr * dr * r(i)**2
      end do
      thl = sn / sd
   end if

   do i = 1, len
      vtemp = (1.0_dp - beta) * ((1.0_dp - thl) * vin(i) + thl * vin_old(i)) &
         + beta * ((1.0_dp - thl) * vout(i) + thl * vout_old(i))
      vin_old(i) = vin(i)
      vout_old(i) = vout(i)
      vin(i) = vtemp
      vout(i) = 0.0_dp
   end do

   return
end subroutine mix_anderson

end module sratom_m

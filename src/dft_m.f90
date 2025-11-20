module dft_m
   use thomas_fermi_m, only: tfapot
   use schroedinger_m, only: lschfb, lschvkbb
   use dirac_m, only: ldiracfb
   use vout_m, only: vout
   implicit none
   private
   public :: sratom
   public :: relatom
   public :: psatom
   public :: psatom_r
contains
!> self-consistent scalar-relativistic all-electron atom
!> calculation using log mesh (non-relativistic when srel=.false.)
subroutine sratom(na, la, ea, fa, rpk, nc, ncv, it, rhoc, rho, &
&           rr, vi, zz, mmax, iexc, etot, ierr, srel)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> mmax  size of log grid
   integer, intent(in) :: mmax
   !> iexc  exchange-correlation function to be used
   integer, intent(in) :: iexc
   !> nc  number of core states
   integer, intent(in) :: nc
   !> ncv  number of core+valence states
   integer, intent(in) :: ncv
   !> na  principal quantum number array, dimension ncv
   integer, intent(in) :: na(ncv)
   !> la  angular-momenta
   integer, intent(in) :: la(ncv)
   !> zz  atomic number
   real(dp), intent(in) :: zz
   !> fa  occupancies
   real(dp), intent(in) :: fa(ncv)
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> srel  .true. for scalar-relativistic, .false. for non-relativistic
   logical, intent(in) :: srel

   !Output variables
   !> it  number of iterations (output)
   integer, intent(out) :: it
   !> ierr  error flag
   integer, intent(out) :: ierr
   !> etot  all-electron total energy (output)
   real(dp), intent(out) :: etot
   !> ea  eigenvalues (output)
   real(dp), intent(out) :: ea(ncv)
   !> rpk  radius of outermost peak of wave function
   real(dp), intent(out) :: rpk(ncv)
   real(dp), intent(out) :: rho(mmax)
   real(dp), intent(out) :: rhoc(mmax)
   !> vi  all-electron potential (output)
   real(dp), intent(out) :: vi(mmax)

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

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)

   do ii = 1, mmax
      vi(ii) = tfapot(rr(ii), zz)
   end do

   ! starting approximation for energies
   sf = 0.0d0
   do ii = 1, ncv
      sf = sf + fa(ii)
      zion = zz + 1.0d0 - sf
      ea(ii) = -0.5d0*(zion/na(ii))**2
      if (ea(ii) > vi(mmax)) ea(ii) = 2.0d0*vi(mmax)
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
         call lschfb(na(ii), la(ii), ierr, et, &
         &                rr, vi, u, up, zz, mmax, mch, srel)
         if (ierr .ne. 0) then
            write (6, '(/a,3i4)') 'sratom123: lschfb convergence ERROR n,l,iter=', &
            &       na(ii), la(ii), it
            stop
         end if

         ! overall convergence criterion based on eps within lschfb
         if (ea(ii) /= et) convg = .false.
         ea(ii) = et

         ! accumulate charge and eigenvalues
         eeig = eeig + fa(ii)*ea(ii)
         rho(:) = rho(:) + fa(ii)*(u(:)/rr(:))**2
         if (ii <= nc) then
            rhoc(:) = rhoc(:) + fa(ii)*(u(:)/rr(:))**2
         end if

         ! find outermost peak of wavefunction
         do jj = mch - 1, 1, -1
            if (up(jj)*up(jj + 1) < 0.0d0) then
               rpk(ii) = rr(jj)
               exit
            end if
         end do

      end do

      if (ierr /= 0) then
         exit
      end if

      ! output potential
      call vout(0, rho, rhoc, vo, vxc, sf - zz, eeel, eexc, &
      &            rr, mmax, iexc)

      etot = eeig + eexc - 0.5d0*eeel

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
            sn = sn + rl*dr*rr(ii)**2
            sd = sd + dr*dr*rr(ii)**2
         end do
         thl = sn/sd
      end if

      do ii = 1, mmax
         vn = (1.0d0 - bl)*((1.0d0 - thl)*vi(ii) + thl*vi1(ii)) &
         &   + bl*((1.0d0 - thl)*vo(ii) + thl*vo1(ii))
         vi1(ii) = vi(ii)
         vo1(ii) = vo(ii)
         vi(ii) = vn
      end do

      if (convg) exit

      if (it == 100 .and. .not. convg) then
         write (6, '(/a)') 'sratom: WARNING failed to converge'
      end if

   end do !it

   if (.not. convg .and. ierr == 0) then
      ierr = 100
   end if

   ! total energy output

   ! output potential for e-e interactions

   call vout(0, rho, rhoc, vo, vxc, sf, eeel, eexc, &
   &          rr, mmax, iexc)

   etot = eeig + eexc - 0.5d0*eeel

   deallocate (u, up)
   deallocate (vo, vi1, vo1, vxc)
   return

end subroutine sratom
!> self-consistent fully-relativistic all-electron atom
!> partially-occupied orbitals are weighted proportionally to 2j+1
!> calculation using log mesh
subroutine relatom(na, la, ea, fa, rpk, nc, ncv, it, rhoc, rho, &
&           rr, vi, zz, mmax, iexc, etot, ierr)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> mmax  size of log grid
   integer, intent(in) :: mmax
   !> iexc  exchange-correlation function to be used
   integer, intent(in) :: iexc
   !> nc  number of core states
   integer, intent(in) :: nc
   !> ncv  number of core+valence states
   integer, intent(in) :: ncv
   !> na  principal quantum number array, dimension ncv
   integer, intent(in) :: na(ncv)
   !> la  angular-momenta
   integer, intent(in) :: la(ncv)
   !> zz  atomic number
   real(dp), intent(in) :: zz
   !> fa  occupancies
   real(dp), intent(in) :: fa(ncv)
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)

   !Output variables
   !> it  number of iterations (output)
   integer, intent(out) :: it
   !> ierr  error flag
   integer, intent(out) :: ierr
   !> etot  all-electron total energy (output)
   real(dp), intent(out) :: etot
   !> ea  eigenvalues (output)
   real(dp), intent(out) :: ea(30, 2)
   !> rpk  radius of outermost peak of wave function
   real(dp), intent(out) :: rpk(30, 2)
   real(dp), intent(out) :: rho(mmax)
   real(dp), intent(out) :: rhoc(mmax)
   !> vi  all-electron potential (output)
   real(dp), intent(out) :: vi(mmax)

   !Local variables
   integer :: nin, mch
   real(dp) :: amesh, al
   real(dp) :: dr, eeel, eexc, et, rl, rl1, sd, sf, sn, eeig
   real(dp) :: thl, vn, zion, fj
   integer :: ii, jj, ll, ikap, kap, mkap
   logical :: convg

   real(dp), allocatable :: uu(:, :), up(:, :)
   real(dp), allocatable :: vo(:), vi1(:), vo1(:), vxc(:)

   ! blend parameter for Anderson iterative potential mixing
   real(dp), parameter :: bl = 0.5d0

   allocate (uu(mmax, 2), up(mmax, 2))
   allocate (vo(mmax), vi1(mmax), vo1(mmax), vxc(mmax))

   ! why all this is necessary is unclear, but it seems to be
   uu(:, :) = 0.d0; up(:, :) = 0.d0; vo(:) = 0.d0; vi1(:) = 0.d0; vo1(:) = 0.d0; vxc(:) = 0.d0
   dr = 0.d0; eeel = 0.d0; eexc = 0.d0; et = 0.d0; rl = 0.d0; rl1 = 0.d0
   sd = 0.d0; sf = 0.d0; sn = 0.d0; eeig = 0.d0; thl = 0.d0; vn = 0.d0; zion = 0.d0
   nin = 0; mch = 0

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)

   do ii = 1, mmax
      vi(ii) = tfapot(rr(ii), zz)
   end do

   ! starting approximation for energies
   sf = 0.0d0
   do ii = 1, ncv
      sf = sf + fa(ii)
      zion = zz + 1.0d0 - sf
      do ikap = 1, 2
         ea(ii, ikap) = -0.5d0*(zion/na(ii))**2
         if (ea(ii, ikap) > vi(mmax)) ea(ii, ikap) = 2.0d0*vi(mmax)
      end do
      if (la(ii) == 0) ea(ii, 2) = 0.0d0
   end do

   ! big self  self-consietency loop

   do it = 1, 100
      convg = .true.

      rhoc(:) = 0.0d0
      rho(:) = 0.0d0

      ! solve for bound states in turn
      eeig = 0.0d0
      do ii = 1, ncv
         ll = la(ii)
         if (ll == 0) then
            mkap = 1
         else
            mkap = 2
         end if
         ! loop on J = ll +/- 1/2
         do ikap = 1, mkap
            if (ikap == 1) then
               kap = -(ll + 1)
               fj = dble(ll + 1)*fa(ii)/dble(2*ll + 1)
            else
               kap = ll
               fj = dble(ll)*fa(ii)/dble(2*ll + 1)
            end if
            et = ea(ii, ikap)
            ierr = 0

            call ldiracfb(na(ii), ll, kap, ierr, et, &
            &                rr, zz, vi, uu, up, mmax, mch)

            if (ierr .ne. 0) then
               write (6, '(/2a,5i4)') 'relatom: ldiracfb convergence ERROR', &
               &           ' n,l,kap,iter,ierr=', na(ii), ll, kap, it, ierr
               stop
            end if

            ! overall convergence criterion based on eps within lschfb
            if (ea(ii, ikap) /= et) convg = .false.
            ea(ii, ikap) = et

            ! accumulate charge and eigenvalues
            eeig = eeig + fj*ea(ii, ikap)
            rho(:) = rho(:) + fj*((uu(:, 1)/rr(:))**2 + (uu(:, 2)/rr(:))**2)
            if (ii <= nc) then
               rhoc(:) = rhoc(:) + fj*((uu(:, 1)/rr(:))**2 + (uu(:, 2)/rr(:))**2)
            end if

            ! find outermost peak of wavefunction
            do jj = mch - 1, 1, -1
               if (up(jj, 1)*up(jj + 1, 1) < 0.0d0) then
                  rpk(ii, ikap) = rr(jj)
                  exit
               end if
            end do
         end do
      end do

      if (ierr /= 0) then
         exit
      end if

      ! output potential
      call vout(0, rho, rhoc, vo, vxc, sf - zz, eeel, eexc, &
      &            rr, mmax, iexc)

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
            sn = sn + rl*dr*rr(ii)**2
            sd = sd + dr*dr*rr(ii)**2
         end do
         thl = sn/sd
      end if

      do ii = 1, mmax
         vn = (1.0d0 - bl)*((1.0d0 - thl)*vi(ii) + thl*vi1(ii)) &
         &   + bl*((1.0d0 - thl)*vo(ii) + thl*vo1(ii))
         vi1(ii) = vi(ii)
         vo1(ii) = vo(ii)
         vi(ii) = vn
      end do

      if (convg) exit

      if (it == 100 .and. .not. convg) then
         write (6, '(/a)') 'relatom: WARNING failed to converge'
      end if

   end do !it

   if (.not. convg .and. ierr == 0) then
      ierr = 100
   end if

   ! total energy output

   ! output potential for e-e interactions

   call vout(0, rho, rhoc, vo, vxc, sf, eeel, eexc, &
   &          rr, mmax, iexc)

   etot = eeig + eexc - 0.5d0*eeel

   deallocate (uu, up)
   deallocate (vo, vi1, vo1, vxc)
   return

end subroutine relatom
!> self-consistent pseudoatom calculation
subroutine psatom(na, la, ea, fat, nv, it, rhoc, rho, &
&           rr, rcmax, mmax, mxprj, iexc, etot, nproj, vpuns, lloc, vkb, evkb, ierr)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> mmax  size of log grid
   integer :: mmax
   !> mxprj  dimension of number of projectors
   integer :: mxprj
   !> iexc  exchange-correlation function to be used
   integer :: iexc
   !> nv  number of valence states
   integer :: nv
   !> lloc  index-1 of local potential
   integer :: lloc
   !> na  principal quantum number array, dimension nv
   integer :: na(nv)
   !> la  angular-momenta
   integer :: la(nv)
   !> nproj  number of VKB projectora to use for each l
   integer :: nproj(5)
   !> rcmax  maximum core radius for psp
   real(dp) :: rcmax
   !> fa  occupancies
   real(dp) :: fat(30, 2)
   !> rr  log radial mesh
   real(dp) :: rr(mmax)
   !> vpuns  unscreened semi-local pseudopotentials (plus differenv vloc if lloc==4)
   real(dp) :: vpuns(mmax, 5)
   !> vkb   Vanderbilt-Kleinman-Bylander projectors
   real(dp) :: vkb(mmax, mxprj, 4)
   !> evkb VKB projector coefficients
   real(dp) :: evkb(mxprj, 4)

   !Output variables
   !> it  number of iterations (output)
   integer :: it
   !> etot  pseudoatom total energy (output)
   real(dp) :: etot
   !> ea  eigenvalues (input starting guess, output)
   real(dp) :: ea(nv)
   real(dp) :: rho(mmax)
   real(dp) :: rhoc(mmax)
   real(dp) :: vi(mmax)

   !Local variables
   integer :: nin, mch
   real(dp) :: amesh, al
   real(dp) :: dr, eeel, eexc, et, emin, emax, rl, rl1, sd, sn, sls, eeig
   real(dp) :: thl, vn, zval, dfa
   real(dp) :: fa(30)
   integer :: ii, jj, l1, ierr, icx, nprj
   logical :: convg

   real(dp), allocatable :: uu(:), up(:)
   real(dp), allocatable :: vo(:), vi1(:), vo1(:), vxc(:), vp(:)
   real(dp), allocatable :: vtot(:)

   ! blend parameter for Anderson iterative potential mixing
   real(dp), parameter :: bl = 0.5d0

   allocate (uu(mmax), up(mmax))
   allocate (vo(mmax), vi1(mmax), vo1(mmax), vxc(mmax), vp(mmax))
   allocate (vtot(mmax))

   ! this seems necessary in sratom, so I might as well do it  here too
   ! don't ask why!
   uu(:) = 0.d0; up(:) = 0.d0
   vo(:) = 0.d0; vi1(:) = 0.d0; vo1(:) = 0.d0; vxc(:) = 0.d0; vp(:) = 0.d0
   vtot(:) = 0.d0
   dr = 0.d0; eeel = 0.d0; eexc = 0.d0; et = 0.d0; emin = 0.d0; emax = 0.d0;
   rl = 0.d0; rl1 = 0.d0; sd = 0.d0; sn = 0.d0; sls = 0.d0; eeig = 0.d0;
   thl = 0.d0; vn = 0.d0; zval = 0.d0;
   nin = 0; mch = 0

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)

   nprj = 0
   do l1 = 1, 4
      nprj = max(nprj, nproj(l1))
   end do

   ! total valence charge, initial config
   zval = 0.0d0
   do ii = 1, nv
      fa(ii) = fat(ii, 1)
      zval = zval + fa(ii)
   end do

   ! test if new states are added
   icx = 50
   !do ii=1,nv
   !  if(fat(ii,2)>fat(ii,1)) icx=100
   !end do

   ! screening potential for pseudocharge
   ! input rho is assumed to be valence rho from all-electron calculation
   ! which shouldn''t be a bad approximation for a starting screened
   ! potential

   call vout(1, rho, rhoc, vi, vxc, zval, eeel, eexc, rr, mmax, iexc)

   ! big self  self-consietency loop
   do it = 1, 100

      ! convergence is only to be considered if we are fully in  the target
      ! configurtion

      if (it > 2) then
         convg = .true.
      else
         convg = .false.
      end if

      ! solve for bound states in turn

      eeig = 0.0d0

      vtot(:) = vpuns(:, lloc + 1) + vi(:)

      rho(:) = 0.0d0

      do ii = 1, nv

         if (fa(ii) == 0.0d0) then
            cycle
         end if
         et = ea(ii)
         ierr = 0
         l1 = la(ii) + 1

         emin = 1.10d0*ea(ii)
         emax = 0.90d0*ea(ii)
         call lschvkbb(na(ii), la(ii), nproj(l1), ierr, et, emin, emax, &
         &                 rr, vtot, vkb(1, 1, l1), evkb(1, l1), uu, up, mmax, mch)

         if (ierr .ne. 0) then
            write (6, '(/a,3i4)') 'psatom: WARNING lschvkbb convergence error n,l,iter=', &
            &       na(ii), la(ii), it
            write (6, '(a,1p,4e14.6)') 'ea(ii),et,emin,emax', ea(ii), et, emin, emax
            exit
         end if

         ! overall convergence criterion based on eps within lschf*
         if (ea(ii) /= et) convg = .false.
         ea(ii) = et

         ! accumulate charge and eigenvalues
         eeig = eeig + fa(ii)*ea(ii)
         rho(:) = rho(:) + fa(ii)*(uu(:)/rr(:))**2

      end do !ii

      if (ierr /= 0) exit

      ! total valence charge
      zval = 0.0d0
      do ii = 1, nv
         zval = zval + fa(ii)
      end do

      ! output potential
      call vout(1, rho, rhoc, vo, vxc, zval, eeel, eexc, rr, mmax, iexc)

      ! generate next iteration using d. g. anderson''s
      ! method
      thl = 0.0d0
      if (it > icx + 1) then
         sn = 0.0d0
         sd = 0.0d0
         do ii = 1, mmax
            rl = vo(ii) - vi(ii)
            rl1 = vo1(ii) - vi1(ii)
            dr = rl - rl1
            sn = sn + rl*dr*rr(ii)**2
            sd = sd + dr*dr*rr(ii)**2
         end do
         thl = sn/sd
      end if

      do ii = 1, mmax
         vn = (1.0d0 - bl)*((1.0d0 - thl)*vi(ii) + thl*vi1(ii)) &
         &   + bl*((1.0d0 - thl)*vo(ii) + thl*vo1(ii))
         vi1(ii) = vi(ii)
         vo1(ii) = vo(ii)
         vi(ii) = vn
      end do

      if (convg) exit

      if (it == 400 .and. .not. convg) then
         write (6, '(/a)') 'psatom: potential failed to converge'
         ierr = 101
      end if

      ! switch from reference to new configuration
      ! new or increased occupancy added in second stage
      do ii = 1, nv
         dfa = 0.02d0*(fat(ii, 2) - fat(ii, 1))
         if (it <= 50) then
            fa(ii) = 0.02d0*(50 - it)*fat(ii, 1) + 0.02d0*it*fat(ii, 2)
         else
            fa(ii) = fat(ii, 2)
         end if
      end do

   end do !it

   ! total energy output

   ! output potential for e-e interactions

   call vout(1, rho, rhoc, vo, vxc, zval, eeel, eexc, &
   &          rr, mmax, iexc)

   etot = eeig + eexc - 0.5d0*eeel

   deallocate (uu, up)
   deallocate (vo, vi1, vo1, vxc, vp)
   deallocate (vtot)
   return

end subroutine psatom
!> self-consistent pseudoatom calculation
subroutine psatom_r(na, la, ea, fat, nv, it, rhoc, rho, &
&           rr, rcmax, mmax, mxprj, iexc, etot, nproj, vpuns, lloc, vkb, evkb, ierr)

   !> rpk  radius of outermost peak of wave function

   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> mmax  size of log grid
   integer, intent(in) :: mmax
   !> mxprj  dimension of number of projectors
   integer, intent(in) :: mxprj
   !> iexc  exchange-correlation function to be used
   integer, intent(in) :: iexc
   !> nv  number of valence states
   integer, intent(in) :: nv
   !> lloc  index-1 of local potential
   integer, intent(in) :: lloc
   !> na  principal quantum number array, dimension nv
   integer, intent(in) :: na(30)
   !> la  angular-momenta
   integer, intent(in) :: la(30)
   !> nproj  number of VKB projectora to use for each l
   integer, intent(in) :: nproj(5)
   !> rcmax  maximum core radius for psp
   real(dp), intent(in) :: rcmax
   !> fa  occupancies
   real(dp), intent(in) :: fat(30, 2)
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> vpuns  unscreened semi-local pseudopotentials (plus differenv vloc if lloc==4)
   real(dp), intent(in) :: vpuns(mmax, 5)
   !> vkb   Vanderbilt-Kleinman-Bylander projectors
   real(dp), intent(in) :: vkb(mmax, mxprj, 4, 2)
   !> evkb VKB projector coefficients
   real(dp), intent(in) :: evkb(mxprj, 4, 2)

   !Output variables
   !> it  number of iterations (output)
   integer, intent(out) :: it
   !> etot  pseudoatom total energy (output)
   real(dp), intent(out) :: etot
   !> ea  eigenvalues (input starting guess, output)
   real(dp), intent(out) :: ea(30, 2)
   real(dp), intent(out) :: rho(mmax)
   real(dp), intent(out) :: rhoc(mmax)

   !Local variables
   integer :: nin, mch
   real(dp) :: amesh, al
   real(dp) :: dr, eeel, eexc, et, emin, emax, rl, rl1, sd, sn, sls, eeig
   real(dp) :: thl, vn, zval, dfa, fj
   real(dp) :: fa(30)
   real(dp) :: vi(mmax)
   integer :: ii, jj, l1, ierr, icx, nprj
   integer :: ikap, kap, ll, mkap
   logical :: convg

   real(dp), allocatable :: uu(:), up(:)
   real(dp), allocatable :: vo(:), vi1(:), vo1(:), vxc(:), vp(:)
   real(dp), allocatable :: vtot(:)

   ! blend parameter for Anderson iterative potential mixing
   real(dp), parameter :: bl = 0.5d0

   allocate (uu(mmax), up(mmax))
   allocate (vo(mmax), vi1(mmax), vo1(mmax), vxc(mmax), vp(mmax))
   allocate (vtot(mmax))

   ! this seems necessary in sratom, so I might as well do it  here too
   ! don't ask why!
   uu(:) = 0.d0; up(:) = 0.d0
   vo(:) = 0.d0; vi1(:) = 0.d0; vo1(:) = 0.d0; vxc(:) = 0.d0; vp(:) = 0.d0
   vtot(:) = 0.d0
   dr = 0.d0; eeel = 0.d0; eexc = 0.d0; et = 0.d0; emin = 0.d0; emax = 0.d0;
   rl = 0.d0; rl1 = 0.d0; sd = 0.d0; sn = 0.d0; sls = 0.d0; eeig = 0.d0;
   thl = 0.d0; vn = 0.d0; zval = 0.d0;
   nin = 0; mch = 0

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)

   nprj = 0
   do l1 = 1, 4
      nprj = max(nprj, nproj(l1))
   end do

   ! total valence charge, initial config
   zval = 0.0d0
   do ii = 1, nv
      fa(ii) = fat(ii, 1)
      zval = zval + fa(ii)
   end do

   ! test if new states are added
   icx = 50

   ! screening potential for pseudocharge
   ! input rho is assumed to be valence rho from all-electron calculation
   ! which shouldn''t be a bad approximation for a starting screened
   ! potential

   call vout(1, rho, rhoc, vi, vxc, zval, eeel, eexc, rr, mmax, iexc)

   ! big self  self-consietency loop
   do it = 1, 100

      ! add well do deal with initially unbound states which may arise from
      ! discrepancy between all-electron valence charge used to initiate screening
      ! and self-consistent pseudocharge
      ! start well at 0.5 Ha at infinity and scale down with iterations

      ! convergence is only to be considered if we are fully in  the target
      ! configurtion
      if (it > 2) then
         convg = .true.
      else
         convg = .false.
      end if

      ! solve for bound states in turn

      eeig = 0.0d0

      vtot(:) = vpuns(:, lloc + 1) + vi(:)

      rho(:) = 0.0d0

      do ii = 1, nv

         if (fa(ii) == 0.0d0) then
            cycle
         end if

         ierr = 0
         l1 = la(ii) + 1
         ll = la(ii)
         if (ll == 0) then
            mkap = 1
         else
            mkap = 2
         end if
         ! loop on J = ll +/- 1/2
         do ikap = 1, mkap
            if (ikap == 1) then
               kap = -(ll + 1)
               fj = dble(ll + 1)*fa(ii)/dble(2*ll + 1)
            else
               kap = ll
               fj = dble(ll)*fa(ii)/dble(2*ll + 1)
            end if
            et = ea(ii, ikap)

            emin = 1.10d0*ea(ii, ikap)
            emax = 0.90d0*ea(ii, ikap)
            call lschvkbb(na(ii), la(ii), nproj(l1), ierr, et, emin, emax, &
            &                   rr, vtot, vkb(1, 1, l1, ikap), evkb(1, l1, ikap), uu, up, mmax, mch)

            if (ierr .ne. 0) then
               write (6, '(/a,3i4)') &
               &              'psatom_r: WARNING lschvkbb convergence error n,l,kap,iter=', &
               &               na(ii), la(ii), kap, it
               write (6, '(a,1p,4e14.6)') 'ea(ii),et,emin,emax', ea(ii, ikap), et, emin, emax
               exit
            end if

            !   overall convergence criterion based on eps within lschf*
            if (ea(ii, ikap) /= et) convg = .false.
            ea(ii, ikap) = et

            ! accumulate charge and eigenvalues
            eeig = eeig + fj*ea(ii, ikap)
            rho(:) = rho(:) + fj*(uu(:)/rr(:))**2
         end do !ikap
      end do !ii

      if (ierr /= 0) exit

      ! total valence charge
      zval = 0.0d0
      do ii = 1, nv
         zval = zval + fa(ii)
      end do

      ! output potential
      call vout(1, rho, rhoc, vo, vxc, zval, eeel, eexc, rr, mmax, iexc)

      ! generate next iteration using d. g. anderson''s
      ! method
      thl = 0.0d0
      if (it > icx + 1) then
         sn = 0.0d0
         sd = 0.0d0
         do ii = 1, mmax
            rl = vo(ii) - vi(ii)
            rl1 = vo1(ii) - vi1(ii)
            dr = rl - rl1
            sn = sn + rl*dr*rr(ii)**2
            sd = sd + dr*dr*rr(ii)**2
         end do
         thl = sn/sd
      end if

      do ii = 1, mmax
         vn = (1.0d0 - bl)*((1.0d0 - thl)*vi(ii) + thl*vi1(ii)) &
         &   + bl*((1.0d0 - thl)*vo(ii) + thl*vo1(ii))
         vi1(ii) = vi(ii)
         vo1(ii) = vo(ii)
         vi(ii) = vn
      end do

      if (convg) exit

      if (it == 100 .and. .not. convg) then
         write (6, '(/a)') 'psatom_r: WARNING potential failed to converge'
         ierr = 101
      end if

      ! switch from reference to new configuration
      ! new or increased occupancy added in second stage
      do ii = 1, nv
         dfa = 0.02d0*(fat(ii, 2) - fat(ii, 1))
         if (it <= 50) then
            fa(ii) = 0.02d0*(50 - it)*fat(ii, 1) + 0.02d0*it*fat(ii, 2)
         else
            fa(ii) = fat(ii, 2)
         end if
      end do

   end do !it

   ! total energy output
   ! output potential for e-e interactions

   call vout(1, rho, rhoc, vo, vxc, zval, eeel, eexc, &
   &          rr, mmax, iexc)

   etot = eeig + eexc - 0.5d0*eeel

   deallocate (uu, up)
   deallocate (vo, vi1, vo1, vxc, vp)
   deallocate (vtot)
   return

end subroutine psatom_r
end module dft_m

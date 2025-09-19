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
program oncvpsp
   !
   ! Creates and tests optimized norm-conserving Vanderbilt or Kleinman-Bylander
   ! pseudopotentials based on D. R. Hamann, Phys. Rev. B 88, 085117 (2013)
   ! and references therein.
   !
   !   D. R. Hamann
   !   Mat-Sim Research LLC
   !   P.O. Box 742
   !   Murray Hill, NJ 07974
   !   USA
   !
   !   Developed from original "gncpp" code of March 8,1987
   !
   !   Output format for ABINIT pspcod=8 and upf format for quantumespresso
   !
   use precision_m, only: dp
   use psmlout_m, only: psmlout
   use read_input_m, only: read_input
   use lsch_m, only: lschfb, lschvkbb
   use write_output_hdf5_m, only: write_output_hdf5
   implicit none
   ! Constants
   !> Scalar-relativistic flag
   logical, parameter :: srel = .true.
   !> "Dimension of number of projectors"
   integer, parameter :: mxprj = 5
   !> Dimension of number of states
   integer, parameter :: mxstate = 30
   !> Dimension of number of angular momenta
   integer, parameter :: mxnl = 6
   !> Dimension of number of test configurations
   integer, parameter :: mxcnf = 5
   !> Convergence / tolerance parameter
   real(dp), parameter :: eps = 1.0d-8
   !> Logarithmic mesh spacing
   real(dp), parameter :: amesh = 1.006d0
   !> Logarithmic mesh upper bound
   real(dp), parameter :: rrmax = 45.0d0

   ! Input parameters
   ! [Atom and reference configuration]
   !> Atomic symbol
   character(len=2) :: atsym
   !> Atomic number
   real(dp) :: zz
   !> Number of core states
   integer :: nc
   !> Number of valence states
   integer :: nv
   !> Exchange-correlation functional code
   integer :: iexc
   !> Output pseudopotential file format
   character(len=4) :: psfile
   !> Principal quantum numbers for reference states
   integer :: na(mxstate)
   !> Angular quantum numbers for reference states
   integer :: la(mxstate)
   !> Occupation numbers for reference states
   real(dp) :: fa(mxstate)
   ! [Pseudopotential and optimization]
   !> Maximum angular momentum
   integer :: lmax
   !> Core radii for pseudopotentials
   real(dp) :: rc(mxnl)
   !> Energies at which to construct the pseudopotentials
   real(dp) :: ep(mxnl)
   !> Number of matching constraints for each pseudopotential
   integer :: ncon(mxnl)
   !> Number of basis functions for each pseudopotential
   integer :: nbas(mxnl)
   !> Maximum wave number for each pseudopotential
   real(dp) :: qcut(mxnl)
   ! [Local potential]
   !> Angular momentum channel used for local potential
   integer :: lloc
   !>
   integer :: lpopt
   !> Local potential offset at origin
   real(dp) :: dvloc0
   ! [Vanderbilt-Kleinman-Bylander projectors]
   !> Number of projectors for each angular momentum
   integer :: nproj(mxnl)
   !>
   real(dp) :: debl(mxnl)
   ! [Model core charge]
   !> Model core charge type
   integer :: icmod
   !> Amplitude factor for model core charge
   real(dp) :: fcfact
   !> Width factor for model core charge
   real(dp) :: rcfact
   !> Minimum value for fcfact coarse grid (icmod>=4)
   real(dp) :: fcfact_min
   !> Maximum value for fcfact coarse grid (icmod>=4)
   real(dp) :: fcfact_max
   !> Step size for fcfact coarse grid (icmod>=4)
   real(dp) :: fcfact_step
   !> Minimum value for rcfact coarse grid (icmod>=4)
   real(dp) :: rcfact_min
   !> Maximum value for rcfact coarse grid (icmod>=4)
   real(dp) :: rcfact_max
   !> Step size for rcfact coarse grid (icmod>=4)
   real(dp) :: rcfact_step
   ! [Log derivative analysis]
   !> Minimum energy for logarithmic derivative analysis
   real(dp) :: epsh1
   !> Maximum energy for logarithmic derivative analysis
   real(dp) :: epsh2
   !> Energy spacing for logarithmic derivative analysis
   real(dp) :: depsh
   !> Radius at which logarithmic derivative is evaluated
   real(dp) :: rxpsh
   ! [Output grid]
   !> Linear mesh spacing
   real(dp) :: drl
   !> Maximum linear mesh value
   real(dp) :: rlmax
   ! [Test configurations]
   !> Number of test configurations
   integer :: ncnf
   !> Number of valence states in test configurations
   integer :: nvcnf(mxcnf)
   !> Principal quantum numbers for test configurations
   integer :: nacnf(mxstate, mxcnf)
   !> Angular quantum numbers for test configurations
   integer :: lacnf(mxstate, mxcnf)
   !> Occupation numbers for test configurations
   real(dp) :: facnf(mxstate, mxcnf)

   ! Local variables
   !> Loop index
   integer :: ii
   !> Loop index
   integer :: jj
   !> Loop index
   integer :: kk
   !> Angular momentum quantum number
   integer :: ll
   !> Angular momentum index (l + 1)
   integer :: l1
   !> Error code
   integer :: ierr
   !> I/O status code
   integer :: ios
   !> Line number in input file
   integer :: inline
   !> log(amesh)
   real(dp) :: al
   !> Number of logarithmic mesh points
   integer :: mmax
   !> First logarithmic mesh value
   real(dp) :: rr1
   !> Number of linear mesh points
   integer :: nrl
   !>
   real(dp) :: ea(mxstate)
   !> Logarithmic mesh point of maximum rc
   integer :: irps
   !> Iteration count for all-electron atom solver
   integer :: it
   !> Logarithmic mesh point at which inward and outward solutions are matched
   integer :: mch
   !> "Index of maximum rc including that of vlocal"
   ! Used in: run_vkb
   integer :: nlim
   !> "Index of point at which rho_core is matched"
   ! Used in: psmlout
   integer :: irct
   !> Projector index
   integer :: iprj
   !> Principal quantum number for each projector
   ! Indexed by (projector, angular momentum)
   integer, allocatable :: npa(:, :)
   !> Index of logarithmic mesh point nearest input rc for each l
   integer :: irc(mxnl)
   !>
   integer :: nodes(4)

   !> Electron-electron energy
   real(dp) :: eeel
   !>
   real(dp) :: eeig
   !> Correction to eigenvalue sum from exchange-correlation for total energy calculation
   real(dp) :: eexc
   !>
   real(dp) :: emax
   !>
   real(dp) :: et
   !>
   real(dp) :: emin
   !>
   real(dp) :: rcmax
   !>
   real(dp) :: rct
   !> Ionic charge (sum of core occupation, zz - zval)
   real(dp) :: zion
   !> Valence charge (sum of valence occupation, zz - zion)
   real(dp) :: zval
   !> Total all-electron energy
   real(dp) :: etot
   !>
   real(dp) :: qmsbf(mxnl)
   !>
   real(dp) :: rc0(mxnl)
   !>
   real(dp) :: rpk(mxstate)
   !>
   real(dp) :: epstot
   !>
   real(dp), allocatable :: evkb(:, :)
   !>
   real(dp), allocatable :: cvgplt(:, :, :, :)
   !>
   real(dp), allocatable :: qq(:, :)
   !>
   real(dp), allocatable :: rr(:)
   !> Total charge density of all-electron atom
   real(dp), allocatable :: rho(:)
   !> Charge density of core states in all-electron atom
   real(dp), allocatable :: rhoc(:)
   !>
   real(dp), allocatable :: rhot(:)
   !> Wavefunction r * ψ(r) on logarithmic mesh
   real(dp), allocatable :: uu(:)
   !> Wavefunction derivative d(r * ψ(r))/dr on logarithmic mesh
   real(dp), allocatable :: up(:)
   !>
   real(dp), allocatable :: vp(:, :)
   !> Total all-electron potential
   real(dp), allocatable :: vfull(:)
   !>
   real(dp), allocatable :: vkb(:, :, :)
   !>
   real(dp), allocatable :: pswf(:, :, :)
   !>
   real(dp), allocatable :: vwell(:)
   !>
   real(dp), allocatable :: vpuns(:, :)
   !>
   real(dp), allocatable :: vo(:)
   !>
   real(dp), allocatable :: vxc(:)
   !>
   real(dp), allocatable :: rhomod(:, :)
   !>
   real(dp), allocatable :: rhoae(:, :)
   !>
   real(dp), allocatable :: rhops(:, :)
   !>
   real(dp), allocatable :: rhotae(:)
   !> Pseudo wave functions
   real(dp), allocatable :: uupsa(:, :)
   !>
   real(dp), allocatable :: epa(:, :)
   !>
   real(dp), allocatable :: fpa(:, :)
   !>
   real(dp), allocatable :: uua(:, :)
   !>
   real(dp), allocatable :: upa(:, :)
   !>
   real(dp), allocatable :: vr(:, :, :)
   !> Command line argument
   character(len=1024) :: arg
   !> Input file name
   character(len=1024) :: infile = ''
   !> HDF5 output file name
   character(len=1024) :: h5file = ''
   !> Unit for input file
   integer :: inunit

   write (6, '(a/a//)') &
   &      'ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
   &      'scalar-relativistic version 4.0.1 03/01/2019'

   write (6, '(a/a/a//)') &
   &      'While it is not required under the terms of the GNU GPL, it is',&
   &      'suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)', &
   &      'in any publication utilizing these pseudopotentials.'

   do ii = 1, command_argument_count()
      call get_command_argument(ii, arg)
      select case (arg)
       case ('-i', '--input')
         if (ii + 1 > command_argument_count()) then
            write (6, '(a)') 'Error: --input requires a filename argument'
            stop 1
         end if
         call get_command_argument(ii + 1, infile)
       case ('--hdf5-output')
         if (ii + 1 > command_argument_count()) then
            write (6, '(a)') 'Error: --hdf5-output requires a filename argument'
            stop 1
         end if
         call get_command_argument(ii + 1, h5file)
       case default
         ! Ignore unknown arguments for now
      end select
   end do

   if (trim(infile) == '') then
      inunit = 5
   else
      open (newunit=inunit, file=trim(infile), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write (6, '(a,a)') 'Error: Could not open input file ', trim(infile)
         stop 1
      end if
   end if

   call read_input(inunit, inline, atsym, zz, nc, nv, iexc, psfile, na, la, fa, lmax, rc, ep, &
   &                       ncon, nbas, qcut, lloc, lpopt, dvloc0, nproj, debl, icmod, fcfact, &
   &                       rcfact, fcfact_min, fcfact_max, fcfact_step, rcfact_min, &
   &                       rcfact_max, rcfact_step, epsh1, epsh2, depsh, rxpsh, rlmax, drl, &
   &                       ncnf, nvcnf, nacnf, lacnf, facnf)

   nvcnf(1) = nv
   do ii = 1, nc + nv
      nacnf(ii, 1) = na(ii)
      lacnf(ii, 1) = la(ii)
      facnf(ii, 1) = fa(ii)
   end do

   do jj = 2, ncnf + 1
      do ii = 1, nc
         nacnf(ii, jj) = na(ii)
         lacnf(ii, jj) = la(ii)
         facnf(ii, jj) = fa(ii)
      end do
   end do

   if (icmod == 0) then
      fcfact = 0.0d0
   end if

   call check_data(atsym, zz, fcfact, rcfact, epsh1, epsh2, depsh, rlmax, drl, fa, facnf, &
   &                rc, ep, qcut, debl, nc, nv, iexc, lmax, lloc, lpopt, icmod, &
   &                ncnf, na, la, nvcnf, nacnf, lacnf, ncon, nbas, nproj, psfile)

   nrl = int((rlmax / drl) - 0.5d0) + 1

   !PWSCF wants an even number of mesh pointe
   if (mod(nrl, 2) /= 0) nrl = nrl + 1

   al = log(amesh)
   rr1 = 0.0005d0 / zz
   rr1 = min(rr1, 0.0005d0 / 10)
   mmax = int(log(rrmax / rr1) / al)

   !calculate zion for output
   zion = zz
   do ii = 1, nc
      zion = zion - fa(ii)
   end do

   allocate (rr(mmax))
   allocate (rho(mmax), rhoc(mmax), rhot(mmax))
   allocate (uu(mmax), up(mmax), uupsa(mmax, 30))
   allocate (evkb(mxprj, 4), cvgplt(2, 7, mxprj, 4), qq(mxprj, mxprj))
   allocate (vp(mmax, 5), vfull(mmax), vkb(mmax, mxprj, 4), pswf(mmax, mxprj, 4))
   allocate (vwell(mmax))
   allocate (vpuns(mmax, 5))
   allocate (vo(mmax), vxc(mmax))
   allocate (rhoae(mmax, nv), rhops(mmax, nv), rhotae(mmax))
   allocate (npa(mxprj, 6))
   allocate (epa(mxprj, 6), fpa(mxprj, 6))
   allocate (uua(mmax, mxprj), upa(mmax, mxprj))
   allocate (vr(mmax, mxprj, 6))

   vr(:, :, :) = 0.0d0
   vp(:, :) = 0.0d0
   vkb(:, :, :) = 0.0d0
   epa(:, :) = 0.0d0

   do ii = 1, mmax
      rr(ii) = rr1 * exp(al * (ii - 1))
   end do

   !
   ! full potential atom solution
   !
   call sratom(na, la, ea, fa, rpk, nc, nc + nv, it, rhoc, rho, &
   &              rr, vfull, zz, mmax, iexc, etot, ierr, srel)
   !
   !

   ! Drop digits beyond 5 decimals for input rcs before making any use of them
   do l1 = 1, max(lmax + 1, lloc + 1)
      jj = int(rc(l1) * 10.0d5)
      rc(l1) = jj / 10.0d5
   end do

   rcmax = 0.0d0
   do l1 = 1, lmax + 1
      rcmax = dmax1(rcmax, rc(l1))
   end do
   do l1 = 1, lmax + 1
      if (rc(l1) == 0.0d0) then
         rc(l1) = rcmax
      end if
   end do

   nproj(lloc + 1) = 0
   rc0(:) = rc(:)

   ! output printing (echos input data, with all-electron eigenvalues added)

   write (6, '(a)') '# ATOM AND REFERENCE CONFIGURATION'
   write (6, '(a)') '# atsym  z   nc   nv     iexc    psfile'
   write (6, '(a,a,f6.2,2i5,i8,2a)') '  ', trim(atsym), zz, nc, nv, iexc, &
   &      '      ', psfile
   write (6, '(a/a)') '#', '#   n    l    f        energy (Ha)'
   do ii = 1, nc + nv
      write (6, '(2i5,f8.2,1pe18.7)') na(ii), la(ii), fa(ii), ea(ii)
   end do

   write (6, '(a/a/a)') '#', '# PSEUDOPOTENTIAL AND OPTIMIZATION', '# lmax'
   write (6, '(i5)') lmax
   write (6, '(a/a)') '#', '#   l,   rc,      ep,       ncon, nbas, qcut'
   do l1 = 1, lmax + 1
      write (6, '(i5,2f10.5,2i5,f10.5)') l1 - 1, rc(l1), ep(l1), ncon(l1),&
      &        nbas(l1), qcut(l1)
   end do

   write (6, '(a/a/a,a)') '#', '# LOCAL POTENTIAL', '# lloc, lpopt,  rc(5),', &
   &      '   dvloc0'
   write (6, '(2i5,f10.5,a,f10.5)') lloc, lpopt, rc(5), '   ', dvloc0

   write (6, '(a/a/a)') '#', '# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', &
   &      '# l, nproj, debl'
   do l1 = 1, lmax + 1
      write (6, '(2i5,f10.5)') l1 - 1, nproj(l1), debl(l1)
   end do

   write (6, '(a/a/a)') '#', '# MODEL CORE CHARGE', &
   &      '# icmod, fcfact, rcfact'
   write (6, '(i5,2f10.5)') icmod, fcfact, rcfact

   write (6, '(a/a/a)') '#', '# LOG DERIVATIVE ANALYSIS', &
   &      '# epsh1, epsh2, depsh'
   write (6, '(3f8.2)') epsh1, epsh2, depsh

   write (6, '(a/a/a)') '#', '# OUTPUT GRID', '# rlmax, drl'
   write (6, '(2f8.4)') rlmax, drl

   write (6, '(a/a/a)') '#', '# TEST CONFIGURATIONS', '# ncnf'
   write (6, '(i5)') ncnf
   write (6, '(a/a)') '# nvcnf', '#   n    l    f'
   do jj = 2, ncnf + 1
      write (6, '(i5)') nvcnf(jj)
      do ii = nc + 1, nc + nvcnf(jj)
         write (6, '(2i5,f8.2)') nacnf(ii, jj), lacnf(ii, jj), facnf(ii, jj)
      end do
      write (6, '(a)') '#'
   end do

   write (6, '(//a)') 'Reference configufation results'
   write (6, '(a,i6)') '  iterations', it
   if (it >= 100) then
      write (6, '(a)') 'oncvpsp: ERROR all-electron reference atom not converged'
      stop
   end if
   write (6, '(a,1p,d18.8)') '  all-electron total energy (Ha)', etot

   !find log mesh point nearest input rc
   rcmax = 0.0d0
   irc(:) = 0
   do l1 = 1, max(lmax + 1, lloc + 1)
      rct = rc(l1)
      irc(l1) = 0
      do ii = 2, mmax
         if (rr(ii) > rct) then
            irc(l1) = ii
            rc(l1) = rr(ii)
            exit
         end if
      end do
      rcmax = dmax1(rcmax, rc(l1))
   end do

   !
   cvgplt(:, :, :, :) = 0.0d0
   !
   ! loop to construct pseudopotentials for all angular momenta
   !
   write (6, '(/a/a)') 'Begin loop to  construct optimized pseudo wave functions',&
   &      'and semi-local pseudopoentials for all angular momenta'

   !temporarily set this to 1 so that the pseudo wave function needed for the
   !local potential will be generated.  Reset after run_vkb.
   nproj(lloc + 1) = 1
   do l1 = 1, lmax + 1
      ll = l1 - 1
      uu(:) = 0.0d0; qq(:, :) = 0.0d0
      iprj = 0

      !get principal quantum number for the highest core state for this l
      npa(1, l1) = l1
      do kk = 1, nc
         if (la(kk) == l1 - 1) npa(1, l1) = na(kk) + 1
      end do  !kk

      !get all-electron bound states for projectors
      if (nv /= 0) then
         do kk = nc + 1, nc + nv
            if (la(kk) == l1 - 1) then
               iprj = iprj + 1
               et = ea(kk)
               call lschfb(na(kk), la(kk), ierr, et, &
               &                      rr, vfull, uu, up, zz, mmax, mch, srel)
               if (ierr /= 0) then
                  write (6, '(/a,3i4)') 'oncvpsp-387: lschfb convergence ERROR n,l,iter=', &
                  &           na(ii), la(ii), it
                  stop
               end if
               epa(iprj, l1) = ea(kk)
               npa(iprj, l1) = na(kk)
               uua(:, iprj) = uu(:)
               upa(:, iprj) = up(:)
            end if  !la(kk)==l1-1
            if (iprj == nproj(l1)) exit
         end do  !kk
      end if  !nv/=0

      !get all-electron well states for projectors
      !if there were no valence states, use ep from input data for 1st well state
      !otherwise shift up by input debl
      if (iprj == 0) epa(1, l1) = ep(l1)
      if (iprj < nproj(l1)) then
         do kk = 1, nproj(l1) - iprj
            iprj = iprj + 1
            if (iprj > 1 .and. debl(l1) <= 0.0d0) then
               write (6, '(a,f8.3,a/a)') 'oncvpsp: ERROR debl =', debl, 'for l=', &
               &              ' ERROR not allowed with 2 or more scattering states', &
               &              'program will stop'
               stop
            end if
            if (iprj > 1) then
               epa(iprj, l1) = epa(iprj - 1, l1) + debl(l1)
               npa(iprj, l1) = npa(iprj - 1, l1) + 1
            end if

            call wellstate(npa(iprj, l1), ll, irc(l1), epa(iprj, l1), rr, &
            &                     vfull, uu, up, zz, mmax, mch, srel)
            uua(:, iprj) = uu(:)
            upa(:, iprj) = up(:)
         end do  !kk
      end if  !iprj<nproj(l1)

      do iprj = 1, nproj(l1)

         !calculate relativistic correction to potential to force projectors to 0 at rc
         call vrel(ll, epa(iprj, l1), rr, vfull, vr(1, iprj, l1), uua(1, iprj), upa(1, iprj), &
         &              zz, mmax, irc(l1), srel)

      end do

      !get all-electron overlap matrix
      do jj = 1, nproj(l1)
         do ii = 1, jj
            call fpovlp(uua(1, ii), uua(1, jj), irc(l1), ll, zz, qq(ii, jj), rr, srel)
            qq(jj, ii) = qq(ii, jj)
         end do
      end do

      call run_optimize(epa(1, l1), ll, mmax, mxprj, rr, uua, qq, &
      &                    irc(l1), qcut(l1), qmsbf(l1), ncon(l1), nbas(l1), nproj(l1), &
      &                    pswf(1, 1, l1), vp(1, l1), vkb(1, 1, l1), vfull, cvgplt(1, 1, 1, l1))

   end do  !l1

   ! construct Vanderbilt / Kleinman-Bylander projectors

   write (6, '(/a,a)') 'Construct Vanderbilt / Kleinmman-Bylander projectors'

   call run_vkb(lmax, lloc, lpopt, dvloc0, irc, nproj, rr, mmax, mxprj, pswf, vfull, vp, &
   &             evkb, vkb, nlim, vr)

   !restore this to its proper value
   nproj(lloc + 1) = 0

   deallocate (uua, upa)

   ! accumulate charge and eigenvalues
   ! pseudo wave functions are calculated with VKB projectors for
   ! maximum consistency of unscreening
   ! get all-electron and pseudopotential valence-state by valence-state
   ! charge densities

   ! null charge and eigenvalue accumulators
   uupsa(:, :) = 0.0d0
   eeig = 0.0d0
   zval = 0.0d0
   rho(:) = 0.0d0
   nodes(:) = 0
   rhotae(:) = 0.0d0
   irps = 0
   do kk = 1, nv

      et = ea(nc + kk)
      ll = la(nc + kk)
      l1 = ll + 1
      call lschfb(na(nc + kk), ll, ierr, et, &
      &                rr, vfull, uu, up, zz, mmax, mch, srel)
      if (ierr /= 0) then
         write (6, '(/a,3i4)') 'oncvpsp-461: lschfb convergence ERROR n,l,iter=', &
         &     na(ii), la(ii), it
         stop
      end if

      rhoae(:, kk) = (uu(:) / rr(:))**2

      rhotae(:) = rhotae(:) + fa(nc + kk) * rhoae(:, kk)

      emax = 0.75d0 * et
      emin = 1.25d0 * et

      call lschvkbb(ll + nodes(l1) + 1, ll, nproj(l1), ierr, et, emin, emax, &
      &                rr, vp(1, lloc + 1), vkb(1, 1, l1), evkb(1, l1), &
      &                uu, up, mmax, mch)

      if (ierr /= 0) then
         write (6, '(a,3i4)') 'oncvpsp: lschvkbb ERROR', ll + nodes(l1) + 1, ll, ierr
         flush (6)
         stop
      end if

      ! save valence pseudo wave functions for upfout
      uupsa(:, kk) = uu(:)

      rhops(:, kk) = (uu(:) / rr(:))**2
      rho(:) = rho(:) + fa(nc + kk) * rhops(:, kk)
      eeig = eeig + fa(nc + kk) * et

      zval = zval + fa(nc + kk)
      nodes(l1) = nodes(l1) + 1
      irps = max(irps, irc(l1))
   end do  !kk

   allocate (rhomod(mmax, 5))

   rhomod(:, :) = 0.0d0

   ! construct model core charge based on monotonic polynomial fit
   ! or Teter function fit

   if (icmod == 1) then
      call modcore(rhops, rho, rhoc, rhoae, rhotae, rhomod, &
      &               fcfact, irps, mmax, rr, nc, nv, la, zion, iexc)

   else if (icmod == 2) then
      call modcore2(rhops, rho, rhoc, rhoae, rhotae, rhomod, &
      &               fcfact, mmax, rr, nc, nv, la, zion, iexc)

   else if (icmod >= 3) then
      call modcore3(icmod, rhops, rho, rhoc, rhoae, rhotae, rhomod, &
      &               fcfact, fcfact_min, fcfact_max, fcfact_step, rcfact, rcfact_min, rcfact_max, rcfact_step, &
      &               mmax, rr, nc, nv, la, zion, iexc)

   end if

   ! screening potential for pseudocharge

   call vout(1, rho, rhomod(1, 1), vo, vxc, zval, eeel, eexc, rr, mmax, iexc)

   ! total energy output

   epstot = eeig + eexc - 0.5d0 * eeel
   write (6, '(/a,f12.6/)') 'Pseudoatom total energy', epstot

   call run_diag(lmax, npa, epa, lloc, irc, &
   &                    vkb, evkb, nproj, rr, vfull, vp, zz, mmax, mxprj, srel)

   call run_ghosts(lmax, la, ea, nc, nv, lloc, irc, qmsbf, &
   &                    vkb, evkb, nproj, rr, vp, mmax, mxprj)

   ! unscreen semi-local potentials

   do l1 = 1, max(lmax + 1, lloc + 1)
      vpuns(:, l1) = vp(:, l1) - vo(:)
   end do

   !fix unscreening error due to greater range of all-electron charge
   do ii = mmax, 1, -1
      if (rho(ii) == 0.0d0) then
         do l1 = 1, max(lmax + 1, lloc + 1)
            vpuns(ii, l1) = -zion / rr(ii)
         end do
      else
         exit
      end if
   end do

   ! loop over reference plus test atom configurations
   !if(.false.) then

   !ncnf=0
   rhot(:) = rho(:)
   do jj = 1, ncnf + 1

      write (6, '(/a,i2)') 'Test configuration', jj - 1

      ! charge density is initialized to that of reference configuration

      rhot(:) = rho(:)

      call run_config(jj, nacnf, lacnf, facnf, nc, nvcnf, rhot, rhomod, rr, zz, &
      &                  mmax, mxprj, iexc, ea, etot, epstot, nproj, vpuns, &
      &                  lloc, vkb, evkb, srel)

   end do  !jj

   call run_plot(lmax, npa, epa, lloc, irc, &
   &                    vkb, evkb, nproj, rr, vfull, vp, vpuns, zz, mmax, mxprj, drl, nrl, &
   &                    rho, rhoc, rhomod, srel, cvgplt)

   if (trim(h5file) /= '') then
      call write_output_hdf5(trim(h5file), lmax, npa, epa, lloc, irc, &
                             vkb, evkb, nproj, rr, vfull, vp, vpuns, zz, mmax, mxprj, drl, nrl, &
                             rho, rhoc, rhomod, srel, cvgplt)
   end if

   call run_phsft(lmax, lloc, nproj, epa, epsh1, epsh2, depsh, vkb, evkb, &
   &               rr, vfull, vp, zz, mmax, mxprj, irc, rxpsh, srel)

   call gnu_script(epa, evkb, lmax, lloc, mxprj, nproj)

   if (trim(psfile) == 'psp8' .or. trim(psfile) == 'both') then

      call linout(lmax, lloc, rc, vkb, evkb, nproj, rr, vpuns, rho, rhomod, &
      &             rhotae, rhoc, zz, zion, mmax, mxprj, iexc, icmod, nrl, drl, atsym, &
      &             na, la, ncon, nbas, nvcnf, nacnf, lacnf, nc, nv, lpopt, ncnf, &
      &             fa, rc0, ep, qcut, debl, facnf, dvloc0, fcfact, rcfact, &
      &             epsh1, epsh2, depsh, rlmax, psfile)
   end if

   if (trim(psfile) == 'upf' .or. trim(psfile) == 'both') then
      call upfout(lmax, lloc, rc, vkb, evkb, nproj, rr, vpuns, rho, rhomod, &
      &             zz, zion, mmax, mxprj, iexc, icmod, nrl, drl, atsym, epstot, &
      &             na, la, ncon, nbas, nvcnf, nacnf, lacnf, nc, nv, lpopt, ncnf, &
      &             fa, rc0, ep, qcut, debl, facnf, dvloc0, fcfact, rcfact, &
      &             epsh1, epsh2, depsh, rlmax, psfile, uupsa, ea)
   end if

   if (trim(psfile) == 'psml' .or. trim(psfile) == 'both') then
      !
      ! Write info for PSML format
      !
      print *, 'calling psmlout'
      call psmlout(lmax, lloc, rc, vkb, evkb, nproj, rr, vpuns, rho, rhomod, &
      &             irct, srel, &
      &             zz, zion, mmax, iexc, icmod, drl, atsym, &
      &             na, la, ncon, nbas, nvcnf, nacnf, lacnf, nc, nv, lpopt, ncnf, &
      &             fa, ep, qcut, debl, facnf, dvloc0, fcfact, &
      &             epsh1, epsh2, depsh, rlmax, psfile)
   end if

   stop
end program oncvpsp

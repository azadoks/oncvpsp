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
   use iso_fortran_env, only: stdout => output_unit
   use precision_m, only: dp
   use psmlout_m, only: psmlout
   use read_input_m, only: read_input
   use lsch_m, only: lschfb, lschvkbb
   use model_core_charge_m, only: &
      construct_model_core_charge_1, construct_model_core_charge_2, &
      construct_model_core_charge_3, construct_model_core_charge_4, &
      construct_model_core_charge_5, construct_model_core_charge_6, &
      N_RHOMOD_DERIVATIVES => N_DERIVATIVES
   implicit none
   ! Constants
   !> Scalar-relativistic flag
   logical, parameter :: SREL = .true.
   !> "Dimension of number of projectors"
   integer, parameter :: MXPRJ = 5
   !> Dimension of number of states
   integer, parameter :: MXSTATE = 30
   !> Dimension of number of angular momenta
   integer, parameter :: MXNL = 6
   !> Dimension of number of test configurations
   integer, parameter :: MXCNF = 5
   !> Convergence / tolerance parameter
   real(dp), parameter :: EPS = 1.0e-8_dp
   !> Logarithmic mesh spacing
   real(dp), parameter :: AMESH = 1.006_dp
   !> Logarithmic mesh upper bound
   real(dp), parameter :: RRMAX = 45.0_dp

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
   integer :: na(MXSTATE)
   !> Angular quantum numbers for reference states
   integer :: la(MXSTATE)
   !> Occupation numbers for reference states
   real(dp) :: fa(MXSTATE)
   ! [Pseudopotential and optimization]
   !> Maximum angular momentum
   integer :: lmax
   !> Core radii for pseudopotentials
   real(dp) :: rc(MXNL)
   !> Energies at which to construct the pseudopotentials
   real(dp) :: ep(MXNL)
   !> Number of matching constraints for each pseudopotential
   integer :: ncon(MXNL)
   !> Number of basis functions for each pseudopotential
   integer :: nbas(MXNL)
   !> Maximum wave number for each pseudopotential
   real(dp) :: qcut(MXNL)
   ! [Local potential]
   !> Angular momentum channel used for local potential
   integer :: lloc
   !>
   integer :: lpopt
   !> Local potential offset at origin
   real(dp) :: dvloc0
   ! [Vanderbilt-Kleinman-Bylander projectors]
   !> Number of projectors for each angular momentum
   integer :: nproj(MXNL)
   !>
   real(dp) :: debl(MXNL)
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
   integer :: nvcnf(MXCNF)
   !> Principal quantum numbers for test configurations
   integer :: nacnf(MXSTATE, MXCNF)
   !> Angular quantum numbers for test configurations
   integer :: lacnf(MXSTATE, MXCNF)
   !> Occupation numbers for test configurations
   real(dp) :: facnf(MXSTATE, MXCNF)

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
   !> log(AMESH)
   real(dp) :: al
   !> Number of logarithmic mesh points
   integer :: mmax
   !> First logarithmic mesh value
   real(dp) :: rr1
   !> Number of linear mesh points
   integer :: nrl
   !>
   real(dp) :: ea(MXSTATE)
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
   integer :: irc(MXNL)
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
   real(dp) :: qmsbf(MXNL)
   !>
   real(dp) :: rc0(MXNL)
   !>
   real(dp) :: rpk(MXSTATE)
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
   !> Unit for input file
   integer :: inunit
   !> Error message
   character(len=1024) :: error_msg

   write (stdout, '(a/a//)') &
      'ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
      'scalar-relativistic version 4.0.1 03/01/2019'

   write (stdout, '(a/a/a//)') &
      'While it is not required under the terms of the GNU GPL, it is',&
      'suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)', &
      'in any publication utilizing these pseudopotentials.'

   do ii = 1, command_argument_count()
      call get_command_argument(ii, arg)
      select case (arg)
       case ('-i', '--input')
         if (ii + 1 > command_argument_count()) then
            write (stdout, '(a)') 'Error: --input requires a filename argument'
            stop 1
         end if
         call get_command_argument(ii + 1, infile)
       case default
         ! Ignore unknown arguments for now
      end select
   end do

   if (trim(infile) == '') then
      inunit = 5
   else
      open (newunit=inunit, file=trim(infile), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write (stdout, '(a,a)') 'Error: Could not open input file ', trim(infile)
         stop 1
      end if
   end if

   call read_input(inunit, inline, atsym, zz, nc, nv, iexc, psfile, na, la, fa, lmax, rc, ep, &
                   ncon, nbas, qcut, lloc, lpopt, dvloc0, nproj, debl, icmod, fcfact, &
                   rcfact, fcfact_min, fcfact_max, fcfact_step, rcfact_min, &
                   rcfact_max, rcfact_step, epsh1, epsh2, depsh, rxpsh, rlmax, drl, &
                   ncnf, nvcnf, nacnf, lacnf, facnf)

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
      fcfact = 0.0_dp
   end if

   call check_data(atsym, zz, fcfact, rcfact, epsh1, epsh2, depsh, rlmax, drl, fa, facnf, &
                   rc, ep, qcut, debl, nc, nv, iexc, lmax, lloc, lpopt, icmod, &
                   ncnf, na, la, nvcnf, nacnf, lacnf, ncon, nbas, nproj, psfile)

   nrl = int((rlmax / drl) - 0.5_dp) + 1

   !PWSCF wants an even number of mesh pointe
   if (mod(nrl, 2) /= 0) nrl = nrl + 1

   al = log(AMESH)
   rr1 = 0.0005_dp / zz
   rr1 = min(rr1, 0.0005_dp / 10)
   mmax = int(log(RRMAX / rr1) / al)

   !calculate zion for output
   zion = zz
   do ii = 1, nc
      zion = zion - fa(ii)
   end do

   allocate (rr(mmax))
   allocate (rho(mmax), rhoc(mmax), rhot(mmax))
   allocate (uu(mmax), up(mmax), uupsa(mmax, 30))
   allocate (evkb(MXPRJ, 4), cvgplt(2, 7, MXPRJ, 4), qq(MXPRJ, MXPRJ))
   allocate (vp(mmax, 5), vfull(mmax), vkb(mmax, MXPRJ, 4), pswf(mmax, MXPRJ, 4))
   allocate (vwell(mmax))
   allocate (vpuns(mmax, 5))
   allocate (vo(mmax), vxc(mmax))
   allocate (rhoae(mmax, nv), rhops(mmax, nv), rhotae(mmax))
   allocate (npa(MXPRJ, 6))
   allocate (epa(MXPRJ, 6), fpa(MXPRJ, 6))
   allocate (uua(mmax, MXPRJ), upa(mmax, MXPRJ))
   allocate (vr(mmax, MXPRJ, 6))

   vr(:, :, :) = 0.0_dp
   vp(:, :) = 0.0_dp
   vkb(:, :, :) = 0.0_dp
   epa(:, :) = 0.0_dp

   do ii = 1, mmax
      rr(ii) = rr1 * exp(al * (ii - 1))
   end do

   !
   ! full potential atom solution
   !
   call sratom(na, la, ea, fa, rpk, nc, nc + nv, it, rhoc, rho, rr, vfull, zz, mmax, iexc, etot, ierr, SREL)
   !
   !

   ! Drop digits beyond 5 decimals for input rcs before making any use of them
   do l1 = 1, max(lmax + 1, lloc + 1)
      jj = int(rc(l1) * 10.0e5_dp)
      rc(l1) = jj / 10.0e5_dp
   end do

   rcmax = 0.0_dp
   do l1 = 1, lmax + 1
      rcmax = dmax1(rcmax, rc(l1))
   end do
   do l1 = 1, lmax + 1
      if (rc(l1) == 0.0_dp) then
         rc(l1) = rcmax
      end if
   end do

   nproj(lloc + 1) = 0
   rc0(:) = rc(:)

   ! output printing (echos input data, with all-electron eigenvalues added)

   write (stdout, '(a)') '# ATOM AND REFERENCE CONFIGURATION'
   write (stdout, '(a)') '# atsym  z   nc   nv     iexc    psfile'
   write (stdout, '(a,a,f6.2,2i5,i8,2a)') '  ', trim(atsym), zz, nc, nv, iexc, &
      '      ', psfile
   write (stdout, '(a/a)') '#', '#   n    l    f        energy (Ha)'
   do ii = 1, nc + nv
      write (stdout, '(2i5,f8.2,1pe18.7)') na(ii), la(ii), fa(ii), ea(ii)
   end do

   write (stdout, '(a/a/a)') '#', '# PSEUDOPOTENTIAL AND OPTIMIZATION', '# lmax'
   write (stdout, '(i5)') lmax
   write (stdout, '(a/a)') '#', '#   l,   rc,      ep,       ncon, nbas, qcut'
   do l1 = 1, lmax + 1
      write (stdout, '(i5,2f10.5,2i5,f10.5)') l1 - 1, rc(l1), ep(l1), ncon(l1),&
         nbas(l1), qcut(l1)
   end do

   write (stdout, '(a/a/a,a)') '#', '# LOCAL POTENTIAL', '# lloc, lpopt,  rc(5),', &
      '   dvloc0'
   write (stdout, '(2i5,f10.5,a,f10.5)') lloc, lpopt, rc(5), '   ', dvloc0

   write (stdout, '(a/a/a)') '#', '# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', &
      '# l, nproj, debl'
   do l1 = 1, lmax + 1
      write (stdout, '(2i5,f10.5)') l1 - 1, nproj(l1), debl(l1)
   end do

   write (stdout, '(a/a/a)') '#', '# MODEL CORE CHARGE', &
      '# icmod, fcfact, rcfact'
   write (stdout, '(i5,2f10.5)') icmod, fcfact, rcfact

   write (stdout, '(a/a/a)') '#', '# LOG DERIVATIVE ANALYSIS', &
      '# epsh1, epsh2, depsh'
   write (stdout, '(3f8.2)') epsh1, epsh2, depsh

   write (stdout, '(a/a/a)') '#', '# OUTPUT GRID', '# rlmax, drl'
   write (stdout, '(2f8.4)') rlmax, drl

   write (stdout, '(a/a/a)') '#', '# TEST CONFIGURATIONS', '# ncnf'
   write (stdout, '(i5)') ncnf
   write (stdout, '(a/a)') '# nvcnf', '#   n    l    f'
   do jj = 2, ncnf + 1
      write (stdout, '(i5)') nvcnf(jj)
      do ii = nc + 1, nc + nvcnf(jj)
         write (stdout, '(2i5,f8.2)') nacnf(ii, jj), lacnf(ii, jj), facnf(ii, jj)
      end do
      write (stdout, '(a)') '#'
   end do

   write (stdout, '(//a)') 'Reference configufation results'
   write (stdout, '(a,i6)') '  iterations', it
   if (it >= 100) then
      write (error_msg, '(a)') 'oncvpsp: ERROR all-electron reference atom not converged'
      error stop error_msg
   end if
   write (stdout, '(a,1p,d18.8)') '  all-electron total energy (Ha)', etot

   !find log mesh point nearest input rc
   rcmax = 0.0_dp
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
   cvgplt(:, :, :, :) = 0.0_dp
   !
   ! loop to construct pseudopotentials for all angular momenta
   !
   write (stdout, '(/a/a)') 'Begin loop to  construct optimized pseudo wave functions',&
      'and semi-local pseudopoentials for all angular momenta'

   !temporarily set this to 1 so that the pseudo wave function needed for the
   !local potential will be generated.  Reset after run_vkb.
   nproj(lloc + 1) = 1
   do l1 = 1, lmax + 1
      ll = l1 - 1
      uu(:) = 0.0_dp;
      qq(:, :) = 0.0_dp
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
                           rr, vfull, uu, up, zz, mmax, mch, SREL)
               if (ierr /= 0) then
                  write (error_msg, '(/a,3i4)') 'oncvpsp-387: lschfb convergence ERROR n,l,iter=', na(ii), la(ii), it
                  error stop error_msg
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
            if (iprj > 1 .and. debl(l1) <= 0.0_dp) then
               write (error_msg, '(a,f8.3,a/a)') 'oncvpsp: ERROR debl =', debl, 'for l=', &
                  ' ERROR not allowed with 2 or more scattering states', &
                  'program will stop'
               error stop error_msg
            end if
            if (iprj > 1) then
               epa(iprj, l1) = epa(iprj - 1, l1) + debl(l1)
               npa(iprj, l1) = npa(iprj - 1, l1) + 1
            end if

            call wellstate(npa(iprj, l1), ll, irc(l1), epa(iprj, l1), rr, vfull, uu, up, zz, mmax, mch, SREL)
            uua(:, iprj) = uu(:)
            upa(:, iprj) = up(:)
         end do  !kk
      end if  !iprj<nproj(l1)

      do iprj = 1, nproj(l1)

         !calculate relativistic correction to potential to force projectors to 0 at rc
         call vrel(ll, epa(iprj, l1), rr, vfull, vr(1, iprj, l1), uua(1, iprj), upa(1, iprj), zz, mmax, irc(l1), SREL)

      end do

      !get all-electron overlap matrix
      do jj = 1, nproj(l1)
         do ii = 1, jj
            call fpovlp(uua(1, ii), uua(1, jj), irc(l1), ll, zz, qq(ii, jj), rr, SREL)
            qq(jj, ii) = qq(ii, jj)
         end do
      end do

      call run_optimize(epa(1, l1), ll, mmax, MXPRJ, rr, uua, qq, irc(l1), qcut(l1), qmsbf(l1), ncon(l1), nbas(l1), &
                        nproj(l1), pswf(1, 1, l1), vp(1, l1), vkb(1, 1, l1), vfull, cvgplt(1, 1, 1, l1))

   end do  !l1

   ! construct Vanderbilt / Kleinman-Bylander projectors

   write (stdout, '(/a,a)') 'Construct Vanderbilt / Kleinmman-Bylander projectors'

   call run_vkb(lmax, lloc, lpopt, dvloc0, irc, nproj, rr, mmax, MXPRJ, pswf, vfull, vp, &
                evkb, vkb, nlim, vr)

   !restore this to its proper value
   nproj(lloc + 1) = 0

   deallocate (uua, upa)

   ! accumulate charge and eigenvalues
   ! pseudo wave functions are calculated with VKB projectors for
   ! maximum consistency of unscreening
   ! get all-electron and pseudopotential valence-state by valence-state
   ! charge densities

   ! null charge and eigenvalue accumulators
   uupsa(:, :) = 0.0_dp
   eeig = 0.0_dp
   zval = 0.0_dp
   rho(:) = 0.0_dp
   nodes(:) = 0
   rhotae(:) = 0.0_dp
   irps = 0
   do kk = 1, nv
      et = ea(nc + kk)
      ll = la(nc + kk)
      l1 = ll + 1
      call lschfb(na(nc + kk), ll, ierr, et, &
                  rr, vfull, uu, up, zz, mmax, mch, SREL)
      if (ierr /= 0) then
         write (error_msg, '(/a,3i4)') 'oncvpsp: lschfb convergence ERROR n,l,iter=', na(ii), la(ii), it
         error stop error_msg
      end if

      rhoae(:, kk) = (uu(:) / rr(:))**2

      rhotae(:) = rhotae(:) + fa(nc + kk) * rhoae(:, kk)

      emax = 0.75d0 * et
      emin = 1.25d0 * et

      call lschvkbb(ll + nodes(l1) + 1, ll, nproj(l1), ierr, et, emin, emax, &
                    rr, vp(1, lloc + 1), vkb(1, 1, l1), evkb(1, l1), &
                    uu, up, mmax, mch)

      if (ierr /= 0) then
         write (error_msg, '(a,3i4)') 'oncvpsp: lschvkbb ERROR', ll + nodes(l1) + 1, ll, ierr
         error stop error_msg
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

   allocate (rhomod(mmax, N_RHOMOD_DERIVATIVES))
   select case (icmod)
    case (0)
      rhomod = 0.0_dp
    case (1)
      call construct_model_core_charge_1(rhops, rho, rhoc, rhoae, rhotae, rhomod, &
                                         fcfact, irps, mmax, rr, nc, nv, la, zion, iexc)
    case (2)
      call construct_model_core_charge_2(rhops, rho, rhoc, rhoae, rhotae, rhomod, &
                                         fcfact, mmax, rr, nc, nv, la, zion, iexc)
    case (3)
      call construct_model_core_charge_3(icmod, rhops, rho, rhoc, rhoae, rhotae, rhomod, &
                                         fcfact, rcfact, mmax, rr, nc, nv, la, zion, iexc)
    case (4)
      call construct_model_core_charge_4(icmod, rhops, rho, rhoc, rhoae, rhotae, rhomod, &
                                         fcfact_min, fcfact_max, fcfact_step, &
                                         rcfact_min, rcfact_max, rcfact_step, &
                                         mmax, rr, nc, nv, la, zion, iexc)
    case (5)
      call construct_model_core_charge_5(icmod, rhops, rho, rhoc, rhoae, rhotae, rhomod, &
                                         fcfact_min, fcfact_max, fcfact_step, &
                                         rcfact_min, rcfact_max, rcfact_step, &
                                         mmax, rr, nc, nv, la, zion, iexc)
    case default
      write (error_msg, '(a,i4)') 'oncvpsp: ERROR icmod =', icmod, ' not implemented'
      error stop error_msg
   end select

   ! screening potential for pseudocharge
   call vout(1, rho, rhomod(:, 1), vo, vxc, zval, eeel, eexc, rr, mmax, iexc)

   ! total energy output
   epstot = eeig + eexc - 0.5_dp * eeel
   write (stdout, '(/a,f12.6/)') 'Pseudoatom total energy', epstot

   call run_diag(lmax, npa, epa, lloc, irc, vkb, evkb, nproj, rr, vfull, vp, zz, mmax, MXPRJ, SREL)
   call run_ghosts(lmax, la, ea, nc, nv, lloc, irc, qmsbf, vkb, evkb, nproj, rr, vp, mmax, MXPRJ)

   ! unscreen semi-local potentials
   do l1 = 1, max(lmax + 1, lloc + 1)
      vpuns(:, l1) = vp(:, l1) - vo(:)
   end do

   !fix unscreening error due to greater range of all-electron charge
   do ii = mmax, 1, -1
      if (rho(ii) == 0.0_dp) then
         do l1 = 1, max(lmax + 1, lloc + 1)
            vpuns(ii, l1) = -zion / rr(ii)
         end do
      else
         exit
      end if
   end do

   ! loop over reference plus test atom configurations
   rhot(:) = rho(:)
   do jj = 1, ncnf + 1
      write (stdout, '(/a,i2)') 'Test configuration', jj - 1
      ! charge density is initialized to that of reference configuration
      rhot(:) = rho(:)
      call run_config(jj, nacnf, lacnf, facnf, nc, nvcnf, rhot, rhomod, rr, zz, &
                      mmax, MXPRJ, iexc, ea, etot, epstot, nproj, vpuns, &
                      lloc, vkb, evkb, SREL)
   end do  ! jj

   call run_plot(lmax, npa, epa, lloc, irc, &
                 vkb, evkb, nproj, rr, vfull, vp, vpuns, zz, mmax, MXPRJ, drl, nrl, &
                 rho, rhoc, rhomod, SREL, cvgplt)

   call run_phsft(lmax, lloc, nproj, epa, epsh1, epsh2, depsh, vkb, evkb, &
                  rr, vfull, vp, zz, mmax, MXPRJ, irc, rxpsh, SREL)

   call gnu_script(epa, evkb, lmax, lloc, MXPRJ, nproj)

   select case (trim(psfile))
    case ('psp8')
      call linout(lmax, lloc, rc, vkb, evkb, nproj, rr, vpuns, rho, rhomod, &
                  rhotae, rhoc, zz, zion, mmax, MXPRJ, iexc, icmod, nrl, drl, atsym, &
                  na, la, ncon, nbas, nvcnf, nacnf, lacnf, nc, nv, lpopt, ncnf, &
                  fa, rc0, ep, qcut, debl, facnf, dvloc0, fcfact, rcfact, &
                  epsh1, epsh2, depsh, rlmax, psfile)
    case ('upf')
      call upfout(lmax, lloc, rc, vkb, evkb, nproj, rr, vpuns, rho, rhomod, &
                  zz, zion, mmax, MXPRJ, iexc, icmod, nrl, drl, atsym, epstot, &
                  na, la, ncon, nbas, nvcnf, nacnf, lacnf, nc, nv, lpopt, ncnf, &
                  fa, rc0, ep, qcut, debl, facnf, dvloc0, fcfact, rcfact, &
                  epsh1, epsh2, depsh, rlmax, psfile, uupsa, ea)
    case ('psml')
      call psmlout(lmax, lloc, rc, vkb, evkb, nproj, rr, vpuns, rho, rhomod, &
                   irct, SREL, &
                   zz, zion, mmax, iexc, icmod, drl, atsym, &
                   na, la, ncon, nbas, nvcnf, nacnf, lacnf, nc, nv, lpopt, ncnf, &
                   fa, ep, qcut, debl, facnf, dvloc0, fcfact, &
                   epsh1, epsh2, depsh, rlmax, psfile)
    case ('both')
      call linout(lmax, lloc, rc, vkb, evkb, nproj, rr, vpuns, rho, rhomod, &
                  rhotae, rhoc, zz, zion, mmax, MXPRJ, iexc, icmod, nrl, drl, atsym, &
                  na, la, ncon, nbas, nvcnf, nacnf, lacnf, nc, nv, lpopt, ncnf, &
                  fa, rc0, ep, qcut, debl, facnf, dvloc0, fcfact, rcfact, &
                  epsh1, epsh2, depsh, rlmax, psfile)
      call upfout(lmax, lloc, rc, vkb, evkb, nproj, rr, vpuns, rho, rhomod, &
                  zz, zion, mmax, MXPRJ, iexc, icmod, nrl, drl, atsym, epstot,&
                  na, la, ncon, nbas, nvcnf, nacnf, lacnf, nc, nv, lpopt, ncnf,&
                  fa, rc0, ep, qcut,debl,facnf,dvloc0 ,fcfact ,rcfact,&
                  epsh1 ,epsh2 ,depsh ,rlmax ,psfile,uupsa ,ea)
      call psmlout(lmax,lloc ,rc ,vkb ,evkb ,nproj ,rr ,vpuns ,rho ,rhomod ,&
                   irct,SREL ,&
                   zz,zion ,mmax ,iexc ,icmod ,drl ,atsym ,&
                   na ,la ,ncon ,nbas ,nvcnf ,nacnf ,lacnf ,nc ,nv ,lpopt ,ncnf ,&
                   fa ,ep ,qcut ,debl ,facnf ,dvloc0 ,fcfact ,&
                   epsh1 ,epsh2 ,depsh ,rlmax ,psfile)
    case default
      write (error_msg, '(a,a,a)') 'oncvpsp: ERROR psfile =', trim(psfile), ' not recognized'
      error stop error_msg
   end select
   stop
end program oncvpsp

module model_core_charge_m
   use iso_fortran_env, only: stdout => output_unit
   use precision_m, only: dp
   use constants_m, only: pi, twopi
   use interpolate_m, only: interpolate
   implicit none
   private
   public :: construct_polynomial_model_core_charge
   public :: construct_model_core_charge_2
   public :: construct_model_core_charge_3
   public :: construct_model_core_charge_4
   public :: construct_model_core_charge_5
   public :: construct_model_core_charge_6

   !> Number of model core charge derivatives to calculate
   integer, parameter, public :: N_DERIVATIVES = 4
   !> convergence criterion
   real(dp), parameter :: EPS = 1.0e-7_dp
contains

!> Constructs a model core charge density by fitting a monotonically-decaying 5th-order
!> polynomial which matches the all-electron core charge density and its first four
!> derivatives at the point where the all-electron core charge density crosses the
!> pseudo valence charge density scaled by a factor of `fcfact`
!> (i.e. at `r` where `rhoc(r) = fcfact * rhotps(r)`).
subroutine construct_polynomial_model_core_charge(rhops, rhotps, rhoc, rhoae, rhotae, &
                                                  rhomod, fcfact, irps, mmax, rr, nc, nv, &
                                                  la, zion, iexc)
   implicit none
   ! Constants
   !> Dimension of polynomial arrays
   integer, parameter :: POLY_DIM = N_DERIVATIVES + 1
   !> Maximum number of iterations for polynomial fitting
   integer, parameter :: MAX_ITER = 50

   ! Input variables
   !> Dimension of log grid
   integer, intent(in) :: mmax
   !> Number of core states
   integer, intent(in) :: nc
   !> Number of valence states
   integer, intent(in) :: nv
   !> State-by-state pseudocharge density
   real(dp), intent(in) :: rhops(mmax, nv)
   !> Total pseudocharge density
   real(dp), intent(in) :: rhotps(mmax)
   !> Total all-electron core charge density
   real(dp), intent(in) :: rhoc(mmax)
   !> State-by-state all-electron valence charge density
   real(dp), intent(in) :: rhoae(mmax, nv)
   !> Total all-electron valence charge density
   real(dp), intent(in) :: rhotae(mmax)
   !> Prefactor for model amplitude (multiplies crossover value)
   real(dp), intent(in) :: fcfact
   !> r index of maximum rc
   integer, intent(in) :: irps
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Angular-momenta of the states
   integer, intent(in) :: la(nv + nc)
   !> Ion charge (atomic charge - valence charge)
   real(dp), intent(in) :: zion
   !> Exchange-correlation function to be used
   integer, intent(in) :: iexc

   ! Output variables
   !> Model core density and derivatives
   real(dp), intent(out) :: rhomod(mmax, POLY_DIM)

   ! Local variables
   !> Total potential
   real(dp), allocatable :: vo(:)
   !> Exchange-correlation potentials
   real(dp), allocatable :: vxcae(:), vxcpsp(:)
   !> Exchange-correlation potential derivatives
   real(dp), allocatable :: dvxcae(:, :), dvxcps(:, :)
   !> Second derivatives of exchange-correlation energy
   real(dp), allocatable :: d2excae(:, :), d2excps(:, :)
   !> All-electron core charge density and its derivatives at the crossover point
   real(dp) :: fmatch(POLY_DIM)
   !> Index of core-valence crossover point
   integer :: ircc
   !> Radial coordinate of scaled crossover point
   real(dp) :: rcc
   !> Index of maximum radius for d2Exc calculation
   integer :: ir_d2exc_cut
   !> Polynomial constants
   real(dp) :: constm(POLY_DIM, POLY_DIM)
   !> Polynomial coefficients
   real(dp) :: polym(POLY_DIM, POLY_DIM)
   !> Powers of rcc
   real(dp) :: xpow(2 * POLY_DIM)
   !> Right-hand side vector for linear equation solver
   real(dp) :: aco(POLY_DIM)
   !> Work array for linear equation solver
   real(dp) :: work(POLY_DIM, POLY_DIM)
   !> Pivot indices for linear equation solver
   integer :: ipiv(POLY_DIM)
   !> Maximum derivative found in polynomial fitting
   real(dp) :: dermax
   !> Result of polynomial evaluation
   real(dp) :: psum
   !> ???
   real(dp) :: a0
   !> ???
   real(dp) :: a0min
   !> ???
   real(dp) :: a0max
   !> RMSE in d2Exc
   real(dp) :: d2exc_rmse
   !> Grid spacing in log scale
   real(dp) :: al
   !> Loop indices
   integer :: ii, jj, kk
   !> Iteration index
   integer :: iter
   !> Info flag for linear equation solver
   integer :: info

   allocate(vo(mmax), vxcae(mmax), vxcpsp(mmax))
   allocate(dvxcae(mmax, nv), dvxcps(mmax, nv))
   allocate(d2excae(nv, nv), d2excps(nc, nv))

   rhomod = 0.0_dp
   d2excae = 0.0_dp
   d2excps = 0.0_dp
   ircc = 0

   ! Find the core -- scaled-valence crossover point
   do ii = mmax, 1, -1
      if (rhoc(ii) > fcfact * rhotps(ii)) then
         ircc = ii
         rcc = rr(ircc)
         exit
      end if
   end do
   if (ircc == 0) error stop 'model_core_charge_m: ERROR core--scaled-valence charge crossover not found'
   if (ircc < 20) error stop 'model_core_charge_m: ERROR core--scaled-valence charge crossover too close to origin'

   write (stdout, '(/a/a)') 'Model core correction analysis', '  based on d2Exc/dn_idn_j'

   ! Set limit for d2Exc calculation to radius beyond which rhomod=rhoc by construction
   ir_d2exc_cut = max(irps, ircc)

   ! Compute d2Exc for all-electron atom
   call der2exc(rhotae, rhoc, rhoae, rr, d2excae, zion, iexc, nc, nv, la, ir_d2exc_cut, mmax)
   write (stdout, '(/a/)') 'd2excae - all-electron derivatives'
   do kk = 1, nv
      write (stdout, '(1p,4d16.6)') (d2excae(kk, jj), jj=1, nv)
   end do

   ! Compute d2exc and d2exc_rmse for pseudo atom with no core correction
   call der2exc(rhotps, rhomod(:, 1), rhops, rr, d2excps, zion, iexc, nc, nv, la, ir_d2exc_cut, mmax)
   d2exc_rmse = rmse(d2excps, d2excae, nv)
   write (stdout, '(/a/)') 'd2excps - pseudofunction derivatives with no core correction'
   do kk = 1, nv
      write (stdout, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (stdout, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2exc_rmse

   ! Core charge density for Louie-Froyen-Cohen correction
   rhomod(:, 1) = rhoc(:)
   ! Find the reference values (fmatch) of rhoc and its derivatives at the crossover point
   fmatch(1) = rhoc(ircc)
   ! Compute derivatives of core charge using successive 7-point numerical differentation
   al = 0.010_dp * log(rr(101) / rr(1)) ! grid spacing in log scale
   do jj = 2, 5
      do ii = ircc - 25 + 3 * jj, mmax - 3
         rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) &
                           + 9.0_dp * rhomod(ii - 2, jj - 1) &
                           - 45.0_dp * rhomod(ii - 1, jj - 1) &
                           + 45.0_dp * rhomod(ii + 1, jj - 1) &
                           - 9.0_dp * rhomod(ii + 2, jj - 1) &
                           + rhomod(ii + 3, jj - 1) &
                           ) / (60.0_dp * al * rr(ii))
      end do ! ii
      fmatch(jj) = rhomod(ircc, jj)
   end do ! jj

   ! Initialize polynomial constants
   constm(1, :) = 1.0_dp
   do ii = 1, POLY_DIM
      constm(2, ii) = dble(ii) + 2
   end do
   do ii = 3, POLY_DIM
      do jj = 1, POLY_DIM
         constm(ii, jj) = constm(ii - 1, jj) * (constm(2, jj) - dble(ii) + 2)
      end do
   end do
   ! Initialize powers of rcc
   xpow(1) = 0.0_dp  ! rcc**-inf
   xpow(2) = 1.0_dp  ! rcc**0
   do ii = 3, 2 * POLY_DIM
      xpow(ii) = xpow(ii - 1) * rcc
   end do
   ! Initialize polynomial matrix
   do jj = 1, 5
      do ii = 1, 5
         polym(ii, jj) = constm(ii, jj) * xpow(jj - ii + POLY_DIM)
      end do
   end do

   ! Iteratively fit monotonic polynomial to core charge and its derivatives at crossover point
   a0 = fmatch(1)
   a0min = a0
   a0max = 0.0_dp
   do iter = 1, MAX_ITER
      ! Prepare work array for linear equation solver
      work(:, :) = polym(:, :)
      ! Prepare right-hand side of Ax = b
      aco(1) = fmatch(1) - a0
      aco(2:POLY_DIM) = fmatch(2:POLY_DIM)
      ! Solve linear equations for polynomial coefficients
      call dgesv(POLY_DIM, 1, work, POLY_DIM, ipiv, aco, POLY_DIM, info)
      if (info < 0) error stop 'model_core_charge_m: ERROR dgesv illegal argument'
      if (info > 0) error stop 'model_core_charge_m: ERROR dgesv matrix is singular'
      ! Find maximum derivative
      ! N.B. xpow ends up back as it started at the end of the loop!
      dermax = -huge(1.0_dp)
      do kk = 1, ircc
         do ii = 3, 2 * POLY_DIM
            xpow(ii) = rr(kk) * xpow(ii - 1)
         end do
         psum = 0.0_dp
         do jj = 1, POLY_DIM
            psum = psum + aco(jj) * constm(2, jj) * xpow(jj - 2 + POLY_DIM)
         end do
         dermax = max(dermax, psum)
      end do
      ! Test maximum derivative and adjust a0
      ! Perform interval-halving search after acheiving monotonicity
      ! to get a barely monotonic model core charge
      if (dermax > 0.0_dp) then
         a0min = a0
      else
         a0max = a0
      end if
      ! Success when maximum derivative is bracketed just above zero
      if (dermax > 0.0_dp .and. dermax < EPS * abs(fmatch(2))) exit
      ! Otherwise perform interval halving and try again
      if (abs(a0max) < 1.0e-12_dp) then
         ! a0max not yet set, so just increase a0
         a0 = 2.5_dp * a0
      else
         a0 = 0.5_dp * (a0max + a0min)
      end if
   end do ! iter
   if (iter >= MAX_ITER) write (stdout, '(/a/)') 'WARNING - modcore not converged'

   ! Evaluate the fitted polynomial and its derivatives on the grid points up to ircc
   do kk = 1, ircc - 1
      do ii = 3, 2 * POLY_DIM
         xpow(ii) = rr(kk) * xpow(ii - 1)
      end do
      do ii = 1, POLY_DIM
         if (ii == 1) then
            psum = a0
         else
            psum = 0.0_dp
         end if
         do jj = 1, POLY_DIM
            psum = psum + constm(ii, jj) * xpow(jj - ii + POLY_DIM) * aco(jj)
         end do
         rhomod(kk, ii) = psum
      end do
   end do

   ! Test the fitted model core charge by computing d2Exc
   call der2exc(rhotps, rhomod(:, 1), rhops, rr, d2excps, zion, iexc, nc, nv, la, ir_d2exc_cut, mmax)
   d2exc_rmse = rmse(d2excps, d2excae, nv)
   write (stdout, '(/a/)') 'Polynomial model core charge'
   write (stdout, '(/a/)') 'd2excps - pseudofunction derivatives with core correction'
   do kk = 1, nv
      write (stdout, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (stdout, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2exc_rmse

   deallocate(vo, vxcae, vxcpsp)
   deallocate(dvxcae, dvxcps)
   deallocate(d2excae, d2excps)
end subroutine construct_polynomial_model_core_charge

!> Constructs a model core charge density using a Teter function, where the
!> Teter parameters are are adjusted to match the all-electron core charge density
!> and its first four derivatives at the point where the all-electron core charge
!> density crosses the pseudo valence charge density scaled by a factor of `fcfact`
!> (i.e. at `r` where `rhoc(r) = fcfact * rhotps(r)`).
subroutine construct_model_core_charge_2(rhops, rhotps, rhoc, rhoae, rhotae, &
                                         rhomod, fcfact, mmax, rr, nc, nv, &
                                         la, zion, iexc)
   implicit none
   ! Constants
   ! Input variables
   !> Dimension of log grid
   integer, intent(in) :: mmax
   !> Number of core states
   integer, intent(in) :: nc
   !> Number of valence states
   integer, intent(in) :: nv
   !> State-by-state pseudocharge density
   real(dp), intent(in) :: rhops(mmax, nv)
   !> Total pseudocharge density
   real(dp), intent(in) :: rhotps(mmax)
   !> Total all-electron core charge density
   real(dp), intent(in) :: rhoc(mmax)
   !> State-by-state all-electron valence charge density
   real(dp), intent(in) :: rhoae(mmax, nv)
   !> Total all-electron valence charge density
   real(dp), intent(in) :: rhotae(mmax)
   !> Prefactor for model amplitude (multiplies crossover value)
   real(dp), intent(in) :: fcfact
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Angular-momenta of the states
   integer, intent(in) :: la(nv + nc)
   !> Ion charge (atomic charge - valence charge)
   real(dp), intent(in) :: zion
   !> Exchange-correlation function to be used
   integer, intent(in) :: iexc

   ! Output variables
   !> Model core density and derivatives
   real(dp), intent(out) :: rhomod(mmax, N_DERIVATIVES + 1)

   ! Local variables
   !> Total potential
   real(dp), allocatable :: vo(:)
   !> Exchange-correlation potentials
   real(dp), allocatable :: vxcae(:), vxcpsp(:)
   !> Exchange-correlation potential derivatives
   real(dp), allocatable :: dvxcae(:, :), dvxcps(:, :)
   !> Second derivatives of exchange-correlation energy
   real(dp), allocatable :: d2excae(:, :), d2excps(:, :)
   !> All-electron core charge density and its derivatives at the crossover point
   real(dp) :: fmatch(N_DERIVATIVES + 1)
   !> Index of core-valence crossover point
   integer :: ircc
   !> Radial coordinate of scaled crossover point
   real(dp) :: rcc
   !> Radial coordinate of true crossover point
   real(dp) :: rmatch
   !> Value of all-electron core charge at crossover point
   real(dp) :: rhocmatch
   !> Index of maximum radius for d2Exc calculation
   integer :: ir_d2exc_cut
   !> RMSE in d2Exc
   real(dp) :: d2exc_rmse
   !> Grid spacing in log scale
   real(dp) :: al
   !> Loop indices
   integer :: ii, jj, kk
   !> ???
   real(dp) :: ymatch, x0max, x0min, xx, yy, ytrial, a0, b0
   !> "ad-hoc smoothing variables"
   real(dp) :: drint, rtst, rint(20), fint(20)
   !> "ad-hoc smoothing variables"
   integer :: iint

   allocate(vo(mmax), vxcae(mmax), vxcpsp(mmax))
   allocate(dvxcae(mmax, nv), dvxcps(mmax, nv))
   allocate(d2excae(nv, nv), d2excps(nc, nv))

   rhomod = 0.0_dp
   d2excae = 0.0_dp
   d2excps = 0.0_dp

   write (stdout, '(/a/a)') 'Model core correction analysis', '  based on d2Exc/dn_idn_j'

   ! Find the core -- valence crossover point
   ! "This is needed for compatability with icmod=3"
   do ii = mmax, 1, -1
      if (rhoc(ii) > rhotps(ii)) then
         ircc = ii
         rmatch = rr(ircc)
         rhocmatch = rhoc(ircc)
         exit
      end if
   end do
   if (ircc == 0) error stop 'model_core_charge_m: ERROR core-valence charge crossover not found'
   if (ircc < 20) error stop 'model_core_charge_m: ERROR core-valence charge crossover too close to origin'
   ! Find the core -- scaled-valence crossover point
   do ii = mmax, 1, -1
      if (rhoc(ii) > fcfact * rhotps(ii)) then
         ircc = ii
         rcc = rr(ircc)
         exit
      end if
   end do
   if (ircc == 0) error stop 'model_core_charge_m: ERROR core--scaled-valence charge crossover not found'
   if (ircc < 20) error stop 'model_core_charge_m: ERROR core--scaled-valence charge crossover too close to origin'

   ! Set limit for d2Exc calculation to the end of the radial grid
   ir_d2exc_cut = mmax

   ! Compute d2Exc for all-electron atom
   call der2exc(rhotae, rhoc, rhoae, rr, d2excae, zion, iexc, nc, nv, la, ir_d2exc_cut, mmax)
   write (stdout, '(/a/)') 'd2excae - all-electron derivatives'
   do kk = 1, nv
      write (stdout, '(1p,4d16.6)') (d2excae(kk, jj), jj=1, nv)
   end do

   ! Compute d2exc and d2exc_rmse for pseudo atom with no core correction
   call der2exc(rhotps, rhomod(:, 1), rhops, rr, d2excps, zion, iexc, nc, nv, la, ir_d2exc_cut, mmax)
   d2exc_rmse = rmse(d2excps, d2excae, nv)
   write (stdout, '(/a/)') 'd2excps - pseudofunction derivatives with no core correction'
   do kk = 1, nv
      write (stdout, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (stdout, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2exc_rmse

   ! Core charge density for Louie-Froyen-Cohen correction
   rhomod(:, 1) = rhoc(:)
   ! Find the reference values (fmatch) of rhoc and its derivatives at the crossover point
   fmatch(1) = rhoc(ircc)

   !7-point numerical first derivative
   jj = 2; ii = ircc
   rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) + 9.0_dp * rhomod(ii - 2, jj - 1)&
                     - 45.0_dp * rhomod(ii - 1, jj - 1) + 45.0_dp * rhomod(ii + 1, jj - 1)&
                     - 9.0_dp * rhomod(ii + 2, jj - 1) + rhomod(ii + 3, jj - 1))&
      / (60.0_dp * al * rr(ii))
   fmatch(jj) = rhomod(ircc, jj)

   ! Fit Teter function to value and slope at ircc
   ymatch = rr(ircc) * fmatch(2) / fmatch(1)

   ! interval-halving search for dimensionless match point
   x0max = 1.480_dp
   x0min = 0.00_dp
   do jj = 1, 50
      xx = 0.50_dp * (x0max + x0min)
      ytrial = xx * teter_deriv_1(xx) / teter(xx)
      if (abs(ytrial - ymatch) < eps) exit
      if (ytrial < ymatch) then
         x0max = xx
      else
         x0min = xx
      end if
   end do

   b0 = xx / rr(ircc)
   a0 = fmatch(1) / yy

   write (stdout, '(/a)') 'amplitude prefactor, scale prefactor'
   write (stdout, '(2f10.4)') a0 / rhocmatch, 1.00_dp / (b0 * rmatch)

   rhomod(:, :) = 0.00_dp
   do ii = 1, mmax
      xx = b0 * rr(ii)
      if (xx < 3.00_dp) then
         rhomod(ii, 1) = a0 * teter(xx)
         rhomod(ii, 2) = a0 * teter_deriv_1(xx) * b0
         rhomod(ii, 3) = a0 * teter_deriv_2(xx) * b0**2
      end if
   end do

   ! Test the fitted model core charge by computing d2Exc
   call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, zion, iexc, nc, nv, la, ir_d2exc_cut, mmax)
   d2exc_rmse = rmse(d2excps, d2excae, nv)
   write (stdout, '(/a/)') 'd2excps - pseudofunction derivatives with core correction'
   do kk = 1, nv
      write (stdout, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (stdout, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2exc_rmse

   ! Compute higher derivatives of the model core charge
   ! ad-hoc treatment of numerical noise near origin
   ! set up a mesh on which 2nd derivative will have a stable polynomial representation
   ! assumes interpolate remains 7th order
   drint = 0.05_dp * rr(ircc)
   rtst = 0.5_dp * drint
   do jj = 1, 4
      do ii = 1, ircc
         if (rr(ii) > rtst) then
            rint(4 + jj) = rr(ii)
            rint(5 - jj) = -rr(ii)
            fint(4 + jj) = rhomod(ii, 3)
            fint(5 - jj) = rhomod(ii, 3)
            iint = ii
            rtst = rr(ii) + drint
            exit
         end if
      end do
   end do
   call interpolate(rint, fint, 8, rr, rhomod(1, 3), iint - 1)

   jj = 4
   do ii = 3 * jj - 8, mmax - 3
      rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) + 9.0_dp * rhomod(ii - 2, jj - 1)&
         &     - 45.0_dp * rhomod(ii - 1, jj - 1) + 45.0_dp * rhomod(ii + 1, jj - 1)&
         &     - 9.0_dp * rhomod(ii + 2, jj - 1) + rhomod(ii + 3, jj - 1))&
         &     / (60.0_dp * al * rr(ii))
   end do
   ! set up a mesh on which 3rd derivative will have a stable polynomial representation
   ! assumes interpolate remains 7th order
   drint = 0.95_dp * drint
   rtst = 0.5_dp * drint
   do jj = 1, 4
      do ii = 1, ircc
         if (rr(ii) > rtst) then
            rint(4 + jj) = rr(ii)
            rint(5 - jj) = -rr(ii)
            fint(4 + jj) = rhomod(ii, 4)
            fint(5 - jj) = -rhomod(ii, 4)
            iint = ii
            rtst = rr(ii) + drint
            exit
         end if
      end do
   end do
   call interpolate(rint, fint, 8, rr, rhomod(1, 4), iint - 1)

   jj = 5
   do ii = 3 * jj - 8, mmax - 3
      rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) + 9.0_dp * rhomod(ii - 2, jj - 1)&
         &     - 45.0_dp * rhomod(ii - 1, jj - 1) + 45.0_dp * rhomod(ii + 1, jj - 1)&
         &     - 9.0_dp * rhomod(ii + 2, jj - 1) + rhomod(ii + 3, jj - 1))&
         &     / (60.0_dp * al * rr(ii))
   end do
   ! set up a mesh on which 4th derivative will have a stable polynomial representation
   ! assumes interpolate remains 7th order
   drint = 0.95_dp * drint
   rtst = 0.5_dp * drint
   do jj = 1, 4
      do ii = 1, ircc
         if (rr(ii) > rtst) then
            rint(4 + jj) = rr(ii)
            rint(5 - jj) = -rr(ii)
            fint(4 + jj) = rhomod(ii, 5)
            fint(5 - jj) = rhomod(ii, 5)
            iint = ii
            rtst = rr(ii) + drint
            exit
         end if
      end do
   end do
   call interpolate(rint, fint, 8, rr, rhomod(1, 5), iint - 1)

   deallocate(vo, vxcae, vxcpsp)
   deallocate(dvxcae, dvxcps)
   deallocate(d2excae, d2excps)
   return
end subroutine construct_model_core_charge_2

subroutine construct_model_core_charge_3(rhops, rhotps, rhoc, rhoae, rhotae, &
                                         rhomod, fcfact, rcfact, mmax, rr, nc, nv, &
                                         la, zion, iexc)
   implicit none
   !Input variables
   !> Dimension of log grid
   integer, intent(in) :: mmax
   !> Number of core states
   integer, intent(in) :: nc
   !> Number of valence states
   integer, intent(in) :: nv
   !> State-by-state pseudocharge density
   real(dp), intent(in) :: rhops(mmax, nv)
   !> Total pseudocharge density
   real(dp), intent(in) :: rhotps(mmax)
   !> Total all-electron core charge density
   real(dp), intent(in) :: rhoc(mmax)
   !> State-by-state all-electron valence charge density
   real(dp), intent(in) :: rhoae(mmax, nv)
   !> Total all-electron valence charge density
   real(dp), intent(in) :: rhotae(mmax)
   !> ???
   real(dp), intent(in) :: fcfact
   !> ???
   real(dp), intent(in) :: rcfact
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Angular-momenta of the states
   integer, intent(in) :: la(nv + nc)
   !> Ion charge (atomic charge - valence charge)
   real(dp), intent(in) :: zion
   !> Exchange-correlation function to be used
   integer, intent(in) :: iexc

   !Output variables
   real(dp) :: rhomod(mmax, N_DERIVATIVES + 1)

   !Local variables
   !> Total potential
   real(dp), allocatable :: vo(:)
   !> Exchange-correlation potentials
   real(dp), allocatable :: vxcae(:), vxcpsp(:)
   !> Exchange-correlation potential derivatives
   real(dp), allocatable :: dvxcae(:, :), dvxcps(:, :)
   !> Second derivatives of exchange-correlation energy
   real(dp), allocatable :: d2excae(:, :), d2excps(:, :)
   !> Index of AE core -- PS valence crossover point
   integer :: i_ae_core_ps_val_match
   !> Radius at which AE core and PS valence cross
   real(dp) :: r_ae_core_ps_val_match
   !> Value of AE core charge at AE core -- PS valence crossover point
   real(dp) :: rhoc_ae_core_ps_val_match
   !> Index of AE core -- model core crossover point
   integer :: ir_ae_model_match
   !> Radius at which AE core and model core cross
   real(dp) :: r_ae_model_match
   !> Teter model core charge scaling parameters
   real(dp) :: teter_params(2)
   !> Match search radial cutoff
   real(dp) :: r_search_max
   !> Loop indices
   integer :: i, j, k
   character(len=1024) :: error_msg

   allocate(vo(mmax), vxcae(mmax), vxcpsp(mmax))
   allocate(dvxcae(mmax, nv), dvxcps(mmax, nv))
   allocate(d2excae(nv, nv), d2excps(nc, nv))

   rhomod = 0.0_dp
   d2excae = 0.0_dp
   d2excps = 0.0_dp

   write (stdout, '(/a/a)') 'Model core correction analysis', '  based on d2Exc/dn_idn_j'

   ! Compute d2Exc for all-electron atom
   call der2exc(rhotae, rhoc, rhoae, rr, d2excae, zion, iexc, nc, nv, la, mmax, mmax)
   write (stdout, '(/a/)') 'd2excae - all-electron derivatives'
   do k = 1, nv
      write (stdout, '(1p,4d16.6)') (d2excae(k, j), j=1, nv)
   end do

   ! Compute d2exc and d2exc_rmse for pseudo atom with no core correction
   call der2exc(rhotps, rhomod(:, 1), rhops, rr, d2excps, zion, iexc, nc, nv, la, mmax, mmax)
   write (stdout, '(/a/)') 'd2excps - pseudofunction derivatives with no core correction'
   do k = 1, nv
      write (stdout, '(1p,4d16.6)') (d2excps(k, j), j=1, nv)
   end do
   write (stdout, '(/a,1p,e16.6)') 'rms 2nd derivative error', rmse(d2excps, d2excae, nv)

   ! Find the core -- valence crossover point
   do i = mmax, 1, -1
      if (rhoc(i) > rhotps(i)) then
         i_ae_core_ps_val_match = i
         r_ae_core_ps_val_match = rr(i_ae_core_ps_val_match)
         rhoc_ae_core_ps_val_match = rhoc(i_ae_core_ps_val_match)
         exit
      end if
   end do
   if (i_ae_core_ps_val_match == 0) then
      error stop 'model_core_charge_m: ERROR core-valence charge crossover not found'
   end if
   if (i_ae_core_ps_val_match < 20) then
      error stop 'model_core_charge_m: ERROR core-valence charge crossover too close to origin'
   end if
   write (stdout, '(/a,2f10.4)') 'rmatch, rhocmatch', r_ae_core_ps_val_match, rhoc_ae_core_ps_val_match

   teter_params(1) = fcfact * rhoc_ae_core_ps_val_match
   teter_params(2) = rcfact * r_ae_core_ps_val_match

   ! Find crossing radius between Teter model core charge and all-electron core charge
   ! Search inward starting from 1.5 times scaled core -- valence crossover radius
   ir_ae_model_match = 0
   r_ae_model_match = 0.0_dp
   r_search_max = 1.5_dp * teter_params(2)
   do i = mmax, 1, -1
      ! A.Z.: I think this is how it _should_ be:
      ! if ((rr(i) < r_search_max) .and. teter_params(1) * teter(rr(i) / teter_params(2)) < rhoc(i)) then
      ! A.Z.: But to reproduce previous results, I leave it as it was:
      if ((rr(i) < r_search_max) .and. teter_params(1) * teter(0.0_dp) <= rhoc(i)) then
         ir_ae_model_match = i
         r_ae_model_match = rr(i)
         exit
      end if
   end do
   if (ir_ae_model_match == 0) then
      error stop 'model_core_charge_m: ERROR AE core charge -- model core charge crossover not found'
   end if

   ! Compute the model core charge and its first two derivatives
   rhomod = 0.0_dp
   do i = 1, mmax
      rhomod(i, 1) = teter_params(1) * teter(rr(i) / teter_params(2))
      if (i < ir_ae_model_match) then
         rhomod(i, 2) = teter_params(1) * teter_deriv_1(rr(i) / teter_params(2)) / teter_params(2)
         rhomod(i, 3) = teter_params(1) * teter_deriv_2(rr(i) / teter_params(2)) / teter_params(2)**2
      end if
   end do

   ! Test the fitted model core charge by computing d2Exc
   call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, zion, iexc, nc, nv, la, mmax, mmax)
   write (stdout, '(/a/)') 'd2excps - pseudofunction derivatives with core correction'
   do k = 1, nv
      write (stdout, '(1p,4d16.6)') (d2excps(k, j), j=1, nv)
   end do
   write (stdout, '(/a,1p,e16.6)') 'rms 2nd derivative error', rmse(d2excps, d2excae, nv)

   ! Blend the Teter function tail into the all-electron core charge density
   call model_core_charge_finite_diff_345()

   deallocate(vo, vxcae, vxcpsp)
   deallocate(dvxcae, dvxcps)
   deallocate(d2excae, d2excps)
   return
end subroutine construct_model_core_charge_3

subroutine construct_model_core_charge_4()
end subroutine construct_model_core_charge_4

subroutine construct_model_core_charge_5()
end subroutine construct_model_core_charge_5

subroutine construct_model_core_charge_6()
end subroutine construct_model_core_charge_6

!> Computes the 2nd derivative of the contribution to the exchange-
!> correlation energy from the region where the valence pseudo wave functions
!> differ from the all-electron wave functions
subroutine der2exc(rhotot, rhoc, rhops, rr, d2exc, zion, iexc, nc, nv, la, ircut, mmax)
   use precision_m, only: dp
   implicit none

   ! Input variables
   !> Dimension of log grid
   integer, intent(in) :: mmax
   !> Number of core states
   integer, intent(in) :: nc
   !> Number of valence states
   integer, intent(in) :: nv
   !> Total valence charge, all-electron of pseudo
   real(dp), intent(in) :: rhotot(mmax)
   !> Core charge, all-electron or model
   real(dp), intent(in) :: rhoc(mmax)
   !> Valence state-by-state charge (one-electron)
   real(dp), intent(in) :: rhops(mmax, nv)
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Ion charge
   real(dp), intent(in) :: zion
   !> Exchange-correlation type
   integer, intent(in) :: iexc
   !> Array of angular momenta for all-electron atom
   integer, intent(in) :: la(nv + nc)
   !> Index of cutoff radius for 2nd derivative calculation
   integer, intent(in) :: ircut

   ! Output variables
   real(dp), intent(out) :: d2exc(nv, nv)

   ! Local variables
   real(dp) :: hh, eeel, eexc
   real(dp), allocatable :: vo(:), vxct(:), rhot(:), dvxc(:, :)
   integer :: jj, kk, l1

   allocate (vo(mmax), vxct(mmax), rhot(mmax), dvxc(mmax, nv))

   hh = 0.150_dp

   do kk = 1, nv
      dvxc(:, kk) = 0.0_dp

      rhot(:) = rhotot(:) - 2.0_dp * hh * rhops(:, kk)
      call vout(1, rhot, rhoc, vo, vxct, zion, eeel, eexc, rr, mmax, iexc)
      dvxc(:, kk) = dvxc(:, kk) + (2.0_dp / (24.0_dp * hh)) * vxct(:)

      rhot(:) = rhotot(:) - hh * rhops(:, kk)
      call vout(1, rhot, rhoc, vo, vxct, zion, eeel, eexc, rr, mmax, iexc)
      dvxc(:, kk) = dvxc(:, kk) + (-16.0_dp / (24.0_dp * hh)) * vxct(:)

      rhot(:) = rhotot(:) + hh * rhops(:, kk)
      call vout(1, rhot, rhoc, vo, vxct, zion, eeel, eexc, rr, mmax, iexc)
      dvxc(:, kk) = dvxc(:, kk) + (16.0_dp / (24.0_dp * hh)) * vxct(:)

      rhot(:) = rhotot(:) - 2.0_dp * hh * rhops(:, kk)
      call vout(1, rhot, rhoc, vo, vxct, zion, eeel, eexc, rr, mmax, iexc)
      dvxc(:, kk) = dvxc(:, kk) + (-2.0_dp / (24.0_dp * hh)) * vxct(:)
   end do  !kk

   ! compute Exc 2nd-derivative wrt occupation numbers matrix
   do kk = 1, nv
      l1 = la(nc + kk) + 1
      rhot(:) = rhops(:, kk) * rr(:)**2
      do jj = 1, nv
         call vpinteg(rhot, dvxc(1, jj), ircut, 2 * l1, d2exc(kk, jj), rr)
      end do  !jj
   end do  !kk

   deallocate (vo, vxct, rhot, dvxc)
end subroutine der2exc

!> Compute the root mean square error between two 2D arrays:
!> $\sqrt{\frac{1}{n^2}\sum_{i=1,j=1}^{i=n,j=n} (x_{ij} - y_{ij})^2}$
function rmse(x, y, n) result(error)
   use precision_m, only: dp
   implicit none
   !> Square 2D array size
   integer, intent(in) :: n
   !> First 2D array
   real(dp), intent(in) :: x(n, n)
   !> Second 2D array
   real(dp), intent(in) :: y(n, n)
   !> Root mean square error between the two arrays
   real(dp) :: error
   !> Loop indices
   integer :: i, j
   error = 0.0_dp
   do i = 1, n
      do j = 1, n
         error = error + (x(j,i) - y(j,i))**2
      end do  ! i
   end do  ! j
   error = sqrt(error / dble(n**2))
   return
end function rmse

!> NAME
!> gg1cc
!>
!> FUNCTION
!> yy=$(\frac{\sin(2\pi xx)}{(2\pi xx)(1-4xx^2)(1-xx^2)})^2$
!>
!> COPYRIGHT
!> Copyright (C) 1998-2014 ABINIT group (XG, DCA, MM)
!> This file is distributed under the terms of the
!> GNU General Public License, see ~abinit/COPYING
!> or http://www.gnu.org/copyleft/gpl.txt .
!> For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!>
!> INPUTS
!>  xx= abscisse to which yy is calculated
!>
!> OUTPUT
!>  yy=teter(xx)
function teter(xx) result(yy)
   implicit none

   !Arguments ------------------------------------
   real(dp), intent(in) :: xx

   !Return value -------------------------------------------
   real(dp) :: yy

   !Local variables -------------------------------------------
   ! The cs are coefficients for Taylor expansion of the analytic form near xx=0, 1/2, and 1.
   real(dp), parameter :: c21 = 4.0_dp / 9.0_dp
   real(dp), parameter :: c22 = -40.0_dp / 27.0_dp
   real(dp), parameter :: c23 = 20.0_dp / 3.0_dp - 16.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c24 = -4160.0_dp / 243.0_dp + 160.0_dp * pi**2 / 81.0_dp
   real(dp), parameter :: c31 = 1.0_dp / 36.0_dp
   real(dp), parameter :: c32 = -25.0_dp / 108.0_dp
   real(dp), parameter :: c33 = 485.0_dp / 432.0_dp - pi**2 / 27.0_dp
   real(dp), parameter :: c34 = -4055.0_dp / 972.0_dp + 25.0_dp * pi**2 / 81.0_dp
   real(dp) :: sox
   real(dp) :: zz

   ! Take care of difficult limits near x=0, 1/2, and 1
   if (xx > 3.0_dp) then
      ! Cut off beyond 3/gcut=xcccrc
      yy = 0.0_dp
   else if (abs(xx) <= 1.0e-3_dp) then
      zz = 2.0_dp * pi * xx
      sox = 1.0_dp - zz**2 / 6.0_dp + zz**4 / 120.0_dp + zz**6 / 5040.0_dp
      yy = (sox / ((1.0_dp - 4.0_dp * xx**2) * (1.0_dp - xx**2)))**2
   else if (abs(xx - 0.5_dp) <= 1.0e-4_dp) then
      !  (this limit and next are more troublesome for numerical cancellation)
      yy = c21 + (xx - 0.5_dp) * (c22 + (xx - 0.5_dp) * (c23 + (xx - 0.5_dp) * c24))
   else if (abs(xx - 1.0_dp) <= 1.0e-4_dp) then
      yy = c31 + (xx - 1.0_dp) * (c32 + (xx - 1.0_dp) * (c33 + (xx - 1.0_dp) * c34))
   else
      ! The following is the square of the Fourier transform of a
      ! function built out of two spherical bessel functions in G
      ! space and cut off absolutely beyond gcut
      yy = (sin(2.0_dp * pi * xx) / ((2.0_dp * pi * xx) * (1.0_dp - 4.0_dp * xx**2) * (1.0_dp - xx**2)))**2
   end if
   return
end function teter

!> NAME
!> gp1cc
!>
!> FUNCTION
!> Derivative of teter(xx) wrt xx.
!>
!> COPYRIGHT
!> Copyright (C) 1998-2014 ABINIT group (XG, DCA, MM)
!> This file is distributed under the terms of the
!> GNU General Public License, see ~abinit/COPYING
!> or http://www.gnu.org/copyleft/gpl.txt .
!> For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!>
!> INPUTS
!>  xx=abscisse to which yy is calculated
!>
!> OUTPUT
!>  yy=derivative of teter(xx) wrt xx.
!>
!> NOTES
!> $ phi(x) = \frac{\sin(2\pi x)}{(2\pi x)(1-4x^2)(1-x^2)}$
!> $ teter(x)= phi(x)^2$
!> $ gp(x)= 2 * phi(x) * phi''(x)$
!> $ phi''(x)=\frac{\cos(2\pi x)-(1-15x^2+20x^4) phi(x)}{x(1-4x^2)(1-x^2)}$
function teter_deriv_1(xx) result(yy)
   implicit none

   ! Arguments ------------------------------------
   real(dp), intent(in) :: xx

   ! Return value -------------------------------------------
   real(dp) :: yy

   ! Local variables -------------------------------------------
   real(dp), parameter :: c11 = 20.0_dp - 8.0_dp * pi**2 / 3.0_dp
   real(dp), parameter :: c12 = 268.0_dp - 160.0_dp / 3.0_dp * pi**2 + 128.0_dp / 45.0_dp * pi**4
   real(dp), parameter :: c21 = -40.0_dp / 27.0_dp
   real(dp), parameter :: c22 = 40.0_dp / 3.0_dp - 32.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c23 = -4160.0_dp / 81.0_dp + 160.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c24 = 157712.0_dp / 729.0_dp - 320.0_dp * pi**2 / 9.0_dp + 512.0_dp * pi**4 / 405.0_dp
   real(dp), parameter :: c25 = -452200.0_dp / 729.0_dp + 83200.0_dp * pi**2 / 729.0_dp - 1280.0_dp * pi**4 / 243.0_dp
   real(dp), parameter :: c31 = -25.0_dp / 108.0_dp
   real(dp), parameter :: c32 = 485.0_dp / 216.0_dp - 2.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c33 = -4055.0_dp / 324.0_dp + 25.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c34 = 616697.0_dp / 11664.0_dp - 485.0_dp * pi**2 / 81.0_dp + 32.0_dp * pi**4 / 405.0_dp
   real(dp), parameter :: c35 = -2933875.0_dp / 15552.0_dp + 20275.0_dp * pi**2 / 729.0_dp - 200.0_dp * pi**4 / 243.0_dp
   real(dp), parameter :: invtwopi = 1.0_dp / twopi
   real(dp) :: denom, phi, phip

   if (xx > 3.0_dp) then
      !  Cut off beyond 3/gcut=3*xcccrc
      yy = 0.0_dp
   else if (xx > 1.0010_dp) then
      ! The part that follows will be repeated later, but written in this way,
      ! only one "if" condition is tested in most of the cases (1.001 < x < 3.0)
      denom = 1.0_dp / (xx * (1.0_dp - 4.0_dp * xx**2) * (1.0_dp - xx**2))
      phi = denom * sin(twopi * xx) * invtwopi
      phip = denom * (cos(twopi * xx) - (1.0_dp - xx**2 * (15.0_dp - xx**2 * 20.0_dp)) * phi)
      yy = 2.0_dp * phi * phip
      ! Handle limits where denominator vanishes
   else if (abs(xx) < 1.0e-3_dp) then
      yy = xx * (c11 + xx**2 * c12)
   else if (abs(xx - 0.5_dp) <= 1.0e-3_dp) then
      yy = c21 + (xx - 0.5_dp) * (c22 + (xx - 0.5_dp) * (c23 + (xx - 0.5_dp) * (c24 + (xx - 0.5_dp) * c25)))
   else if (abs(xx - 1.0_dp) <= 1.0e-3_dp) then
      yy = c31 + (xx - 1.0_dp) * (c32 + (xx - 1.0_dp) * (c33 + (xx - 1.0_dp) * (c34 + (xx - 1.0_dp) * c35)))
   else
      ! Here is the repeated part ...
      denom = 1.0_dp / (xx * (1.0_dp - 4.0_dp * xx**2) * (1.0_dp - xx**2))
      phi = denom * sin(twopi * xx) * invtwopi
      phip = denom * (cos(twopi * xx) - (1.0_dp - xx**2 * (15.0_dp - xx**2 * 20.0_dp)) * phi)
      yy = 2.0_dp * phi * phip
   end if
   return
end function teter_deriv_1

!> NAME
!> gpp1cc
!>
!> FUNCTION
!> Second derivative of teter wrt xx.
!>
!> COPYRIGHT
!> Copyright (C) 1998-2014 ABINIT group (XG, DCA, MM, DRH)
!> This file is distributed under the terms of the
!> GNU General Public License, see ~abinit/COPYING
!> or http://www.gnu.org/copyleft/gpl.txt .
!> For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!>
!> INPUTS
!>  xx= abscisse to which yy is calculated
!>
!> OUTPUT
!>  yy=second derivative of teter wrt xx.
function teter_deriv_2(xx) result(yy)
   implicit none

   ! Arguments ------------------------------------
   real(dp), intent(in) :: xx

   ! Return value -------------------------------------------
   real(dp) :: yy

   ! Local variables -------------------------------------------
   real(dp), parameter :: c2 = 40.0_dp / 3.0_dp - 32.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c3 = -8320.0_dp / 81.0_dp + 320.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c4 = 157712.0_dp / 243.0_dp - 320.0_dp * pi**2 / 3.0_dp + 512.0_dp * pi**4 / 135.0_dp
   real(dp), parameter :: c5 = -18088.d2 / 729.0_dp + 3328.d2 * pi**2 / 729.0_dp - 5120.0_dp * pi**4 / 243.0_dp
   real(dp), parameter :: c6 = 485.0_dp / 216.0_dp - 2.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c7 = -4055.0_dp / 162.0_dp + 50.0_dp * pi**2 / 27.0_dp
   real(dp), parameter :: c8 = 616697.0_dp / 3888.0_dp - 485.0_dp * pi**2 / 27.0_dp + 32.0_dp * pi**4 / 135.0_dp
   real(dp), parameter :: c9 = -2933875.0_dp / 3888.0_dp + 81100.0_dp * pi**2 / 729.0_dp - 800.0_dp * pi**4 / 243.0_dp
   !Series expansion coefficients around xx=0 from Mathematica
   real(dp), parameter :: cpp0 = 20.0_dp - 8.0_dp * pi**2 / 3.0_dp
   real(dp), parameter :: cpp2 = 804.0_dp - 160.0_dp * pi**2 + (128.0_dp * pi**4) / 15.0_dp
   real(dp), parameter :: cpp4 = 11400.0_dp - 2680.0_dp * pi**2 + (640.0_dp * pi**4) / 3.0_dp &
      - (128.0_dp * pi**6) / 21.0_dp
   real(dp), parameter :: cpp6 = 8.0_dp * (27967275.0_dp - 7182000.0_dp * pi**2 &
                                           + 675360.0_dp * pi**4 &
                                           - 28800.0_dp * pi**6 + 512.0_dp * pi**8) / 2025.0_dp

   real(dp) :: t1, t10, t100, t11, t12, t120, t121, t122, t127, t138, t14, t140, t15, t152
   real(dp) :: t157, t16, t160, t17, t174, t175, t18, t19, t2, t20, t21, t23, t24, t3, t31, t33
   real(dp) :: t34, t4, t41, t42, t44, t45, t46, t5, t54, t55, t56, t57, t6, t62, t64, t65, t7
   real(dp) :: t72, t78, t79, t8, t85, t9, t93

   if (xx > 3.0_dp) then
      !  Cut off beyond 3/gcut=3*xcccrc
      yy = 0.0_dp
      !  Take care of difficult limits near xx=0, 1/2, and 1
   else if (abs(xx) <= 1.0e-3_dp) then
      yy = cpp0 + cpp2 * xx**2 + cpp4 * xx**4 + cpp6 * xx**6
   else if (abs(xx - 0.5_dp) <= 1.0e-4_dp) then
      !  (this limit and next are more troublesome for numerical cancellation)
      yy = c2 + (xx - 0.5_dp) * (c3 + (xx - 0.5_dp) * (c4 + (xx - 0.5_dp) * c5))
   else if (abs(xx - 1.0_dp) <= 1.0e-4_dp) then
      yy = c6 + (xx - 1.0_dp) * (c7 + (xx - 1.0_dp) * (c8 + (xx - 1.0_dp) * c9))
   else
      !  Should fix up this Maple fortran later
      t1 = xx**2
      t2 = 1 / t1
      t3 = 1 / pi
      t4 = 2 * xx
      t5 = t4 - 1
      t6 = t5**2
      t7 = 1 / t6
      t8 = t4 + 1
      t9 = t8**2
      t10 = 1 / t9
      t11 = xx - 1
      t12 = t11**2
      t14 = 1 / t12 / t11
      t15 = xx + 1
      t16 = t15**2
      t17 = 1 / t16
      t18 = pi * xx
      t19 = sin(t18)
      t20 = cos(t18)
      t21 = t20**2
      t23 = t19 * t21 * t20
      t24 = t17 * t23
      t31 = t19**2
      t33 = t31 * t19 * t20
      t34 = t17 * t33
      t41 = pi**2
      t42 = 1 / t41
      t44 = 1 / t16 / t15
      t45 = t31 * t21
      t46 = t44 * t45
      t54 = 1 / t1 / xx
      t55 = 1 / t12
      t56 = t55 * t46
      t57 = t10 * t56
      t62 = t9**2
      t64 = t17 * t45
      t65 = t55 * t64
      t72 = 1 / t9 / t8
      t78 = t14 * t64
      t79 = t10 * t78
      t85 = t12**2
      t93 = t21**2
      t100 = t31**2
      t120 = 1 / t6 / t5
      t121 = t55 * t34
      t122 = t10 * t121
      t127 = t16**2
      t138 = t6**2
      t140 = t10 * t65
      t152 = t72 * t65
      t157 = t7 * t140
      t160 = t1**2
      t174 = t55 * t24
      t175 = t10 * t174
      yy = 8 * t2 * t3 * t7 * t10 * t14 * t34 + 8 * t2 * t42 * t7 * t10 * t14 * t46&
         &   - 8 * t2 * t3 * t7 * t10 * t14 * t24 + 8 * t2 * t3 * t7 * t10 * t55 * t44 * t33 +&
         &   6 * t2 * t42 * t7 * t10 * t55 / t127 * t45 + 24 * t2 * t42 / t138 * t140 +&
         &   16 * t54 * t42 * t120 * t140 + 16 * t2 * t3 * t120 * t122 + 16 * t2&
         &   * t42 * t7 * t72 * t78 - 8 * t2 * t3 * t7 * t10 * t55 * t44 * t23 - 8 * t54 * t3 * t7 * t175&
         &   + 2 * t2 * t7 * t10 * t55 * t17 * t100 + 2 * t2 * t7 * t10 * t55 * t17 * t93 +&
         &   8 * t54 * t42 * t7 * t79 + 16 * t2 * t42 * t7 * t72 * t56 + 6 * t2 * t42 * t7 * t10 / t85&
         &   * t64 + 24 * t2 * t42 * t7 / t62 * t65 + 8 * t54 * t42 * t7 * t57 -&
         &   16 * t2 * t3 * t7 * t72 * t174 + 8 * t54 * t3 * t7 * t122 - 16 * t2 * t3 * t120 * t175&
         &   + 16 * t2 * t42 * t120 * t79 + 16 * t2 * t42 * t120 * t57 + 16 * t54 * t42 * t7 * t152 +&
         &   32 * t2 * t42 * t120 * t152 + 16 * t2 * t3 * t7 * t72 * t121 - 12 * t2 * t157 +&
         &   6 / t160 * t42 * t157
   end if
   return
end function teter_deriv_2

subroutine optimize_teter_factors_d2exc_rmse(x)
   implicit none
   ! Input/Output variables
   !> Parameters to be optimized: rhom(r) = x(1) * teter(r / x(2))
   real(dp), intent(inout) :: x(2)


end subroutine optimize_teter_factors_d2exc_rmse


end module model_core_charge_m

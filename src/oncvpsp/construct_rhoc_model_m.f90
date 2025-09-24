module construct_rhoc_model_m
   use precision_m, only: dp
   use teter_m, only: teter, teter_deriv_1, teter_deriv_2
   implicit none
   private
   public :: construct_rhoc_model_1
   public :: construct_rhoc_model_2
   public :: construct_rhoc_model_3
contains

subroutine construct_rhoc_model_1(rhops, rhotps, rhoc, rhoae, rhotae, rhomod, fcfact, &
                                  irps, mmax, rr, nc, nv, la, zion, iexc)
   ! Constants
   real(dp), parameter :: EPS = 1.0d-7

   !Input variables
   integer :: nv, nc, iexc, irps, mmax
   integer :: la(30)
   real(dp) :: rhoae(mmax, nv), rhops(mmax, nv), rhotae(mmax)
   real(dp) :: rhotps(mmax), rhoc(mmax), rr(mmax)
   real(dp) :: zion, fcfact

   !Output variables
   real(dp) :: rhomod(mmax, 5)

   !Local variables
   real(dp) :: a0, al, a0min, a0max, dermax, psum
   real(dp) :: d2mdiff
   real(dp) :: xx
   real(dp) :: aco(5), polym(5, 5), work(5, 5), constm(5, 5), xpow(9), fmatch(5)
   real(dp) :: dr(mmax)
   real(dp), allocatable :: vxcae(:), vxcpsp(:), vo(:), d2excae(:, :), d2excps(:, :)
   real(dp), allocatable :: dvxcae(:, :), dvxcps(:, :), vxct(:)
   integer :: ii, ircc, irmod, iter, jj, kk
   integer :: ipvt(5)

   allocate (vxcae(mmax), vxcpsp(mmax), vo(mmax))
   allocate (dvxcae(mmax, nv), dvxcps(mmax, nv), vxct(mmax))
   allocate (d2excae(nv, nv), d2excps(nv, nv))

   d2excae(:, :) = 0.0d0
   d2excps(:, :) = 0.0d0
   ircc = 0

   do ii = mmax, 1, -1
      if (rhoc(ii) > fcfact * rhotps(ii)) then
         ircc = ii
         exit
      end if
   end do

   if (ircc == 0) then
      write (6, '(/a)') 'rhomod: ERROR ircc (core-valence charge crossover) &
      &        not found'
      stop
   end if

   !set limit for d2Exc calculation to radius beyond which rhomod=rhoc
   irmod = max(irps, ircc)

   write (6, '(/a/a)') 'Model core correction analysis',&
      '  based on d2Exc/dn_idn_j'

   ! get derivatives of all-electron xc energy

   call der2exc(rhotae, rhoc, rhoae, rr, d2excae, d2excps, d2mdiff, &
                zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excae - all-electron derivatives'
   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excae(kk, jj), jj=1, nv)
   end do

   ! set model charge to zero
   rhomod(:, :) = 0.0d0

   ! compute d2excps with no core correction
   rhomod(:, :) = 0.0d0

   call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
                zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excps - pseudofunction derivatives with no core correction'
   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (6, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2mdiff

   al = 0.01d0 * dlog(rr(101) / rr(1))
   dr(:) = al * rr(:)

   ! core charge density for Louie-Froyen-Cohen correction

   rhomod(:, 1) = rhoc(:)

   xx = rr(ircc)
   fmatch(1) = rhoc(ircc)

   ! core charge derivatives
   do jj = 2, 5

      !7-point numerical first derivatives applied successively
      call fin_diff_1st_6(rhomod(:, jj - 1), dr, ircc - 25 + 3 * jj, mmax - 3, mmax, rhomod(:, jj))
      ! do ii = ircc - 25 + 3 * jj, mmax - 3
      !    rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) + 9.d0 * rhomod(ii - 2, jj - 1)&
      !       &     - 45.d0 * rhomod(ii - 1, jj - 1) + 45.d0 * rhomod(ii + 1, jj - 1)&
      !       &     - 9.d0 * rhomod(ii + 2, jj - 1) + rhomod(ii + 3, jj - 1))&
      !       &     / (60.d0 * al * rr(ii))
      ! end do
      fmatch(jj) = rhomod(ircc, jj)
   end do

   ! constants for polynomial matrix
   do jj = 1, 5
      constm(1, jj) = 1.0d0
   end do
   do jj = 1, 5
      constm(2, jj) = jj + 2
   end do
   do ii = 3, 5
      do jj = 1, 5
         constm(ii, jj) = constm(ii - 1, jj) * (constm(2, jj) - ii + 2)
      end do
   end do

   ! powers of xx=rr(ircc)

   xpow(1) = 0.0d0
   xpow(2) = 1.0d0
   do ii = 3, 9
      xpow(ii) = xx * xpow(ii - 1)
   end do

   ! polynomial matrix
   do jj = 1, 5
      do ii = 1, 5
         polym(ii, jj) = constm(ii, jj) * xpow(jj - ii + 5)
      end do
   end do

   ! begin iteration for monotonic polynomial match

   a0 = fmatch(1)
   a0min = a0
   a0max = 0.0d0

   do iter = 1, 50

      work(:, :) = polym(:, :)

      ! fill target value vector
      aco(1) = fmatch(1) - a0
      do jj = 2, 5
         aco(jj) = fmatch(jj)
      end do

      ! solve linear equations for coefficients
      call dgesv(5, 1, work, 5, ipvt, aco, 5, kk)
      if (kk /= 0) then
         if (kk > 0) write (6, '(a,i4)') &
         &      'modcore:ERROR stop - singular polym matrix', kk
         if (kk < 0) write (6, '(a,i4)') &
         &      'modcore:ERROR stop - dgesv input error', kk
         stop
      end if

      ! find maximum derivative
      ! note that xpow end up properly re-written for continued iteation
      dermax = -1.0d20
      do kk = 1, ircc
         do ii = 3, 9
            xpow(ii) = rr(kk) * xpow(ii - 1)
         end do
         ii = 2
         psum = 0.0d0
         do jj = 1, 5
            psum = psum + constm(ii, jj) * xpow(jj - ii + 5) * aco(jj)
         end do
         dermax = dmax1(dermax, psum)
      end do

      ! test maximum derivative and adjust a0
      ! interval-halving search after achieving monotonicity to get
      ! barely monotonic model
      if (dermax > 0.0d0) then
         a0min = a0
      else
         a0max = a0
      end if

      ! success when maximum derivative bracketed just above zero
      if (dermax > 0.0d0 .and. dermax < EPS * dabs(fmatch(2))) exit

      if (a0max == 0.0d0) then
         a0 = 2.5d0 * a0
      else
         a0 = 0.5d0 * (a0max + a0min)
      end if

   end do  !iterr
   if (iter > 50) write (6, '(/a/)') 'WARNING - modcore not conveged'

   ! fill in model and derivatives with fitted polynomial
   do kk = 1, ircc - 1
      do ii = 3, 9
         xpow(ii) = rr(kk) * xpow(ii - 1)
      end do
      do ii = 1, 5
         if (ii == 1) then
            psum = a0
         else
            psum = 0.0d0
         end if
         do jj = 1, 5
            psum = psum + constm(ii, jj) * xpow(jj - ii + 5) * aco(jj)
         end do
         rhomod(kk, ii) = psum
      end do
   end do

   !test model
   call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
                zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'Polynomial model core charge'
   write (6, '(/a/)') 'd2excps - pseudofunction derivatives with core correction'

   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (6, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2mdiff

   deallocate (vxcae, vxcpsp, vo)
   deallocate (dvxcae, dvxcps)
   deallocate (d2excae, d2excps)

   return
end subroutine construct_rhoc_model_1

subroutine construct_rhoc_model_2(rhops, rhotps, rhoc, rhoae, rhotae, rhomod, fcfact, &
                                  mmax, rr, nc, nv, la, zion, iexc)
   ! Constants
   real(dp), parameter :: eps = 1.0d-7

   ! Input variables
   integer :: nv, nc, iexc, mmax
   integer :: la(30)
   real(dp) :: rhoae(mmax, nv), rhops(mmax, nv), rhotae(mmax)
   real(dp) :: rhotps(mmax), rhoc(mmax), rr(mmax)
   real(dp) :: zion, fcfact

   ! Output variables
   real(dp) :: rhomod(mmax, 5)

   ! Local variables
   real(dp) :: al
   real(dp) :: d2mdiff, rmatch, rhocmatch
   real(dp) :: xx
   real(dp) :: x0max, x0min, a0, b0, ymatch, ytrial
   real(dp) :: fmatch(5)
   real(dp) :: drint, rtst, rint(20), fint(20)  !ad-hoc smoothing variables
   real(dp) :: dr(mmax)
   real(dp), allocatable :: vxcae(:), vxcpsp(:), vo(:), d2excae(:, :), d2excps(:, :)
   real(dp), allocatable :: dvxcae(:, :), dvxcps(:, :), vxct(:)
   integer :: ii, ircc, irmod, jj, kk
   integer :: iint  !ad-hoc smoothing variables

   allocate (vxcae(mmax), vxcpsp(mmax), vo(mmax))
   allocate (dvxcae(mmax, nv), dvxcps(mmax, nv), vxct(mmax))
   allocate (d2excae(nv, nv), d2excps(nv, nv))

   d2excae(:, :) = 0.0d0
   d2excps(:, :) = 0.0d0

   ! find valence pseudocharge - core charge crossover
   ! this is needed for compatability with icmod=3 option
   ircc = 0
   do ii = mmax, 1, -1
      if (rhoc(ii) > rhotps(ii)) then
         ircc = ii
         rmatch = rr(ircc)
         rhocmatch = rhoc(ircc)
         exit
      end if
   end do

   ! find scaled valence pseudocharge - core charge crossover
   ircc = 0
   do ii = mmax, 1, -1
      if (rhoc(ii) > fcfact * rhotps(ii)) then
         ircc = ii
         exit
      end if
   end do

   if (ircc == 0) then
      write (6, '(/a)') 'modcore2: ERROR ircc (core-valence charge crossover) not found'
      stop
   end if

   !set limit for d2Exc calculation to radius beyond which rhomod=rhoc
   irmod = mmax

   write (6, '(/a/a)') 'Model core correction analysis',&
      '  based on d2Exc/dn_idn_j'

   ! get derivatives of all-electron xc energy

   call der2exc(rhotae, rhoc, rhoae, rr, d2excae, d2excps, d2mdiff, &
                zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excae - all-electron derivatives'
   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excae(kk, jj), jj=1, nv)
   end do

   ! set model charge to zero
   rhomod(:, :) = 0.0d0

   ! compute d2excps with no core correction
   rhomod(:, :) = 0.0d0

   call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
                zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excps - pseudofunction derivatives with no core correction'
   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (6, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2mdiff

   al = 0.01d0 * dlog(rr(101) / rr(1))
   dr(:) = al * rr(:)

   ! core charge density for Louie-Froyen-Cohen correction

   rhomod(:, 1) = rhoc(:)

   xx = rr(ircc)
   fmatch(1) = rhoc(ircc)

   ! core charge derivatives

   !7-point numerical first derivative
   ii = ircc
   jj = 2
   rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) + 9.d0 * rhomod(ii - 2, jj - 1)&
   &     - 45.d0 * rhomod(ii - 1, jj - 1) + 45.d0 * rhomod(ii + 1, jj - 1)&
   &     - 9.d0 * rhomod(ii + 2, jj - 1) + rhomod(ii + 3, jj - 1))&
   &     / (60.d0 * al * rr(ii))
   fmatch(jj) = rhomod(ircc, jj)

   ! Fit Teter function to value and slope at ircc

   ymatch = rr(ircc) * fmatch(2) / fmatch(1)

   ! interval-halving search for dimensionless match point

   x0max = 1.48d0
   x0min = 0.0d0

   do jj = 1, 50

      xx = 0.5d0 * (x0max + x0min)

      ! call gg1cc(yy, xx)
      ! call gp1cc(dy, xx)
      ! ytrial = xx * dy / yy
      ytrial = xx * teter_deriv_1(xx) / teter(xx)

      if (abs(ytrial - ymatch) < eps) exit

      if (ytrial < ymatch) then
         x0max = xx
      else
         x0min = xx
      end if
   end do

   b0 = xx / rr(ircc)
   a0 = fmatch(1) / teter(xx)

   write (6, '(/a)') 'amplitude prefactor, scale prefactor'
   write (6, '(2f10.4)') a0 / rhocmatch, 1.0d0 / (b0 * rmatch)

   rhomod(:, :) = 0.0d0
   do ii = 1, mmax
      xx = b0 * rr(ii)
      if (xx < 3.0d0) then
         ! call gg1cc(yy, xx)
         rhomod(ii, 1) = a0 * teter(xx)
         ! call gp1cc(yy, xx)
         rhomod(ii, 2) = a0 * teter_deriv_1(xx) * b0
         ! call gpp1cc(yy, xx)
         rhomod(ii, 3) = a0 * teter_deriv_2(xx) * b0**2
      end if
   end do

   ! test model
   call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
                zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excps - pseudofunction derivatives with core correction'

   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (6, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2mdiff

   !7-point numerical first derivatives applied successively

   !ad-hoc treatment of numerical noise near origin
   !set up a mesh on which 2nd derivative will have a stable
   !polynomial representation
   !assumes dpnint remains 7th order
   drint = 0.05d0 * rr(ircc)
   rtst = 0.5d0 * drint
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

   call dpnint(rint, fint, 8, rr, rhomod(1, 3), iint - 1)

   jj = 4
   call fin_diff_1st_6(rhomod(:, jj - 1), dr, 3 * jj - 8, mmax - 3, mmax, rhomod(:, jj))
   ! do ii = 3 * jj - 8, mmax - 3
   !    rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) + 9.d0 * rhomod(ii - 2, jj - 1)&
   !       &     - 45.d0 * rhomod(ii - 1, jj - 1) + 45.d0 * rhomod(ii + 1, jj - 1)&
   !       &     - 9.d0 * rhomod(ii + 2, jj - 1) + rhomod(ii + 3, jj - 1))&
   !       &     / (60.d0 * al * rr(ii))
   ! end do

   !set up a mesh on which 3rd derivative will have a stable
   drint = 0.95d0 * drint
   rtst = 0.5d0 * drint
   rtst = 0.5d0 * drint
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

   call dpnint(rint, fint, 8, rr, rhomod(1, 4), iint - 1)

   jj = 5
   call fin_diff_1st_6(rhomod(:, jj - 1), dr, 3 * jj - 8, mmax - 3, mmax, rhomod(:, jj))
   ! do ii = 3 * jj - 8, mmax - 3
   !    rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) + 9.d0 * rhomod(ii - 2, jj - 1)&
   !       &     - 45.d0 * rhomod(ii - 1, jj - 1) + 45.d0 * rhomod(ii + 1, jj - 1)&
   !       &     - 9.d0 * rhomod(ii + 2, jj - 1) + rhomod(ii + 3, jj - 1))&
   !       &     / (60.d0 * al * rr(ii))
   ! end do

   !set up a mesh on which 4th derivative will have a stable
   drint = 0.95d0 * drint
   rtst = 0.5d0 * drint
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

   call dpnint(rint, fint, 8, rr, rhomod(1, 5), iint - 1)

   deallocate (vxcae, vxcpsp, vo)
   deallocate (dvxcae, dvxcps)
   deallocate (d2excae, d2excps)

   return
end subroutine construct_rhoc_model_2

subroutine construct_rhoc_model_3(rhops, rhotps, rhoc, rhoae, rhotae, rhomod, fcfact, rcfact, &
                                  mmax, rr, nc, nv, la, zion, iexc)

   ! Constants
   real(dp), parameter :: EPS = 1.0d-7
   real(dp), parameter :: BLEND = 2.0d0

   ! Input variables
   integer :: nv
   integer :: nc
   integer :: iexc
   integer :: mmax
   integer :: la(30)
   real(dp) :: rhoae(mmax, nv)
   real(dp) :: rhops(mmax, nv)
   real(dp) :: rhotae(mmax)
   real(dp) :: rhotps(mmax)
   real(dp) :: rhoc(mmax)
   real(dp) :: rr(mmax)
   real(dp) :: zion, fcfact, rcfact

   ! Output variables
   real(dp) :: rhomod(mmax, 5)

   ! Local variables
   real(dp) :: al
   real(dp) :: dr(mmax)
   real(dp) :: d2diff, rmatch, rhocmatch, r0
   real(dp) :: gg, yy
   real(dp) :: drint, rtst, rint(20), fint(20)  !ad-hoc smoothing variables
   real(dp), allocatable :: vxcae(:), vxcpsp(:), vo(:), d2excae(:, :), d2excps(:, :)
   real(dp), allocatable :: dvxcae(:, :), dvxcps(:, :), vxct(:)
   integer :: ii, ircc, ircross, irmod, jj, kk
   integer :: iint  !ad-hoc smoothing variables

   !2-dimensional Nelder-Mead variables
   real(dp) :: xx(2, 3)

   allocate (vxcae(mmax), vxcpsp(mmax), vo(mmax))
   allocate (dvxcae(mmax, nv), dvxcps(mmax, nv), vxct(mmax))
   allocate (d2excae(nv, nv), d2excps(nv, nv))

   d2excae(:, :) = 0.0d0
   d2excps(:, :) = 0.0d0

   !set limit for d2Exc calculation to radius beyond which rhomod=rhoc
   irmod = mmax

   write (6, '(/a/a)') 'Model core correction analysis', '  based on d2Exc/dn_idn_j'

   call der2exc(rhotae, rhoc, rhoae, rr, d2excae, d2excps, d2diff, &
                zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excae - all-electron derivatives'
   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excae(kk, jj), jj=1, nv)
   end do

   ! set model charge to zero
   rhomod(:, :) = 0.0d0

   ! compute d2excps with no core correction
   rhomod(:, :) = 0.0d0

   call der2exc(rhotps, rhomod(:, 1), rhops, rr, d2excps, d2excae, d2diff, &
                zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excps - pseudofunction derivatives with no core correction'
   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (6, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2diff

   ! find valence pseudocharge - core charge crossover
   ircc = 0
   do ii = mmax, 1, -1
      if (rhoc(ii) > rhotps(ii)) then
         ircc = ii
         rmatch = rr(ircc)
         rhocmatch = rhoc(ircc)
         exit
      end if
   end do

   if (ircc == 0) then
      write (6, '(/a)') 'construct_rhoc_model_3: ERROR ircc (core-valence charge crossover) &
      &        not found'
      stop
   else
      write (6, '(/a,2f10.4)') 'rmatch, rhocmatch', rmatch, rhocmatch
   end if

   gg = 0.d0
   yy = 0.d0

   xx(1, 1) = fcfact * rhocmatch
   xx(2, 1) = rcfact * rmatch
   write (6, '(/a/)') 'Teter function model core charge with specified parameters'
   write (6, '(a)') 'amplitude prefactor, scale prefactor'
   write (6, '(2f10.4)') xx(1, 1) / rhocmatch, xx(2, 1) / rmatch

   r0 = 1.5d0 * xx(2, 1)
   do ii = mmax, 1, -1
      if (rr(ii) < r0) then
         ! call gg1cc(gg, yy)
         gg = teter(yy)
         if (xx(1, 1) * gg < rhoc(ii)) then
            ircross = ii
            exit
         end if
      end if
   end do

   ! blend the Teter function tail into the all-electron rhoc
   ! first two derivatives are filled in analytically before the blend starts
   rhomod(:, :) = 0.0d0
   do ii = 1, mmax
      yy = rr(ii) / xx(2, 1)
      ! call gg1cc(gg, yy)
      gg = teter(yy)

      rhomod(ii, 1) = xx(1, 1) * gg
      if (ii < ircross) then
         ! call gp1cc(gg, yy)
         rhomod(ii, 2) = xx(1, 1) * teter_deriv_1(yy) / xx(2, 1)
         ! call gpp1cc(gg, yy)
         rhomod(ii, 3) = xx(1, 1) * teter_deriv_2(yy) / xx(2, 1)**2
      end if
   end do

   call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2diff, &
                zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excps - pseudofunction derivatives with core correction'
   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (6, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2diff

   !7-point numerical first derivatives applied successively
   !skip non-blended section for 1st and 2nd derivatives

   al = 0.01d0 * dlog(rr(101) / rr(1))
   dr(:) = al * rr(:)

   do jj = 2, 3
      call fin_diff_1st_6(rhomod(:, jj - 1), dr, ircross - 6, mmax - 3, mmax, rhomod(:, jj))
   end do

   !ad-hoc treatment of numerical noise near origin
   !set up a mesh on which 2nd derivative will have a stable
   !polynomial representation
   !assumes dpnint remains 7th order
   drint = 0.02d0 * rr(ircross)
   rtst = 0.5d0 * drint
   do jj = 1, 4
      do ii = 1, ircross
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

   call dpnint(rint, fint, 8, rr, rhomod(1, 3), iint - 1)

   jj = 4
   call fin_diff_1st_6(rhomod(:, jj - 1), dr, 3 * jj - 8, mmax - 3, mmax, rhomod(:, jj))

   !set up a mesh on which 3rd derivative will have a stable
   drint = 0.95d0 * drint
   rtst = 0.5d0 * drint
   do jj = 1, 4
      do ii = 1, ircross
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

   call dpnint(rint, fint, 8, rr, rhomod(1, 4), iint - 1)

   jj = 5
   call fin_diff_1st_6(rhomod(:, jj - 1), dr, 3 * jj - 8, mmax - 3, mmax, rhomod(:, jj))

   !set up a mesh on which 4th derivative will have a stable
   drint = 0.95d0 * drint
   rtst = 0.5d0 * drint
   do jj = 1, 4
      do ii = 1, ircross
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

   call dpnint(rint, fint, 8, rr, rhomod(1, 5), iint - 1)

   deallocate (vxcae, vxcpsp, vo)
   deallocate (dvxcae, dvxcps)
   deallocate (d2excae, d2excps)
   return

end subroutine construct_rhoc_model_3

!> First-order central finite differences using a 7-point stencil (accuracy=6)
!> Note: does not special-case boundaries
subroutine fin_diff_1st_6(f, dx, first, last, n, y)
   integer, intent(in) :: first, last, n
   real(dp), intent(in) :: f(n)
   real(dp), intent(in) :: dx(n)
   real(dp), intent(out) :: y(n)
   integer :: i

   if (last - first < 6) error stop 'Error in fin_diff_1st_6: array too small'

   do i = first, last
      y(i) = (-f(i - 3) + 9.d0 * f(i - 2) - 45.d0 * f(i - 1) &
              + 45.d0 * f(i + 1) - 9.d0 * f(i + 2) + f(i + 3)) &
         / (60.d0 * dx(i))
   end do
   return
end subroutine fin_diff_1st_6

end module construct_rhoc_model_m

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

   call der2exc(rhotps, rhomod(:, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
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
   call der2exc(rhotps, rhomod(:, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
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

   call der2exc(rhotps, rhomod(:, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
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
   call der2exc(rhotps, rhomod(:, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
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

   call der2exc(rhotps, rhomod(:, 1), rhops, rr, d2excps, d2excae, d2diff, &
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

subroutine construct_rhoc_model_4(rhops, rhotps, rhoc, rhoae, rhotae, rhomod, &
                                  mmax, rr, nc, nv, la, zion, iexc)
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
   real(dp), allocatable :: vxcae(:), vxcpsp(:), vo(:), d2exc_true(:, :)
   real(dp), allocatable :: dvxcae(:, :), dvxcps(:, :), vxct(:)

   allocate (vxcae(mmax), vxcpsp(mmax), vo(mmax))
   allocate (dvxcae(mmax, nv), dvxcps(mmax, nv), vxct(mmax))
   allocate (d2exc_true(nv, nv))

   d2excae(:, :) = 0.0d0
   d2excps(:, :) = 0.0d0

   ! Set limit for d2Exc calculation to radius beyond which rhomod=rhoc
   irmod = mmax

   call der2exc(rhotae, rhoc, rhoae, rr, d2excae, d2excps, d2diff, &
                zion, iexc, nc, nv, la, irmod, mmax)

   ! Compute d2excps with no core correction
   rhomod(:, :) = 0.0d0
   call der2exc(rhotps, rhomod(:, 1), rhops, rr, d2excps, d2excae, d2diff, &
                zion, iexc, nc, nv, la, irmod, mmax)

   ! Find valence pseudocharge - core charge crossover
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
      write (6, '(/a)') 'construct_rhoc_model_4: ERROR ircc (core-valence charge crossover) &
      &        not found'
      stop
   end if

   ! Coarse grid for initial guess
   fcfact_guess = 0.0_dp
   rcfact_guess = 0.0_dp

   ! Nelder-Mead optimization of fcfact and rcfact w.r.t. RMS error in d2Exc
   call teter_nelder_mead(rhops, rhotps, rhoc, rhomod(:,1), d2exc_true, mmax, rr, nc, nv, la, zion, iexc, ircut, &
                          rhocmatch, rmatch, fcfact_guess, rcfact_guess, atol, &
                          fcfact_optimal, rcfact_optimal, obj_optimal)

   call construct_rhoc_model_3(rhops, rhotps, rhoc, rhoae, rhotae, rhomod, fcfact_optimal, rcfact_optimal, &
                               mmax, rr, nc, nv, la, zion, iexc)

   return
end subroutine construct_rhoc_model_4

!> Optimize the Teter function metaparameters `fcfact` and `rcfact` using
!> the Nelder-Mead simplex algorithm and the objective function `objective`
subroutine teter_nelder_mead(rhops, rhotps, rhoc, rhomod, &
                             d2exc_true, mmax, rr, nc, nv, la, &
                             zion, iexc, ircut, rhocmatch, rmatch, &
                             fcfact_guess, rcfact_guess, atol, &
                             fcfact_optimal, rcfact_optimal, obj_optimal)
   ! Constants
   !> Nelder-Mead convergence criterion
   real(dp), parameter :: EPS_NM = 1.0e-4_dp
   !> Nelder-Mead reflection parameter
   real(dp), parameter :: ALPHA_NM = 1.0_dp
   !> Nelder-Mead expansion parameter
   real(dp), parameter :: GAMMA_NM = 2.0_dp
   !> Nelder-Mead contraction parameter
   real(dp), parameter :: RHO_NM = -0.5_dp
   !> Nelder-Mead reduction parameter
   real(dp), parameter :: SIGMA_NM = 0.5_dp
   !> Maximum number of iterations
   integer, parameter :: MAXIT_NM = 100

   ! Input variables
   integer :: mmax
   integer :: nc
   integer :: nv
   real(dp) :: rhops(mmax, nv)
   real(dp) :: rhotps(mmax)
   real(dp) :: rhoc(mmax)
   real(dp) :: rr(mmax)
   integer :: la(nc + nv)
   real(dp) :: zion
   integer :: iexc
   integer :: ircut
   real(dp) :: rhocmatch
   real(dp) :: rmatch
   real(dp) :: fcfact_guess
   real(dp) :: rcfact_guess
   real(dp) :: d2exc_true(nv, nv)
   !> Absolute convergence tolerance
   real(dp) :: atol

   ! Output variables
   real(dp) :: rhomod(mmax)
   real(dp) :: fcfact_optimal
   real(dp) :: rcfact_optimal
   real(dp) :: obj_optimal

   ! Local variables
   integer :: ii, jj, kk
   logical :: converged
   !> Simplex
   real(dp) :: xx(2, 3), ff(3)
   !> Centroid
   real(dp) :: x0(2)
   !> Temporary (sorting)
   real(dp) :: xt(2), ft
   !> Reflection
   real(dp) :: xr(2), fr
   !> Expansion
   real(dp) :: xe(2), fe
   !> Contraction
   real(dp) :: xc(2), fc

   converged = .false.

   interface
      subroutine objective(xt_o, rhot_o, rhov_o, rhoc_o, rhomod_o, rr_o, &
                           d2exc_true_o, d2exc_model_o, zion_o, la_o, iexc_o, &
                           ircut_o, nc_o, nv_o, mmax_o, y_o)
         use precision_m, only: dp
         implicit none
         ! Input variables
         real(dp), intent(in) :: xt_o(2)
         integer, intent(in) :: mmax_o
         integer, intent(in) :: nc_o
         integer, intent(in) :: nv_o
         real(dp), intent(in) :: rhot_o(mmax)
         real(dp), intent(in) :: rhov_o(mmax, nv)
         real(dp), intent(in) :: rhoc_o(mmax)
         real(dp), intent(in out) :: rhomod_o(mmax)
         real(dp), intent(in) :: d2exc_true_o(nv, nv)
         real(dp), intent(in out) :: d2exc_model_o(nv, nv)
         integer, intent(in) :: la_o(nc + nv)
         real(dp), intent(in) :: rr_o(mmax)
         real(dp), intent(in) :: zion_o
         integer, intent(in) :: iexc_o
         integer, intent(in) :: ircut_o
         ! Output variable
         real(dp), intent(out) :: y_o
      end subroutine objective
   end interface

   ! Set up initial simplex
   xx(1, 1) = (fcfact_guess - 0.125d0) * rhocmatch
   xx(2, 1) = (rcfact_guess - 0.025d0) * rmatch
   xx(1, 2) = (fcfact_guess + 0.25d0) * rhocmatch
   xx(2, 2) = rcfact_guess * rmatch
   xx(1, 3) = fcfact_guess * rhocmatch
   xx(2, 3) = (rcfact_guess + 0.05d0) * rmatch
   do kk = 1, 3
      call objective(xx(:, kk), rhotps, rhops, rhoc, rhomod, rr, &
                     d2exc_true, d2exc_model, zion, la, iexc, &
                     ircut, nc, nv, mmax, ff(kk))
   end do

   ! Nelder-Mead iteration loop
   write (6, '(/a)') 'Nelder-Mead iteration'
   do kk = 1, MAXIT_NM
      ! (1) Order (dumb bubble sort)
      do jj = 1, 3
         if (ff(3) < ff(2)) then
            ft = ff(2); xt(:) = xx(:, 2)
            ff(2) = ff(3); xx(:, 2) = xx(:, 3)
            ff(3) = ft; xx(:, 3) = xt(:)
         end if
         if (ff(2) < ff(1)) then
            ft = ff(1); xt(:) = xx(:, 1)
            ff(1) = ff(2); xx(:, 1) = xx(:, 2)
            ff(2) = ft; xx(:, 2) = xt(:)
         end if
      end do

      ! Stopping criterion
      if (abs(ff(3) - ff(1)) < atol) then
         converged = .true.
         write (6, '(a,i4,a)') ' converged in', kk - 1, ' steps'
         write (6, '(a)') 'amplitude prefactor, scale prefactor'
         write (6, '(2f10.4)') xx(1, 1) / rhocmatch, xx(2, 1) / rmatch
         exit
      end if

      ! (2) Centroid
      x0(:) = 0.5d0 * (xx(:, 1) + xx(:, 2))

      ! (3) Reflection
      xr(:) = x0(:) + ALPHA_NM * (x0(:) - xx(:, 3))
      call objective(xr, rhotps, rhops, rhoc, rhomod, rr, &
                     d2exc_true, d2exc_model, zion, la, iexc, &
                     ircut, nc, nv, mmax, fr)
      if (ff(1) <= fr .and. fr < ff(2)) then
         ff(3) = fr; xx(:, 3) = xr(:)
         cycle  ! kk loop
      end if

      ! (4) Expansion
      if (fr < ff(1)) then
         xe(:) = x0(:) + GAMMA_NM * (x0(:) - xx(:, 3))
         call objective(xe, rhotps, rhops, rhoc, rhomod, rr, &
                        d2exc_true, d2exc_model, zion, la, iexc, &
                        ircut, nc, nv, mmax, fe)
         if (fe < fr) then
            ff(3) = fe; xx(:, 3) = xe(:)
            cycle  ! kk
         else
            ff(3) = fr; xx(:, 3) = xr(:)
            cycle  ! kk
         end if
      end if

      ! (5) Contraction
      xc(:) = x0(:) + RHO_NM * (x0(:) - xx(:, 3))
      call objective(xc, rhotps, rhops, rhoc, rhomod, rr, &
                     d2exc_true, d2exc_model, zion, la, iexc, &
                     ircut, nc, nv, mmax, fc)
      if (fc < ff(3)) then
         ff(3) = fc; xx(:, 3) = xc(:)
         cycle  ! kk
      end if

      ! (6) Reduction
      do jj = 2, 3
         xx(:, jj) = xx(:, 1) + SIGMA_NM * (xx(:, jj) - xx(:, 1))
         call objective(xx(:,jj), rhotps, rhops, rhoc, rhomod, rr, &
                        d2exc_true, d2exc_model, zion, la, iexc, &
                        ircut, nc, nv, mmax, ff(jj))
      end do  ! jj
   end do  ! kk

   if (.not. converged) then
      write (6, '(a,i4,a)') ' WARNING: not fully converged in', MAXIT_NM ,' steps'
      write (6, '(a)') 'amplitude prefactor, scale prefactor'
      write (6, '(2f10.4)') xx(1, 1) / rhocmatch, xx(2, 1) / rmatch
      exit
   end if

end subroutine teter_nelder_mead

!> Objective function for Nelder-Mead optimization of Teter parameters
!> which computes the RMSE between the model and reference d2Exc matrices
subroutine rmse_d2exc_objective(xt, rhot, rhov, rhoc, rhomod, rr, &
                                d2exc_true, d2exc_model, zion, la, iexc, &
                                ircut, nc, nv, mmax, y)
   ! Input variables
   real(dp), intent(in) :: xt(2)
   integer, intent(in) :: mmax
   integer, intent(in) :: nc
   integer, intent(in) :: nv
   real(dp), intent(in) :: rhot(mmax)
   real(dp), intent(in) :: rhov(mmax, nv)
   real(dp), intent(in) :: rhoc(mmax)
   real(dp), intent(in out) :: rhomod(mmax)
   real(dp), intent(in) :: d2exc_true(nv, nv)
   real(dp), intent(in out) :: d2exc_model(nv, nv)
   integer, intent(in) :: la(nc + nv)
   real(dp), intent(in) :: rr(mmax)
   real(dp), intent(in) :: zion
   integer, intent(in) :: iexc
   integer, intent(in) :: ircut
   ! Output variable
   real(dp), intent(out) :: y
   ! Local variables
   integer :: i

   ! Construct model core charge
   do i = 1, mmax
      rhomod(i) = xt(1) * teter(rr(i) / xt(2))
   end do
   ! Compute d2Exc RMSE
   call der2exc(rhot, rhomod, rhov, rr, d2exc_model, d2exc_true, y, zion, iexc, nc, nv, la, ircut, mmax)
   deallocate(rhomod)

   return
end subroutine rmse_d2exc_objective

!> Computes the 2nd derivative of the contribution to the exchange-
!> correlation energy from the region where the valence pseudo wave functions
!> differ from the all-electron wave functions
subroutine der2exc(rhotot, rhoc, rho, rr, d2exc, d2ref, d2mdiff, &
                   zion, iexc, nc, nv, la, ircut, mmax)
   use precision_m, only: dp
   ! Input variables
   !> Number of core states
   integer, intent(in) :: nc
   !> Number of valence states
   integer, intent(in) :: nv
   !> Size of logarithmic radial grid
   integer, intent(in) :: mmax
   !> Total valence charge
   real(dp), intent(in) :: rhotot(mmax)
   !> Total core charge
   real(dp), intent(in) :: rhoc(mmax)
   !> Valence state-by-state charge
   real(dp), intent(in) :: rho(mmax, nv)
   !> Logarithmic radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Reference 2nd-derivative matrix
   real(dp), intent(in) :: d2ref(nv, nv)
   !> Ionic charge
   real(dp), intent(in) :: zion
   !> Array of angular momentum quantum numbers for all-electron atom
   integer, intent(in) :: la(nv + nc)
   !> Exchange-correlation type
   integer, intent(in) :: iexc
   !> Maximum radius for which AE and PS charges differ
   integer, intent(in) :: ircut

   ! Output variables
   !> 2nd-derivative matrix
   real(dp), intent(out) :: d2exc(nv, nv)
   !> RMSE b/w d2exc and d2ref
   real(dp), intent(out) :: d2mdiff

   ! Local variables
   real(dp) :: hh
   real(dp) :: eeel
   real(dp) :: eexc
   real(dp) :: ss
   real(dp), allocatable :: vo(:)
   real(dp), allocatable :: vxct(:)
   real(dp), allocatable :: rhot(:)
   real(dp), allocatable :: dvxc(:, :)
   integer :: jj
   integer :: kk
   integer :: l1

   allocate (vo(mmax), vxct(mmax), rhot(mmax), dvxc(mmax, nv))

   hh = 0.15d0

   do kk = 1, nv
      dvxc(:, kk) = 0.0d0

      rhot(:) = rhotot(:) - 2.d0 * hh * rho(:, kk)
      call vout(1, rhot, rhoc, vo, vxct, zion, eeel, eexc, rr, mmax, iexc)
      dvxc(:, kk) = dvxc(:, kk) + (2.0d0 / (24.0d0 * hh)) * vxct(:)

      rhot(:) = rhotot(:) - hh * rho(:, kk)
      call vout(1, rhot, rhoc, vo, vxct, zion, eeel, eexc, rr, mmax, iexc)
      dvxc(:, kk) = dvxc(:, kk) + (-16.0d0 / (24.0d0 * hh)) * vxct(:)

      rhot(:) = rhotot(:) + hh * rho(:, kk)
      call vout(1, rhot, rhoc, vo, vxct, zion, eeel, eexc, rr, mmax, iexc)
      dvxc(:, kk) = dvxc(:, kk) + (16.0d0 / (24.0d0 * hh)) * vxct(:)

      rhot(:) = rhotot(:) - 2.d0 * hh * rho(:, kk)
      call vout(1, rhot, rhoc, vo, vxct, zion, eeel, eexc, rr, mmax, iexc)
      dvxc(:, kk) = dvxc(:, kk) + (-2.0d0 / (24.0d0 * hh)) * vxct(:)

   end do  !kk

   ! compute Exc 2nd-derivative wrt occupation numbers matrix
   do kk = 1, nv
      l1 = la(nc + kk) + 1
      rhot(:) = rho(:, kk) * rr(:)**2
      do jj = 1, nv
         call vpinteg(rhot, dvxc(:, jj), ircut, 2 * l1, d2exc(kk, jj), rr)
      end do  !jj
   end do  !kk

   ss = 0.0d0
   do kk = 1, nv
      do jj = 1, nv
         ss = ss + (d2exc(jj, kk) - d2ref(jj, kk))**2
      end do
   end do
   d2mdiff = sqrt(ss / dble(nv**2))

   deallocate (vo, vxct, rhot, dvxc)
end subroutine der2exc

!> Fermi-Dirac function $f_{FD}(x; \mu, \sigma) = \frac{1}{1 + \exp(\frac{x - \mu}{\sigma})}$
function fermi_dirac(x, mu, sigma) result(f)
   ! Input variables
   real(dp), intent(in) :: x, mu, sigma
   ! Output variables
   real(dp) :: f

   f = 1.d0 / (1.d0 + exp((x - mu) / sigma))
   return
end function fermi_dirac

function compute_combined_metric(d2mdiff, iminus) result(metric)
   ! Constants
   real(dp), parameter :: MU = -3.0_dp
   real(dp), parameter :: SIGMA = 1.0_dp
   real(dp), parameter :: IMINUS_SHIFT = 1.0d-1
   ! Input variables
   real(dp), intent(in) :: d2mdiff, iminus
   ! Output variables
   real(dp) :: metric

   metric = log10(d2mdiff) * fermi_dirac(log10(iminus + IMINUS_SHIFT), MU, SIGMA)
   return
end function compute_combined_metric

!> Comupte $I^- = \int_0^\infty dr \max(0, \rho_{mod}(r) - \rho_c(r))$
function compute_iminus(rhoc, rhomod, rr, mmax) result(iminus)
   ! Input variables
   integer, intent(in) :: mmax
   real(dp), intent(in) :: rhoc(mmax), rhomod(mmax), rr(mmax)
   ! Output variables
   real(dp) :: iminus
   ! Local variables
   integer :: ii
   real(dp) :: integrand(mmax)

   do ii = 1, mmax
      integrand(ii) = max(0.d0, rr(ii)**2 * rhomod(ii) - rr(ii)**2 * rhoc(ii))
   end do
   iminus = nonuniform_trapezoidal(integrand, rr, mmax)
   return
end function compute_iminus

!> Non-uniform trapezoidal integration
function nonuniform_trapezoidal(fx, xx, nn) result(integral)
   ! Input variables
   integer, intent(in) :: nn
   real(dp), intent(in) :: fx(nn), xx(nn)
   ! Output variables
   real(dp) :: integral
   ! Local variables
   integer :: ii

   integral = 0.d0
   do ii = 1, nn - 1
      integral = integral + 0.5d0 * (fx(ii) + fx(ii + 1)) * (xx(ii + 1) - xx(ii))
   end do
   return
end function nonuniform_trapezoidal

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

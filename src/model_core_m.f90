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
module model_core_m
   use dpnint_m, only: dpnint
   use der2exc_m, only: der2exc
   use teter_m, only: gg1cc, gp1cc, gpp1cc
   implicit none
   private
   public :: modcore
   public :: modcore2
   public :: modcore3
contains
!> Creates monotonic polynomial model core charge matching all-electron
!> core charge and 4 derivatives at "crossover" radius.
!> Polynomial is 8th-order with no linear term.
!> Performs analysis and based on "hardness" criterion described in
!> Teter, Phys. Rev. B 48, 5031 (1993) , Appendix, as
subroutine modcore(icmod, rhops, rhotps, rhoc, rhoae, rhotae, rhomod, &
&                   fcfact, rcfact, irps, mmax, rr, nc, nv, la, zion, iexc)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> icmod  3 coefficient optimizaion, 4 for specivied fcfact and rfact
   integer, intent(in) :: icmod
   !> nv  number of valence states
   integer, intent(in) :: nv
   !> nc  number of core states
   integer, intent(in) :: nc
   !> iexc  exchange-correlation function to be used
   integer, intent(in) :: iexc
   !> irps  rr index of maximum rc
   integer, intent(in) :: irps
   !> mmax  dimension of log grid
   integer, intent(in) :: mmax
   !> la  angular-momenta
   integer, intent(in) :: la(30)
   !> rhoae  state-by-state all-electron valence charge density
   real(dp), intent(in) :: rhoae(mmax, nv)
   !> rhops  state-by-state pseudocharge density
   real(dp), intent(in) :: rhops(mmax, nv)
   !> rhotae  total all-electron valence charge density
   real(dp), intent(in) :: rhotae(mmax)
   !> rhotps  total pseudocharge density
   real(dp), intent(in) :: rhotps(mmax)
   !> rhoc  core-charge density
   real(dp), intent(in) :: rhoc(mmax)
   !> rr log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> zion  ion charge
   real(dp), intent(in) :: zion
   !> fcfact  prefactor for model amplitude (multiplies crossover value)
   real(dp), intent(in) :: fcfact
   !> rcfact  prefactor for model scale (multiplies crossover radius)
   real(dp), intent(in) :: rcfact

   !Output variables
   !> rhomod  model core density and 4 derivatives
   real(dp), intent(out) :: rhomod(mmax, 5)

   !convergence criterion
   real(dp), parameter :: eps = 1.0d-7

   !Local variables
   real(dp) :: a0, al, et, yy, gg, a0min, a0max, dermax, psum, sf, eeel, eexc
   real(dp) :: d2mdiff, rmatch, rhocmatch
   real(dp) :: emin, emax, sls, ss, hh, rgg, xx
   real(dp) :: aco(5), polym(5, 5), work(5, 5), constm(5, 5), xpow(9), fmatch(5)
   real(dp), allocatable :: vxcae(:), vxcpsp(:), vo(:), d2excae(:, :), d2excps(:, :)
   real(dp), allocatable :: dvxcae(:, :), dvxcps(:, :), vxct(:)
   integer :: ii, ierr, ircc, irmod, iter, jj, kk, ll, l1, mch
   integer :: ipvt(5), nodes(4)

   allocate (vxcae(mmax), vxcpsp(mmax), vo(mmax))
   allocate (dvxcae(mmax, nv), dvxcps(mmax, nv), vxct(mmax))
   allocate (d2excae(nv, nv), d2excps(nv, nv))

   d2excae(:, :) = 0.0d0
   d2excps(:, :) = 0.0d0
   ircc = 0

   do ii = mmax, 1, -1
      if (rhoc(ii) .gt. fcfact*rhotps(ii)) then
         ircc = ii
         exit
      end if
   end do

   if (ircc .eq. 0) then
      write (6, '(/a)') 'rhomod: ERROR ircc (core-valence charge crossover) &
      &        not found'
      stop
   end if

   !set limit for d2Exc calculation to radius beyond which rhomod=rhoc
   irmod = max(irps, ircc)

   write (6, '(/a/a)') 'Model core correction analysis',&
   &                  '  based on d2Exc/dn_idn_j'

   ! get derivatives of all-electron xc energy

   call der2exc(rhotae, rhoc, rhoae, rr, d2excae, d2excps, d2mdiff, &
   &                   zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excae - all-electron derivatives'
   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excae(kk, jj), jj=1, nv)
   end do

   ! set model charge to zero
   rhomod(:, :) = 0.0d0

   ! compute d2excps with no core correction
   rhomod(:, :) = 0.0d0

   call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
   &                   zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excps - pseudofunction derivatives with no core correction'
   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (6, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2mdiff

   al = 0.01d0*dlog(rr(101)/rr(1))

   ! core charge density for Louie-Froyen-Cohen correction

   rhomod(:, 1) = rhoc(:)

   xx = rr(ircc)
   fmatch(1) = rhoc(ircc)

   ! core charge derivatives
   do jj = 2, 5

      !7-point numerical first derivatives applied successively
      do ii = ircc - 25 + 3*jj, mmax - 3
         rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) + 9.d0*rhomod(ii - 2, jj - 1)&
         &     - 45.d0*rhomod(ii - 1, jj - 1) + 45.d0*rhomod(ii + 1, jj - 1)&
         &     - 9.d0*rhomod(ii + 2, jj - 1) + rhomod(ii + 3, jj - 1))&
         &     /(60.d0*al*rr(ii))
      end do
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
         constm(ii, jj) = constm(ii - 1, jj)*(constm(2, jj) - ii + 2)
      end do
   end do

   ! powers of xx=rr(ircc)

   xpow(1) = 0.0d0
   xpow(2) = 1.0d0
   do ii = 3, 9
      xpow(ii) = xx*xpow(ii - 1)
   end do

   ! polynomial matrix
   do jj = 1, 5
      do ii = 1, 5
         polym(ii, jj) = constm(ii, jj)*xpow(jj - ii + 5)
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
            xpow(ii) = rr(kk)*xpow(ii - 1)
         end do
         ii = 2
         psum = 0.0d0
         do jj = 1, 5
            psum = psum + constm(ii, jj)*xpow(jj - ii + 5)*aco(jj)
         end do
         dermax = dmax1(dermax, psum)
      end do

      ! test maximum derivative and adjust a0
      ! interval-halving search after achieving monotonicity to get
      ! barely monotonic model
      if (dermax .gt. 0.0d0) then
         a0min = a0
      else
         a0max = a0
      end if

      ! success when maximum derivative bracketed just above zero
      if (dermax > 0.0d0 .and. dermax < eps*dabs(fmatch(2))) exit

      if (a0max == 0.0d0) then
         a0 = 2.5d0*a0
      else
         a0 = 0.5d0*(a0max + a0min)
      end if

   end do !iterr
   if (iter > 50) write (6, '(/a/)') 'WARNING - modcore not conveged'

   ! fill in model and derivatives with fitted polynomial
   do kk = 1, ircc - 1
      do ii = 3, 9
         xpow(ii) = rr(kk)*xpow(ii - 1)
      end do
      do ii = 1, 5
         if (ii .eq. 1) then
            psum = a0
         else
            psum = 0.0d0
         end if
         do jj = 1, 5
            psum = psum + constm(ii, jj)*xpow(jj - ii + 5)*aco(jj)
         end do
         rhomod(kk, ii) = psum
      end do
   end do

   !test model
   call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
   &                   zion, iexc, nc, nv, la, irmod, mmax)

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
end subroutine modcore
!> Creates monotonic polynomial model core charge matching all-electron
!> core charge and 4 derivatives at "crossover" radius.
!> Polynomial is 8th-order with no linear term.
!> Performs analysis and based on "hardness" criterion described in
!> Teter, Phys. Rev. B 48, 5031 (1993) , Appendix, as
subroutine modcore2(icmod, rhops, rhotps, rhoc, rhoae, rhotae, rhomod, &
&                   fcfact, rcfact, irps, mmax, rr, nc, nv, la, zion, iexc)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> icmod  3 coefficient optimizaion, 4 for specivied fcfact and rfact
   integer :: icmod
   !> nv  number of valence states
   integer :: nv
   !> nc  number of core states
   integer :: nc
   !> iexc  exchange-correlation function to be used
   integer :: iexc
   !> irps  rr index of maximum rc
   integer :: irps
   !> mmax  dimension of log grid
   integer :: mmax
   !> la  angular-momenta
   integer :: la(30)
   !> rhoae  state-by-state all-electron valence charge density
   real(dp) :: rhoae(mmax, nv)
   !> rhops  state-by-state pseudocharge density
   real(dp) :: rhops(mmax, nv)
   !> rhotae  total all-electron valence charge density
   real(dp) :: rhotae(mmax)
   !> rhotps  total pseudocharge density
   real(dp) :: rhotps(mmax)
   !> rhoc  core-charge density
   real(dp) :: rhoc(mmax)
   !> rr log radial grid
   real(dp) :: rr(mmax)
   !> zion  ion charge
   real(dp) :: zion
   !> fcfact  prefactor for model amplitude (multiplies crossover value)
   real(dp) :: fcfact
   !> rcfact  prefactor for model scale (multiplies crossover radius)
   real(dp) :: rcfact

   !Output variables
   !> rhomod  model core density and 4 derivatives
   real(dp) :: rhomod(mmax, 5)

   !convergence criterion
   real(dp), parameter :: eps = 1.0d-7

   !Local variables
   real(dp) :: al, eeel, eexc
   real(dp) :: d2mdiff, rmatch, rhocmatch
   real(dp) :: xx, yy, dy
   real(dp) :: x0max, x0min, a0, b0, r0, tt, ymatch, ytrial
   real(dp) :: fmatch(5)
   real(dp) :: drint, rtst, rint(20), fint(20) !ad-hoc smoothing variables
   real(dp), allocatable :: vxcae(:), vxcpsp(:), vo(:), d2excae(:, :), d2excps(:, :)
   real(dp), allocatable :: dvxcae(:, :), dvxcps(:, :), vxct(:)
   integer :: ii, ierr, ircc, irmod, iter, jj, kk
   integer :: iint !ad-hoc smoothing variables

   allocate (vxcae(mmax), vxcpsp(mmax), vo(mmax))
   allocate (dvxcae(mmax, nv), dvxcps(mmax, nv), vxct(mmax))
   allocate (d2excae(nv, nv), d2excps(nv, nv))

   d2excae(:, :) = 0.0d0
   d2excps(:, :) = 0.0d0

   ! find valence pseudocharge - core charge crossover
   ! this is needed for compatability with icmod=3 option
   ircc = 0
   do ii = mmax, 1, -1
      if (rhoc(ii) .gt. rhotps(ii)) then
         ircc = ii
         rmatch = rr(ircc)
         rhocmatch = rhoc(ircc)
         exit
      end if
   end do

   ! find scaled valence pseudocharge - core charge crossover
   ircc = 0
   do ii = mmax, 1, -1
      if (rhoc(ii) .gt. fcfact*rhotps(ii)) then
         ircc = ii
         exit
      end if
   end do

   if (ircc .eq. 0) then
      write (6, '(/a)') 'modcore2: ERROR ircc (core-valence charge crossover) &
      &        not found'
      stop
   end if

   !set limit for d2Exc calculation to radius beyond which rhomod=rhoc
   irmod = mmax

   write (6, '(/a/a)') 'Model core correction analysis',&
   &                  '  based on d2Exc/dn_idn_j'

   ! get derivatives of all-electron xc energy

   call der2exc(rhotae, rhoc, rhoae, rr, d2excae, d2excps, d2mdiff, &
   &                   zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excae - all-electron derivatives'
   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excae(kk, jj), jj=1, nv)
   end do

   ! set model charge to zero
   rhomod(:, :) = 0.0d0

   ! compute d2excps with no core correction
   rhomod(:, :) = 0.0d0

   call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
   &                   zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excps - pseudofunction derivatives with no core correction'
   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (6, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2mdiff

   al = 0.01d0*dlog(rr(101)/rr(1))

   ! core charge density for Louie-Froyen-Cohen correction

   rhomod(:, 1) = rhoc(:)

   xx = rr(ircc)
   fmatch(1) = rhoc(ircc)

   ! core charge derivatives

   !7-point numerical first derivative
   jj = 2; ii = ircc
   rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) + 9.d0*rhomod(ii - 2, jj - 1)&
   &     - 45.d0*rhomod(ii - 1, jj - 1) + 45.d0*rhomod(ii + 1, jj - 1)&
   &     - 9.d0*rhomod(ii + 2, jj - 1) + rhomod(ii + 3, jj - 1))&
   &     /(60.d0*al*rr(ii))
   fmatch(jj) = rhomod(ircc, jj)

   ! Fit Teter function to value and slope at ircc

   ymatch = rr(ircc)*fmatch(2)/fmatch(1)

   ! interval-halving search for dimensionless match point

   x0max = 1.48d0
   x0min = 0.0d0

   do jj = 1, 50

      xx = 0.5d0*(x0max + x0min)

      call gg1cc(yy, xx)
      call gp1cc(dy, xx)
      ytrial = xx*dy/yy

      if (abs(ytrial - ymatch) < eps) exit

      if (ytrial < ymatch) then
         x0max = xx
      else
         x0min = xx
      end if
   end do

   b0 = xx/rr(ircc)
   a0 = fmatch(1)/yy

   write (6, '(/a)') 'amplitude prefactor, scale prefactor'
   write (6, '(2f10.4)') a0/rhocmatch, 1.0d0/(b0*rmatch)

   rhomod(:, :) = 0.0d0
   do ii = 1, mmax
      xx = b0*rr(ii)
      if (xx < 3.0d0) then
         call gg1cc(yy, xx)
         rhomod(ii, 1) = a0*yy
         call gp1cc(yy, xx)
         rhomod(ii, 2) = a0*yy*b0
         call gpp1cc(yy, xx)
         rhomod(ii, 3) = a0*yy*b0**2
      end if
   end do

   ! test model
   call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
   &                   zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excps - pseudofunction derivatives with core correction'

   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (6, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2mdiff

   !7-point numerical first derivatives applied successively

   ! do jj=4,5
   !  do ii=3*jj-8,mmax-3
   !    if(rhomod(ii,1)==0.0d0) exit
   !     rhomod(ii,jj)=(-rhomod(ii-3,jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
   !&     -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
   !&     -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
   !&     /(60.d0*al*rr(ii))
   !  end do
   ! end do

   !ad-hoc treatment of numerical noise near origin
   !set up a mesh on which 2nd derivative will have a stable
   !polynomial representation
   !assumes dpnint remains 7th order
   drint = 0.05d0*rr(ircc)
   rtst = 0.5d0*drint
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
   do ii = 3*jj - 8, mmax - 3
      rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) + 9.d0*rhomod(ii - 2, jj - 1)&
      &     - 45.d0*rhomod(ii - 1, jj - 1) + 45.d0*rhomod(ii + 1, jj - 1)&
      &     - 9.d0*rhomod(ii + 2, jj - 1) + rhomod(ii + 3, jj - 1))&
      &     /(60.d0*al*rr(ii))
   end do

   !set up a mesh on which 3rd derivative will have a stable
   drint = 0.95d0*drint
   rtst = 0.5d0*drint
   rtst = 0.5d0*drint
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
   do ii = 3*jj - 8, mmax - 3
      rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) + 9.d0*rhomod(ii - 2, jj - 1)&
      &     - 45.d0*rhomod(ii - 1, jj - 1) + 45.d0*rhomod(ii + 1, jj - 1)&
      &     - 9.d0*rhomod(ii + 2, jj - 1) + rhomod(ii + 3, jj - 1))&
      &     /(60.d0*al*rr(ii))
   end do

   !set up a mesh on which 4th derivative will have a stable
   drint = 0.95d0*drint
   rtst = 0.5d0*drint
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
end subroutine modcore2
!> Creates monotonic polynomial model core charge matching all-electron
!> core charge and 4 derivatives at "crossover" radius.
!> Polynomial is 8th-order with no linear term.
!> Performs analysis and based on "hardness" criterion described in
!> Teter, Phys. Rev. B 48, 5031 (1993) , Appendix, as
subroutine modcore3(icmod, rhops, rhotps, rhoc, rhoae, rhotae, rhomod, &
&                   fcfact, rcfact, irps, mmax, rr, nc, nv, la, zion, iexc)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> icmod  3 coefficient optimizaion, 4 for specivied fcfact and rfact
   integer, intent(in) :: icmod
   !> nv  number of valence states
   integer, intent(in) :: nv
   !> nc  number of core states
   integer, intent(in) :: nc
   !> iexc  exchange-correlation function to be used
   integer, intent(in) :: iexc
   !> irps  rr index of maximum rc
   integer, intent(in) :: irps
   !> mmax  dimension of log grid
   integer, intent(in) :: mmax
   !> la  angular-momenta
   integer, intent(in) :: la(30)
   !> rhoae  state-by-state all-electron valence charge density
   real(dp), intent(in) :: rhoae(mmax, nv)
   !> rhops  state-by-state pseudocharge density
   real(dp), intent(in) :: rhops(mmax, nv)
   !> rhotae  total all-electron valence charge density
   real(dp), intent(in) :: rhotae(mmax)
   !> rhotps  total pseudocharge density
   real(dp), intent(in) :: rhotps(mmax)
   !> rhoc  core-charge density
   real(dp), intent(in) :: rhoc(mmax)
   !> rr log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> zion  ion charge
   real(dp), intent(in) :: zion
   !> fcfact  prefactor for model amplitude (multiplies crossover value)
   real(dp), intent(in out) :: fcfact
   !> rcfact  prefactor for model scale (multiplies crossover radius)
   real(dp), intent(in) :: rcfact

   !Output variables
   !> rhomod  model core density and 4 derivatives
   real(dp), intent(out) :: rhomod(mmax, 5)

   !convergence criterion
   real(dp), parameter :: eps = 1.0d-7
   real(dp), parameter :: blend = 2.0d0

   !Local variables
   real(dp) :: al, eeel, eexc
   real(dp) :: d2mdiff, rmatch, rhocmatch, r0, rcross
   real(dp) :: gg, tt, yy
   real(dp) :: drint, rtst, rint(20), fint(20) !ad-hoc smoothing variables
   real(dp), allocatable :: vxcae(:), vxcpsp(:), vo(:), d2excae(:, :), d2excps(:, :)
   real(dp), allocatable :: dvxcae(:, :), dvxcps(:, :), vxct(:)
   integer :: ii, ierr, ircc, ircross, irmod, iter, jj, kk
   integer :: iint !ad-hoc smoothing variables

   !2-dimensional Nelder-Mead variables
   real(dp), parameter :: alpha_nm = 1.0d0
   real(dp), parameter :: gamma_nm = 2.0d0
   real(dp), parameter :: rho_nm = -0.5d0
   real(dp), parameter :: sigma_nm = 0.5d0
   real(dp) :: xx(2, 3), ff(3), xt(2), ft, x0(2), xr(2), fr, xe(2), fe, xc(2), fc
   real(dp) :: d2ref, fta(10), d2min

   allocate (vxcae(mmax), vxcpsp(mmax), vo(mmax))
   allocate (dvxcae(mmax, nv), dvxcps(mmax, nv), vxct(mmax))
   allocate (d2excae(nv, nv), d2excps(nv, nv))

   d2excae(:, :) = 0.0d0
   d2excps(:, :) = 0.0d0

   !set limit for d2Exc calculation to radius beyond which rhomod=rhoc
   irmod = mmax

   write (6, '(/a/a)') 'Model core correction analysis',&
   &                  '  based on d2Exc/dn_idn_j'

   call der2exc(rhotae, rhoc, rhoae, rr, d2excae, d2excps, d2mdiff, &
   &                   zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excae - all-electron derivatives'
   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excae(kk, jj), jj=1, nv)
   end do

   ! set model charge to zero
   rhomod(:, :) = 0.0d0

   ! compute d2excps with no core correction
   rhomod(:, :) = 0.0d0

   call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
   &                   zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excps - pseudofunction derivatives with no core correction'
   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (6, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2mdiff

   d2ref = d2mdiff

   ! find valence pseudocharge - core charge crossover
   ircc = 0
   do ii = mmax, 1, -1
      if (rhoc(ii) .gt. rhotps(ii)) then
         ircc = ii
         rmatch = rr(ircc)
         rhocmatch = rhoc(ircc)
         exit
      end if
   end do

   if (ircc .eq. 0) then
      write (6, '(/a)') 'modcore3: ERROR ircc (core-valence charge crossover) &
      &        not found'
      stop
   else
      write (6, '(/a,2f10.4)') 'rmatch, rhocmatch', rmatch, rhocmatch
   end if

   gg = 0.d0
   yy = 0.d0

   !option for optimization
   if (icmod == 4) then

      fcfact = 1.0d0

      !Coarse-grid search for minimum

      write (6, '(/a)') 'Coarse scan for minimum error'
      write (6, '(a)') '  matrix elements: rms 2nd-derivative errors (mHa)'
      write (6, '(a)') '  column index : amplitude prefactor to rhocmatch'
      write (6, '(a)') '  row index : scale prefactor to rmatch'

      d2min = 10.0d0
      write (6, '(/7x,10f7.3/)') (1.5d0 + 0.5d0*(jj - 1), jj=1, 10)
      do kk = 1, 10
         xt(2) = (1.0d0 + 0.1d0*(kk - 1))*rmatch
         do jj = 1, 10
            xt(1) = (1.5d0 + 0.5d0*(jj - 1))*rhocmatch

            r0 = 1.5d0*xt(2)
            do ii = mmax, 1, -1
               if (rr(ii) < r0) then
                  call gg1cc(gg, yy)
                  if (xt(1)*gg < rhoc(ii)) then
                     rcross = rr(ii)
                     exit
                  end if
               end if
            end do

            do ii = 1, mmax
               yy = rr(ii)/xt(2)
               call gg1cc(gg, yy)
               rhomod(ii, 1) = xt(1)*gg
            end do

            call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
            &                     zion, iexc, nc, nv, la, irmod, mmax)
            fta(jj) = 1.0d3*d2mdiff

            if (d2mdiff < d2min) then
               xx(:, 1) = xt(:)
               d2min = d2mdiff
            end if

         end do
         write (6, '(f5.1,f9.3,9f7.3)') 1.0d0 + 0.1d0*(kk - 1), (fta(jj), jj=1, 10)
      end do

      !Initial Nelder-Mead simplex from coarse-search minimum

      xx(1, 2) = xx(1, 1) + 0.25d0*rhocmatch
      xx(2, 2) = xx(2, 1)
      xx(1, 3) = xx(1, 1)
      xx(2, 3) = xx(2, 1) + 0.05d0*rmatch
      xx(1, 1) = xx(1, 1) - 0.125d0*rhocmatch
      xx(2, 1) = xx(2, 1) - 0.025d0*rmatch

      !Fill function values for initial simplex

      do kk = 1, 3
         xt(:) = xx(:, kk)

         r0 = 1.5d0*xt(2)
         do ii = mmax, 1, -1
            if (rr(ii) < r0) then
               call gg1cc(gg, yy)
               if (xt(1)*gg < rhoc(ii)) then
                  rcross = rr(ii)
                  exit
               end if
            end if
         end do

         do ii = 1, mmax
            yy = rr(ii)/xt(2)
            call gg1cc(gg, yy)
            tt = (rr(ii) - r0 - 2.0d0*rcross)/(r0 - rcross)
            rhomod(ii, 1) = xt(1)*gg
         end do

         call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
         &                    zion, iexc, nc, nv, la, irmod, mmax)
         ff(kk) = d2mdiff

         !  write(6,'(i4,a,2f10.4,1p,e14.4)') kk,'  xt',xt(1),xt(2),ff(kk)
      end do

      !Nelder-Mead iteration loop
      write (6, '(/a)') 'Nelder-Mead iteration'
      do kk = 1, 101

         !(1) Order (dumb bubble sort)
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

         !stopping criterion
         if (ff(3) - ff(1) < 1.0d-4*d2ref) then
            write (6, '(a,i4,a)') ' converged in', kk - 1, ' steps'
            write (6, '(a)') 'amplitude prefactor, scale prefactor'
            write (6, '(2f10.4)') xx(1, 1)/rhocmatch, xx(2, 1)/rmatch
            exit
         else if (kk == 101) then
            write (6, '(a,i4,a)') ' WARNING: not fully converged in 100 steps'
            write (6, '(a)') 'amplitude prefactor, scale prefactor'
            write (6, '(2f10.4)') xx(1, 1)/rhocmatch, xx(2, 1)/rmatch
            exit
         end if

         !(2) Centroid
         x0(:) = 0.5d0*(xx(:, 1) + xx(:, 2))

         !(3) Reflection
         xr(:) = x0(:) + alpha_nm*(x0(:) - xx(:, 3))

         r0 = 1.5d0*xr(2)
         do ii = mmax, 1, -1
            if (rr(ii) < r0) then
               call gg1cc(gg, yy)
               if (xr(1)*gg < rhoc(ii)) then
                  rcross = rr(ii)
                  exit
               end if
            end if
         end do

         do ii = 1, mmax
            yy = rr(ii)/xr(2)
            call gg1cc(gg, yy)
            tt = (rr(ii) - r0 - 2.0d0*rcross)/(r0 - rcross)
            rhomod(ii, 1) = xr(1)*gg
         end do

         call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
         &                    zion, iexc, nc, nv, la, irmod, mmax)
         fr = d2mdiff
         ! write(6,'(i4,a,2f10.4,1p,e14.4)') kk,'  xr',xr(1),xr(2),fr

         if (ff(1) <= fr .and. fr < ff(2)) then
            ff(3) = fr; xx(:, 3) = xr(:)
            cycle !kk loop
         end if

         !(4) Expansion
         if (fr < ff(1)) then
            xe(:) = x0(:) + gamma_nm*(x0(:) - xx(:, 3))

            r0 = 1.5d0*xe(2)
            do ii = mmax, 1, -1
               if (rr(ii) < r0) then
                  call gg1cc(gg, yy)
                  if (xe(1)*gg < rhoc(ii)) then
                     rcross = rr(ii)
                     exit
                  end if
               end if
            end do

            do ii = 1, mmax
               yy = rr(ii)/xe(2)
               call gg1cc(gg, yy)
               tt = (rr(ii) - r0 - 2.0d0*rcross)/(r0 - rcross)
               rhomod(ii, 1) = xe(1)*gg
            end do

            call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
            &                     zion, iexc, nc, nv, la, irmod, mmax)
            fe = d2mdiff
            ! write(6,'(i4,a,2f10.4,1p,e14.4)') kk,'  xe',xe(1),xe(2),fe

            if (fe < fr) then
               ff(3) = fe; xx(:, 3) = xe(:)
               cycle !kk
            else
               ff(3) = fr; xx(:, 3) = xr(:)
               cycle !kk
            end if
         end if

         !(5) Contraction
         xc(:) = x0(:) + rho_nm*(x0(:) - xx(:, 3))

         r0 = 1.5d0*xc(2)
         do ii = mmax, 1, -1
            if (rr(ii) < r0) then
               call gg1cc(gg, yy)
               if (xc(1)*gg < rhoc(ii)) then
                  rcross = rr(ii)
                  exit
               end if
            end if
         end do

         do ii = 1, mmax
            yy = rr(ii)/xc(2)
            call gg1cc(gg, yy)
            tt = (rr(ii) - r0 - 2.0d0*rcross)/(r0 - rcross)
            rhomod(ii, 1) = xc(1)*gg
         end do

         call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
         &                    zion, iexc, nc, nv, la, irmod, mmax)
         fc = d2mdiff
         !  write(6,'(i4,a,2f10.4,1p,e14.4)') kk,'  xc',xc(1),xc(2),fc
         if (fc < ff(3)) then
            ff(3) = fc; xx(:, 3) = xc(:)
            cycle !kk
         end if

         !(6) Reduction
         do jj = 2, 3
            xx(:, jj) = xx(:, 1) + sigma_nm*(xx(:, jj) - xx(:, 1))

            r0 = 1.5d0*xx(2, jj)
            do ii = mmax, 1, -1
               if (rr(ii) < r0) then
                  call gg1cc(gg, yy)
                  if (xx(1, jj)*gg < rhoc(ii)) then
                     rcross = rr(ii)
                     exit
                  end if
               end if
            end do

            do ii = 1, mmax
               yy = rr(ii)/xx(2, jj)
               call gg1cc(gg, yy)
               tt = (rr(ii) - r0 - 2.0d0*rcross)/(r0 - rcross)
               rhomod(ii, 1) = xx(1, jj)*gg
            end do

            call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
            &                     zion, iexc, nc, nv, la, irmod, mmax)
            ff(jj) = d2mdiff
            !  write(6,'(i4,a,2f10.4,1p,e14.4)') kk,' xrd',xx(1,jj),xx(2,jj),ff(jj)
         end do !jj

      end do !kk

      write (6, '(/a)') 'Optimized Teter model core charge'

      !option for specifying prefactors as input
   else if (icmod == 3) then
      xx(1, 1) = fcfact*rhocmatch
      xx(2, 1) = rcfact*rmatch

      write (6, '(/a/)') &
      &       'Teter function model core charge with specified parameters'
   end if

   r0 = 1.5d0*xx(2, 1)
   do ii = mmax, 1, -1
      if (rr(ii) < r0) then
         call gg1cc(gg, yy)
         if (xx(1, 1)*gg < rhoc(ii)) then
            rcross = rr(ii)
            ircross = ii
            exit
         end if
      end if
   end do

   ! blend the Teter function tail into the all-electron rhoc
   ! first two derivatives are filled in analytically before the blend starts
   rhomod(:, :) = 0.0d0
   do ii = 1, mmax
      yy = rr(ii)/xx(2, 1)
      call gg1cc(gg, yy)
      tt = (rr(ii) - r0 - 2.0d0*rr(ircross))/(r0 - rr(ircross))

      rhomod(ii, 1) = xx(1, 1)*gg
      if (ii < ircross) then
         call gp1cc(gg, yy)
         rhomod(ii, 2) = xx(1, 1)*gg/xx(2, 1)
         call gpp1cc(gg, yy)
         rhomod(ii, 3) = xx(1, 1)*gg/xx(2, 1)**2
      end if
   end do

   call der2exc(rhotps, rhomod(1, 1), rhops, rr, d2excps, d2excae, d2mdiff, &
   &                   zion, iexc, nc, nv, la, irmod, mmax)

   write (6, '(/a/)') 'd2excps - pseudofunction derivatives with core correction'
   do kk = 1, nv
      write (6, '(1p,4d16.6)') (d2excps(kk, jj), jj=1, nv)
   end do
   write (6, '(/a,1p,e16.6)') 'rms 2nd derivative error', d2mdiff

   !7-point numerical first derivatives applied successively
   !skip non-blended section for 1st and 2nd derivatives

   al = 0.01d0*dlog(rr(101)/rr(1))

   do jj = 2, 3
      do ii = ircross - 6, mmax - 3
         rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) + 9.d0*rhomod(ii - 2, jj - 1)&
         &     - 45.d0*rhomod(ii - 1, jj - 1) + 45.d0*rhomod(ii + 1, jj - 1)&
         &     - 9.d0*rhomod(ii + 2, jj - 1) + rhomod(ii + 3, jj - 1))&
         &     /(60.d0*al*rr(ii))
      end do
   end do

   !ad-hoc treatment of numerical noise near origin
   !set up a mesh on which 2nd derivative will have a stable
   !polynomial representation
   !assumes dpnint remains 7th order
   drint = 0.02d0*rr(ircross)
   rtst = 0.5d0*drint
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
   do ii = 3*jj - 8, mmax - 3
      rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) + 9.d0*rhomod(ii - 2, jj - 1)&
      &     - 45.d0*rhomod(ii - 1, jj - 1) + 45.d0*rhomod(ii + 1, jj - 1)&
      &     - 9.d0*rhomod(ii + 2, jj - 1) + rhomod(ii + 3, jj - 1))&
      &     /(60.d0*al*rr(ii))
   end do

   !set up a mesh on which 3rd derivative will have a stable
   drint = 0.95d0*drint
   rtst = 0.5d0*drint
   rtst = 0.5d0*drint
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
   do ii = 3*jj - 8, mmax - 3
      rhomod(ii, jj) = (-rhomod(ii - 3, jj - 1) + 9.d0*rhomod(ii - 2, jj - 1)&
      &     - 45.d0*rhomod(ii - 1, jj - 1) + 45.d0*rhomod(ii + 1, jj - 1)&
      &     - 9.d0*rhomod(ii + 2, jj - 1) + rhomod(ii + 3, jj - 1))&
      &     /(60.d0*al*rr(ii))
   end do

   !set up a mesh on which 4th derivative will have a stable
   drint = 0.95d0*drint
   rtst = 0.5d0*drint
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
end subroutine modcore3
end module model_core_m

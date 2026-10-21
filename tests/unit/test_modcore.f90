! Unit tests of the pieces of the model core charges at continuous radii
! (modcore_m; M-1, M-3, M-4, M-5 of doc/smoothness_plan.md), written before
! the code, on the program's default mesh (amesh=1.006, rr1=5e-5, rmax=45):
!   teter      the Teter function g and its derivatives to the 4th, with no
!              cutoff: against the closed form and gg1cc, gp1cc, gpp1cc away
!              from the removable singularities (1e-13 relative), against
!              their series constants at x = 0, 1/2, 1 (1e-12), each
!              derivative against an 8th-order difference of the one below
!              (1e-11 of (2 pi)**k) from 0 to 6, g > 0 beyond x = 3; and
!              smooth: the 4th derivative on 41-point windows of width 0.04
!              from 0 to 3.5 fits a polynomial of degree 8 to 1e-11
!   crossover  the outermost root of ln(rhoc) - ln(fcfact rhotps) for
!              analytic densities, one crossing and three, to 1e-12 bohr;
!              smooth in fcfact (41 values, degree 6, 1e-12); no crossing
!              gives ierr=1
!   poly       the icmod 1 polynomial matches the value and 4 derivatives at
!              rstar (1e-11 relative), its slope over r**2 is <= 1e-12 of
!              |rho'|/rstar**2 on [0, rstar] and rises above 0 with a0
!              lowered by 1e-8 (at 0, if that is where it touches 0), a0
!              agrees with a bisection on a fine grid to 1e-10, and is
!              smooth in rstar (41 values over one mesh step, 1e-11)
!   scale      the argument of the Teter function at a given logarithmic
!              slope to 1e-13; a slope out of range gives ierr=1
! Stops with a non-zero code if any check failed.
 program test_modcore
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
 use modcore_m, only: teter_derivs, core_crossover, poly_core_fit, &
&                     teter_scale
 implicit none

!Local variables
 real(dp), parameter :: pi = 3.14159265358979323846_dp
 real(dp), parameter :: rr1 = 5.0e-5_dp, rmax = 45.0_dp, amesh = 1.006_dp
 real(dp), allocatable :: rr(:), rhoc(:), rhotps(:)
 real(dp) :: al, gd(0:4), gl, xx, err, rstar, ref, ts(41), ys(41)
 real(dp) :: resid(41), fmatch(5), a0, aco(5), a0ref
 real(dp) :: xs(5), x0s(5)
 integer :: nfail = 0, mmax, ii, jj, kk, ierr

 al = log(amesh)
 mmax = int(log(rmax/rr1)/al)
 allocate(rr(mmax), rhoc(mmax), rhotps(mmax))
 do ii = 1, mmax
   rr(ii) = rr1*exp(al*(ii - 1))
 end do

! teter: the closed form and the legacy routines, away from 0, 1/2, 1
 xs = [0.13_dp, 0.37_dp, 0.71_dp, 1.23_dp, 2.61_dp]
 err = 0.0_dp
 do ii = 1, size(xs)
   call teter_derivs(xs(ii), 2, gd(0:2))
   err = max(err, abs(gd(0) - gdirect(xs(ii)))/gdirect(xs(ii)))
   call gg1cc(gl, xs(ii))
   err = max(err, abs(gd(0) - gl)/gl)
   call gp1cc(gl, xs(ii))
   err = max(err, abs(gd(1) - gl)/abs(gl))
   call gpp1cc(gl, xs(ii))
   err = max(err, abs(gd(2) - gl)/abs(gl))
 end do
 call check('teter closed form', err, 1.0e-13_dp)

! teter: the series constants of gg1cc and gpp1cc
 call teter_derivs(0.0_dp, 4, gd)
 err = max(abs(gd(0) - 1.0_dp), abs(gd(1)), abs(gd(3))/100.0_dp, &
&          abs(gd(2) - (20.0_dp - 8.0_dp*pi**2/3.0_dp))/10.0_dp, &
&          abs(gd(4) - 2.0_dp*(804.0_dp - 160.0_dp*pi**2 &
&                              + 128.0_dp*pi**4/15.0_dp))/1.0e3_dp)
 call check('teter series at 0', err, 1.0e-12_dp)
 call teter_derivs(0.5_dp, 4, gd)
 err = max(abs(gd(0) - 4.0_dp/9.0_dp), abs(gd(1) + 40.0_dp/27.0_dp), &
&          abs(gd(2) - (40.0_dp/3.0_dp - 32.0_dp*pi**2/27.0_dp))/10.0_dp, &
&          abs(gd(3) - (-8320.0_dp/81.0_dp + 320.0_dp*pi**2/27.0_dp)) &
&          /100.0_dp, &
&          abs(gd(4) - 2.0_dp*(157712.0_dp/243.0_dp - 320.0_dp*pi**2/3.0_dp &
&                              + 512.0_dp*pi**4/135.0_dp))/1.0e3_dp)
 call check('teter series at 1/2', err, 1.0e-12_dp)
 call teter_derivs(1.0_dp, 4, gd)
 err = max(abs(gd(0) - 1.0_dp/36.0_dp), abs(gd(1) + 25.0_dp/108.0_dp), &
&          abs(gd(2) - (485.0_dp/216.0_dp - 2.0_dp*pi**2/27.0_dp))/10.0_dp, &
&          abs(gd(3) - (-4055.0_dp/162.0_dp + 50.0_dp*pi**2/27.0_dp)) &
&          /100.0_dp, &
&          abs(gd(4) - 2.0_dp*(616697.0_dp/3888.0_dp - 485.0_dp*pi**2/27.0_dp &
&                              + 32.0_dp*pi**4/135.0_dp))/1.0e3_dp)
 call check('teter series at 1', err, 1.0e-12_dp)

! teter: each derivative against the difference of the one below
 err = 0.0_dp
 do ii = 0, 600
   xx = 0.01_dp*ii + 0.003_dp
   do kk = 1, 4
     call teter_derivs(xx, kk, gd(0:kk))
     err = max(err, abs(gd(kk) - fdiff(xx, kk - 1))/(2.0_dp*pi)**kk)
   end do
 end do
 call check('teter derivatives', err, 1.0e-11_dp)

! teter: no cutoff
 err = 0.0_dp
 do ii = 1, 4
   xx = 3.0_dp + 0.37_dp*ii
   call teter_derivs(xx, 0, gd(0:0))
   if(.not. gd(0) > 0.0_dp) err = 1.0_dp
   err = max(err, abs(gd(0) - gdirect(xx)))
 end do
 call check('teter beyond 3', err, 1.0e-15_dp)

! teter: smooth on windows of width 0.04
 err = 0.0_dp
 do jj = 0, 174
   do ii = 1, 41
     ts(ii) = 2.0_dp*(ii - 1)/40.0_dp - 1.0_dp
     xx = 0.02_dp*jj + 0.02_dp*(ts(ii) + 1.0_dp)
     call teter_derivs(xx, 4, gd)
     ys(ii) = gd(4)
   end do
   call fit_residual(ts, ys, 8, resid)
   err = max(err, maxval(abs(resid))/(2.0_dp*pi)**4)
 end do
 call check('teter smooth', err, 1.0e-11_dp)

! crossover: rhoc = 40 exp(-4r), rhotps = 2 r**2 exp(-r), one crossing
 rhoc(:) = 40.0_dp*exp(-4.0_dp*rr)
 rhotps(:) = 2.0_dp*rr**2*exp(-rr)
 call core_crossover(rr, rhoc, rhotps, mmax, 0.25_dp, rstar, ierr)
 ref = outer_root(0.25_dp, 0.0_dp)
 call check('crossover one', abs(rstar - ref) + ierr, 1.0e-12_dp)
 do ii = 1, 41
   ts(ii) = 2.0_dp*(ii - 1)/40.0_dp - 1.0_dp
   call core_crossover(rr, rhoc, rhotps, mmax, 0.25_dp + 0.01_dp*ts(ii), &
&                      ys(ii), ierr)
 end do
 call fit_residual(ts, ys, 6, resid)
 call check('crossover smooth', maxval(abs(resid))/maxval(abs(ys)), &
&           1.0e-12_dp)

! crossover: a shell of core charge at 2 bohr, three crossings
 rhoc(:) = rhoc(:) + exp(-(rr - 2.0_dp)**2/0.08_dp)
 call core_crossover(rr, rhoc, rhotps, mmax, 0.25_dp, rstar, ierr)
 ref = outer_root(0.25_dp, 1.0_dp)
 call check('crossover outermost', abs(rstar - ref) + ierr, 1.0e-12_dp)
 if(ref < 1.9_dp) call check('crossover reference', 1.0_dp, 0.0_dp)

! crossover: none
 call core_crossover(rr, 0.0_dp*rhoc, rhotps, mmax, 0.25_dp, rstar, ierr)
 call check('crossover none', real(abs(ierr - 1), dp), 0.0_dp)

! poly: the core 40 exp(-4r) matched at 1.1 bohr, and 3/(1+r**4) at 0.8
 do jj = 1, 2
   if(jj == 1) then
     rstar = 1.1_dp
     fmatch = 40.0_dp*exp(-4.0_dp*rstar)*[1.0_dp, -4.0_dp, 16.0_dp, &
&                                          -64.0_dp, 256.0_dp]
   else
     rstar = 0.8_dp
     fmatch = lorentz(rstar)
   end if
   call poly_core_fit(rstar, fmatch, a0, aco)
   err = 0.0_dp
   do kk = 0, 4
     err = max(err, abs(pval(a0, aco, rstar, kk) - fmatch(kk + 1)) &
&                   /abs(fmatch(kk + 1)))
   end do
   call check('poly match', err, 1.0e-11_dp)
   call check('poly monotone', max(0.0_dp, max_slope(rstar, a0, aco)) &
&             *rstar**2/abs(fmatch(2)), 1.0e-12_dp)
   call coefficients(rstar, fmatch, a0 - 1.0e-8_dp*abs(a0), aco)
   if(.not. max_slope(rstar, a0 - 1.0e-8_dp*abs(a0), aco) > 0.0_dp) then
     call check('poly barely monotone', 1.0_dp, 0.0_dp)
   end if
   a0ref = bisect_a0(rstar, fmatch)
   call check('poly a0', abs(a0 - a0ref)/abs(a0ref), 1.0e-10_dp)
 end do
 do ii = 1, 41
   ts(ii) = 2.0_dp*(ii - 1)/40.0_dp - 1.0_dp
   rstar = 1.1_dp*amesh**(0.5_dp*(ts(ii) + 1.0_dp))
   fmatch = 40.0_dp*exp(-4.0_dp*rstar)*[1.0_dp, -4.0_dp, 16.0_dp, &
&                                        -64.0_dp, 256.0_dp]
   call poly_core_fit(rstar, fmatch, ys(ii), aco)
 end do
 call fit_residual(ts, ys, 6, resid)
 call check('poly smooth', maxval(abs(resid))/maxval(abs(ys)), 1.0e-11_dp)

! scale: x g'(x)/g(x) at x0, back to x0
 x0s = [0.2_dp, 0.45_dp, 0.7_dp, 1.1_dp, 1.4_dp]
 err = 0.0_dp
 do ii = 1, size(x0s)
   call teter_scale(log_slope(x0s(ii)), xx, ierr)
   err = max(err, abs(xx - x0s(ii)) + ierr)
 end do
 call check('scale', err, 1.0e-13_dp)
 call teter_scale(0.5_dp, xx, ierr)
 call check('scale out of range', real(abs(ierr - 1), dp), 0.0_dp)

 if(nfail > 0) then
   write(stderr, '(i0,a)') nfail, ' check(s) failed'
   error stop 1
 end if
 write(stderr, '(a)') 'modcore: all checks passed'

 contains

 real(dp) function gdirect(xq)

!! the Teter function in closed form, away from x = 0, 1/2, 1

 real(dp), intent(in) :: xq

 gdirect = (sin(2.0_dp*pi*xq)/(2.0_dp*pi*xq*(1.0_dp - 4.0_dp*xq**2) &
&          *(1.0_dp - xq**2)))**2
 end function gdirect

 real(dp) function log_slope(xq)

!! x g'(x)/g(x) of the Teter function in closed form, away from 1/2 and 1

 real(dp), intent(in) :: xq

 log_slope = 2.0_dp*xq*(2.0_dp*pi/tan(2.0_dp*pi*xq) - 1.0_dp/xq &
&            + 8.0_dp*xq/(1.0_dp - 4.0_dp*xq**2) &
&            + 2.0_dp*xq/(1.0_dp - xq**2))
 end function log_slope

 real(dp) function fdiff(xq, kk)

!! the 8th-order centred difference of the kk-th derivative at xq; for
!! xq < 4 h, of its even or odd extension to x < 0 (g is even)

 real(dp), intent(in) :: xq
 integer, intent(in) :: kk
 real(dp), parameter :: hh = 1.0e-3_dp
 real(dp), parameter :: cw(4) = [4.0_dp/5.0_dp, -1.0_dp/5.0_dp, &
&                                4.0_dp/105.0_dp, -1.0_dp/280.0_dp]
 real(dp) :: gp(0:4), gm(0:4), xm
 integer :: mm

 fdiff = 0.0_dp
 do mm = 1, 4
   call teter_derivs(xq + mm*hh, kk, gp(0:kk))
   xm = xq - mm*hh
   call teter_derivs(abs(xm), kk, gm(0:kk))
   if(xm < 0.0_dp .and. mod(kk, 2) == 1) gm(kk) = -gm(kk)
   fdiff = fdiff + cw(mm)*(gp(kk) - gm(kk))
 end do
 fdiff = fdiff/hh
 end function fdiff

 real(dp) function lnratio(rq, fc, shell)

!! ln(rhoc) - ln(fc rhotps) of the analytic densities

 real(dp), intent(in) :: rq, fc, shell

 lnratio = log(40.0_dp*exp(-4.0_dp*rq) &
&              + shell*exp(-(rq - 2.0_dp)**2/0.08_dp)) &
&          - log(fc*2.0_dp*rq**2*exp(-rq))
 end function lnratio

 real(dp) function outer_root(fc, shell)

!! the outermost root of lnratio below 10 bohr, by bisection to the last bit

 real(dp), intent(in) :: fc, shell
 real(dp) :: ra, rb, rm
 integer :: it

 ra = 0.0_dp
 do it = 10000, 1, -1
   if(lnratio(1.0e-3_dp*it, fc, shell) > 0.0_dp) then
     ra = 1.0e-3_dp*it
     exit
   end if
 end do
 rb = ra + 1.0e-3_dp
 do it = 1, 200
   rm = 0.5_dp*(ra + rb)
   if(rm <= ra .or. rm >= rb) exit
   if(lnratio(rm, fc, shell) > 0.0_dp) then
     ra = rm
   else
     rb = rm
   end if
 end do
 outer_root = 0.5_dp*(ra + rb)
 end function outer_root

 function lorentz(rq) result(ff)

!! 3/(1+r**4) and its first 4 derivatives, from (1+r**4) f = 3

 real(dp), intent(in) :: rq
 real(dp) :: ff(5), dd(0:4)
 integer :: kk, jj

! the derivatives of 1 + r**4
 dd = [1.0_dp + rq**4, 4.0_dp*rq**3, 12.0_dp*rq**2, 24.0_dp*rq, 24.0_dp]
 ff(1) = 3.0_dp/dd(0)
 do kk = 1, 4
   ff(kk + 1) = 0.0_dp
   do jj = 1, kk
     ff(kk + 1) = ff(kk + 1) - binom(kk, jj)*dd(jj)*ff(kk - jj + 1)
   end do
   ff(kk + 1) = ff(kk + 1)/dd(0)
 end do
 end function lorentz

 real(dp) function binom(nn, kk)

 integer, intent(in) :: nn, kk
 integer :: jj

 binom = 1.0_dp
 do jj = 1, kk
   binom = binom*(nn - jj + 1)/jj
 end do
 end function binom

 real(dp) function pval(a0q, acoq, rq, kk)

!! the kk-th derivative of a0 + sum_j aco(j) r**(j+2)

 real(dp), intent(in) :: a0q, acoq(5), rq
 integer, intent(in) :: kk
 real(dp) :: fac
 integer :: jj, mm

 pval = 0.0_dp
 if(kk == 0) pval = a0q
 do jj = 1, 5
   fac = 1.0_dp
   do mm = 0, kk - 1
     fac = fac*(jj + 2 - mm)
   end do
   pval = pval + acoq(jj)*fac*rq**(jj + 2 - kk)
 end do
 end function pval

 subroutine coefficients(rq, fm, a0q, acoq)

!! the coefficients that match fm at rq for the constant a0q (dgesv)

 real(dp), intent(in) :: rq, fm(5), a0q
 real(dp), intent(out) :: acoq(5)
 real(dp) :: amat(5, 5), unit(5)
 integer :: ipvt(5), info, jj, kk
 external :: dgesv

 do jj = 1, 5
   unit(:) = 0.0_dp
   unit(jj) = 1.0_dp
   do kk = 0, 4
     amat(kk + 1, jj) = pval(0.0_dp, unit, rq, kk)
   end do
 end do
 acoq(:) = fm(:)
 acoq(1) = acoq(1) - a0q
 call dgesv(5, 1, amat, 5, ipvt, acoq, 5, info)
 end subroutine coefficients

 real(dp) function max_slope(rq, a0q, acoq)

!! the largest p'/r**2 on [0, rq], on 200001 points (3 aco(1) at 0): p'
!! vanishes as r**2 at 0, where the polynomial may be barely monotone

 real(dp), intent(in) :: rq, a0q, acoq(5)
 real(dp) :: rj
 integer :: jj

 max_slope = 3.0_dp*acoq(1)
 do jj = 1, 200000
   rj = rq*jj/200000.0_dp
   max_slope = max(max_slope, pval(a0q, acoq, rj, 1)/rj**2)
 end do
 end function max_slope

 real(dp) function bisect_a0(rq, fm)

!! the least a0 whose polynomial has p' <= 0 on the points of max_slope

 real(dp), intent(in) :: rq, fm(5)
 real(dp) :: alo, ahi, am, acoq(5)
 integer :: it

 alo = fm(1)
 ahi = fm(1)
 do it = 1, 200
   call coefficients(rq, fm, ahi, acoq)
   if(max_slope(rq, ahi, acoq) <= 0.0_dp) exit
   ahi = 2.0_dp*ahi
 end do
 do it = 1, 200
   call coefficients(rq, fm, alo, acoq)
   if(max_slope(rq, alo, acoq) > 0.0_dp) exit
   alo = 0.5_dp*alo
 end do
 do it = 1, 80
   am = 0.5_dp*(alo + ahi)
   if(am <= alo .or. am >= ahi) exit
   call coefficients(rq, fm, am, acoq)
   if(max_slope(rq, am, acoq) > 0.0_dp) then
     alo = am
   else
     ahi = am
   end if
 end do
 bisect_a0 = ahi
 end function bisect_a0

 subroutine fit_residual(tt, yq, ndeg, res)

!! residual of the least-squares polynomial of degree ndeg, by dgels

 real(dp), intent(in) :: tt(41), yq(41)
 integer, intent(in) :: ndeg
 real(dp), intent(out) :: res(41)
 real(dp) :: aa(41, ndeg + 1), bb(41), work(2048)
 integer :: jj, kk, info
 external :: dgels

 do jj = 1, 41
   do kk = 0, ndeg
     aa(jj, kk + 1) = tt(jj)**kk
   end do
 end do
 bb(:) = yq(:)
 call dgels('N', 41, ndeg + 1, 1, aa, 41, bb, 41, work, size(work), info)
 do jj = 1, 41
   res(jj) = yq(jj)
   do kk = 0, ndeg
     res(jj) = res(jj) - bb(kk + 1)*tt(jj)**kk
   end do
 end do
 end subroutine fit_residual

 subroutine check(label, errq, tol)

!! errq against tol, printed and counted

 character(len=*), intent(in) :: label
 real(dp), intent(in) :: errq, tol

 write(stderr, '(a,a,es9.2)') label, ': error ', errq
 if(.not. (errq <= tol)) then
   write(stderr, '(a,a,es10.3,a,es9.2)') 'FAILED: ', label, errq, ' > ', tol
   nfail = nfail + 1
 end if
 end subroutine check

 end program test_modcore

! Unit tests of the polynomial local potential at a continuous radius
! (vploc_m; R-7 of doc/smoothness_plan.md), written before the code, on
! the program's default mesh (amesh=1.006, rr1=5e-5, rmax=45), for
! V(r) = -2/r + exp(-r) (its derivatives of every order O(1) at rc):
!   derivatives  V, V', V'', V''' at rc = 1.7 (lm_derivs) against the exact
!                ones to 1e-13, 1e-12, 1e-10, 1e-8
!   matching     the polynomial of each lpopt, differentiated analytically,
!                reproduces the matched V, V', V'' (and V''' for lpopt 3-5)
!                at rc to 1e-12 relative, with dvloc0 = 0; with dvloc0 =
!                0.1 the value at r = 0 rises by 0.1 and the match at rc
!                stays
!   smooth       the coefficients of lpopt 5 over 41 values of rc between
!                1.7 and 1.7 amesh (one mesh step) fit a polynomial of
!                degree 6 to 1e-11
! Stops with a non-zero code if any check failed.
 program test_vploc
 use constants_m, only: dp
 use unit_check_m, only: check, check_done
 use vploc_m, only: vploc_t, vploc_fit, vploc_eval
 use lapack_m, only: dgels
 implicit none

!Local variables
 real(dp), parameter :: rr1 = 5.0e-5_dp, rmax = 45.0_dp, amesh = 1.006_dp
 real(dp), parameter :: tols(0:3) = [1.0e-13_dp, 1.0e-12_dp, 1.0e-10_dp, &
&                                     1.0e-8_dp]
 real(dp), allocatable :: rr(:), vv(:)
 type(vploc_t) :: cf
 real(dp) :: al, rc, ts(41), cs(41, 4), resid(41), scale
 integer :: mmax, ii, kk, lpopt

 al = log(amesh)
 mmax = int(log(rmax/rr1)/al)
 allocate(rr(mmax), vv(mmax))
 do ii = 1, mmax
   rr(ii) = rr1*exp(al*(ii - 1))
 end do
 vv(:) = -2.0_dp/rr + exp(-rr)

! derivatives
 rc = 1.7_dp
 call vploc_fit(rr, vv, mmax, rc, 0.0_dp, 5, cf)
 call check('derivative 0', abs(cf%vder(0) - vex(rc, 0)), tols(0))
 call check('derivative 1', abs(cf%vder(1) - vex(rc, 1)), tols(1))
 call check('derivative 2', abs(cf%vder(2) - vex(rc, 2)), tols(2))
 call check('derivative 3', abs(cf%vder(3) - vex(rc, 3)), tols(3))

! matching: the polynomial's derivatives at rc, analytically
 do lpopt = 1, 5
   call vploc_fit(rr, vv, mmax, rc, 0.0_dp, lpopt, cf)
   do kk = 0, merge(3, 2, lpopt >= 3)
     scale = max(abs(cf%vder(kk)), 1.0_dp)
     call check('match derivative', abs(pder(cf, rc, kk) - cf%vder(kk))/scale, &
&               1.0e-12_dp)
   end do
   call vploc_fit(rr, vv, mmax, rc, 0.1_dp, lpopt, cf)
   call check('dvloc0 at 0', &
&             abs(vploc_eval(cf, 0.0_dp) - vpoly0(cf) - 0.1_dp), &
&             1.0e-12_dp)
   call check('dvloc0 at rc', abs(vploc_eval(cf, rc) - cf%vder(0)), 1.0e-12_dp)
 end do

! smooth
 do ii = 1, 41
   ts(ii) = 2.0_dp*(ii - 1)/40.0_dp - 1.0_dp
   rc = 1.7_dp*amesh**(0.5_dp*(ts(ii) + 1.0_dp))
   call vploc_fit(rr, vv, mmax, rc, 0.0_dp, 5, cf)
   cs(ii, :) = [cf%aco, cf%bco, cf%cco, cf%dco]
 end do
 do kk = 1, 4
   call fit_residual(ts, cs(:, kk), 6, resid)
   call check('smooth coefficient', maxval(abs(resid))/maxval(abs(cs(:, kk))), &
&             1.0e-11_dp)
 end do

 call check_done('vploc')

 contains

 real(dp) function vex(rq, kk)

!! the k-th derivative of -2/r + exp(-r)

 real(dp), intent(in) :: rq
 integer, intent(in) :: kk
 real(dp) :: fact
 integer :: jj

 fact = 1.0_dp
 do jj = 1, kk
   fact = fact*jj
 end do
 vex = -2.0_dp*(-1.0_dp)**kk*fact/rq**(kk + 1) + (-1.0_dp)**kk*exp(-rq)
 end function vex

 real(dp) function pder(cfq, rq, kk)

!! the k-th derivative at rq of the matched polynomial alone,
!! aco + bco r**p1 + cco r**p2 + dco r**p3 with the powers of lpopt

 type(vploc_t), intent(in) :: cfq
 real(dp), intent(in) :: rq
 integer, intent(in) :: kk
 integer :: pw(3), jj, mm
 real(dp) :: cc(3), fac

 select case(cfq%lpopt)
 case(1)
   pw = [2, 4, 0]
 case(2)
   pw = [4, 6, 0]
 case(3)
   pw = [4, 5, 6]
 case(4)
   pw = [4, 6, 8]
 case default
   pw = [2, 4, 6]
 end select
 cc = [cfq%bco, cfq%cco, cfq%dco]
 pder = 0.0_dp
 if(kk == 0) pder = cfq%aco
 do jj = 1, 3
   if(pw(jj) < kk) cycle
   fac = 1.0_dp
   do mm = 0, kk - 1
     fac = fac*(pw(jj) - mm)
   end do
   pder = pder + cc(jj)*fac*rq**(pw(jj) - kk)
 end do
 end function pder

 real(dp) function vpoly0(cfq)

!! the matched polynomial alone at r=0 (its constant term)

 type(vploc_t), intent(in) :: cfq

 vpoly0 = cfq%aco
 end function vpoly0

 subroutine fit_residual(tt, ys, ndeg, res)

!! residual of the least-squares polynomial of degree ndeg, by dgels

 real(dp), intent(in) :: tt(41), ys(41)
 integer, intent(in) :: ndeg
 real(dp), intent(out) :: res(41)
 real(dp) :: aa(41, ndeg + 1), bb(41), work(2048)
 integer :: jj, kk, info

 do jj = 1, 41
   do kk = 0, ndeg
     aa(jj, kk + 1) = tt(jj)**kk
   end do
 end do
 bb(:) = ys(:)
 call dgels('N', 41, ndeg + 1, 1, aa, 41, bb, 41, work, size(work), info)
 do jj = 1, 41
   res(jj) = ys(jj)
   do kk = 0, ndeg
     res(jj) = res(jj) - bb(kk + 1)*tt(jj)**kk
   end do
 end do
 end subroutine fit_residual

 end program test_vploc

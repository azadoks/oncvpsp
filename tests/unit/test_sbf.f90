! Unit tests of the double-precision spherical Bessel functions (sbf_m,
! sbf_values; Q-10 of doc/smoothness_plan.md), written before the code. Errors
! are against the quad-precision sbfq_values, relative to the envelope
! sqrt(j_l**2 + j_(l+1)**2) (no zero of j_l is a zero of j_(l+1), so it
! is |j_l| at small x and the amplitude 1/x at large x; at x=0, where
! it is 0 for l>0, the error is absolute), on x=0 and 2001 points
! log-spaced in [1e-8, 50]:
!   trig        sbf_trig, j_l for l=0..3 from x, sin x and cos x, within
!               1e-15 (5 eps): the closed forms lose about 3 eps/x**2,
!               45 eps/x**4 and 1575 eps/x**6 for l=1, 2, 3 at small x
!   trig_abs    sbf_trig given sin x and cos x with an absolute error of
!               1e-14 (eresid rotates them from one q to the next) within
!               1e-13: below its series threshold the result must not
!               depend on them
!   rising      sbf_trig_rising, r*j_l(q*r) on the same points as radii
!               with q=1, within 1e-15
!   sbf_values        every order sbf_values returns, for nm=1..8, within
!               eps*max(10, x), and none left unset (the array is filled
!               with a NaN first): the downward recurrence of nm>=5 gains
!               about 0.6 eps of round-off at each of the x oscillatory
!               orders it passes; with a start at nm+1 below x=1/15 its
!               top orders were off by 1e-9
! Stops with a non-zero code if any check failed.
 program test_sbf
 use constants_m, only: dp, qp
 use unit_check_m, only: check, check_done
 use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
&  ieee_is_finite
 use sbf_m, only: sbf_trig, sbf_trig_rising, sbf_values
 use sbf_q_m, only: sbfq_values
 implicit none

!Local variables
 integer, parameter :: npt = 2002, lmax = 7
 real(dp) :: xs(npt), env(0:lmax, npt), jq(0:lmax, npt)
 real(dp) :: sb(lmax + 1), err, jl, rjl(npt)
 real(qp) :: sbq(lmax + 2)
 integer :: ii, ll, nm, nunset

 xs(1) = 0.0_dp
 do ii = 2, npt
   xs(ii) = 1.0e-8_dp*(50.0_dp/1.0e-8_dp)**(dble(ii - 2)/dble(npt - 2))
 end do
 do ii = 1, npt
   call sbfq_values(lmax + 2, real(xs(ii), qp), sbq)
   do ll = 0, lmax
     jq(ll, ii) = real(sbq(ll + 1), dp)
     env(ll, ii) = max(real(sqrt(sbq(ll + 1)**2 + sbq(ll + 2)**2), dp), &
&                      tiny(1.0_dp))
   end do
 end do

! trig
 err = 0.0_dp
 do ii = 1, npt
   do ll = 0, 3
     jl = sbf_trig(ll, xs(ii), sin(xs(ii)), cos(xs(ii)))
     err = max(err, abs(jl - jq(ll, ii))/env(ll, ii))
   end do
 end do
 call check('trig', err, 1.0e-15_dp)

! trig_abs
 err = 0.0_dp
 do ii = 1, npt
   do ll = 0, 3
     jl = sbf_trig(ll, xs(ii), sin(xs(ii)) + 1.0e-14_dp, &
&                  cos(xs(ii)) - 1.0e-14_dp)
     err = max(err, abs(jl - jq(ll, ii))/env(ll, ii))
   end do
 end do
 call check('trig_abs', err, 1.0e-13_dp)

! rising
 err = 0.0_dp
 do ll = 0, 3
   call sbf_trig_rising(ll, npt, 1.0_dp, xs, sin(xs), cos(xs), rjl)
   do ii = 2, npt
     err = max(err, abs(rjl(ii)/xs(ii) - jq(ll, ii))/env(ll, ii))
   end do
   err = max(err, abs(rjl(1)))
 end do
 call check('rising', err, 1.0e-15_dp)

! sbf_values
 err = 0.0_dp
 nunset = 0
 do nm = 1, lmax + 1
   do ii = 1, npt
     sb(:) = ieee_value(1.0_dp, ieee_quiet_nan)
     call sbf_values(nm, xs(ii), sb)
     do ll = 0, nm - 1
       if(.not. ieee_is_finite(sb(ll + 1))) then
         nunset = nunset + 1
       else
         err = max(err, abs(sb(ll + 1) - jq(ll, ii))/env(ll, ii) &
&                       /max(10.0_dp, xs(ii)))
       end if
     end do
   end do
 end do
 call check('sbf_values', err, epsilon(1.0_dp))
 call check('sbf_values orders unset', dble(nunset), 0.0_dp)

 call check_done('sbf')

 end program test_sbf

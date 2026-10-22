! Unit tests of the free tail of the well states (free_tail_m; phase 8 of
! doc/smoothness_plan.md), for a free particle (V = 0, l = 0) from rc = 2:
!   solution  the tail of sin(k r) at ep = k**2/2 = 0.5 is sin(k r) times
!             the window to 1e-10 of its largest value (Runge-Kutta)
!   kinetic   tail_kinetic, t t'/2 at rc + int t'**2/2, equals int
!             t (-t''/2) of the tail of sin(k r) to 1e-9
!   leakage   the window's share of E_resid above qcut = 5: the transforms
!             of the tails with L and 2 L, F(q) = int t(r) sin(q r)/q dr,
!             differ above qcut by (1/pi) int q**4 |dF|**2 dq below 1e-13
!             of the tail's own (1/pi) int q**4 |F|**2 dq there (from its
!             start at rc, the same in both, so dF is the window's alone)
!   growing   the same for sinh(kappa r) at ep = -kappa**2/2 = -0.125,
!             with L and 1.1 L (at 2 L the tail reaches exp(34), and the
!             transforms' round-off would dominate)
!   width     L has no kink at ep = 0: its second difference there over
!             h**2 is the same for h = 1e-4 and 2e-4 to 1e-2
! Stops with a non-zero code if any check failed.
 program test_free_tail
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
 use free_tail_m, only: free_tail, free_tail_width, tail_kinetic
 use ugrid_m, only: ug_weights
 implicit none

!Local variables
 real(dp), parameter :: pi = 3.14159265358979323846_dp
 real(dp), parameter :: rc = 2.0_dp, qcut = 5.0_dp, qmax = 20.0_dp
 integer, parameter :: mmax = 2000
 real(dp), allocatable :: ulin(:), dlin(:), u2(:), d2(:), wr(:)
 real(dp) :: rr(mmax), vv(mmax), uu(mmax), up(mmax)
 real(dp) :: al, hh, kk, kp, ee, wl, r0, err, rad, xx
 real(dp) :: dm, dp_, tk, tk2
 integer :: nfail = 0, nr, n2, ii

 al = log(1.006_dp)
 do ii = 1, mmax
   rr(ii) = 5.0e-4_dp*exp(al*dble(ii - 1))
 end do
 vv(:) = 0.0_dp
 hh = 2.0_dp*pi/(75.0_dp*qmax)

! solution
 kk = 1.0_dp
 ee = 0.5_dp*kk**2
 uu(:) = sin(kk*rr)
 up(:) = al*rr*kk*cos(kk*rr)
 call free_tail(0, 0, -1, ee, rr, vv, uu, up, mmax, rc, qcut, hh, nr, &
&               ulin, dlin)
 call free_tail_width(qcut, ee, wl)
 r0 = rc + 6.0_dp*wl
 err = 0.0_dp
 do ii = 1, nr
   rad = rc + dble(ii - 1)*hh
   err = max(err, abs(ulin(ii) - sin(kk*rad)*0.5_dp*erfc((rad - r0)/wl)))
 end do
 call check('solution sin(k r)', err, 1.0e-10_dp)

! kinetic
 allocate(wr(nr))
 call ug_weights(nr, hh, wr)
 tk = tail_kinetic(0, rc, hh, nr, ulin, dlin)
 tk2 = 0.0_dp
 do ii = 1, nr
   rad = rc + dble(ii - 1)*hh
   xx = (rad - r0)/wl
! t'' of sin(k r) w(r), with w' = -exp(-x**2)/(L sqrt(pi)), w'' = -2 x w'/L
   tk2 = tk2 + wr(ii)*ulin(ii)*(-0.5_dp)*( &
&        -kk**2*sin(kk*rad)*0.5_dp*erfc(xx) &
&        - 2.0_dp*kk*cos(kk*rad)*exp(-xx**2)/(wl*sqrt(pi)) &
&        + sin(kk*rad)*2.0_dp*xx*exp(-xx**2)/(wl**2*sqrt(pi)))
 end do
 deallocate(wr)
 write(stderr, '(a,2es22.14)') 'kinetic: tail_kinetic, int t(-t''''/2)', &
&  tk, tk2
 call check('kinetic', abs(tk - tk2)/abs(tk2), 1.0e-9_dp)

! leakage
 call free_tail(0, 0, -1, ee, rr, vv, uu, up, mmax, rc, qcut, hh, n2, &
&               u2, d2, 2.0_dp)
 call check('leakage sin(k r)', leak(), 1.0e-13_dp)

! growing
 kp = 0.5_dp
 ee = -0.5_dp*kp**2
 uu(:) = sinh(kp*rr)
 up(:) = al*rr*kp*cosh(kp*rr)
 call free_tail(0, 0, -1, ee, rr, vv, uu, up, mmax, rc, qcut, hh, nr, &
&               ulin, dlin)
 call free_tail(0, 0, -1, ee, rr, vv, uu, up, mmax, rc, qcut, hh, n2, &
&               u2, d2, 1.1_dp)
 call check('leakage sinh(kappa r)', leak(), 1.0e-13_dp)

! width
 dm = second(1.0e-4_dp)
 dp_ = second(2.0e-4_dp)
 write(stderr, '(a,2es14.6)') 'width: d2L/dep2 at 0 with h, 2 h', dm, dp_
 call check('width smooth at ep = 0', abs(dm - dp_)/abs(dm), 1.0e-2_dp)

 if(nfail > 0) then
   write(stderr, '(i0,a)') nfail, ' check(s) failed'
   error stop 1
 end if
 write(stderr, '(a)') 'all free_tail checks passed'

 contains

 real(dp) function second(hd)

!! the second difference of L at ep = 0 over hd**2

 real(dp), intent(in) :: hd
 real(dp) :: lm, l0, lp

 call free_tail_width(qcut, -hd, lm)
 call free_tail_width(qcut, 0.0_dp, l0)
 call free_tail_width(qcut, hd, lp)
 second = (lm - 2.0_dp*l0 + lp)/hd**2
 end function second

 real(dp) function leak()

!! (1/pi) int_qcut^qmax q**4 |F1 - F2|**2 dq over (1/pi) int_qcut^qmax
!! q**4 |F1|**2 dq, F1 and F2 the transforms of ulin (nr nodes) and u2
!! (n2 nodes)

 integer, parameter :: nq = 4001
 real(dp) :: wq(nq), wa(nr), wb(n2), qq, dq, fa, fb, above, total
 integer :: iq, jj

 call ug_weights(nr, hh, wa)
 call ug_weights(n2, hh, wb)
 dq = qmax/dble(nq - 1)
 call ug_weights(nq, dq, wq)
 above = 0.0_dp
 total = 0.0_dp
 do iq = 2, nq
   qq = dble(iq - 1)*dq
   fa = 0.0_dp
   do jj = 1, nr
     fa = fa + wa(jj)*ulin(jj)*sin(qq*(rc + dble(jj - 1)*hh))
   end do
   fb = 0.0_dp
   do jj = 1, n2
     fb = fb + wb(jj)*u2(jj)*sin(qq*(rc + dble(jj - 1)*hh))
   end do
   fa = fa/qq
   fb = fb/qq
   if(qq >= qcut) then
     total = total + wq(iq)*qq**4*fa**2
     above = above + wq(iq)*qq**4*(fa - fb)**2
   end if
 end do
 write(stderr, '(a,2es12.4)') 'leakage: window''s, tail''s above qcut', &
&  above, total
 leak = above/total
 end function leak

 subroutine check(name, val, tol)
 character(len=*), intent(in) :: name
 real(dp), intent(in) :: val, tol
 if(.not. (val <= tol)) then
   write(stderr, '(3a,es12.4,a,es12.4)') 'FAILED: ', name, ': ', val, &
&    ' > ', tol
   nfail = nfail + 1
 else
   write(stderr, '(3a,es12.4)') 'ok: ', name, ': ', val
 end if
 end subroutine check

 end program test_free_tail

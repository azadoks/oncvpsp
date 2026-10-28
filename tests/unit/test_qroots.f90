! Unit tests of the basis wave vectors of run_optimize (qroots; Q-5 of
! doc/smoothness_plan.md): the roots of g(q) = ulgd*j_l(q rc) - q j_l'(q rc)
! that qroots returns for l=0..3, rc=1.7, 8 roots (6 matching ones, at
! qroot(2) and qroot(4:8)), scanned in steps of 0.02 as run_optimize does:
!   roots     at ulgd = -5, -1, 0.5*l/rc - 0.3, 20 and 65 (a Pb well state):
!             increasing, and each matching root within 1e-13 relative
!             of the root of g, by a Newton step |g/g'|/q (g itself can be
!             small at the root: ulgd near 0 puts it near a zero of j_l')
!   pole      ulgd = -1e12 and 1e12 (a node of the wave function at rc):
!             the same roots within 1e-9 relative, as g/ulgd is smooth in
!             1/ulgd
!   lcross    ulgd scanned through l/rc, where the lowest matching root
!             of the first branch of j_l (below its first zero) goes to
!             q=0 and leaves: at ulgd = l/rc + s, s from -0.05 to 0.05 in
!             steps of 1e-4, no root may change between neighbouring
!             points by more than 0.02 (the scan step); a root that leaves
!             shifts every index by one, a step of the size of the
!             spacing of the roots (~pi/rc). This is the jump of the
!             basis where the wave function's R=u/r has the log
!             derivative of r**l at rc (for l=0 the peak of R: rc=2 in
!             69_Tm_ae and 70_Yb_ae is beyond it, at ulgd - l/rc = 1.1 and
!             0.97; the other inputs are at least 0.65 below), and it
!             fails today: run alone with the argument lcross, as the
!             test unit_qroots_lcross that is expected to fail
!             (WILL_FAIL) until qroots follows the root through q=0
! Stops with a non-zero code if any check failed.
 program test_qroots
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: stderr => error_unit
 use unit_check_m, only: check, check_done
 use qroots_m, only: qroots
 use sbf_rc_der_m, only: sbf_rc_der
 implicit none
 integer, parameter :: iout=6

!Local variables
 integer :: ifail
 integer, parameter :: nroot = 8, nscan = 1001
 integer, parameter :: imatch(6) = [2, 4, 5, 6, 7, 8]
 real(dp), parameter :: rc = 1.7_dp, dltq = 0.02_dp, qmax = 200.0_dp
 real(dp) :: qroot(nroot), qa(nroot), qb(nroot), qlast(nroot), ulgds(5)
 real(dp) :: sbfd(5), err, ss, ulgd, jump, sjump, gg, dg, qq
 integer :: ll, ii, kk
 character(len=16) :: arg
 logical :: lcross

 arg = ''
 if(command_argument_count() > 0) call get_command_argument(1, arg)
 lcross = trim(arg) == 'lcross'

 if(.not. lcross) then
! roots
   err = 0.0_dp
   do ll = 0, 3
     ulgds = [-5.0_dp, -1.0_dp, 0.5_dp*ll/rc - 0.3_dp, 20.0_dp, 65.0_dp]
     do kk = 1, size(ulgds)
       call qroots(ll, rc, ulgds(kk), nroot, dltq, qmax, qroot, iout, ifail)
       if(ifail /= 0) error stop 1
       do ii = 2, nroot
         if(qroot(ii) <= qroot(ii - 1)) err = huge(1.0_dp)
       end do
       do ii = 1, size(imatch)
         qq = qroot(imatch(ii))
         call sbf_rc_der(ll, qq, rc, sbfd)
         gg = ulgds(kk)*sbfd(1) - sbfd(2)
         dg = (ulgds(kk)*rc*sbfd(2) - sbfd(2) - rc*sbfd(3))/qq
         err = max(err, abs(gg/dg)/qq)
       end do
     end do
   end do
   call check('roots', err, 1.0e-13_dp)

! pole
   err = 0.0_dp
   do ll = 0, 3
     call qroots(ll, rc, -1.0e12_dp, nroot, dltq, qmax, qa, iout, ifail)
     if(ifail /= 0) error stop 1
     call qroots(ll, rc, 1.0e12_dp, nroot, dltq, qmax, qb, iout, ifail)
     if(ifail /= 0) error stop 1
     err = max(err, maxval(abs(qa - qb)/qb))
   end do
   call check('pole', err, 1.0e-9_dp)
 end if

! lcross
 if(lcross) then
   do ll = 0, 3
     jump = 0.0_dp
     sjump = 0.0_dp
     do ii = 1, nscan
       ss = -0.05_dp + 1.0e-4_dp*(ii - 1)
       ulgd = dble(ll)/rc + ss
       call qroots(ll, rc, ulgd, nroot, dltq, qmax, qroot, iout, ifail)
       if(ifail /= 0) error stop 1
       if(ii > 1) then
         if(maxval(abs(qroot - qlast)) > jump) then
           jump = maxval(abs(qroot - qlast))
           sjump = ss
         end if
       end if
       qlast = qroot
     end do
     write(stderr, '(a,i0,a,es10.3,a,es10.3)') '  lcross l=', ll, &
&      ': largest step ', jump, ' at s=', sjump
     call check('lcross l='//achar(48 + ll), jump, dltq)
   end do
 end if

 call check_done('qroots')

 end program test_qroots

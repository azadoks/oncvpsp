! Unit tests of the quad-precision Bessel overlaps and the symmetric
! eigen-solver (sbf_q_m, qlinalg_m; doc/smoothness_plan.md, R-5), written
! before the code:
!   values     j_0(x) = sin(x)/x and j_1(x) = sin(x)/x**2 - cos(x)/x from
!              sbfq_values at x = 0.3, 2.5, 17 to 1e-32; the small-argument
!              series at x = 1e-4 against sin(x)/x and the series of j_1
!   overlap    sbfq_overlap against sbf_overlap (double) at l = 0..3 and
!              several (q1, q2, rc), to 1e-14 relative; and the exact
!              int_0^rc r**2 j_0(q r)**2 dr = (rc/2 - sin(2 q rc)/(4 q))/q**2
!              to 1e-32 relative
!   rc_der     sbfq_rc_der against sbf_rc_der (double) to 1e-13 relative,
!              and its first derivative against a quad finite difference of
!              sbfq_values to 1e-20
!   jacobi     a 6x6 symmetric matrix with known spectrum (Q diag(w) Q^T
!              from a Householder Q): eigenvalues to 1e-32 absolute,
!              eigenvectors orthonormal and reproducing the matrix to 1e-32
!   basis      the overlap matrix of the 8 basis functions of 14_Si l=2
!              (wave vectors of qroots) has its smallest eigenvalue near
!              4.8e-9 and positive, and V diag(w) V^T reproduces it to 1e-32
! Stops with a non-zero code if any check failed.
 program test_quadbasis
 use constants_m, only: dp, qp
 use, intrinsic :: iso_fortran_env, only: stderr => error_unit
 use unit_check_m, only: check, check_done, nfail
 use sbf_q_m, only: sbfq_values, sbfq_overlap, sbfq_rc_der
 use qlinalg_m, only: qjacobi
 use sbf_overlap_m, only: sbf_overlap
 use sbf_rc_der_m, only: sbf_rc_der
 implicit none

!Local variables
 real(qp), parameter :: xs(3) = [0.3_qp, 2.5_qp, 17.0_qp]
 real(dp), parameter :: dqs(4) = [0.0_dp, 1.0e-7_dp, 0.2_dp, 3.0_dp]
 real(qp) :: sb(8), xq, sq, dq, der(5), hq, fp, fm, aq(6, 6), wq(6), &
&            vq(6, 6), uu(6), ss(8, 8), ws(8), vs(8, 8), err
 real(qp), parameter :: qroot(8) = [0.6479769_qp, 1.295954_qp, 2.003985_qp, &
&   2.712016_qp, 3.936362_qp, 5.129723_qp, 6.310968_qp, 7.486094_qp]
 real(dp) :: sd, derd(5), q1, q2, rc
 integer :: ii, jj, ll, kk

! values
 do ii = 1, 3
   xq = xs(ii)
   call sbfq_values(8, xq, sb)
   call check('values j_0', abs(sb(1) - sin(xq)/xq), 1.0e-32_qp)
   call check('values j_1', abs(sb(2) - (sin(xq)/xq**2 - cos(xq)/xq)), &
&             1.0e-32_qp)
 end do
 xq = 1.0e-4_qp
 call sbfq_values(8, xq, sb)
 call check('values j_0 small', abs(sb(1) - sin(xq)/xq), 1.0e-32_qp)
 call check('values j_1 small', abs(sb(2) - (xq/3.0_qp - xq**3/30.0_qp &
&           + xq**5/840.0_qp - xq**7/45360.0_qp)), 1.0e-32_qp)

! overlap
 do ll = 0, 3
   do ii = 1, 4
     rc = 2.7_dp
     q1 = 0.4_dp*ii + 0.3_dp*ll
     q2 = q1 + dqs(ii)
     call sbf_overlap(ll, q1, q2, rc, sd)
     call sbfq_overlap(ll, real(q1, qp), real(q2, qp), real(rc, qp), sq)
     call check('overlap vs double', real(abs(real(sq, dp) - sd)/abs(sd), qp), &
&               1.0e-14_qp)
   end do
 end do
 xq = 3.7_qp
 call sbfq_overlap(0, xq, xq, 2.7_qp, sq)
 err = (1.35_qp - sin(2.0_qp*xq*2.7_qp)/(4.0_qp*xq))/xq**2
 call check('overlap j_0 exact', abs(sq - err)/abs(err), 1.0e-32_qp)

! rc_der
 do ll = 0, 3
   q1 = 2.1_dp + 0.7_dp*ll
   call sbf_rc_der(ll, q1, 2.7_dp, derd)
   call sbfq_rc_der(ll, real(q1, qp), 2.7_qp, der)
   do kk = 1, 5
     call check('rc_der vs double', real(abs(real(der(kk), dp) - derd(kk)) &
&               /maxval(abs(derd)), qp), 1.0e-13_qp)
   end do
   hq = 1.0e-10_qp
   call sbfq_values(ll + 1, real(q1, qp)*(2.7_qp + hq), sb)
   fp = sb(ll + 1)
   call sbfq_values(ll + 1, real(q1, qp)*(2.7_qp - hq), sb)
   fm = sb(ll + 1)
   call check('rc_der first derivative', abs(der(2) - (fp - fm)/(2.0_qp*hq)), &
&             1.0e-19_qp)
 end do

! jacobi: Q = I - 2 u u^T/(u^T u), A = Q diag(w) Q^T
 do ii = 1, 6
   uu(ii) = real(ii, qp)**2 - 3.0_qp
   wq(ii) = real(ii, qp)*0.5_qp - 1.7_qp
 end do
 do ii = 1, 6
   do jj = 1, 6
     vq(ii, jj) = -2.0_qp*uu(ii)*uu(jj)/sum(uu**2)
   end do
   vq(ii, ii) = vq(ii, ii) + 1.0_qp
 end do
 do ii = 1, 6
   do jj = 1, 6
     aq(ii, jj) = sum(vq(ii, :)*wq(:)*vq(jj, :))
   end do
 end do
 call qjacobi(6, aq, ws(1:6), vs(1:6, 1:6))
 do ii = 1, 6
   call check('jacobi eigenvalue', abs(ws(ii) - wq(ii)), 1.0e-32_qp)
 end do
 err = 0.0_qp
 do ii = 1, 6
   do jj = 1, 6
     err = max(err, abs(sum(vs(1:6, ii)*vs(1:6, jj)) &
&              - merge(1.0_qp, 0.0_qp, ii == jj)))
     err = max(err, abs(sum(vs(ii, 1:6)*ws(1:6)*vs(jj, 1:6)) - aq(ii, jj)))
   end do
 end do
 call check('jacobi orthonormal and reproducing', err, 1.0e-32_qp)

! basis
 do ii = 1, 8
   do jj = 1, 8
     call sbfq_overlap(2, qroot(ii), qroot(jj), 2.7030231215135725_qp, &
&                       ss(ii, jj))
   end do
 end do
 call qjacobi(8, ss, ws, vs)
 write(stderr, '(a,2es12.3)') 'basis: smallest and largest eigenvalue', &
&  ws(1), ws(8)
 if(.not. (ws(1) > 1.0e-9_qp .and. ws(1) < 1.0e-8_qp)) then
   write(stderr, '(a)') 'FAILED: basis: smallest eigenvalue not near 4.8e-9'
   nfail = nfail + 1
 end if
 err = 0.0_qp
 do ii = 1, 8
   do jj = 1, 8
     err = max(err, abs(sum(vs(ii, :)*ws(:)*vs(jj, :)) - ss(ii, jj)))
   end do
 end do
 call check('basis reproduced', err, 1.0e-32_qp)

 call check_done('quadbasis')

 end program test_quadbasis

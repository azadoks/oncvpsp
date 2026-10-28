! Unit tests of the continuation of the phases (phsft_m, phsft_pair) across
! a resonance much narrower than the energy step: hydrogen-like, zz=3, l=1,
! non-relativistic, whose 2p level at -1.125 Ha is a pole of r phi'/phi
! about 1e-6 Ha wide at 4 and 5 bohr, scanned from -0.5 Ha down by 0.2 Ha.
! The all-electron phases must differ from the Pruefer angle
! (phsft_node_phase) by the same multiple of pi at every energy, to 1e-12
! (before the Pruefer angle's check, by pi more below the level); the
! pseudopotential's, the same potential as the local one with no
! projectors, which are continued by halving the steps, must equal them to
! 1e-6 (the two solvers differ by 4e-7).
! Stops with a non-zero code if any check failed.
 program test_phsft
 use constants_m, only: dp
 use unit_check_m, only: check, check_done
 use pieces_m, only: piecewise_t, pw_alloc
 use phsft_m, only: phsft_t, phsft_pair, phsft_ae_init, phsft_node_phase, &
&                   phsft_free
 use, intrinsic :: iso_fortran_env, only: output_unit
 use adams_m, only: adams_t
 implicit none

!Local variables
!> the Adams-Moulton steps, of the default order
 type(adams_t) :: adams
 integer, parameter :: mmax = 950, ll = 1, nq = 2, npsh = 7
 real(dp), parameter :: zz = 3.0_dp, rr1 = 0.0005_dp, amesh = 1.012_dp
 real(dp), parameter :: epsh2 = -0.5_dp, depsh = 0.2_dp, eal = -0.5_dp
 real(dp) :: rr(mmax), vv(mmax), rq(nq), evkb(1)
 real(dp) :: pshae(npsh, nq), pshps(npsh, nq), pr(npsh, nq)
 real(dp) :: err_ae, err_ps, ee
 type(piecewise_t) :: vlocpw, vkbpw(1)
 type(phsft_t) :: ps
 integer :: ii, iq, ifail

 do ii = 1, mmax
   rr(ii) = rr1*amesh**(ii - 1)
 end do
 vv(:) = -zz/rr(:)
 call pw_alloc(vlocpw, mmax, 1, [0.0_dp, rr(mmax)])
 vlocpw%vals(:, 1) = vv
 evkb(:) = 0.0_dp
 rq = [4.0_dp, 5.0_dp]

 call phsft_pair(0, 0.0_dp, ll, ll, 0, rr, vv, zz, vv, vlocpw, vkbpw, &
&                evkb, mmax, adams, nq, rq, npsh, epsh2, depsh, eal, pshae, &
&                pshps, output_unit, ifail)
 call check('phsft_pair: ifail', real(ifail, dp), 0.0_dp)

 call phsft_ae_init(0, 0.0_dp, ll, ll, rr, vv, zz, mmax, adams, nq, rq, ps, &
&                   output_unit, ifail)
 call check('phsft_ae_init: ifail', real(ifail, dp), 0.0_dp)
 do iq = 1, nq
   do ii = 1, npsh
     ee = epsh2 - depsh*real(ii - 1, dp)
     call phsft_node_phase(ps, iq, ee, pr(ii, iq), iout=output_unit, &
&                          ifail=ifail)
     call check('phsft_node_phase: ifail', real(ifail, dp), 0.0_dp)
   end do
 end do
 call phsft_free(ps)

 do iq = 1, nq
   err_ae = maxval(abs(pshae(:, iq) - pr(:, iq) - (pshae(1, iq) - pr(1, iq))))
   err_ps = maxval(abs(pshps(:, iq) - pshae(:, iq)))
   call check('all-electron against the Pruefer angle', err_ae, 1.0e-12_dp)
   call check('pseudopotential against all-electron', err_ps, 1.0e-6_dp)
 end do

 call check_done('test_phsft')

 end program test_phsft

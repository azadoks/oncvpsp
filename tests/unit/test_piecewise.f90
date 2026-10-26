! Unit tests of the operations on functions built from pieces (T6, phase 3
! of doc/smoothness_plan.md), on the program's default mesh (amesh=1.006,
! rr1=5e-5, rmax=45). The derivatives at rc are unit_ode_rc_der's:
!   kink        f = r**2 exp(-r) inside a radius k near 1.5 bohr, plus
!               (r-k)**2 r exp(-r) outside (C1 at k): its integral to 6
!               bohr, split at k and each piece integrated as its own smooth
!               function (lm_integral of the piece's extension), is exact to
!               1e-13 wherever k sits between mesh points; integrated as one
!               sampled array it is off by 5e-12 to 1e-9 depending on where
!               (printed, not checked)
!   wronskian   hydrogen 1s and 2s solved by boundstate: int_0^a u1 u2 dr
!               (lm_integral) against the Wronskian at a alone,
!               [u1 u2' - u1' u2](a) / (2 (E1 - E2)), with u, u' at a from
!               ode_rc_der, to 1e-9 (the solver's error in u)
! Stops with a non-zero code if any check failed.
 program test_piecewise
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
 use logmesh_m, only: lm_integral
 use ode_rc_der_m, only: ode_rc_der
 use boundstate_m, only: boundstate
 implicit none

!Local variables
 integer :: ifail
 real(dp), parameter :: rr1 = 5.0e-5_dp, rmax = 45.0_dp, amesh = 1.006_dp
 real(dp), parameter :: rb = 6.0_dp
 real(dp), parameter :: offsets(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, &
&                                      0.999_dp]
 real(dp), parameter :: radii(3) = [0.7_dp, 1.2345_dp, 2.718281828459045_dp]
 real(dp), allocatable :: rr(:), vv(:), f1(:), f2(:), fk(:), uu(:, :, :), &
&                         up(:, :, :)
 real(dp) :: al, rk, s1, s2, sk, exact, err, ee(2), ord1(5), ord2(5), wr, ss
 integer :: nfail = 0, mmax, ii, io, ierr, mch, ir

 al = log(amesh)
 mmax = int(log(rmax/rr1)/al)
 allocate(rr(mmax), vv(mmax), f1(mmax), f2(mmax), fk(mmax), &
&         uu(mmax, 2, 2), up(mmax, 2, 2))
 do ii = 1, mmax
   rr(ii) = rr1*exp(al*(ii - 1))
 end do

! kink: the inside piece and the outside piece's extension, both smooth
 f1(:) = rr**2*exp(-rr)
 do io = 1, size(offsets)
   rk = rr(1 + int(log(1.5_dp/rr1)/al))*amesh**offsets(io)
   f2(:) = f1(:) + (rr - rk)**2*rr*exp(-rr)
   fk(:) = merge(f1, f2, rr <= rk)
   call lm_integral(rr, f1, mmax, 2, 0.0_dp, rk, s1)
   call lm_integral(rr, f2, mmax, 1, rk, rb, s2)
   call lm_integral(rr, fk, mmax, 2, 0.0_dp, rb, sk)
   exact = prim1(rk) + prim2(rb, rk) - prim2(rk, rk)
   err = abs(s1 + s2 - exact)/exact
   write(stderr, '(a,f6.3,a,es9.2,a,es9.2)') 'kink at offset', offsets(io), &
&    ': split', err, '  one array', abs(sk - exact)/exact
   if(.not. (err <= 1.0e-13_dp)) then
     write(stderr, '(a,es10.3)') 'FAILED: kink: split integral off by ', err
     nfail = nfail + 1
   end if
 end do

! wronskian: hydrogen 1s and 2s, nr
 vv(:) = -1.0_dp/rr
 uu(:, :, :) = 0.0_dp
 up(:, :, :) = 0.0_dp
 do ii = 1, 2
   ee(ii) = -0.5_dp/ii**2
   mch = 0
   call boundstate(0, ii, 0, -1, ierr, ee(ii), rr, 1.0_dp, vv, uu(:, :, ii), &
&                  up(:, :, ii), mmax, mch, ifail)
   if(ifail /= 0) error stop 1
   if(ierr /= 0) then
     write(stderr, '(a,i0,a,i0)') 'FAILED: wronskian: boundstate n ', ii, &
&      ' ierr ', ierr
     nfail = nfail + 1
   end if
 end do
 f1(:) = uu(:, 1, 1)*uu(:, 1, 2)
 do ir = 1, size(radii)
   call ode_rc_der(0, 0, -1, ee(1), rr, vv, uu(:, 1, 1), up(:, 1, 1), mmax, &
&                  radii(ir), ord1)
   call ode_rc_der(0, 0, -1, ee(2), rr, vv, uu(:, 1, 2), up(:, 1, 2), mmax, &
&                  radii(ir), ord2)
! u = r (u/r), u' = u/r + r (u/r)'
   wr = radii(ir)*ord1(1)*(ord2(1) + radii(ir)*ord2(2)) &
&     - (ord1(1) + radii(ir)*ord1(2))*radii(ir)*ord2(1)
   wr = wr/(2.0_dp*(ee(1) - ee(2)))
   call lm_integral(rr, f1, mmax, 2, 0.0_dp, radii(ir), ss)
   err = abs(ss - wr)/maxval(abs(f1))
   write(stderr, '(a,f7.4,a,es21.13,a,es9.2)') 'wronskian at', radii(ir), &
&    ': integral', ss, '  error', err
   if(.not. (err <= 1.0e-9_dp)) then
     write(stderr, '(a,es10.3)') 'FAILED: wronskian: off by ', err
     nfail = nfail + 1
   end if
 end do

 if(nfail > 0) then
   write(stderr, '(i0,a)') nfail, ' check(s) failed'
   error stop 1
 end if
 write(stderr, '(a)') 'piecewise: all checks passed'

 contains

 function prim1(rq)

!! int_0^rq r**2 exp(-r) dr

 real(dp), intent(in) :: rq
 real(dp) :: prim1

 prim1 = 2.0_dp - (rq**2 + 2.0_dp*rq + 2.0_dp)*exp(-rq)
 end function prim1

 function prim2(rq, rk)

!! an antiderivative of the outside piece r**2 exp(-r) + (r-rk)**2 r exp(-r),
!! its second term expanded in powers of r: r**3 - 2 rk r**2 + rk**2 r

 real(dp), intent(in) :: rq, rk
 real(dp) :: prim2

 prim2 = -(rq**2 + 2.0_dp*rq + 2.0_dp)*exp(-rq) &
&        - (rq**3 + 3.0_dp*rq**2 + 6.0_dp*rq + 6.0_dp)*exp(-rq) &
&        + 2.0_dp*rk*(rq**2 + 2.0_dp*rq + 2.0_dp)*exp(-rq) &
&        - rk**2*(rq + 1.0_dp)*exp(-rq)
 end function prim2

 end program test_piecewise

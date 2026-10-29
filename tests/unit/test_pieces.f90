! Unit tests of the piecewise functions on the log mesh (pieces_m; R-6 of
! doc/smoothness_plan.md), written before the code, on the program's
! default mesh (amesh=1.006, rr1=5e-5, rmax=45):
!   kink      f = r exp(-r) on [0, a], r exp(-r) + (r-a)**2 on (a, b], 0
!             beyond (C1 at a, a jump at b), a = 1.3, b = 2.7, as two
!             pieces each held on the whole mesh: int_0^b r**2 f dr by
!             pw_integral, and int_0^b f**2 dr by pw_product_integral,
!             against the closed forms to 1e-13; pw_node_values gives the
!             values of the piece of each node's segment, 0 beyond b, to
!             8 epsilon (exact at -O2; -ffast-math moves an ulp)
!   moving    the integrals of the same f with a and b stepped by a tenth
!             of a mesh step: smooth to 1e-13 (the closed forms move with
!             them)
!   segments  pw_segment at the bounds and beyond
! Stops with a non-zero code if any check failed.
 program test_pieces
 use constants_m, only: dp, qp
 use unit_check_m, only: check, check_done
 use pieces_m, only: piecewise_t, pw_alloc, pw_free, pw_integral, &
&                    pw_product_integral, pw_node_values, pw_segment
 implicit none

!Local variables
 real(dp), parameter :: rr1 = 5.0e-5_dp, rmax = 45.0_dp, amesh = 1.006_dp
 real(dp), allocatable :: rr(:), r2(:), fnode(:)
 type(piecewise_t) :: pw
 real(dp) :: al, aa, bb, ss, exact, err, bound(0:2)
 integer :: mmax, ii, io, kk, ifail

 al = log(amesh)
 mmax = int(log(rmax/rr1)/al)
 allocate(rr(mmax), r2(mmax), fnode(mmax))
 do ii = 1, mmax
   rr(ii) = rr1*exp(al*(ii - 1))
 end do
 r2(:) = rr**2

! kink
 aa = 1.3_dp
 bb = 2.7_dp
 call build(aa, bb)
 call pw_integral(rr, mmax, 3, r2, pw, ss)
 exact = i_r2(aa, bb)
 call check('kink: int r**2 f', abs(ss - exact)/abs(exact), 1.0e-13_dp)
 call pw_product_integral(rr, mmax, 2, pw, pw, ss, ifail)
 exact = i_ff(aa, bb)
 call check('kink: int f**2', abs(ss - exact)/abs(exact), 1.0e-13_dp)
 call pw_node_values(rr, mmax, pw, fnode)
 err = 0.0_dp
 do ii = 1, mmax
   if(rr(ii) <= aa) then
     err = max(err, abs(fnode(ii) - rr(ii)*exp(-rr(ii))))
   else if(rr(ii) <= bb) then
     err = max(err, abs(fnode(ii) - (rr(ii)*exp(-rr(ii)) + (rr(ii) - aa)**2)))
   else
     err = max(err, abs(fnode(ii)))
   end if
 end do
 call check('kink: node values', err, 8*epsilon(1.0_dp))

! moving
 err = 0.0_dp
 do io = 0, 10
   call pw_free(pw)
   call build(aa*amesh**(0.1_dp*io), bb*amesh**(0.1_dp*io))
   call pw_integral(rr, mmax, 3, r2, pw, ss)
   exact = i_r2(pw%bound(1), pw%bound(2))
   err = max(err, abs(ss - exact)/abs(exact))
   call pw_product_integral(rr, mmax, 2, pw, pw, ss, ifail)
   exact = i_ff(pw%bound(1), pw%bound(2))
   err = max(err, abs(ss - exact)/abs(exact))
 end do
 call check('moving bounds', err, 1.0e-13_dp)

! segments
 kk = 0
 if(pw_segment(pw, 0.5_dp*pw%bound(1)) /= 1) kk = kk + 1
 if(pw_segment(pw, pw%bound(1)) /= 1) kk = kk + 1
 if(pw_segment(pw, pw%bound(1)*1.000001_dp) /= 2) kk = kk + 1
 if(pw_segment(pw, pw%bound(2)) /= 2) kk = kk + 1
 if(pw_segment(pw, pw%bound(2)*1.000001_dp) /= 0) kk = kk + 1
 call check('segments', dble(kk), 0.0_dp)

 call check_done('pieces')

 contains

 subroutine build(a1, b1)

!! the two pieces of f with the bounds a1, b1, each on the whole mesh

 real(dp), intent(in) :: a1, b1

 bound(0) = 0.0_dp
 bound(1) = a1
 bound(2) = b1
 call pw_alloc(pw, mmax, 2, bound, ifail)
 pw%vals(:, 1) = rr*exp(-rr)
 pw%vals(:, 2) = rr*exp(-rr) + (rr - a1)**2
 end subroutine build

 real(dp) function i_r2(a1, b1)

!! int_0^b1 r**2 f dr: int_0^b1 r**3 e**-r + int_a1^b1 r**2 (r-a1)**2

 real(dp), intent(in) :: a1, b1
 real(qp) :: a, b, s

 a = a1
 b = b1
 s = 6.0_qp - (b**3 + 3.0_qp*b**2 + 6.0_qp*b + 6.0_qp)*exp(-b)
 s = s + ((b**5 - a**5)/5.0_qp - 2.0_qp*a*(b**4 - a**4)/4.0_qp &
&         + a**2*(b**3 - a**3)/3.0_qp)
 i_r2 = real(s, dp)
 end function i_r2

 real(dp) function i_ff(a1, b1)

!! int_0^b1 f**2 dr: int_0^b1 r**2 e**-2r + int_a1^b1 [2 r e**-r (r-a1)**2
!! + (r-a1)**4]; with s = r-a1, r e**-r (r-a1)**2 = e**-a1 (s**3 + a1 s**2)
!! e**-s

 real(dp), intent(in) :: a1, b1
 real(qp) :: a, b, c, s, t3, t2

 a = a1
 b = b1
 c = b - a
 s = 0.25_qp - exp(-2.0_qp*b)*(0.5_qp*b**2 + 0.5_qp*b + 0.25_qp)
 t3 = 6.0_qp - exp(-c)*(c**3 + 3.0_qp*c**2 + 6.0_qp*c + 6.0_qp)
 t2 = 2.0_qp - exp(-c)*(c**2 + 2.0_qp*c + 2.0_qp)
 s = s + 2.0_qp*exp(-a)*(t3 + a*t2) + c**5/5.0_qp
 i_ff = real(s, dp)
 end function i_ff

 end program test_pieces

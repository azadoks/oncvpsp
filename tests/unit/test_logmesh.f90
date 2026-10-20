! Unit tests of the operations at any radius on the log mesh (logmesh_m) and
! of the closed-form Bessel overlaps (sbf_overlap), written before the code
! (phase 1 of doc/smoothness_plan.md). Each check is on a function whose
! answer is known in closed form, on the program's default mesh (amesh=1.006,
! rr1=5e-5, rmax=45) unless it says otherwise:
!   interp_poly      a polynomial of degree lm_blend in ln(r) is reproduced
!   interp_node      at a mesh point, the mesh value is returned exactly
!   interp_order     r**2 exp(-r): error O(al**(lm_blend+1)), and small
!   interp_c1        rough data: the interpolant's slope is continuous
!                    across mesh points (a local Lagrange stencil that
!                    switches at the mesh points is not)
!   derivs           r**2 exp(-r): value and four derivatives at radii off
!                    the mesh
!   cumint           r**2 exp(-r) and r**4 exp(-r) integrated from 0 to each
!                    mesh point
!   integral         between radii off the mesh
!   integral_c1      rough data: d/drb of the integral is continuous across
!                    mesh points
!   root             a root off the mesh, to round-off; no bracket, ierr=1
!   sbf_overlap      against a 200-point Gauss-Legendre rule, for l=0..4,
!                    including q1=q2 and q1 within 1e-9 of q2
! Stops with a non-zero code after all checks if any failed.
 program test_logmesh
 use, intrinsic :: iso_fortran_env, only: dp => real64, i8 => int64, &
&                                       stderr => error_unit
 use logmesh_m, only: lm_blend, lm_interp, lm_derivs, lm_cumint, &
&                     lm_integral, lm_root
 implicit none
 external :: sbf_overlap, sbf8

!Local variables
 real(dp), parameter :: rr1 = 5.0e-5_dp, rmax = 45.0_dp
 integer :: nfail = 0
 integer(i8) :: seed = 20261020_i8

 call check_interp_poly()
 call check_interp_node()
 call check_interp_order()
 call check_interp_c1()
 call check_derivs()
 call check_cumint()
 call check_integral()
 call check_integral_c1()
 call check_root()
 call check_sbf_overlap()

 if(nfail > 0) then
   write(stderr, '(i0,a)') nfail, ' check(s) failed'
   error stop 1
 end if
 write(stderr, '(a)') 'logmesh: all checks passed'

 contains

 subroutine make_mesh(amesh, rr, mmax)

!! the log mesh of the programs: rr(i)=rr1*amesh**(i-1) up to rmax
!!
!! @param amesh  mesh ratio
!! @param rr  the mesh
!! @param mmax  its size

 real(dp), intent(in) :: amesh
 real(dp), allocatable, intent(out) :: rr(:)
 integer, intent(out) :: mmax
 real(dp) :: al
 integer :: ii

 al = log(amesh)
 mmax = int(log(rmax/rr1)/al)
 allocate(rr(mmax))
 do ii = 1, mmax
   rr(ii) = rr1*exp(al*(ii - 1))
 end do
 end subroutine make_mesh

 function uniform()

!! pseudo-random number in [-1,1), the same on every compiler

 real(dp) :: uniform
 seed = modulo(16807_i8*seed, 2147483647_i8)
 uniform = 2.0_dp*real(seed, dp)/2147483647.0_dp - 1.0_dp
 end function uniform

 subroutine check(name, err, tol)

!! counts a failure if err > tol, and reports it
!!
!! @param name  what was checked
!! @param err  its error
!! @param tol  the tolerance

 character(len=*), intent(in) :: name
 real(dp), intent(in) :: err, tol

 if(.not. (err <= tol)) then
   write(stderr, '(3a,es10.3,a,es10.3)') 'FAILED: ', name, ': error ', err, &
&    ' > ', tol
   nfail = nfail + 1
 end if
 end subroutine check

 elemental real(dp) function texp(kk, rr)

!! the k-th derivative of r**2 exp(-r)
!!
!! @param kk  order, 0..4
!! @param rr  radius

 integer, intent(in) :: kk
 real(dp), intent(in) :: rr

 select case(kk)
 case(0); texp = rr**2
 case(1); texp = 2.0_dp*rr - rr**2
 case(2); texp = 2.0_dp - 4.0_dp*rr + rr**2
 case(3); texp = -6.0_dp + 6.0_dp*rr - rr**2
 case default; texp = 12.0_dp - 8.0_dp*rr + rr**2
 end select
 texp = texp*exp(-rr)
 end function texp

 subroutine check_interp_poly()

!! a polynomial of degree lm_blend in ln(r) is reproduced

 real(dp), allocatable :: rr(:), ff(:)
 real(dp) :: rq, fq, exact, err, fmax
 integer :: mmax, ii

 call make_mesh(1.006_dp, rr, mmax)
 allocate(ff(mmax))
 do ii = 1, mmax
   ff(ii) = poly(log(rr(ii)))
 end do
 fmax = maxval(abs(ff))
 err = 0.0_dp
 do ii = 0, 49
   rq = 0.01_dp*2000.0_dp**(ii/49.0_dp)*1.0001234_dp
   call lm_interp(rr, ff, mmax, rq, fq)
   exact = poly(log(rq))
   err = max(err, abs(fq - exact))
 end do
 call check('interp_poly', err, 1.0e-12_dp*fmax)
 end subroutine check_interp_poly

 real(dp) function poly(xx)

!! sum_k (x/4)**k/(k+1), k=0..lm_blend
!!
!! @param xx  ln(r)

 real(dp), intent(in) :: xx
 integer :: kk

 poly = 0.0_dp
 do kk = 0, lm_blend
   poly = poly + (0.25_dp*xx)**kk/(kk + 1)
 end do
 end function poly

 subroutine check_interp_node()

!! at a mesh point, the mesh value is returned exactly

 real(dp), allocatable :: rr(:), ff(:)
 real(dp) :: fq, err
 integer :: mmax, ii

 call make_mesh(1.006_dp, rr, mmax)
 allocate(ff(mmax))
 do ii = 1, mmax
   ff(ii) = uniform()
 end do
 err = 0.0_dp
 do ii = 1, mmax, 37
   call lm_interp(rr, ff, mmax, rr(ii), fq)
   err = max(err, abs(fq - ff(ii)))
 end do
 call lm_interp(rr, ff, mmax, rr(mmax), fq)
 err = max(err, abs(fq - ff(mmax)))
 call check('interp_node', err, 0.0_dp)
 end subroutine check_interp_node

 subroutine check_interp_order()

!! r**2 exp(-r): the error falls as al**(lm_blend+1) on coarse meshes
!! (al = 0.1, 0.05, on r in [0.3,3] where they resolve it; finer ones are
!! at round-off already), and is round-off on the default mesh on [0.3,12]

 real(dp), parameter :: amesh(3) = [1.1051709180756477_dp, &
&                                   1.0512710963760241_dp, 1.006_dp]
 real(dp), parameter :: rtop(3) = [10.0_dp, 10.0_dp, 40.0_dp]
 real(dp) :: err(3)
 real(dp), allocatable :: rr(:), ff(:)
 real(dp) :: rq, fq
 integer :: mmax, ii, im

 do im = 1, 3
   call make_mesh(amesh(im), rr, mmax)
   allocate(ff(mmax))
   ff(:) = texp(0, rr)
   err(im) = 0.0_dp
   do ii = 0, 199
     rq = 0.3_dp*rtop(im)**(ii/199.0_dp)*1.000777_dp
     call lm_interp(rr, ff, mmax, rq, fq)
     err(im) = max(err(im), abs(fq - texp(0, rq)))
   end do
   deallocate(rr, ff)
 end do
 call check('interp_order al 0.1/0.05', 2.0_dp**lm_blend*err(2), err(1))
 call check('interp_order 1.006', err(3), 1.0e-13_dp)
 end subroutine check_interp_order

 subroutine check_interp_c1()

!! rough data: one-sided slopes at mesh points agree

 real(dp), allocatable :: rr(:), ff(:)
 real(dp) :: al, dx, fl, f0, fr, sl, sr, err, smax
 integer :: mmax, ii

 call make_mesh(1.006_dp, rr, mmax)
 allocate(ff(mmax))
 do ii = 1, mmax
   ff(ii) = uniform()
 end do
 al = log(rr(2)/rr(1))
 dx = 1.0e-4_dp*al
 err = 0.0_dp
 smax = 1.0_dp/al
 do ii = mmax/2, mmax/2 + 20
   call lm_interp(rr, ff, mmax, rr(ii)*exp(-dx), fl)
   call lm_interp(rr, ff, mmax, rr(ii), f0)
   call lm_interp(rr, ff, mmax, rr(ii)*exp(dx), fr)
   sl = (f0 - fl)/dx
   sr = (fr - f0)/dx
   smax = max(smax, abs(sl), abs(sr))
   err = max(err, abs(sr - sl))
   if(f0 /= ff(ii)) err = huge(err)
 end do
 call check('interp_c1', err, 1.0e-2_dp*smax)
 end subroutine check_interp_c1

 subroutine check_derivs()

!! r**2 exp(-r): value and four derivatives at radii off the mesh; the
!! tolerance of the k-th is tol(k)/r**k, the round-off of its stencil

 real(dp), parameter :: tol(0:4) = [1.0e-13_dp, 1.0e-11_dp, 1.0e-9_dp, &
&                                   1.0e-7_dp, 1.0e-5_dp]
 real(dp), parameter :: radii(4) = [0.7_dp, 1.3_dp, 2.718281828459045_dp, &
&                                    5.5_dp]
 real(dp), allocatable :: rr(:), ff(:)
 real(dp) :: dd(5), err(0:4)
 character(len=16) :: name
 integer :: mmax, ii, kk

 call make_mesh(1.006_dp, rr, mmax)
 allocate(ff(mmax))
 ff(:) = texp(0, rr)
 err(:) = 0.0_dp
 do ii = 1, 4
   call lm_derivs(rr, ff, mmax, radii(ii), 4, dd)
   do kk = 0, 4
     err(kk) = max(err(kk), abs(dd(kk + 1) - texp(kk, radii(ii))) &
&                          *radii(ii)**kk)
   end do
 end do
 do kk = 0, 4
   write(name, '(a,i0)') 'derivs order ', kk
   call check(trim(name), err(kk), tol(kk))
 end do
 end subroutine check_derivs

 subroutine check_cumint()

!! r**2 exp(-r) (mm=2) and r**4 exp(-r) (mm=4) from 0 to each mesh point

 real(dp), allocatable :: rr(:), ff(:), cum(:)
 real(dp) :: err2, err4, xx
 integer :: mmax, ii

 call make_mesh(1.006_dp, rr, mmax)
 allocate(ff(mmax), cum(mmax))
 ff(:) = texp(0, rr)
 call lm_cumint(rr, ff, mmax, 2, cum)
 err2 = 0.0_dp
 do ii = 1, mmax
   xx = rr(ii)
   err2 = max(err2, abs(cum(ii) - (2.0_dp - exp(-xx)*(xx**2 + 2.0_dp*xx &
&                                                      + 2.0_dp))))
 end do
 ff(:) = rr**4*exp(-rr)
 call lm_cumint(rr, ff, mmax, 4, cum)
 err4 = 0.0_dp
 do ii = 1, mmax
   xx = rr(ii)
   err4 = max(err4, abs(cum(ii) - (24.0_dp - exp(-xx)*(xx**4 + 4.0_dp*xx**3 &
&             + 12.0_dp*xx**2 + 24.0_dp*xx + 24.0_dp)))/12.0_dp)
 end do
 call check('cumint r**2 exp(-r)', err2, 1.0e-13_dp)
 call check('cumint r**4 exp(-r)', err4, 1.0e-13_dp)
 end subroutine check_cumint

 subroutine check_integral()

!! r**2 exp(-r) between radii off the mesh, and from 0

 real(dp), allocatable :: rr(:), ff(:)
 real(dp) :: ss, err, ra, rb
 integer :: mmax

 call make_mesh(1.006_dp, rr, mmax)
 allocate(ff(mmax))
 ff(:) = texp(0, rr)
 ra = 0.123456_dp
 rb = 3.14159_dp
 call lm_integral(rr, ff, mmax, 2, ra, rb, ss)
 err = abs(ss - (prim(rb) - prim(ra)))
 call lm_integral(rr, ff, mmax, 2, 0.0_dp, rb, ss)
 err = max(err, abs(ss - prim(rb)))
 call lm_integral(rr, ff, mmax, 2, rb, ra, ss)
 err = max(err, abs(ss + (prim(rb) - prim(ra))))
 call check('integral', err, 1.0e-13_dp)
 end subroutine check_integral

 real(dp) function prim(xx)

!! int_0^x r**2 exp(-r) dr
!!
!! @param xx  upper limit

 real(dp), intent(in) :: xx
 prim = 2.0_dp - exp(-xx)*(xx**2 + 2.0_dp*xx + 2.0_dp)
 end function prim

 subroutine check_integral_c1()

!! rough integrand: one-sided derivatives of the integral in its upper limit
!! agree at mesh points

 real(dp), allocatable :: rr(:), ff(:)
 real(dp) :: al, dx, sl, s0, sr, dl, dr, err, dmax
 integer :: mmax, ii

 call make_mesh(1.006_dp, rr, mmax)
 allocate(ff(mmax))
 do ii = 1, mmax
   ff(ii) = 1.0_dp + 0.5_dp*uniform()
 end do
 al = log(rr(2)/rr(1))
 dx = 1.0e-4_dp*al
 err = 0.0_dp
 dmax = 0.0_dp
 do ii = mmax/2, mmax/2 + 20
   call lm_integral(rr, ff, mmax, 0, 0.0_dp, rr(ii)*exp(-dx), sl)
   call lm_integral(rr, ff, mmax, 0, 0.0_dp, rr(ii), s0)
   call lm_integral(rr, ff, mmax, 0, 0.0_dp, rr(ii)*exp(dx), sr)
   dl = (s0 - sl)/(rr(ii)*(1.0_dp - exp(-dx)))
   dr = (sr - s0)/(rr(ii)*(exp(dx) - 1.0_dp))
   dmax = max(dmax, abs(dl), abs(dr))
   err = max(err, abs(dr - dl))
! the derivative is the integrand there, which lies in [0.5,1.5]
   if(min(dl, dr) < 0.25_dp .or. max(dl, dr) > 1.75_dp) err = huge(err)
 end do
 call check('integral_c1', err, 1.0e-2_dp*dmax)
 end subroutine check_integral_c1

 subroutine check_root()

!! ln(r/1.2345678) (linear in ln r) and r**2 exp(-r) - 0.3, roots off the
!! mesh; a bracket without a sign change gives ierr=1

 real(dp), allocatable :: rr(:), ff(:)
 real(dp) :: rq, exact, err
 integer :: mmax, ierr, ii

 call make_mesh(1.006_dp, rr, mmax)
 allocate(ff(mmax))
 ff(:) = log(rr/1.2345678_dp)
 call lm_root(rr, ff, mmax, 1.0_dp, 2.0_dp, rq, ierr)
 err = abs(rq/1.2345678_dp - 1.0_dp)
 if(ierr /= 0) err = huge(err)
 call check('root ln(r/r0)', err, 1.0e-14_dp)

 ff(:) = texp(0, rr) - 0.3_dp
 exact = 0.8_dp
 do ii = 1, 50
   exact = exact - (texp(0, exact) - 0.3_dp)/texp(1, exact)
 end do
 call lm_root(rr, ff, mmax, 0.5_dp, 1.5_dp, rq, ierr)
 err = abs(rq - exact)
 if(ierr /= 0) err = huge(err)
 call check('root r**2 exp(-r) - 0.3', err, 1.0e-12_dp)

 call lm_root(rr, ff, mmax, 1.5_dp, 2.5_dp, rq, ierr)
 call check('root without a bracket', real(abs(ierr - 1), dp), 0.0_dp)
 end subroutine check_root

 subroutine check_sbf_overlap()

!! int_0^rc r**2 j_l(q1 r) j_l(q2 r) dr against Gauss-Legendre on [0,rc],
!! relative to sqrt(S(q1,q1) S(q2,q2))

 integer, parameter :: ngl = 200, ncase = 6
 real(dp), parameter :: rc = 2.7182818_dp
 real(dp) :: qpair(2, ncase), xg(ngl), wg(ngl), ss, s11, s22, ref, err
 integer :: ll, ic

 qpair(:, 1) = [1.3_dp, 2.9_dp]
 qpair(:, 2) = [4.1_dp, 4.1_dp]
 qpair(:, 3) = [4.1_dp, 4.1_dp*(1.0_dp + 1.0e-9_dp)]
 qpair(:, 4) = [0.2_dp, 9.7_dp]
 qpair(:, 5) = [12.3_dp, 15.9_dp]
 qpair(:, 6) = [0.05_dp, 0.07_dp]
 call gauss_legendre(ngl, xg, wg)
 xg(:) = 0.5_dp*rc*(xg + 1.0_dp)
 wg(:) = 0.5_dp*rc*wg
 err = 0.0_dp
 do ll = 0, 4
   do ic = 1, ncase
     ref = gl_overlap(ll, qpair(1, ic), qpair(2, ic), xg, wg, ngl)
     s11 = gl_overlap(ll, qpair(1, ic), qpair(1, ic), xg, wg, ngl)
     s22 = gl_overlap(ll, qpair(2, ic), qpair(2, ic), xg, wg, ngl)
     call sbf_overlap(ll, qpair(1, ic), qpair(2, ic), rc, ss)
     err = max(err, abs(ss - ref)/sqrt(s11*s22))
   end do
 end do
 call check('sbf_overlap', err, 1.0e-13_dp)
 end subroutine check_sbf_overlap

 real(dp) function gl_overlap(ll, q1, q2, xg, wg, ngl)

!! int_0^rc r**2 j_l(q1 r) j_l(q2 r) dr by a Gauss-Legendre rule on [0,rc]
!!
!! @param ll  angular momentum
!! @param q1,q2  wave vectors
!! @param xg,wg  nodes and weights of the rule on [0,rc]
!! @param ngl  number of nodes

 integer, intent(in) :: ll, ngl
 real(dp), intent(in) :: q1, q2, xg(ngl), wg(ngl)
 real(dp) :: sb1(8), sb2(8)
 integer :: ii

 gl_overlap = 0.0_dp
 do ii = 1, ngl
   call sbf8(8, q1*xg(ii), sb1)
   call sbf8(8, q2*xg(ii), sb2)
   gl_overlap = gl_overlap + wg(ii)*xg(ii)**2*sb1(ll + 1)*sb2(ll + 1)
 end do
 end function gl_overlap

 subroutine gauss_legendre(nn, xx, ww)

!! nodes and weights of the nn-point Gauss-Legendre rule on [-1,1]
!!
!! @param nn  number of nodes
!! @param xx  nodes
!! @param ww  weights

 integer, intent(in) :: nn
 real(dp), intent(out) :: xx(nn), ww(nn)
 real(dp) :: p0, p1, p2, dp1, zz, dz
 integer :: ii, jj, it

 do ii = 1, nn
   zz = cos(acos(-1.0_dp)*(ii - 0.25_dp)/(nn + 0.5_dp))
   do it = 1, 100
     p0 = 1.0_dp
     p1 = zz
     do jj = 2, nn
       p2 = ((2*jj - 1)*zz*p1 - (jj - 1)*p0)/jj
       p0 = p1
       p1 = p2
     end do
     dp1 = nn*(zz*p1 - p0)/(zz**2 - 1.0_dp)
     dz = p1/dp1
     zz = zz - dz
     if(abs(dz) < 1.0e-16_dp) exit
   end do
   xx(ii) = zz
   ww(ii) = 2.0_dp/((1.0_dp - zz**2)*dp1**2)
 end do
 end subroutine gauss_legendre

 end program test_logmesh

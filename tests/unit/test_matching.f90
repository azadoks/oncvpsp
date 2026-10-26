! Unit tests of the bound-state solvers' smoothness in the potential (S-8
! of doc/smoothness_plan.md): hydrogen-like -Z/r with Z stepped through 41
! values in [1, 1.1], which moves each state's turning point across about
! 16 mesh points. The eigenvalue over Z**2 (-1/(2n**2) plus the mesh's
! discretization error, a smooth function of Z) and u at a fixed mesh point
! beyond the turning point are fitted with a polynomial of degree 10 in Z;
! the residual must be below 1e-14 of max(|E|, 0.2) for the eigenvalue and
! 1e-12 of the largest |u| for u. Matched at the turning point's mesh
! index, the eigenvalue had a sawtooth of about 1e-13 relative each time
! that index moved. On the program's default mesh (amesh=1.006, rr1=5e-5,
! rmax=45):
!   lschfb     2s, 3p; non- and scalar-relativistic
!   ldiracfb   2s, 3p (j=3/2)
!   lschvkbb   2s, 3p with no projectors
! Stops with a non-zero code if any check failed.
 program test_matching
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
 use pieces_m, only: piecewise_t, pw_alloc, pw_free
 use ldiracfb_m, only: ldiracfb
 use lschfb_m, only: lschfb
 use lschvkbb_m, only: lschvkbb
 use lapack_m, only: dgels
 implicit none

!Local variables
 integer :: ifail
 integer, parameter :: npts = 41, ndeg = 10
 real(dp), parameter :: rr1 = 5.0e-5_dp, rmax = 45.0_dp, amesh = 1.006_dp
 real(dp), parameter :: tol_e = 1.0e-14_dp, tol_u = 1.0e-12_dp
 real(dp), parameter :: zlo = 1.0_dp, zhi = 1.1_dp
 real(dp), allocatable :: rr(:), vv(:), uu(:, :), up(:, :)
 type(piecewise_t) :: vpw, vkbpw(0)
 real(dp) :: al, zz, ee, emin, emax, evkb(1), es(npts), us(npts), zs(npts)
 real(dp) :: rq
 integer :: nfail = 0, mmax, ii, iz, ierr, mch, isolv, nn, ll, ist
 character(len=40) :: name

 al = log(amesh)
 mmax = int(log(rmax/rr1)/al)
 allocate(rr(mmax), vv(mmax), uu(mmax, 2), up(mmax, 2))
 do ii = 1, mmax
   rr(ii) = rr1*exp(al*(ii - 1))
 end do
 evkb(1) = 0.0_dp

 do isolv = 1, 4
   do ist = 1, 2
     nn = ist + 1
     ll = ist - 1
     do iz = 1, npts
       zz = zlo + (zhi - zlo)*(iz - 1)/real(npts - 1, dp)
       zs(iz) = zz
       vv(:) = -zz/rr
       ee = -0.5_dp*zz**2/nn**2
       uu(:, :) = 0.0_dp
       up(:, :) = 0.0_dp
       mch = 0
       select case(isolv)
       case(1, 2)
         call lschfb(nn, ll, ierr, ee, rr, vv, uu(:, 1), up(:, 1), zz, &
&                    mmax, mch, isolv - 1, ifail)
         if(ifail /= 0) error stop 1
       case(3)
         call ldiracfb(nn, ll, -(ll + 1), ierr, ee, rr, zz, vv, uu, up, &
&                      mmax, mch, ifail)
         if(ifail /= 0) error stop 1
       case(4)
         emin = -zz**2
         emax = 0.0_dp
         call pw_alloc(vpw, mmax, 1, [0.0_dp, rr(mmax)])
         vpw%vals(:, 1) = vv
         call lschvkbb(nn, ll, 0, ierr, ee, emin, emax, rr, vv, vpw, &
&                      vkbpw, evkb, uu(:, 1), up(:, 1), mmax, mch, ifail)
         if(ifail /= 0) error stop 1
         call pw_free(vpw)
       case default
         ierr = -1
       end select
       if(ierr /= 0) then
         write(stderr, '(a,i0,a,i0,a,i0)') 'FAILED: solver ', isolv, &
&          ' n ', nn, ' ierr ', ierr
         nfail = nfail + 1
       end if
       es(iz) = ee/zz**2
! u at the mesh point nearest 3 n**2, beyond the turning point 2 n**2/Z
! for every Z
       rq = 3.0_dp*nn**2
       ii = 1 + nint(log(rq/rr1)/al)
       us(iz) = uu(ii, 1)
     end do
     write(name, '(a,i0,a,i0,a,i0)') 'solver ', isolv, ' n=', nn, ' l=', ll
     call check(trim(name)//' E/Z**2', es, max(maxval(abs(es)), 0.2_dp), &
&               tol_e)
     call check(trim(name)//' u', us, maxval(abs(us)), tol_u)
   end do
 end do

 if(nfail > 0) then
   write(stderr, '(i0,a)') nfail, ' check(s) failed'
   error stop 1
 end if
 write(stderr, '(a)') 'matching: all checks passed'

 contains

 subroutine check(label, ys, scale, tol)

!! the largest residual of the degree-ndeg least-squares fit of ys in zs,
!! relative to scale, against tol

 character(len=*), intent(in) :: label
 real(dp), intent(in) :: ys(npts), scale, tol
 real(dp) :: aa(npts, ndeg + 1), bb(npts), work(64*(ndeg + 1)), tt
 integer :: jj, kk, info

 do jj = 1, npts
   tt = 2.0_dp*(zs(jj) - zlo)/(zhi - zlo) - 1.0_dp
   aa(jj, 1) = 1.0_dp
   aa(jj, 2) = tt
   do kk = 2, ndeg
     aa(jj, kk + 1) = ((2*kk - 1)*tt*aa(jj, kk) - (kk - 1)*aa(jj, kk - 1))/kk
   end do
 end do
 bb(:) = ys(:)
 call dgels('N', npts, ndeg + 1, 1, aa, npts, bb, npts, work, size(work), &
&           info)
 aa(:, :) = 0.0_dp
 do jj = 1, npts
   tt = 2.0_dp*(zs(jj) - zlo)/(zhi - zlo) - 1.0_dp
   aa(jj, 1) = 1.0_dp
   aa(jj, 2) = tt
   do kk = 2, ndeg
     aa(jj, kk + 1) = ((2*kk - 1)*tt*aa(jj, kk) - (kk - 1)*aa(jj, kk - 1))/kk
   end do
 end do
 tt = maxval(abs(ys - matmul(aa, bb(1:ndeg + 1))))/scale
 write(stderr, '(a,a,es9.2)') label, ': fit residual ', tt
 if(.not. (tt <= tol)) then
   write(stderr, '(a,a,es10.3,a,es9.2)') 'FAILED: ', label, tt, ' > ', tol
   nfail = nfail + 1
 end if
 end subroutine check

 end program test_matching

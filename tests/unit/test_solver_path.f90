! Unit tests of the eigen-solvers' independence of their iteration path
! (T5, phase 2 of doc/smoothness_plan.md), written before the code: each
! solver started from a guess below and one above its eigenvalue must return
! the same eigenvalue to 1e-13 relative (of max(|E|, 0.2), the solvers' own
! scale) and the same normalized function to 1e-12 of its largest value. A
! solver that stops at its tolerance with the guess before its last
! correction returns an energy that depends on the path by up to that
! tolerance. On the program's default mesh (amesh=1.006, rr1=5e-5, rmax=45):
!   lschfb     hydrogen 1s, 2p, 3d; non- and scalar-relativistic
!   ldiracfb   hydrogen 1s, 2p (j=1/2, 3/2)
!   lschvkbb   a Gaussian well with one projector per l: 1s, 2p
!   lschvkbbe  the same, at the log derivative of the lschvkbb state
!   lschpsbar  the Gaussian well inside a barrier at 8 bohr: 1s, 2s
! Stops with a non-zero code if any check failed.
 program test_solver_path
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: stderr => error_unit
 use unit_check_m, only: check_done, nfail
 use pieces_m, only: piecewise_t, pw_alloc, pw_free
 use ldiracfb_m, only: ldiracfb
 use lschfb_m, only: lschfb
 use lschpsbar_m, only: lschpsbar
 use lschvkbb_m, only: lschvkbb
 use lschvkbbe_m, only: lschvkbbe
 use phsft_m, only: phsft_mch
 use relativity_m, only: c_light => c_light_default
 implicit none
 integer, parameter :: iout=6

!Local variables
 integer :: ifail
 real(dp), parameter :: rr1 = 5.0e-5_dp, rmax = 45.0_dp, amesh = 1.006_dp
 real(dp), parameter :: tol_e = 1.0e-13_dp, tol_u = 1.0e-12_dp
 real(dp), allocatable :: rr(:), vv(:), vloc(:), uu(:, :), &
&                         up(:, :), ud(:, :, :), udp(:, :, :)
 type(piecewise_t) :: vlocpw, vkbpw(1)
 real(dp) :: al, ee(2), evkb(1), emin, emax, uld, rcut, rm
 integer :: mmax, ii, ip, ierr(2), mch(2), nn, ll, kap, mbar
 character(len=40) :: name

 al = log(amesh)
 mmax = int(log(rmax/rr1)/al)
 allocate(rr(mmax), vv(mmax), vloc(mmax), uu(mmax, 2), &
&         up(mmax, 2), ud(mmax, 2, 2), udp(mmax, 2, 2))
 do ii = 1, mmax
   rr(ii) = rr1*exp(al*(ii - 1))
 end do
 vv(:) = -1.0_dp/rr
 vloc(:) = -3.0_dp*exp(-0.5_dp*rr**2)
 evkb(1) = -1.0_dp
 rcut = 3.0_dp
 call pw_alloc(vlocpw, mmax, 1, [0.0_dp, rr(mmax)])
 vlocpw%vals(:, 1) = vloc

! lschfb: hydrogen, E = -1/(2n**2)
 do ip = 0, 1
   do ll = 0, 2
     nn = ll + 1
     do ii = 1, 2
       ee(ii) = guess(-0.5_dp/nn**2, ii)
       call lschfb(nn, ll, ierr(ii), ee(ii), rr, vv, uu(:, ii), up(:, ii), &
&                  1.0_dp, mmax, mch(ii), ip,c_light, iout, ifail)
       if(ifail /= 0) error stop 1
     end do
     write(name, '(a,i0,a,i0,a,i0)') 'lschfb irel ', ip, ' n ', nn, ' l ', ll
     call compare(name, uu)
   end do
 end do

! ldiracfb: hydrogen, kap = l (j=l-1/2) and -(l+1) (j=l+1/2)
 do kap = -2, 1
   if(kap == 0) cycle
   ll = merge(kap, -kap - 1, kap > 0)
   nn = ll + 1
   do ii = 1, 2
     ee(ii) = guess(-0.5_dp/nn**2, ii)
     call ldiracfb(nn, ll, kap,c_light, ierr(ii), ee(ii), rr, 1.0_dp, vv, &
&                  ud(:, :, ii), udp(:, :, ii), mmax, mch(ii), iout, ifail)
     if(ifail /= 0) error stop 1
     uu(:, ii) = ud(:, 1, ii)
     up(:, ii) = udp(:, 1, ii)
   end do
   write(name, '(a,i0,a,i0)') 'ldiracfb n ', nn, ' kap ', kap
   call compare(name, uu)
 end do

! lschvkbb, then lschvkbbe at its log derivative: a Gaussian well, one
! projector r**(l+1) (1 - (r/rcut)**2)**4 inside rcut per l, one piece
! ending at rcut
 do ll = 0, 1
   nn = ll + 1
   call pw_alloc(vkbpw(1), mmax, 1, [0.0_dp, rcut])
   vkbpw(1)%vals(:, 1) = rr**(ll + 1)*(1.0_dp - (rr/rcut)**2)**4
   do ii = 1, 2
     ee(ii) = guess(merge(-1.3194_dp, -0.3914_dp, ll == 0), ii)
     emin = -3.0_dp
     emax = 0.0_dp
     call lschvkbb(nn, ll, 1, ierr(ii), ee(ii), emin, emax, rr, vloc, &
&                  vlocpw, vkbpw, evkb, uu(:, ii), up(:, ii), mmax, mch(ii), &
&                  iout, ifail)
     if(ifail /= 0) error stop 1
   end do
   write(name, '(a,i0,a,i0)') 'lschvkbb n ', nn, ' l ', ll
   call compare(name, uu)

   uld = up(mch(1), 1)/uu(mch(1), 1)
   rm = rr(mch(1))
   mch(2) = phsft_mch(rr, mmax, rm)
   do ii = 1, 2
     ee(ii) = guess(merge(-1.3194_dp, -0.3914_dp, ll == 0), ii)
     emin = -3.0_dp
     emax = 0.0_dp
     call lschvkbbe(nn, ll, 1, ierr(ii), ee(ii), uld, emin, emax, rr, vloc, &
&                   vlocpw, vkbpw, evkb, uu(:, ii), up(:, ii), mmax, mch(2), &
&                   rm, iout, ifail)
     if(ifail /= 0) error stop 1
   end do
   write(name, '(a,i0,a,i0)') 'lschvkbbe n ', nn, ' l ', ll
   call compare(name, uu)
   call pw_free(vkbpw(1))
 end do

! lschpsbar: the Gaussian well with a barrier at 8 bohr, tht=0
 mbar = 1 + nint(log(8.0_dp/rr1)/al)
 ll = 0
 do nn = 1, 2
   do ii = 1, 2
     ee(ii) = guess(merge(-0.9054_dp, 0.0984_dp, nn == 1), ii)
     emin = 0.0_dp
     emax = 5.0_dp
     call lschpsbar(nn, ll, ierr(ii), ee(ii), emin, emax, rr, vloc, &
&                   uu(:, ii), up(:, ii), mmax, mbar, 0.0_dp)
   end do
   write(name, '(a,i0,a,i0)') 'lschpsbar n ', nn, ' l ', ll
   call compare(name, uu)
 end do

 call check_done('solver_path')

 contains

 function guess(eig, iside)

!! a start guess below (iside=1) or above (2) an eigenvalue, by a tenth of
!! max(|E|, 0.2), well inside the range of the same node count
!!
!! @param eig  the eigenvalue, roughly
!! @param iside  1 below, 2 above

 real(dp), intent(in) :: eig
 integer, intent(in) :: iside
 real(dp) :: guess

 guess = eig + merge(-0.1_dp, 0.1_dp, iside == 1)*max(abs(eig), 0.2_dp)
 end function guess

 subroutine compare(label, uq)

!! counts a failure if either solve failed, or if the two eigenvalues or the
!! two normalized functions differ by more than tol_e and tol_u
!!
!! @param label  what was solved
!! @param uq  the functions of the two solves

 character(len=*), intent(in) :: label
 real(dp), intent(in) :: uq(mmax, 2)
 real(dp) :: de, du, scale

 if(any(ierr /= 0)) then
   write(stderr, '(3a,2i4)') 'FAILED: ', trim(label), ': ierr', ierr
   nfail = nfail + 1
   return
 end if
 de = abs(ee(1) - ee(2))/max(abs(ee(1)), 0.2_dp)
 scale = maxval(abs(uq(:, 1)))
 du = maxval(abs(uq(:, 1) - uq(:, 2)))/scale
 write(stderr, '(a,a,es20.12,a,es9.2,a,es9.2)') trim(label), ': E', ee(1), &
&  '  dE', de, '  du', du
 if(.not. (de <= tol_e)) then
   write(stderr, '(3a,es10.3,a,es10.3)') 'FAILED: ', trim(label), &
&    ': eigenvalues differ by ', de, ' > ', tol_e
   nfail = nfail + 1
 end if
 if(.not. (du <= tol_u)) then
   write(stderr, '(3a,es10.3,a,es10.3)') 'FAILED: ', trim(label), &
&    ': functions differ by ', du, ' > ', tol_u
   nfail = nfail + 1
 end if
 end subroutine compare

 end program test_solver_path

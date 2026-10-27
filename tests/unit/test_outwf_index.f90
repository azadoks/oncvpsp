! Unit tests of the eigenvalue index of the outward solution with
! projectors (outwf_index): the number of eigenvalues below E of the radial
! equation with the projectors on (0, R], u(R) = 0. On the program's default
! mesh (amesh=1.006, rr1=5e-5, rmax=45), a Gaussian well
! V = -3 exp(-r**2/2), one projector beta = r**(l+1) (1 - (r/rc)**2)**4 on
! [0, rc], rc = 3, and R the last point of the grid of the segment beyond
! rc (rc amesh**j) at or below the mesh point nearest 4 bohr, for l = 0, 1
! and coefficients evkb = -4, 4:
!   index    at energies from -3 to 6 Ha in steps of 0.05, the index of
!            outwf_solve against the eigenvalues below E of the same
!            Hamiltonian by finite differences on a uniform grid of (0, R]
!            (800 points; energies within 0.01 Ha of one of those, where its
!            discretization error could change the count, are skipped)
!   origin   among those energies there are some where the node count of
!            the solution without its cutoff (outwf_count from r = 0)
!            differs from the index: a node that entered at r = 0 without
!            an eigenvalue passing, the case the index is for
! Stops with a non-zero code if any check failed.
 program test_outwf_index
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
 use pieces_m, only: piecewise_t, pw_alloc, pw_free
 use outwf_m, only: outwf_t, outwf_setup, outwf_solve, outwf_count, &
&                   outwf_free
 use lapack_m, only: dsyev
 implicit none
 integer, parameter :: iout=6

!Local variables
 integer :: ifail
 real(dp), parameter :: rr1 = 5.0e-5_dp, rmax = 45.0_dp, amesh = 1.006_dp
 real(dp), parameter :: rc = 3.0_dp, rbox = 4.0_dp, de = 0.05_dp
 real(dp), parameter :: emin = -3.0_dp, emax = 6.0_dp, egap = 0.01_dp
 integer, parameter :: nfd = 800, ncase = 4
 integer, parameter :: lcase(ncase) = [0, 0, 1, 1]
 real(dp), parameter :: ecase(ncase) = [-4.0_dp, 4.0_dp, -4.0_dp, 4.0_dp]
 real(dp), allocatable :: rr(:), vloc(:), uu(:), up(:), efd(:)
 type(piecewise_t) :: vlocpw, vkbpw(1)
 type(outwf_t) :: ow
 real(dp) :: al, ee, evkb(1), rlast
 integer :: nfail = 0, mmax, ibox, icase, ll, ie, ne, neig, nall
 integer :: nfd_below, nbad, norigin, nskip

 al = log(amesh)
 mmax = int(log(rmax/rr1)/al)
 allocate(rr(mmax), vloc(mmax), uu(mmax), up(mmax), efd(nfd))
 do ie = 1, mmax
   rr(ie) = rr1*exp(al*(ie - 1))
 end do
 ibox = 1 + nint(log(rbox/rr1)/al)
 rlast = rc*amesh**floor(log(rr(ibox)/rc)/al)
 vloc(:) = -3.0_dp*exp(-0.5_dp*rr**2)
 call pw_alloc(vlocpw, mmax, 1, [0.0_dp, rr(mmax)])
 vlocpw%vals(:, 1) = vloc
 ne = 1 + nint((emax - emin)/de)

 norigin = 0
 do icase = 1, ncase
   ll = lcase(icase)
   evkb(1) = ecase(icase)
   call pw_alloc(vkbpw(1), mmax, 1, [0.0_dp, rc])
   vkbpw(1)%vals(:, 1) = rr**(ll + 1)*(1.0_dp - (rr/rc)**2)**4
   call fd_eigenvalues(ll, evkb(1), rlast, efd)
   call outwf_setup(ll, rr, mmax, vloc, vlocpw, 1, vkbpw, ow, iout, ifail)
   if(ifail /= 0) error stop 1
   nbad = 0
   nskip = 0
   do ie = 1, ne
     ee = emin + de*(ie - 1)
     if(minval(abs(efd - ee)) < egap) then
       nskip = nskip + 1
       cycle
     end if
     call outwf_solve(ow, ee, evkb, ibox, ibox, uu, up, neig, &
&          iout=iout, ifail=ifail)
     if(ifail /= 0) error stop 1
     nall = outwf_count(ow, 1, -1.0_dp, ibox)
     nfd_below = count(efd < ee)
     if(neig /= nfd_below) then
       nbad = nbad + 1
       write(stderr, '(a,i0,a,f6.2,a,f6.2,3(a,i0))') '  l=', ll, &
&        ' evkb=', evkb(1), ' E=', ee, ': index ', neig, ', fd ', &
&        nfd_below, ', nodes ', nall
     end if
     if(nall /= neig) norigin = norigin + 1
   end do
   write(stderr, '(a,i0,a,f6.2,3(a,i0))') 'index l=', ll, ' evkb=', &
&    evkb(1), ': ', ne - nskip, ' energies, ', nbad, ' wrong, skipped ', &
&    nskip
   if(nbad > 0) nfail = nfail + 1
   call outwf_free(ow)
   call pw_free(vkbpw(1))
 end do
 write(stderr, '(a,i0)') 'origin: energies where the node count is not '// &
&  'the index ', norigin
 if(norigin == 0) then
   write(stderr, '(a)') 'FAILED: origin, no node entered at r = 0'
   nfail = nfail + 1
 end if

 call pw_free(vlocpw)
 if(nfail > 0) then
   write(stderr, '(i0,a)') nfail, ' check(s) failed'
   error stop 1
 end if
 write(stderr, '(a)') 'outwf_index: all checks passed'

 contains

 subroutine fd_eigenvalues(lq, eq, rb, ev)

!! the eigenvalues of -u''/2 + (l(l+1)/(2 r**2) + V) u + eq beta <beta|u>
!! on (0, rb], u(0) = u(rb) = 0, by second-order finite differences on nfd
!! interior points (dsyev)

 integer, intent(in) :: lq
 real(dp), intent(in) :: eq, rb
 real(dp), intent(out) :: ev(nfd)
 real(dp), allocatable :: hh(:, :), work(:), rq(:), bq(:)
 real(dp) :: hq
 integer :: i1, info

 allocate(hh(nfd, nfd), work(64*nfd), rq(nfd), bq(nfd))
 hq = rb/(nfd + 1)
 do i1 = 1, nfd
   rq(i1) = hq*i1
 end do
 bq(:) = 0.0_dp
 where(rq < rc) bq = rq**(lq + 1)*(1.0_dp - (rq/rc)**2)**4
 do i1 = 1, nfd
   hh(:, i1) = eq*hq*bq*bq(i1)
   hh(i1, i1) = hh(i1, i1) + 1.0_dp/hq**2 + 0.5_dp*lq*(lq + 1)/rq(i1)**2 &
&    - 3.0_dp*exp(-0.5_dp*rq(i1)**2)
   if(i1 > 1) hh(i1 - 1, i1) = hh(i1 - 1, i1) - 0.5_dp/hq**2
   if(i1 < nfd) hh(i1 + 1, i1) = hh(i1 + 1, i1) - 0.5_dp/hq**2
 end do
 call dsyev('N', 'U', nfd, hh, nfd, ev, work, size(work), info)
 if(info /= 0) error stop 'fd_eigenvalues: dsyev failed'
 deallocate(hh, work, rq, bq)
 end subroutine fd_eigenvalues

 end program test_outwf_index

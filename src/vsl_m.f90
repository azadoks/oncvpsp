!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
module vsl_m
   use sbf_m, only: sbf_basis, sbf_basis_con, sbf_rc_der, qroots
   use wf_rc_der_m, only: wf_rc_der
   use eresid_m, only: eresid
   implicit none
   private
   public :: run_optimize
   public :: optimize
   public :: pspot
   public :: const_basis
contains
!> calls routines to generate optimized pseudo wave function and semi-local
!> pseudopotential and prints diagnostic information on process and
!> convergence performance
subroutine run_optimize(eig, ll, mmax, mxprj, rr, uua, qq,&
&                        irc, qcut, qmsbf, ncon_in, nbas_in, npr, &
&                        psopt, vpsp, vkb, vae, cvgplt)
   implicit none
   integer, parameter :: dp = kind(1.0d0)
   real(dp), parameter :: Ha_eV = 27.21138386d0 ! 1 Hartree, in eV

   !Input variables
   !> ll  angular momentum
   integer, intent(in) :: ll
   !> mmax  dimension of log grid
   integer, intent(in) :: mmax
   !> mxprj  dimension of number of projectors
   integer, intent(in) :: mxprj
   !> irc  index of core radius
   integer, intent(in out) :: irc
   !> ncon_in  number of constraints for matching psuedo and AE wave functions
   integer, intent(in) :: ncon_in
   !> nbas_in  number of basis functions for pseudo wave function
   integer, intent(in) :: nbas_in
   !> np_in   number of projectors = 0,1,2
   integer, intent(in) :: npr
   !> rr  log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> uu  all-electron wave function for first projector
   !> uu2  all-electron wave function for second projector
   real(dp), intent(in) :: uua(mmax, mxprj)
   !> vae  all-electron potential
   real(dp), intent(in) :: vae(mmax)
   !> qq  2x2 matrix of norms and overlaps of uu and uu2
   real(dp), intent(in) :: qq(mxprj, mxprj)
   !> eig  energy at which pseudo wave function is computed
   real(dp), intent(in) :: eig(mxprj)
   !> qcut  q cutoff defining residual energy
   real(dp), intent(in) :: qcut

   !Output variables
   !> qmsbf  maximum q in sbfs for this ll
   real(dp), intent(out) :: qmsbf
   !> psopt  optimized pseudo wave function(s)
   real(dp), intent(out) :: psopt(mmax, mxprj)
   !> vpsp  corresponding pseudopotential
   real(dp), intent(out) :: vpsp(mmax)
   !> vkb  Vanderbilt-Kleinman-Bylander projectors (without local v correction)
   real(dp), intent(out) :: vkb(mmax, mxprj)
   !> cvgplt  energy error vs. cutoff energy for plotting
   real(dp), intent(out) :: cvgplt(2, 7, mxprj)

   !Local variables
   real(dp) :: uord(6)
   real(dp) :: al, rc, ulgd, tt
   real(dp) :: sn, ps0norm, amesh, ro
   real(dp) :: err, lerr, qerr, ehaerr, eeverr
   real(dp) :: sbf(mmax)
   real(dp) :: cons(6)
   real(dp), allocatable :: orbasis(:, :)
   real(dp), allocatable :: orbasis_der(:, :), pswf0_sb(:), pswfnull_sb(:, :)
   real(dp), allocatable :: pswf0_or(:), pswfnull_or(:, :)
   real(dp), allocatable :: work(:)
   integer :: ii, iprj, jj, l1, lmax, nbas, ncon, nconmx
   logical :: found_root

   integer :: nnull, nqout

   real(dp), allocatable :: qroot(:), sbasis(:)
   real(dp), allocatable :: eresid0(:), eresiddot(:, :)
   real(dp), allocatable :: eresidmat(:, :, :), pswfresid_sb(:, :)
   real(dp), allocatable :: qout(:), eresq(:), leresq(:)

   real(dp), allocatable :: pswfopt_sb(:), pswfopt_or(:)
   real(dp) :: eresidmin
   real(dp) :: ekin_anal, ekin_num, qmax

   !parameters whose default values here give well-converged results
   ! increasing dq and/or dr will speed up the code at some cost
   ! in optimization accuracy
   real(dp) :: dq, dr, dqout
   dq = 0.02d0  !linear integration mesh spacing for above
   dr = 0.02d0  !linear integration spacing for Fourier transforms
   !dr=0.002d0  !linear integration spacing for Fourier transforms
   dqout = 0.5000001d0  !linear mesh spacing for final E_resid(q) interpolation
   !increment is to avoid exact commnesurability with dq mesh

   al = 0.01d0*dlog(rr(101)/rr(1))
   rc = rr(irc)
   l1 = ll + 1

   nconmx = ncon_in + npr - 1

   allocate (work(mmax))

   !maximum basis size
   nbas = nbas_in + npr - 1
   allocate (qroot(nbas))

   !calculate derivatives of all-electron wave function at r_c
   !for the first projector to set sbf wave vectors for all projectors

   call wf_rc_der(rr, uua(1, 1), al, rc, irc, mmax, uord)
   ulgd = uord(2)/uord(1)

   qmax = 200.0d0  !very large value should always allow nbas q's to be found

   !select wave vectors for spherical Bessel functions
   call qroots(ll, rc, ulgd, nbas, dq, qmax, qroot)

   qmsbf = qroot(nbas)

   !write(6,'(/a)') 'qroots'
   !write(6,'(6f12.6)') (qroot(ii),ii=1,nbas)

   psopt(:, :) = 0.0d0

   ! loop over projectors
   do iprj = 1, npr

      write (6, '(/a,i4)') 'Calculating optimized projector #', iprj, &
      &        ' for l=', ll

      !basis size and number of constraints for this projector
      nbas = nbas_in + iprj - 1
      ncon = ncon_in + iprj - 1
      nnull = nbas - ncon

      allocate (sbasis(nbas))
      allocate (orbasis(nbas, nbas))
      allocate (orbasis_der(ncon, nbas))
      allocate (pswf0_sb(nbas), pswf0_or(nbas))
      allocate (pswfopt_sb(nbas), pswfopt_or(nbas))
      allocate (pswfresid_sb(nbas, nnull))
      allocate (pswfnull_sb(nbas, nnull))
      allocate (pswfnull_or(nbas, nnull))

      !calculate derivatives of all-electron wave function at r_c
      !to define basis set for current projector and derivative constraints

      call wf_rc_der(rr, uua(1, iprj), al, rc, irc, mmax, uord)

      !q considered infinity for E_resid calculation
      qmax = dq*int(2.0d0*qroot(nbas)/dq)
      qmax = max(qmax, 20.0d0)

      nqout = 2 + 0.75d0*qmax/dqout

      allocate (qout(nqout), eresq(nqout), leresq(nqout), eresid0(nqout))
      allocate (eresiddot(nnull, nqout), eresidmat(nnull, nnull, nqout))

      qout(1) = qcut
      do ii = nqout, 2, -1
         qout(ii) = (ii - 2)*dqout
      end do

      !calculate orthogonal basis and constraint matrix
      call sbf_basis_con(ll, rr, mmax, irc, nbas, qroot, psopt, orbasis, orbasis_der, &
      &                     iprj, mxprj, ncon, ncon_in)

      !load constraint vector for value/derivative matching
      cons(:) = 0.0d0
      do jj = 1, ncon_in
         cons(jj) = uord(jj)
      end do
      !load overlap constraints if present
      if (iprj >= 2) then
         do jj = 1, iprj - 1
            cons(ncon_in + jj) = qq(jj, iprj)
         end do
      end if

      !calculate constrained basis for residual energy minimization
      !satisfying off-diagonal norm conservation
      call const_basis(nbas, ncon, cons, orbasis, orbasis_der, &
      &                   pswf0_or, pswfnull_or, &
      &                   pswf0_sb, pswfnull_sb, ps0norm)

      !calculate eigenvectors, eigenvalues, and inhomogeneous terms for
      !the residual energy for a set of q lower cutoffs

      write (6, '(/a,f10.6)') '    Fraction of norm inside rc', qq(iprj, iprj)

      call eresid(ll, irc, nnull, nbas, mmax, rr, dr, dq, qmax, qroot, &
      &                    uua(1, iprj), pswf0_sb, pswfnull_sb, nqout, qout, &
      &                    eresid0, eresiddot, eresidmat)

      write (6, '(a,f7.2,a,f7.2,a)') '    Optimizing pswf for qcut=', qout(1), &
      &   ' a_B^-1,  ecut=', 0.5d0*qout(1)**2, ' Ha'
      write (6, '(a,f6.2,a,f7.1,a)') '    q_infinity defining residual KE=',&
      &   qmax, '  (E_inf=', 0.5d0*qmax**2, ' Ha)'

      !find the null-basis coefficients which minimize the eresid while
      !satisfying diagonal norm conservation

      call optimize(nnull, nbas, pswf0_sb, pswf0_or, nqout, qout,&
      &                eresid0, eresiddot, eresidmat,&
      &                pswfnull_sb, pswfnull_or, qq(iprj, iprj), ps0norm, eresidmin, &
      &                pswfopt_sb, pswfopt_or, ekin_anal, eresq)

      write (6, '(a,1p,d10.2,a)') '    Residual kinetic energy error=', &
      &        eresidmin, ' Ha'

      ! find the Vanderbilt projectors and optimized wave functions

      call pspot(iprj, ll, rr, irc, mmax, al, nbas, qroot, eig(iprj), uua(1, iprj), &
      &             pswfopt_sb, psopt(1, iprj), vae, work, vkb(1, iprj), ekin_num)

      !semi-local potential
      if (iprj == 1) then
         do ii = 1, irc
            vpsp(ii) = work(ii)
         end do
         do ii = irc + 1, mmax
            vpsp(ii) = vae(ii)
         end do
      end if

      write (6, '(/a)') '    Total kinetic energy consistency test'
      write (6, '(a)') '      Fourier integration compared to d^2/dr^2 integral'
      write (6, '(a,f12.8,a,f12.8,a,f12.8)') '      Fourier', ekin_anal, &
      &    '  r-space=', ekin_num, '  ratio=', ekin_anal/ekin_num
      write (6, '(a)') '    Potential consistency test at r_c'
      write (6, '(a,f12.8,a,f12.8,a,1p,d10.2)') '    "vpsp"=', work(irc), '  vae=', &
      &    vae(irc), '  difference=', vae(irc) - work(irc)

      !interpolate the convergence behavior of the optimized pseudo wave function

      write (6, '(/a)') '    Energy error per electron        Cutoff'
      write (6, '(a)') '         Ha          eV             Ha'
      leresq(:) = log(eresq(:))
      err = 0.01d0
      do jj = 1, 7
         lerr = log(err)
         do ii = 3, nqout
            if (leresq(ii - 1) > lerr .and. leresq(ii) <= lerr) then
               qerr = qout(ii) - dqout*(lerr - leresq(ii))/(leresq(ii - 1) - leresq(ii))
               cvgplt(1, jj, iprj) = 0.5d0*qerr**2
               cvgplt(2, jj, iprj) = err
               if (mod(jj, 2) /= 0) then
                  write (6, '(4x,2f12.5,f12.2)') err, err*Ha_eV, 0.5d0*qerr**2
               end if
            end if
         end do
         err = sqrt(0.1d0)*err
      end do

      deallocate (sbasis)
      deallocate (orbasis)
      deallocate (orbasis_der)
      deallocate (pswf0_sb, pswf0_or)
      deallocate (pswfopt_sb, pswfopt_or)
      deallocate (pswfresid_sb)
      deallocate (pswfnull_sb)
      deallocate (pswfnull_or)
      deallocate (eresiddot, eresidmat)
      deallocate (qout, eresq, leresq, eresid0)

   end do !iprj

   deallocate (work, qroot)

   return
end subroutine run_optimize
!> calculates convergence-optimized pseudo-wave-function coefficients
!> in orthonormal and spherical Bessel function bases
subroutine optimize(nnull, nbas, pswf0_sb, pswf0_or, nqout, qout,&
&                    eresid0, eresiddot, eresidmat,&
&                    pswfnull_sb, pswfnull_or, uunorm, ps0norm, eresidmin,&
&                    pswfopt_sb, pswfopt_or, ekin_anal, eresq)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> nnull  number of unconstrained basis vectors for residual minimization
   integer, intent(in) :: nnull
   !> nbas  number of sbf basis functions
   integer, intent(in) :: nbas
   !> nqout  number of q values in [0,qmax] for which results are to be saved
   integer, intent(in) :: nqout
   !> pswf0_sb(nbas)  sbf basis  coefficients for constrained part of pseudo
   !>                     wave function
   real(dp), intent(in) :: pswf0_sb(nbas)
   !> pswfnull_sb(nbas,nnull)  sbf coefficients of null-space eigenfunctions
   !> pswfnull_sb(nbas,nnull)  or basis coefficients of null-space eigenfunctions
   real(dp), intent(in) :: pswfnull_sb(nbas, nnull)
   !> pswf0_or(nbas)  or basis  coefficients
   real(dp), intent(in) :: pswf0_or(nbas)
   real(dp), intent(in) :: pswfnull_or(nbas, nnull)
   !> eresid0(nqout) set of <pswf0| E_resid |pswf0> matrix elements
   real(dp), intent(in) :: eresid0(nqout)
   !> eresiddot(nnull,nqout) set of <pswfnull| E_resid |pswf0> matrix elements
   real(dp), intent(in) :: eresiddot(nnull, nqout)
   !> eresidmat(nnull,nnull,nqout) set of <pswfnull | E_resid |pswfnull'>
   real(dp), intent(in) :: eresidmat(nnull, nnull, nqout)
   !> qout(nqout)  q values in [0,qmax] for which residuals are to be calculated
   real(dp), intent(in) :: qout(nqout)
   !> eresq(nqout)  E_resid for optimized pswf as a function of cutoff
   real(dp), intent(in out) :: eresq(nqout)
   !> uunorm  all-electron charge inside rc
   real(dp), intent(in) :: uunorm
   !> ps0norm  ps0 charge inside rc
   real(dp), intent(in) :: ps0norm

   !Output variables
   real(dp), intent(out) :: eresidmin
   !> ekin_anal  total pswf kinetic energy calculated analytically from E_resid
   real(dp), intent(out) :: ekin_anal
   !> pswfopt_sb optimized pseudowavefunction sbf coefficients
   real(dp), intent(out) :: pswfopt_sb(nbas)
   !> pswfopt_or optimized pseudowavefunction or basis coefficients
   real(dp), intent(out) :: pswfopt_or(nbas)
   real(dp) :: eres(nqout)

   !Local variables
   real(dp) :: yy, tt, emin, rtsgn, sn
   real(dp) :: emin0, emin1, x1, x1min, x1max
   real(dp), allocatable :: work(:), wmat(:, :), wev(:), wvec(:)
   real(dp), allocatable :: eresidevec_nu(:, :)
   real(dp), allocatable :: pswfopt(:), pswfopt_nu(:)
   real(dp), allocatable :: eresiddot_ev(:), eresideval(:)
   real(dp), parameter :: eps = 1.d-10
   integer :: ii, iter, jj, kk, ll, nn, info
   integer, parameter :: niter = 100
   logical :: converged

   !find eigenvalues and eigenvectors of saved E_resid matrix for qout(1)
   !which is meant to be the cutoff used for optimization

   allocate (work(5*nnull), wmat(nnull, nnull), wev(nnull), wvec(nnull))
   allocate (eresidevec_nu(nnull, nnull), eresideval(nnull))
   allocate (eresiddot_ev(nnull))

   !      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

   wmat(:, :) = eresidmat(:, :, 1)

   call dsyev('V', 'U', nnull, wmat, nnull, wev, work, 5*nnull, info)
   if (info .ne. 0) then
      write (6, '(a,i4)') 'optimize: residual-energy eigenvalue ERROR, info=', info
      stop
   end if

   !convert <ps0| E_resid |psnull> vector from the null-space representation
   !to the eigenvector representation

   wvec(:) = eresiddot(:, 1)
   eresiddot_ev(:) = 0.0d0
   do jj = 1, nnull
      do kk = 1, nnull
         eresiddot_ev(jj) = eresiddot_ev(jj) + wvec(kk)*wmat(kk, jj)
      end do
   end do

   !write(6,'(/a,1p,4d14.4)') 'eresideval',wev(:)
   !write(6,'(/a,1p,4d14.4)') 'eresiddot_ev',eresiddot_ev(:)
   !write(6,'(/a,1p,4d14.4)') 'eresid0',eresid0(1)

   eresideval(:) = wev(:)
   eresidevec_nu(:, :) = wmat(:, :)

   deallocate (work, wmat, wev, wvec)

   !find coefficients of E_resid eigenvectors which must be added to pswf0
   !to minimize E_resid

   !note that E_resid is an exteremely simple function in this representation,
   !consisting of a constant plus as sum of linear and quadratic terms in each of
   !these coefficients

   allocate (pswfopt(nnull), pswfopt_nu(nnull))
   pswfopt(:) = 0.d0
   pswfopt_nu(:) = 0.d0

   !write(6,'(/a,1p,3d16.6)') 'initial emin, final, diff',emin0,emin,emin0-emin

   ! Interval-halving search to minimize E*r

   ! This is all based on Eqs. (14) - (16) in my paper as correccted in the
   ! erratum.  It is obvious from Eq.(14) that to minimize E^r, the signs of x_i
   ! must be the opposite of those of f_i.  Eq.(16) is derived by treating
   ! x_1 as depenent upon x_2...x_N-M through Eq.(15).  We can then turn
   ! around and treat those x's as functions all dependent on x_1.  Since the
   ! signs of x_1 and f_1 must be opposite at the minimum, and since e_1
   ! is the smallest eigenvalue, the denominator of Eq.(16) is always positive,
   ! and the magnitudes of the x_i are all monotonically increasing functions
   ! of x_1.  The limits of |x_1| are zero and D_norm (yy here), so a simple
   ! interval-halving search varying |x_1| to satisfy the norm constraint
   ! (Eq.(15)) should be robust.

   ! Translations between the notation in the paper and here are as follows:
   !   pswfopt(ii) = x_i
   !   eresiddot_ev(ii) = f_i
   !   eresideval(ii) = e_i
   !   D_norm == yy
   !   E^r_00 = eresid0(1)
   !   N-M = nnull

   rtsgn = -sign(1.0d0, eresiddot_ev(1))
   if (ps0norm > uunorm) then
      write (6, '(/a)') 'optimize: ERROR ps0norm > uunorm, code will stop'
      stop
   end if
   yy = sqrt((uunorm - ps0norm))

   x1min = 0.0d0
   x1max = yy
   converged = .false.

   do iter = 1, niter

      x1 = 0.5d0*(x1max + x1min)

      tt = x1**2 - yy**2
      do ii = 2, nnull
         pswfopt(ii) = -eresiddot_ev(ii)/(eresideval(ii) - eresideval(1) &
         &               + abs(eresiddot_ev(1))/x1)
         tt = tt + pswfopt(ii)**2
      end do

      if (x1max - x1min < eps) then
         pswfopt(1) = rtsgn*x1
         converged = .true.
         exit
      end if

      if (tt < 0.0d0) then
         x1min = x1
      else
         x1max = x1
      end if
   end do

   if (converged .eqv. .false.) then
      write (6, '(/a)') 'optimize: ERROR interval-halving search not converged'
      stop
   end if

   emin = eresid0(1)
   do ii = 1, nnull
      emin = emin + (2.0d0*eresiddot_ev(ii)*pswfopt(ii) &
      &            + eresideval(ii)*pswfopt(ii)**2)
   end do
   eresidmin = emin

   ekin_anal = 0.0d0
   !find cutoff-dependence of E_resid for optimized pswf
   !must transform pswfopt to the null basis representation to be
   !consistent with the input dot product and
   allocate (work(nnull))
   do ii = 1, nnull
      pswfopt_nu(ii) = 0.0d0
      do jj = 1, nnull
         pswfopt_nu(ii) = pswfopt_nu(ii) + eresidevec_nu(ii, jj)*pswfopt(jj)
      end do
   end do

   do kk = 1, nqout
      emin = eresid0(kk)
      work(:) = 0.0d0
      do ii = 1, nnull
         emin = emin + 2.0d0*eresiddot(ii, kk)*pswfopt_nu(ii)
         work(:) = work(:) + eresidmat(:, ii, kk)*pswfopt_nu(ii)
      end do
      do ii = 1, nnull
         emin = emin + pswfopt_nu(ii)*work(ii)
      end do
      eresq(kk) = emin
      if (qout(kk) == 0.0d0) ekin_anal = emin
   end do

   !calculate optimized pswf in sbf representation
   pswfopt_sb(:) = pswf0_sb(:)
   pswfopt_or(:) = pswf0_or(:)
   do ii = 1, nnull
      pswfopt_sb(:) = pswfopt_sb(:) + pswfopt_nu(ii)*pswfnull_sb(:, ii)
      pswfopt_or(:) = pswfopt_or(:) + pswfopt_nu(ii)*pswfnull_or(:, ii)
   end do

   deallocate (pswfopt, pswfopt_nu)
   return
end subroutine optimize
!> calculates optimized pseudopotential from coefficients of optimized
!> pseudo wave function and all-electron wave function
subroutine pspot(ipr, ll, rr, irc, mmax, al, nbas, qroot, eig, uu, pswfopt_sb, &
&           psopt, vae, vpsp, vkb, ekin_num)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> ipr 1 for first projector, 2 for second which usually has a node
   integer, intent(in) :: ipr
   !> ll  angular momentum
   integer, intent(in) :: ll
   !> irc  index of rc
   integer, intent(in) :: irc
   !> mmax  last point of rr
   integer, intent(in) :: mmax
   !> nbas  number of sbfs
   integer, intent(in) :: nbas
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> qroot  q's for set of sbf
   real(dp), intent(in) :: qroot(nbas)
   !> uu  rr*all-electron wave function
   real(dp), intent(in) :: uu(mmax)
   !> pswfopt_sb  sbf coefficients for optimized pseudo wave function
   real(dp), intent(in) :: pswfopt_sb(nbas)
   !> vae  all-electron potential
   real(dp), intent(in) :: vae(mmax)
   !> al  log of mesh spacing factor
   real(dp), intent(in) :: al
   !> eig  eigenvalue of all-electron wave function
   real(dp), intent(in) :: eig

   !Output variables
   !> psopt  rr * pseudo wave function on log grid
   real(dp), intent(out) :: psopt(mmax)
   !> vpsp  pseudopotential
   real(dp), intent(out) :: vpsp(mmax)
   !> vkb  (epsilon-T)*phi component of chi projector
   real(dp), intent(out) :: vkb(mmax)
   !> ekin_num  total radial kinetic energy
   real(dp), intent(out) :: ekin_num

   !Local variables
   real(dp), allocatable :: work(:, :)
   real(dp) :: sbfder(5)
   real(dp) :: amesh, tt, ddt, ske, ro
   integer :: ii, jj, kk, ll1
   integer :: ixp(4)

   allocate (work(mmax, 5))
   work(:, :) = 0.0d0

   ll1 = ll + 1
   ixp(1) = 2
   ixp(2) = 4
   ixp(3) = 6
   ixp(4) = 8

   do ii = 1, mmax
      work(ii, 1) = uu(ii)
   end do

   !numerical derivatives of all-electron wave function
   !note that we are taking derivatives of uu here, not uu/rr as for sbf matching
   do jj = 2, 3

      ! 5-point numerical first derivatives applied successively
      !  do ii=1+2*jj,mmax-2*jj
      !     work(ii,jj)=(2.d0*work(ii-2,jj-1)-16.d0*work(ii-1,jj-1)&
      !&     +16.d0*work(ii+1,jj-1)-2.d0*work(ii+2,jj-1))&
      !&     /(24.d0*al*rr(ii))

      ! 7-point numerical first derivatives applied successively
      do ii = 1 + 3*jj, mmax - 3*jj
         work(ii, jj) = (-work(ii - 3, jj - 1) + 9.d0*work(ii - 2, jj - 1)&
         &     - 45.d0*work(ii - 1, jj - 1) + 45.d0*work(ii + 1, jj - 1)&
         &     - 9.d0*work(ii + 2, jj - 1) + work(ii + 3, jj - 1))&
         &     /(60.d0*al*rr(ii))
      end do
   end do

   !fill [0, rc] with rr * analytic pswf and 2nd derivatives
   do ii = 1, irc
      tt = 0.0d0
      ddt = 0.0d0
      do jj = 1, nbas
         call sbf_rc_der(ll, qroot(jj), rr(ii), sbfder)
         tt = tt + pswfopt_sb(jj)*rr(ii)*sbfder(1)
         ddt = ddt + pswfopt_sb(jj)*(rr(ii)*sbfder(3) + 2.0d0*sbfder(2))
      end do
      work(ii, 1) = tt
      work(ii, 3) = ddt
   end do

   do ii = 1, mmax
      psopt(ii) = work(ii, 1)
   end do

   !integration for radial kinetic energy

   amesh = exp(al)
   ro = rr(5)/sqrt(amesh)

   !kinetic operator on psopt including centrifugal term
   work(:, 4) = 0.5d0*(ll*ll1*work(:, 1)/rr(:)**2 - work(:, 3))

   !pseudopotential and node test
   do ii = 1, irc
      if (abs(work(ii, 1)) > 0.0d0) then
         if (ipr <= 1 .and. work(ii, 1)*work(ii + 1, 1) < 0.0d0) then
            write (6, '(a)') ' ERROR pspot:  first pseudo wave function has node, &
            &         program will stop'
            write (6, '(a)') ' ERROR pspot: try changing psp parameters for this l'
            stop
         end if
         vpsp(ii) = -work(ii, 4)/work(ii, 1) + eig
         vkb(ii) = -work(ii, 4) + eig*work(ii, 1)
      else
         vpsp(ii) = 0.0d0
      end if
   end do

   do ii = irc + 1, mmax
      vpsp(ii) = vae(ii)
      vkb(ii) = vae(ii)*work(ii, 1)
   end do

   !kinetic energy integrand
   work(:, 5) = work(:, 1)*work(:, 4)

   ro = rr(5)/sqrt(amesh)
   ske = ((work(5, 5))/rr(5)**ixp(ll1))*ro**ixp(ll1)/ixp(ll1)

   do ii = 5, mmax - 9
      ske = ske + al*rr(ii)*work(ii, 5)
   end do

   ekin_num = ske

   deallocate (work)
   return
end subroutine pspot
!> calculates basis for residual energy minimization
subroutine const_basis(nbas, ncon, cons, orbasis, orbasis_der,&
&                       pswf0_or, pswfnull_or, &
&                       pswf0_sb, pswfnull_sb, ps0norm)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !input variables
   !> ncon  number of constraints
   integer, intent(in) :: ncon
   !> nbas  number of orthogonalized sbf basis functions
   integer, intent(in) :: nbas
   !> orbasis  coefficients of sbfs in orthogonal basis
   real(dp), intent(in) :: orbasis(nbas, nbas)
   !> orbasis_der  values and derivatives of orthogonalized sbf basis
   !>               functions at rc (constraint matrix)
   real(dp), intent(in) :: orbasis_der(ncon, nbas)
   !> cons  constraints (value at rc and derivatives up to 4th)
   real(dp), intent(in) :: cons(ncon)

   !Output variables
   real(dp), intent(out) :: pswf0_or(nbas)
   real(dp), intent(out) :: pswfnull_or(nbas, nbas - ncon)
   !> pswf0_sb  linear combination of sbfs that matches cons at rc
   real(dp), intent(out) :: pswf0_sb(nbas)
   !> pswfnull_sb  nbas-ncon linear combinations of sbfs satisfying value and ncon-1
   !>            derivatives at rc == 0 (null space of constraint matrix)
   !>           These form an orthonormal set, and are all orthogonal to
   !>           the pswf0 combination (which is not normalized)
   real(dp), intent(out) :: pswfnull_sb(nbas, nbas - ncon)
   !> ps0norm  charge contained in ps0 inside rc
   real(dp), intent(out) :: ps0norm

   !Local variables
   real(dp) :: usvd(10, 10)
   real(dp) :: ss(10)
   real(dp) :: work(100)
   real(dp) :: con_tst(10)
   real(dp), allocatable :: amat(:, :)
   real(dp), allocatable :: vt(:, :)
   real(dp) :: sn
   real(dp), parameter :: eps = 1.d-8
   integer :: ii
   integer :: jj
   integer :: kk
   integer :: info

   !DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
   !     $                   WORK, LWORK, INFO )
   allocate (amat(10, nbas), vt(nbas, nbas))

   !save orbasis_der from overwrite
   amat(:, :) = 0.0d0
   amat(1:ncon, :) = orbasis_der(1:ncon, :)

   !singular value decomposition of constraint matrix
   call dgesvd('A', 'A', ncon, nbas, amat, 10, ss, usvd, 10, vt, nbas, &
   &             work, 64, info)

   if (info /= 0) then
      write (6, '(/a,i4)') 'const_basis: ERTOR const_basis: dgesvd info = ', info
      stop
   end if

   !write(6,'(/a,6d12.4)') 'singular values',(ss(ii),ii=1,ncon)

   !write(6,'(/a)') 'usvd'
   !do jj=1,ncon
   ! write(6,'(6f12.6)') (usvd(jj,ii),ii=1,ncon)
   !end do

   !write(6,'(/a)') 'vt'
   !do jj=1,nbas
   ! write(6,'(11f12.6)') (vt(jj,ii),ii=1,nbas)
   !end do

   !operate on constraint vector with usvd^transpose
   work(:) = 0.0d0
   do ii = 1, ncon
      work(ii) = 0.0d0
      do jj = 1, ncon
         work(ii) = work(ii) + usvd(jj, ii)*cons(jj)
      end do
   end do

   !scale work with inverse transpose of singular value "matrix" ss
   do ii = 1, ncon
      if (abs(ss(ii)) > eps) then
         work(ii) = work(ii)/ss(ii)
      end if
   end do

   !form coefficient vector for basic matching pseudo wave function
   !null vectors (with ss=0) can be added to improve kinetic residual
   !convergence without changing match at rc
   do ii = 1, nbas
      pswf0_or(ii) = 0.0d0
      do jj = 1, ncon
         pswf0_or(ii) = pswf0_or(ii) + work(jj)*vt(jj, ii)
      end do
   end do

   !copy coefficient vectors spanning null space of constraint matrix
   !note that vt (transpose) rows become columns of pswfnull
   !write(6,'(/)')
   do jj = 1, nbas - ncon
      sn = 0.0d0
      do ii = 1, nbas
         pswfnull_or(ii, jj) = vt(ncon + jj, ii)
         sn = sn + pswfnull_or(ii, jj)**2
      end do
   end do

   !now, test result
   do ii = 1, ncon
      con_tst(ii) = 0.0d0
      do jj = 1, nbas
         con_tst(ii) = con_tst(ii) + orbasis_der(ii, jj)*pswf0_or(jj)
      end do
   end do

   write (6, '(/a)') '    Constraint consistency test'
   write (6, '(a)') '      pswf0 val/der aewf val/der     diff'
   do ii = 1, ncon
      write (6, '(4x,2f14.8,1p,d12.2)') con_tst(ii), cons(ii), con_tst(ii) - cons(ii)
   end do

   !calculate charge contained in ps0 inside rc
   ps0norm = 0.0d0
   do ii = 1, nbas
      ps0norm = ps0norm + pswf0_or(ii)**2
   end do

   !convert ps0 and psnull in orthogonal representation to sbf represention

   pswf0_sb(:) = 0.0d0
   pswfnull_sb(:, :) = 0.0d0
   do jj = 1, nbas
      do kk = 1, nbas
         pswf0_sb(jj) = pswf0_sb(jj) + orbasis(jj, kk)*pswf0_or(kk)
         do ii = 1, nbas - ncon
            pswfnull_sb(jj, ii) = pswfnull_sb(jj, ii) + orbasis(jj, kk)*pswfnull_or(kk, ii)
         end do
      end do
   end do

   deallocate (amat, vt)
   return
end subroutine const_basis
end module vsl_m

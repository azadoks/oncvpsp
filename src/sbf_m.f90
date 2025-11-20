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
module sbf_m
   implicit none
   private
   public :: sbf8
   public :: sbf_rc_der
   public :: sbf_basis
   public :: sbf_basis_con
   public :: qroots
contains
!> calculates spherical Bessel functions using recursion algorithm
subroutine sbf8(nm, xx, sb_out)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
   real(dp), parameter :: pi2 = 2.0_dp*pi

   !Input variables
   !> nm  maximum angular momentum wanted plus 1
   integer, intent(in) :: nm
   !> xx  argument of spherical bessel function
   real(dp), intent(in) :: xx

   !Output variables
   !> sb_out  output of sbf_l(xx) for l=0,...,nm-1
   real(dp), intent(out) :: sb_out(nm)

   !Local variables
   integer :: nlim, nn
   real(dp) :: cc, fn, sn, ss, xi, xn, xs, xmod
   real(dp), allocatable :: sb(:)

   if (xx <= 1.0d-36) then
      !  zero argument section
      sb_out(:) = 0.0d0
      sb_out(1) = 1.0d0
   else if (xx < 1.0d-3) then
      !  small argument section
      xn = 1.0d0
      xs = 0.5d0*xx**2
      do nn = 1, nm
         sb_out(nn) = xn*(1.0d0 - xs*(1.0d0 - xs/(4*nn + 6))/(2*nn + 1))
         xn = xx*xn/(2*nn + 1)
      end do
      !  trigonometric section (small l values)
   else if (nm == 1) then
      xmod = mod(xx, pi2)
      ss = sin(xmod)
      sb_out(1) = ss/xx
   else if (nm == 2) then
      xmod = mod(xx, pi2)
      ss = sin(xmod)
      cc = cos(xmod)
      sb_out(1) = ss/xx
      sb_out(2) = (ss - xx*cc)/xx**2
   else if (nm == 3) then
      xmod = mod(xx, pi2)
      ss = sin(xmod)
      cc = cos(xmod)
      sb_out(1) = ss/xx
      sb_out(3) = ((3.0d0 - xx**2)*ss - 3.0d0*xx*cc)/xx**3
   else if (nm == 4) then
      xmod = mod(xx, pi2)
      ss = sin(xmod)
      cc = cos(xmod)
      sb_out(1) = ss/xx
      sb_out(3) = ((3.0d0 - xx**2)*ss - 3.0d0*xx*cc)/xx**3
      sb_out(4) = (15.d0*ss - 15.d0*xx*cc - 6.d0*xx**2*ss + xx**3*cc)/xx**4
   else
      !  recursion method (accurate for large l, slow for large arguments)
      if (xx < 1.0d0) then
         nlim = nm + int(15.0d0*xx) + 1
      else
         nlim = nm + int(1.36d0*xx) + 15
      end if
      allocate (sb(nlim + 1))
      nn = nlim
      xi = 1.0d0/xx
      sb(nn + 1) = 0.0d0
      sb(nn) = 1.0d-18
      sn = dble(2*nn - 1)*1.0d-36
      do nn = nlim - 1, 1, -1
         sb(nn) = dble(2*nn + 1)*xi*sb(nn + 1) - sb(nn + 2)
      end do
      do nn = 1, nlim - 1
         sn = sn + dble(2*nn - 1)*sb(nn)*sb(nn)
      end do
      fn = 1.0d0/sqrt(sn)
      sb_out(:) = fn*sb(1:nm)
      deallocate (sb)
   end if
   return
end subroutine sbf8
!> calculate spherical Bessel function values and first four derivatives
subroutine sbf_rc_der(llin, qq, rc, sbfder)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> llin  angular momentum l
   integer, intent(in) :: llin
   !> qq  wave vector
   real(dp), intent(in) :: qq
   !> rc  core radius
   real(dp), intent(in) :: rc

   !Output variables
   !> sbfder  values and first 4 derivatives of j_l(qq*rc)
   real(dp), intent(out) :: sbfder(5)

   !Local variables
   real(dp) :: xx
   real(dp) :: sb_out(10), sbfad(5, 10)
   integer :: ii, ll

   if (llin > 5) then
      write (6, '(a,i4,a)') 'sbf_rc_der: argument ERROR, llin = ', llin,&
      &  ' >5'
      stop
   end if

   ! calculate sbf for needed l values
   xx = qq*rc
   call sbf8(llin + 5, xx, sb_out)

   ! derivatives based on derivative formula applied recursively
   !
   ! (1/z d/dz)^m [z^-n j_n(z)]= (-1)^m z^(-n-m) j_(n+m)(z)
   !

   ! store values
   do ll = llin, llin + 5
      sbfad(1, ll + 1) = sb_out(ll + 1)
   end do

   ! first derivatives
   do ll = llin, llin + 3
      sbfad(2, ll + 1) = ll*sbfad(1, ll + 1)/xx &
      &                - sbfad(1, ll + 2)
   end do

   ! second derivatives
   do ll = llin, llin + 2
      sbfad(3, ll + 1) = -ll*sbfad(1, ll + 1)/xx**2 &
      &              + ll*sbfad(2, ll + 1)/xx&
      &                  - sbfad(2, ll + 2)
   end do

   ! third derivatives
   do ll = llin, llin + 1
      sbfad(4, ll + 1) = 2*ll*sbfad(1, ll + 1)/xx**3 &
      &               - 2*ll*sbfad(2, ll + 1)/xx**2 &
      &                 + ll*sbfad(3, ll + 1)/xx &
      &                    - sbfad(3, ll + 2)
   end do

   ! fourth derivatives
   do ll = llin, llin
      sbfad(5, ll + 1) = -6*ll*sbfad(1, ll + 1)/xx**4 &
      &               + 6*ll*sbfad(2, ll + 1)/xx**3 &
      &               - 3*ll*sbfad(3, ll + 1)/xx**2 &
      &                 + ll*sbfad(4, ll + 1)/xx &
      &                    - sbfad(4, ll + 2)
   end do

   do ii = 1, 5
      sbfder(ii) = sbfad(ii, llin + 1)*qq**(ii - 1)
   end do

   return
end subroutine sbf_rc_der
!> orthonormalize basis functions and derivatives at rc and find all-electron
!> charge inside rc.
subroutine sbf_basis(ll, rr, mmax, irc, nbas, qroot, sbasis, orbasis, orbasis_der, &
&                     nconmx)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> ll  angujlar momentum
   integer, intent(in) :: ll
   !> mmax  number of points in log radial mesh
   integer, intent(in) :: mmax
   integer, intent(in) :: nconmx
   !> irc  index rr such that rr(irc)=rc
   integer, intent(in) :: irc
   !> nbas  number of basis functions
   integer, intent(in) :: nbas
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> qroot  q values for j_l(q*r)
   real(dp), intent(in) :: qroot(nbas)

   !Output variables
   !> sbasis  normalization coefficients for j_l basis set
   real(dp), intent(out) :: sbasis(nbas)
   !> orbasis  matrix for orthonormal basis (rows j_l, columns basis vectors)
   real(dp), intent(out) :: orbasis(nbas, nbas)
   !> orbasis_der  values and derivatives of orthonormal basis set at rc,
   real(dp), intent(out) :: orbasis_der(nconmx, nbas)

   !Local variables
   integer :: ii, jj, ll1, ibas, info
   real(dp) :: al, amesh, ro, rc, sn, xx, tt
   real(dp) :: sb_out(10), sbfder(5), tder(6)
   real(dp), allocatable :: sbfar(:, :), sev(:), work(:), sovlp(:, :)
   real(dp), allocatable :: sovlp_save(:, :), orbasis_ke(:), tor(:)
   logical :: sorted

   ll1 = ll + 1
   rc = rr(irc)
   al = 0.01d0*log(rr(101)/rr(1))
   amesh = exp(al)

   allocate (sbfar(irc, nbas), sev(nbas), work(5*nbas), sovlp(nbas, nbas))
   allocate (sovlp_save(nbas, nbas), orbasis_ke(nbas), tor(nbas))

   do ibas = 1, nbas
      do jj = 1, irc
         xx = qroot(ibas)*rr(jj)
         call sbf8(ll1, xx, sb_out)
         sbfar(jj, ibas) = rr(jj)*sb_out(ll1)
      end do
   end do

   !perform sbf orthonormalization sum for overlap matrix

   sovlp(:, :) = 0.0d0
   ro = rr(1)/sqrt(amesh)
   do ibas = 1, nbas
      do jj = ibas, nbas

         sn = (sbfar(1, ibas)/rr(1)**ll)*(sbfar(1, jj)/rr(1)**ll)&
         &     *ro**(2*ll + 3)/dfloat(2*ll + 3)

         do ii = 1, irc - 3
            sn = sn + al*rr(ii)*sbfar(ii, ibas)*sbfar(ii, jj)
         end do

         sn = sn + al*(23.0d0*rr(irc - 2)*sbfar(irc - 2, ibas)*sbfar(irc - 2, jj)&
         &            + 28.0d0*rr(irc - 1)*sbfar(irc - 1, ibas)*sbfar(irc - 1, jj)&
         &            + 9.0d0*rr(irc)*sbfar(irc, ibas)*sbfar(irc, jj))/24.0d0

         sovlp(ibas, jj) = sn
         sovlp(jj, ibas) = sn
      end do
   end do

   !normalization for j_l basis
   do ibas = 1, nbas
      sbasis(ibas) = 1.0d0/sqrt(sovlp(ibas, ibas))
   end do

   sovlp_save(:, :) = sovlp(:, :)
   !find eigenvalues and eigenvectors of the overlap matrix

   !      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

   call dsyev('V', 'U', nbas, sovlp, nbas, sev, work, 5*nbas, info)
   if (info .ne. 0) then
      write (6, '(a,i4)') 'sbf_basis: ERROR overlap matrix eigenvalue error, info=', info
      stop
   end if

   !write(6,'(/a)') 'sbf overlap matrix eigenvalues'
   !write(6,'(1p,6d12.3)') (sev(ibas),ibas=1,nbas)

   !scale eigenvectors to form orthonormal basis coefficients for sbf's
   !note that we are reversing order so that the leading eigenvector is the
   !most linearly independent
   do ibas = 1, nbas
      if (sev(ibas) > 0.0d0) then
         tt = 1.0d0/sqrt(sev(ibas))
      else
         write (6, '(a,f12.6)') 'sbfbasis: ERROR negative eigenvalue of overlap matrix'
         stop
      end if
      do jj = 1, nbas
         orbasis(jj, nbas - ibas + 1) = tt*sovlp(jj, ibas)
      end do
   end do

   !find rc derivatives of basis funtction
   orbasis_der(:, :) = 0.0d0
   do jj = 1, nbas !sbf loop
      call sbf_rc_der(ll, qroot(jj), rc, sbfder)
      do ibas = 1, nbas !orbasis loop
         do ii = 1, nconmx
            orbasis_der(ii, ibas) = orbasis_der(ii, ibas) + orbasis(jj, ibas)*sbfder(ii)
         end do
      end do
   end do

   !find approximate kinetic energy of orbasis

   orbasis_ke(:) = 0.0d0
   do ibas = 1, nbas
      do jj = 1, nbas
         orbasis_ke(ibas) = orbasis_ke(ibas) + (orbasis(jj, ibas)*qroot(jj))**2
      end do
   end do
   !do  ibas=1,nbas
   !  write(6,*) 'orbasis_ke',ibas,orbasis_ke(ibas)
   !end do

   ! bubble-sort on approximate kinetic energies
   ! (Yes, I know bubble-sort is the least-efficient sorting algorithm.)

   do ii = 1, 100
      sorted = .true.
      do jj = 2, nbas
         if (orbasis_ke(jj - 1) > orbasis_ke(jj)) then

            tt = orbasis_ke(jj)
            tor(:) = orbasis(:, jj)
            tder(:) = orbasis_der(:, jj)

            orbasis_ke(jj) = orbasis_ke(jj - 1)
            orbasis(:, jj) = orbasis(:, jj - 1)
            orbasis_der(:, jj) = orbasis_der(:, jj - 1)

            orbasis_ke(jj - 1) = tt
            orbasis(:, jj - 1) = tor(:)
            orbasis_der(:, jj - 1) = tder(:)

            sorted = .false.
         end if
      end do
      if (sorted) exit
   end do

   !do  ibas=1,nbas
   !  write(6,*) 'sorbasis_ke',ibas,orbasis_ke(ibas)
   !end do

   deallocate (sbfar, sev, work, sovlp)
   deallocate (sovlp_save, orbasis_ke, tor)

   return
end subroutine sbf_basis
!> orthonormalize basis functions and derivatives at rc and form constraint
!> matrix based on derivative to be matched and overlaps with prior
!> optimized wave functions
subroutine sbf_basis_con(ll, rr, mmax, irc, nbas, qroot, psopt, orbasis, orbasis_der, &
&                     iprj, mxprj, ncon, ncon_in)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> ll  angujlar momentum
   integer, intent(in) :: ll
   !> mmax  number of points in log radial mesh
   integer, intent(in) :: mmax
   integer, intent(in) :: ncon
   integer, intent(in) :: ncon_in
   !> irc  index rr such that rr(irc)=rc
   integer, intent(in) :: irc
   !> nbas  number of basis functions
   integer, intent(in) :: nbas
   integer, intent(in) :: iprj
   integer, intent(in) :: mxprj
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> qroot  q values for j_l(q*r)
   real(dp), intent(in) :: qroot(nbas)
   !> psopt  optimized projector wave functions for lower iprj
   real(dp), intent(in) :: psopt(mmax, mxprj)

   !Output variables
   !> orbasis  matrix for orthonormal basis (rows j_l, columns basis vectors)
   real(dp), intent(out) :: orbasis(nbas, nbas)
   !> orbasis_der  values and derivatives of orthonormal basis set at rc, and
   !>  norm constraint vectors from already-computed psopt
   real(dp), intent(out) :: orbasis_der(ncon, nbas)

   !Local variables
   integer :: ii, jj, kk, ll1, ibas, info
   real(dp) :: al, amesh, ro, rc, sn, xx, tt
   real(dp) :: sb_out(10), sbfder(5), tder(6)
   real(dp), allocatable :: sbfar(:, :), sev(:), work(:), sovlp(:, :)
   real(dp), allocatable :: sbf_or(:, :)

   ll1 = ll + 1
   rc = rr(irc)
   al = 0.01d0*log(rr(101)/rr(1))
   amesh = exp(al)

   allocate (sbfar(irc, nbas), sev(nbas), work(5*nbas), sovlp(nbas, nbas))
   allocate (sbf_or(irc, nbas))

   do ibas = 1, nbas
      do jj = 1, irc
         xx = qroot(ibas)*rr(jj)
         call sbf8(ll1, xx, sb_out)
         sbfar(jj, ibas) = rr(jj)*sb_out(ll1)
      end do
   end do

   !perform sbf orthonormalization sum for overlap matrix

   sovlp(:, :) = 0.0d0
   ro = rr(1)/sqrt(amesh)
   do ibas = 1, nbas
      do jj = ibas, nbas

         sn = (sbfar(1, ibas)/rr(1)**ll)*(sbfar(1, jj)/rr(1)**ll)&
         &     *ro**(2*ll + 3)/dfloat(2*ll + 3)

         do ii = 1, irc - 3
            sn = sn + al*rr(ii)*sbfar(ii, ibas)*sbfar(ii, jj)
         end do

         sn = sn + al*(23.0d0*rr(irc - 2)*sbfar(irc - 2, ibas)*sbfar(irc - 2, jj)&
         &            + 28.0d0*rr(irc - 1)*sbfar(irc - 1, ibas)*sbfar(irc - 1, jj)&
         &            + 9.0d0*rr(irc)*sbfar(irc, ibas)*sbfar(irc, jj))/24.0d0

         sovlp(ibas, jj) = sn
         sovlp(jj, ibas) = sn
      end do
   end do

   !find eigenvalues and eigenvectors of the overlap matrix

   !      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

   call dsyev('V', 'U', nbas, sovlp, nbas, sev, work, 5*nbas, info)
   if (info .ne. 0) then
      write (6, '(a,i4)') 'sbf_basis: overlap matrix eigenvalue error, info=', info
      stop
   end if

   !write(6,'(/a)') 'sbf overlap matrix eigenvalues'
   !write(6,'(1p,6d12.3)') (sev(ibas),ibas=1,nbas)

   !scale eigenvectors to form orthonormal basis coefficients for sbf's
   !note that we are reversing order so that the leading eigenvector is the
   !most linearly independent
   do ibas = 1, nbas
      if (sev(ibas) > 0.0d0) then
         tt = 1.0d0/sqrt(sev(ibas))
      else
         write (6, '(a,f12.6)') 'sbfbasis: negative eigenvalue of overlap matrix'
         stop
      end if
      do jj = 1, nbas
         orbasis(jj, nbas - ibas + 1) = tt*sovlp(jj, ibas)
      end do
   end do

   !find rc derivatives of basis funtction
   orbasis_der(:, :) = 0.0d0
   do jj = 1, nbas !sbf loop
      call sbf_rc_der(ll, qroot(jj), rc, sbfder)
      do ibas = 1, nbas !orbasis loop
         do ii = 1, ncon_in
            orbasis_der(ii, ibas) = orbasis_der(ii, ibas) + orbasis(jj, ibas)*sbfder(ii)
         end do
      end do
   end do

   !find orthogonal basis radial functiona
   sbf_or(:, :) = 0.0d0
   do jj = 1, nbas !sbf loop
      do ibas = 1, nbas !orbasis loop
         sbf_or(:, ibas) = sbf_or(:, ibas) + orbasis(jj, ibas)*sbfar(:, jj)
      end do
   end do

   !find new or-basis representation of previous optimized wave functions
   !fill in last rows of orbssis_der constraint matrix for overlap constraint

   if (iprj >= 2) then
      ro = rr(1)/sqrt(amesh)
      do kk = 1, iprj - 1
         do jj = 1, nbas

            sn = (psopt(1, kk)/rr(1)**ll)*(sbf_or(1, jj)/rr(1)**ll)&
            &     *ro**(2*ll + 3)/dfloat(2*ll + 3)

            do ii = 1, irc - 3
               sn = sn + al*rr(ii)*psopt(ii, kk)*sbf_or(ii, jj)
            end do

            sn = sn + al*(23.0d0*rr(irc - 2)*psopt(irc - 2, kk)*sbf_or(irc - 2, jj)&
            &            + 28.0d0*rr(irc - 1)*psopt(irc - 1, kk)*sbf_or(irc - 1, jj)&
            &            + 9.0d0*rr(irc)*psopt(irc, kk)*sbf_or(irc, jj))/24.0d0

            orbasis_der(ncon_in + kk, jj) = sn
         end do !jj
      end do !kk
   end if !iprj>=2

   deallocate (sbfar, sbf_or, sev, work, sovlp)

   return
end subroutine sbf_basis_con
!> calculate q values for spherical Bessel function orthogonal basis inside r_c
!> all having specified log derivative ulgd at rc
!> 1st and 3rd output roots do not correspond to these conditions, rather
!> 0.5 times the first matching root, and the average of the 1st and second
!> matching roots
subroutine qroots(ll, rc, ulgd, nroot, dltq, qmax, qroot)
   implicit none
   integer, parameter :: dp = kind(1.0d0)
   !Input variables
   !> ll  angular momentum
   integer, intent(in) :: ll
   !> nroot  number of spherical Bessel functions matgching ulgd at rc (form
   !>         orthogonal basis)
   integer, intent(in) :: nroot
   !> rc  core radius
   real(dp), intent(in) :: rc
   !> ulgd  log derivaive of all-electron radial wave function at rc
   real(dp), intent(in) :: ulgd
   !> dltq  step for search for qroots satisfying log derivative condition
   real(dp), intent(in) :: dltq
   !> qmax  maximum q for search
   real(dp), intent(in) :: qmax

   !Output variables
   !> qroot  set of q's satisfying dj_l(q*r)/dr / jl(q*r) = ulgd
   real(dp), intent(out) :: qroot(nroot)

   !Local variables
   real(dp), parameter :: eps = 1.0d-12

   real(dp) :: sbfd(5)
   real(dp) al, qq, dlgd, dlgd_last, qhi, qlow, qt
   integer :: ii, jj, ll1, lmax, mmax, iroot, nq
   logical :: found_root

   nq = int(qmax/dltq) + 1

   qroot(:) = 0.0d0
   iroot = 2
   do ii = 1, nq
      qq = dltq*ii
      call sbf_rc_der(ll, qq, rc, sbfd)
      dlgd = ulgd - sbfd(2)/sbfd(1)
      found_root = .false.
      !interval halving to find root
      if (ii > 1 .and. dlgd*dlgd_last < 0.0d0 .and. &
      &         abs(dlgd*dlgd_last) < 1.0d0) then
         if (dlgd > 0.0d0) then
            qhi = qq
            qlow = dltq*(ii - 1)
         else
            qhi = dltq*(ii - 1)
            qlow = qq
         end if

         do jj = 1, 100
            qt = 0.5d0*(qhi + qlow)
            call sbf_rc_der(ll, qt, rc, sbfd)
            dlgd = ulgd - sbfd(2)/sbfd(1)
            if (abs(dlgd) < eps) then
               iroot = iroot + 1
               qroot(iroot) = qt
               found_root = .true.
               exit
            end if
            if (dlgd > 0.0d0) then
               qhi = qt
            else
               qlow = qt
            end if
         end do
         found_root = .true.
         if (.not. found_root) then
            write (6, '(a)') 'qroots: ERROR - failed to find root'
            stop
         end if
      end if
      if (iroot .eq. nroot) exit
      dlgd_last = dlgd
   end do
   if (.not. found_root) then
      write (6, '(a)') 'qroots: ERROR - failed to find nroot roots'
      stop
   end if

   !extra q values for needed flexibility to satisfy constraints
   qroot(1) = 0.5d0*qroot(3)
   qt = 0.5d0*(qroot(3) + qroot(4))
   qroot(2) = qroot(3)
   qroot(3) = qt

   return
end subroutine qroots
end module sbf_m

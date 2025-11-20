
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
module vkb_m
   use vpinteg_m, only: vpinteg
   use vploc_m, only: vploc
   implicit none
   private
   public :: run_vkb
   public :: run_vkb_r
   public :: sr_so_r
contains
!> computes Vanderbilt / Kleinman-Bylander non-local potentials
subroutine run_vkb(lmax, lloc, lpopt, dvloc0, irc, nproj, rr, mmax, mxprj, pswf, &
&                   vfull, vp, evkb, vkb, nlim, vr)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> lmax  maximum angular momentum
   integer, intent(in) :: lmax
   !> lloc  l for local potential (lloc==4 => use linear combination)
   integer, intent(in) :: lloc
   !> lpopt  choice of polynomial for lloc==4
   integer, intent(in) :: lpopt
   !> mmax  size of radial grid
   integer, intent(in) :: mmax
   !> mxprj  dimension of number of projectors
   integer, intent(in) :: mxprj
   !> irc  core radii indices
   integer, intent(in) :: irc(6)
   !> nproj  number of projectors for each l
   integer, intent(in) :: nproj(6)
   !> rr  log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> pswf pseudo wave functions
   real(dp), intent(in) :: pswf(mmax, mxprj, 5)
   !> vfull  all-electron potential
   real(dp), intent(in) :: vfull(mmax)
   !> vp  semi-local pseudopotentials for first projectors
   !>      (vp(:,5) is local potential if lloc=4)
   real(dp), intent(in out) :: vp(mmax, 5)
   !> dvloc0  amplitude at rr==0 to be smoothely added for lloc==4
   real(dp), intent(in) :: dvloc0

   !Input/Output variables
   !> vkb  semi-local "potentials"*pswf (input); VKB projectors (output)
   real(dp), intent(out) :: vkb(mmax, mxprj, 4)
   !> evkb  coefficients of BKB projectors (output)
   real(dp), intent(out) :: evkb(mxprj, 4)
   !> vr  effective scalar-relativisic "potential" calculated in vrel
   real(dp), intent(out) :: vr(mmax, mxprj, 6)
   !> nlim  index of maximum rc including that of vlocal (output)
   integer, intent(out) :: nlim

   !Local variables
   integer :: ii, ipk, jj, kk, ierr, l1, info, np
   integer :: lwork
   real(dp) :: apk, sn, tt, qq12
   real(dp) :: xx, ff
   real(dp) :: bb(mxprj, mxprj), bbev(mxprj), bbi(mxprj, mxprj)
   real(dp) :: sovl(mxprj, mxprj), sovlev(mxprj), smhalf(mxprj, mxprj)
   real(dp) :: sphalf(mxprj, mxprj)
   real(dp) :: bbist(mxprj, mxprj), bbistev(mxprj), bbit(mxprj, mxprj)
   real(dp), allocatable :: vloc(:), vkbt(:, :), vkbst(:, :), work(:)
   logical :: sorted

   evkb(:, :) = 0.0d0

   lwork = 3*mxprj
   allocate (vloc(mmax), vkbt(mmax, mxprj), vkbst(mmax, mxprj), work(lwork))

   ! calculate local potential

   if (lloc == 4) then
      call vploc(rr, vfull, vp, dvloc0, irc(5), mmax, lpopt)
   end if

   vloc(:) = vp(:, lloc + 1)

   ! Vanderbilt / Kleinman-Bylander projector construction
   nlim = irc(lloc + 1)
   do l1 = 1, lmax + 1
      nlim = max(nlim, irc(l1))
      if (l1 == lloc + 1) cycle
      do jj = 1, nproj(l1)
         !incorporates scalar-relativistic "potential" in last 5% of [0:rc] to smoothly
         !force projectors to zero
         do ii = 1, irc(l1)
            xx = 20.0d0*((rr(ii)/rr(irc(l1))) - 1.0d0)
            if (xx > -1.0d0) then
               ff = (1.0d0 - xx**2)**2
            else
               ff = 0.0d0
            end if
            vkb(ii, jj, l1) = vkb(ii, jj, l1) - (vloc(ii) + ff*vr(ii, jj, l1))*pswf(ii, jj, l1)
         end do
         do ii = irc(l1) + 1, mmax
            vkb(ii, jj, l1) = 0.0d0
         end do
      end do
   end do

   ! Vanderbilt B-matrix construction

   do l1 = 1, lmax + 1
      if (l1 == lloc + 1) cycle
      np = nproj(l1)
      bb(:, :) = 0.0d0
      call vpinteg(pswf(1, 1, l1), vkb(1, 1, l1), irc(l1), 2*l1, bb(1, 1), rr)
      evkb(1, l1) = 1.0d0/bb(1, 1)
      ! this is Kleinman-Bylander result if nproj==1

      if (nproj(l1) >= 2) then

         do jj = 1, np
            do ii = 1, np
               call vpinteg(pswf(1, ii, l1), vkb(1, jj, l1), irc(l1), 2*l1, bb(ii, jj), rr)
            end do
         end do

         ! symmetrize exactly
         write (6, '(/a,i4)') 'B matrix Hermiticity error, ll=', l1 - 1
         do jj = 2, np
            do ii = 1, jj - 1
               write (6, '(2i4,1p,e14.4)') ii, jj, bb(ii, jj) - bb(jj, ii)
               tt = 0.5d0*(bb(ii, jj) + bb(jj, ii))
               bb(ii, jj) = tt
               bb(jj, ii) = tt
            end do
         end do

         ! find eigenvalues and eigenvectors of the B matrix

         !      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

         call dsyev('V', 'U', np, bb, mxprj, bbev, work, lwork, info)
         if (info .ne. 0) then
            write (6, '(a,i4)') 'run_vkb: B matrix eigenvalue ERROR, info=', info
            stop
         end if

         ! take linear combinations to form diagonal projectors

         do jj = 1, np
            vkbt(:, jj) = vkb(:, jj, l1)
         end do

         do jj = 1, np
            vkb(:, jj, l1) = 0.0d0
            evkb(jj, l1) = 1.0d0/bbev(jj)
            do ii = 1, np
               vkb(:, jj, l1) = vkb(:, jj, l1) + bb(ii, jj)*vkbt(:, ii)
            end do
         end do

      end if !nproj>=2

      ! normalize projectors
      do jj = 1, np
         call vpinteg(vkb(1, jj, l1), vkb(1, jj, l1), irc(l1), 2*l1, sn, rr)
         if (vkb(irc(l1) - 10, jj, l1) >= 0.0d0) then
            vkb(:, jj, l1) = vkb(:, jj, l1)/sqrt(sn)
         else
            vkb(:, jj, l1) = -vkb(:, jj, l1)/sqrt(sn)
         end if
         evkb(jj, l1) = sn*evkb(jj, l1)
      end do

      ! calculate overlap
      if (nproj(l1) >= 2) then

         bbi(:, :) = 0.0d0
         do jj = 1, np
            sovl(jj, jj) = 1.0d0
            bbi(jj, jj) = evkb(jj, l1)
            do ii = jj + 1, np
               call vpinteg(vkb(1, ii, l1), vkb(1, jj, l1), irc(l1), 2*l1, sn, rr)
               sovl(ii, jj) = sn
               sovl(jj, ii) = sn
            end do
         end do
         call vpinteg(vkb(1, 1, l1), vkb(1, 2, l1), irc(l1), 2*l1, sn, rr)

         ! find eigenvalues and eigenvectors of the overlap matrix

         !      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

         call dsyev('V', 'U', np, sovl, mxprj, sovlev, work, lwork, info)
         if (info .ne. 0) then
            write (6, '(a,i4)') 'run_vkb: sovl matrix ERROR error, info=', info
            stop
         end if

         !   write(6,'(/1p,2e12.4)') (sovlev(ii),ii=1,2)

         ! construct S^(-1/2) AND s^(1/2)

         do jj = 1, np
            tt = sqrt(sovlev(jj))
            do ii = 1, np
               sphalf(ii, jj) = tt*sovl(ii, jj)
               smhalf(ii, jj) = sovl(ii, jj)/tt
            end do
         end do

         ! construct B^(-1)* = S^(1/2)^T B^(-1) S^(1/2)

         bbit(:, :) = 0.0d0
         bbist(:, :) = 0.0d0

         do ii = 1, np
            do jj = 1, np
               do kk = 1, np
                  bbit(ii, jj) = bbit(ii, jj) + bbi(ii, kk)*sphalf(kk, jj)
               end do
            end do
         end do

         do ii = 1, np
            do jj = 1, np
               do kk = 1, np
                  bbist(ii, jj) = bbist(ii, jj) + bbit(kk, jj)*sphalf(kk, ii)
               end do
            end do
         end do

         !      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

         call dsyev('V', 'U', np, bbist, mxprj, bbistev, work, lwork, info)
         if (info .ne. 0) then
            write (6, '(a,i4)') 'run_vkb: B^(-1)* matrix eigenvalue ERROR, info=', info
            stop
         end if

         ! take linear combinations to form orthonormal basis functions

         vkbt(:, :) = vkb(:, :, l1)
         vkbst(:, :) = 0.0d0

         do jj = 1, np
            do ii = 1, np
               vkbst(:, jj) = vkbst(:, jj) + smhalf(ii, jj)*vkbt(:, ii)
            end do
         end do

         ! take linear combinations to form orthonormal projectors

         vkb(:, :, L1) = 0.0d0

         do ii = 1, np
            evkb(ii, l1) = bbistev(ii)
            do jj = 1, np
               vkb(:, ii, l1) = vkb(:, ii, l1) + bbist(jj, ii)*vkbst(:, jj)
            end do
         end do

         ! re-order if necessary so first projector has largest magnitude coefficient
         ! which is consistent with relativistic sr_so_r.f90 output

         do ii = 1, 100
            sorted = .true.
            do jj = 2, np

               if (abs(evkb(jj, l1)) > abs(evkb(jj - 1, l1))) then
                  tt = evkb(jj - 1, l1)
                  vkbt(:, 1) = vkb(:, jj - 1, l1)
                  evkb(jj - 1, l1) = evkb(jj, l1)
                  vkb(:, jj - 1, l1) = vkb(:, jj, l1)
                  evkb(jj, l1) = tt
                  vkb(:, jj, l1) = vkbt(:, 1)
                  sorted = .false.
               end if
            end do
            if (sorted) exit
         end do

         write (6, '(/a,1p,5e12.4)') '  Orthonormal projector coefficients',&
         &        (evkb(jj, l1), jj=1, np)

         ! Set sign of projectors (physically irrelevant) so that they are positive
         ! at their peak (needed for compaisons apparently)

         do jj = 1, np
            apk = 0.0d0
            do ii = 1, mmax
               if (abs(vkb(ii, jj, l1)) > apk) then
                  apk = abs(vkb(ii, jj, l1))
                  ipk = ii
               end if
            end do
            if (vkb(ipk, jj, l1) < 0.0d0) then
               vkb(:, jj, l1) = -vkb(:, jj, l1)
            end if
         end do

      end if !nproj(l1)>=2

      if (nproj(l1) < mxprj) then
         do jj = nproj(l1) + 1, mxprj
            evkb(jj, l1) = 0.0d0
            vkb(:, jj, l1) = 0.0d0
         end do
      end if

   end do  ! l1

   deallocate (vloc, vkbt, vkbst, work)

   return
end subroutine run_vkb
!> computes Vanderbilt / Kleinman-Bylander non-local potentials
subroutine run_vkb_r(lmax, lloc, lpopt, dvloc0, irc, nproj, rr, mmax, mxprj, pswf, vfull, vp, &
&                   evkb, vkb, nlim, vr)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> lmax  maximum angular momentum
   integer, intent(in) :: lmax
   !> lloc  l for local potential (lloc==4 => use linear combination)
   integer, intent(in) :: lloc
   !> lpopt  choice of polynomial for lloc==4
   integer, intent(in) :: lpopt
   !> mmax  size of radial grid
   integer, intent(in) :: mmax
   !> mxprj  dimension of number of projectors
   integer, intent(in) :: mxprj
   !> irc  core radii indices
   integer, intent(in) :: irc(6)
   !> nproj  number of projectors for each lx
   integer, intent(in) :: nproj(6)
   !> rr  log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> pswf pseudo wave functions
   real(dp), intent(in) :: pswf(mmax, mxprj, 4, 2)
   !> vfull  all-electron potential
   real(dp), intent(in) :: vfull(mmax)
   !> vr  effective scalar-relativisic "potential" calculated in vrel
   real(dp), intent(in) :: vr(mmax, mxprj, 6, 2)
   !> dvloc0  amplitude at rr==0 to be smoothely added for lloc==4
   real(dp), intent(in) :: dvloc0

   !Input/Output variables
   !> vp  semi-local pseudopotentials for first projectors
   !>      (vp(:,5) is local potential if lloc=4)
   real(dp), intent(in out) :: vp(mmax, 5, 2)
   !> vkb  semi-local "potentials"*pswf (input); VKB projectors (output)
   real(dp), intent(in out) :: vkb(mmax, mxprj, 4, 2)
   !> evkb  coefficients of BKB projectors (output)
   real(dp), intent(in out) :: evkb(mxprj, 4, 2)
   !> nlim  index of maximum rc including that of vlocal (output)
   integer, intent(in out) :: nlim

   !Local variables
   integer :: ii, ipk, jj, kk, ikap, ll, l1, info, np, kap, mkap
   integer :: lwork
   real(dp) :: apk, sn, tt, qq12, xx, ff
   real(dp) :: bb(mxprj, mxprj), bbev(mxprj), bbi(mxprj, mxprj)
   real(dp) :: sovl(mxprj, mxprj), sovlev(mxprj), smhalf(mxprj, mxprj)
   real(dp) :: sphalf(mxprj, mxprj)
   real(dp) :: bbist(mxprj, mxprj), bbistev(mxprj), bbit(mxprj, mxprj)
   real(dp), allocatable :: vloc(:), vkbt(:, :), vkbst(:, :), work(:)
   logical :: sorted

   lwork = 3*mxprj
   allocate (vloc(mmax), vkbt(mmax, mxprj), vkbst(mmax, mxprj), work(lwork))

   ! calculate local potential
   ! use "scalar-relativistic" weighted average if real ll/=0 is specified
   ! reset semi-local potential to this average (ie., there is no SO for this l)

   l1 = 0
   if (lloc == 4) then
      call vploc(rr, vfull, vp(1, 1, 1), dvloc0, irc(5), mmax, lpopt)
      vloc(:) = vp(:, 5, 1)
      vp(:, 5, 2) = vp(:, 5, 1)
   else if (l1 == 1) then
      vloc(:) = vp(:, 1, 1)
   else
      l1 = lloc + 1
      vloc(:) = (dble(l1)*vp(:, l1, 1) + dble(l1 - 1)*vp(:, l1, 2))/dble(2*l1 - 1)
      vp(:, l1, 1) = vloc(:)
      vp(:, l1, 2) = vloc(:)
   end if

   ! Vanderbilt / Kleinman-Bylander projector construction
   nlim = irc(lloc + 1)
   do l1 = 1, lmax + 1
      if (l1 == lloc + 1) cycle
      ll = l1 - 1
      if (l1 == 1) then
         mkap = 1
      else
         mkap = 2
      end if
      ! loop on J = ll +/- 1/2
      do ikap = 1, mkap
         if (ikap == 1) kap = -(ll + 1)
         if (ikap == 2) kap = ll

         nlim = max(nlim, irc(l1))
         do jj = 1, nproj(l1)
            !incorporates scalar-relativistic "potential" in last 5% of [0:rc] to smoothly
            !force projectors to zero.
            do ii = 1, irc(l1)
               xx = 20.0d0*((rr(ii)/rr(irc(l1))) - 1.0d0)
               if (xx > -1.0d0) then
                  ff = (1.0d0 - xx**2)**2
                  !        ff=0.0d0
               else
                  ff = 0.0d0
               end if

               vkb(ii, jj, l1, ikap) = vkb(ii, jj, l1, ikap) &
               &                         - (vloc(ii) + ff*vr(ii, jj, l1, ikap))*pswf(ii, jj, l1, ikap)
            end do
            do ii = irc(l1) + 1, mmax
               vkb(ii, jj, l1, ikap) = 0.0d0
            end do
         end do
      end do !ikap
   end do !l1

   ! Vanderbilt B-matrix construction

   do l1 = 1, lmax + 1
      if (l1 == lloc + 1) cycle
      ll = l1 - 1
      if (l1 == 1) then
         mkap = 1
      else
         mkap = 2
      end if
      ! loop on J = ll +/- 1/2
      do ikap = 1, mkap
         if (ikap == 1) kap = -(ll + 1)
         if (ikap == 2) kap = ll
         np = nproj(l1)
         bb(:, :) = 0.0d0
         call vpinteg(pswf(1, 1, l1, ikap), vkb(1, 1, l1, ikap), irc(l1), 2*l1, bb(1, 1), rr)
         evkb(1, l1, ikap) = 1.0d0/bb(1, 1)
         ! this is Kleinman-Bylander result if nproj==1

         if (nproj(l1) >= 2) then

            do jj = 1, np
               do ii = 1, np
                  call vpinteg(pswf(1, ii, l1, ikap), vkb(1, jj, l1, ikap), irc(l1), 2*l1, bb(ii, jj), rr)
               end do
            end do

            ! symmetrize exactly
            write (6, '(/a,i3,a,i3)') 'B matrix Hermiticity error, ll=', l1 - 1, &
            &        '  kap=', kap
            do jj = 2, np
               do ii = 1, jj - 1
                  write (6, '(2i4,1p,d14.4)') ii, jj, bb(ii, jj) - bb(jj, ii)
                  tt = 0.5d0*(bb(ii, jj) + bb(jj, ii))
                  bb(ii, jj) = tt
                  bb(jj, ii) = tt
               end do
            end do

            ! find eigenvalues and eigenvectors of the B matrix

            !      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

            call dsyev('V', 'U', np, bb, mxprj, bbev, work, lwork, info)
            if (info .ne. 0) then
               write (6, '(a,i4)') 'run_vkb: B matrix eigenvalue ERROR, info=', info
               stop
            end if

            ! take linear combinations to form diagonal projectors

            do jj = 1, np
               vkbt(:, jj) = vkb(:, jj, l1, ikap)
            end do

            do jj = 1, np
               vkb(:, jj, l1, ikap) = 0.0d0
               evkb(jj, l1, ikap) = 1.0d0/bbev(jj)
               do ii = 1, np
                  vkb(:, jj, l1, ikap) = vkb(:, jj, l1, ikap) + bb(ii, jj)*vkbt(:, ii)
               end do
            end do

         end if !nproj>=2

         ! normalize projectors
         do jj = 1, np
            call vpinteg(vkb(1, jj, l1, ikap), vkb(1, jj, l1, ikap), irc(l1), 2*l1, sn, rr)
            if (vkb(irc(l1) - 10, jj, l1, ikap) >= 0.0d0) then
               vkb(:, jj, l1, ikap) = vkb(:, jj, l1, ikap)/sqrt(sn)
            else
               vkb(:, jj, l1, ikap) = -vkb(:, jj, l1, ikap)/sqrt(sn)
            end if
            evkb(jj, l1, ikap) = sn*evkb(jj, l1, ikap)
         end do

         ! calculate overlap
         if (nproj(l1) >= 2) then

            bbi(:, :) = 0.0d0
            do jj = 1, np
               sovl(jj, jj) = 1.0d0
               bbi(jj, jj) = evkb(jj, l1, ikap)
               do ii = jj + 1, np
                  call vpinteg(vkb(1, ii, l1, ikap), vkb(1, jj, l1, ikap), irc(l1), 2*l1, sn, rr)
                  sovl(ii, jj) = sn
                  sovl(jj, ii) = sn
               end do
            end do
            call vpinteg(vkb(1, 1, l1, ikap), vkb(1, 2, l1, ikap), irc(l1), 2*l1, sn, rr)

            ! find eigenvalues and eigenvectors of the overlap matrix

            !      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

            call dsyev('V', 'U', np, sovl, mxprj, sovlev, work, lwork, info)
            if (info .ne. 0) then
               write (6, '(a,i4)') 'run_vkb: sovl matrix eigenvalue ERROR, info=', info
               stop
            end if

            !   write(6,'(/1p,2e12.4)') (sovlev(ii),ii=1,2)

            ! construct S^(-1/2) AND s^(1/2)

            do jj = 1, np
               tt = sqrt(sovlev(jj))
               do ii = 1, np
                  sphalf(ii, jj) = tt*sovl(ii, jj)
                  smhalf(ii, jj) = sovl(ii, jj)/tt
               end do
            end do

            ! construct B^(-1)* = S^(1/2)^T B^(-1) S^(1/2)

            bbit(:, :) = 0.0d0
            bbist(:, :) = 0.0d0

            do ii = 1, np
               do jj = 1, np
                  do kk = 1, np
                     bbit(ii, jj) = bbit(ii, jj) + bbi(ii, kk)*sphalf(kk, jj)
                  end do
               end do
            end do

            do ii = 1, np
               do jj = 1, np
                  do kk = 1, np
                     bbist(ii, jj) = bbist(ii, jj) + bbit(kk, jj)*sphalf(kk, ii)
                  end do
               end do
            end do

            !      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

            call dsyev('V', 'U', np, bbist, mxprj, bbistev, work, lwork, info)
            if (info .ne. 0) then
               write (6, '(a,i4)') 'run_vkb: B^(-1)* matrix eigenvalue ERROR, info=', info
               stop
            end if

            ! take linear combinations to form orthonormal basis functions

            vkbt(:, :) = vkb(:, :, l1, ikap)
            vkbst(:, :) = 0.0d0

            do jj = 1, np
               do ii = 1, np
                  vkbst(:, jj) = vkbst(:, jj) + smhalf(ii, jj)*vkbt(:, ii)
               end do
            end do

            ! take linear combinations to form orthonormal projectors

            vkb(:, :, L1, ikap) = 0.0d0

            do ii = 1, np
               evkb(ii, l1, ikap) = bbistev(ii)
               do jj = 1, np
                  vkb(:, ii, l1, ikap) = vkb(:, ii, l1, ikap) + bbist(jj, ii)*vkbst(:, jj)
               end do
            end do

            ! re-order if necessary so first projector has largest magnitude coefficient
            ! which is consistent with relativistic sr_so_r.f90 output

            do ii = 1, 100
               sorted = .true.
               do jj = 2, np

                  if (abs(evkb(jj, l1, ikap)) > abs(evkb(jj - 1, l1, ikap))) then
                     tt = evkb(jj - 1, l1, ikap)
                     vkbt(:, 1) = vkb(:, jj - 1, l1, ikap)
                     evkb(jj - 1, l1, ikap) = evkb(jj, l1, ikap)
                     vkb(:, jj - 1, l1, ikap) = vkb(:, jj, l1, ikap)
                     evkb(jj, l1, ikap) = tt
                     vkb(:, jj, l1, ikap) = vkbt(:, 1)
                     sorted = .false.
                  end if
               end do
               if (sorted) exit
            end do

            write (6, '(/a,1p,5e12.4)') '  Orthonormal projector coefficients',&
            &        (evkb(jj, l1, ikap), jj=1, np)

            ! Set sign of projectors (physically irrelevant) so that they are positive
            ! at their peak (needed for compaisons apparently)

            do jj = 1, np
               apk = 0.0d0
               do ii = 1, mmax
                  if (abs(vkb(ii, jj, l1, ikap)) > apk) then
                     apk = abs(vkb(ii, jj, l1, ikap))
                     ipk = ii
                  end if
               end do
               if (vkb(ipk, jj, l1, ikap) < 0.0d0) then
                  vkb(:, jj, l1, ikap) = -vkb(:, jj, l1, ikap)
               end if
            end do

         end if !nproj(l1)>=2

         if (nproj(l1) < mxprj) then
            do jj = nproj(l1) + 1, mxprj
               evkb(jj, l1, ikap) = 0.0d0
               vkb(:, jj, l1, ikap) = 0.0d0
            end do
         end if
      end do !ikap
   end do !l1

   deallocate (vloc, vkbt, vkbst, work)

   return
end subroutine run_vkb_r
!> reformulates non-local potentials based on j = l +/- 1/2 to scalar-
!> relativistic and L dot S projectors
!> uses relationship <L dot S> = (J^2 - L^2 - S^2)/2
!> so L dot S = +/- l/2 for j = l +/- 1/2
subroutine sr_so_r(lmax, irc, nproj, rr, mmax, mxprj, evkb, vkb, &
&                       vsr, esr, vso, eso)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> lmax  maximum angular momentum
   integer, intent(in) :: lmax
   !> mmax  dimension of log grid
   integer, intent(in) :: mmax
   !> mxprj  dimension of number of projectors
   integer, intent(in) :: mxprj
   !> irc  core radii indices
   integer, intent(in) :: irc(6)
   !> nproj  number of projectors for each l
   integer, intent(in) :: nproj(6)
   !> rr  log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> vkb  vkb projectors
   real(dp), intent(in) :: vkb(mmax, mxprj, 4, 2)
   !> evkb  coefficients of BKB projectors
   real(dp), intent(in) :: evkb(mxprj, 4, 2)

   !Output variables
   !> vsr  normalized scalar projectors
   real(dp), intent(out) :: vsr(mmax, 2*mxprj, 4)
   !> esr  energy  coefficients of vscal
   real(dp), intent(out) :: esr(2*mxprj, 4)
   !> vso  normalized spin-orbig projectors
   real(dp), intent(out) :: vso(mmax, 2*mxprj, 4)
   !> esol  energy  coefficients of vso
   real(dp), intent(out) :: eso(2*mxprj, 4)

   !Local variables
   integer :: ii, jj, kk, ierr, ik1, ik2, ip1, ip2, ipk, ll, l1, info, nn
   real(dp) :: apk, sn, tt, qq1mxprj
   real(dp) :: sovl(2*mxprj, 2*mxprj), sovlev(2*mxprj), ascl(2*mxprj, 2*mxprj), aso(2*mxprj, 2*mxprj)
   real(dp) :: asclst(2*mxprj, 2*mxprj), wsclst(2*mxprj), asost(2*mxprj, 2*mxprj), wsost(2*mxprj)
   real(dp) :: asclt(2*mxprj, 2*mxprj), asot(2*mxprj, 2*mxprj)
   real(dp) :: sphalf(2*mxprj, 2*mxprj), smhalf(2*mxprj, 2*mxprj)
   real(dp) :: fscl(mxprj), fso(mxprj), work(10*mxprj)
   real(dp), allocatable :: vkbt(:, :), vkbst(:, :)
   logical :: sorted

   allocate (vkbt(mmax, 2*mxprj), vkbst(mmax, 2*mxprj))

   do l1 = 1, lmax + 1
      ll = l1 - 1

      if (ll == 0) then
         vsr(:, :, l1) = 0.0d0
         vso(:, :, l1) = 0.0d0
         esr(:, l1) = 0.0d0
         eso(:, l1) = 0.0d0
         if (nproj(l1) >= 1) then
            do ii = 1, nproj(l1)
               vsr(:, ii, l1) = vkb(:, ii, l1, 1)
               esr(ii, l1) = evkb(ii, l1, 1)
            end do
         end if
         cycle
      end if

      nn = 2*nproj(l1)

      fscl(1) = (ll + 1)/dble(2*ll + 1)
      fscl(2) = ll/dble(2*ll + 1)
      fso(1) = 2/dble(2*ll + 1)
      fso(2) = -2/dble(2*ll + 1)

      ! construct overlap matrix and diagonal energy matrices

      sovl(:, :) = 0.0d0
      ascl(:, :) = 0.0d0
      aso(:, :) = 0.0d0
      vkbt(:, :) = 0.0d0

      do ik1 = 1, 2
         do ip1 = 1, nproj(l1)
            ii = ip1 + (ik1 - 1)*nproj(l1)

            ascl(ii, ii) = fscl(ik1)*evkb(ip1, l1, ik1)
            aso(ii, ii) = fso(ik1)*evkb(ip1, l1, ik1)

            vkbt(:, ii) = vkb(:, ip1, l1, ik1)

            do ik2 = 1, 2
               do ip2 = 1, nproj(l1)
                  jj = ip2 + (ik2 - 1)*nproj(l1)

                  call vpinteg(vkb(1, ip1, l1, ik1), vkb(1, ip2, l1, ik2), irc(l1), 2*l1, &
                  &                   sovl(ii, jj), rr)

               end do
            end do
         end do
      end do

      call dsyev('V', 'U', nn, sovl, 2*mxprj, sovlev, work, 10*mxprj, info)

      if (info .ne. 0) then
         write (6, '(a,i4)') 'sr_so_r: S matrix eigenvalue ERROR, info=', info
         stop
      end if

      ! construct S^(-1/2) AND s^(1/2)

      do jj = 1, nn
         tt = sqrt(sovlev(jj))
         do ii = 1, nn
            sphalf(ii, jj) = tt*sovl(ii, jj)
            smhalf(ii, jj) = sovl(ii, jj)/tt
         end do
      end do

      ! take linear combinations to form orthonormal basis functions

      vkbst(:, :) = 0.0d0

      do jj = 1, nn
         do ii = 1, nn
            vkbst(:, jj) = vkbst(:, jj) + smhalf(ii, jj)*vkbt(:, ii)
         end do
      end do

      ! construct A^(-1)* = S^(1/2)^T A^(-1) S^(1/2)

      asclt(:, :) = 0.0d0
      asclst(:, :) = 0.0d0
      asot(:, :) = 0.0d0
      asost(:, :) = 0.0d0

      do ii = 1, nn
         do jj = 1, nn
            do kk = 1, nn
               asclt(ii, jj) = asclt(ii, jj) + ascl(ii, kk)*sphalf(kk, jj)
               asot(ii, jj) = asot(ii, jj) + aso(ii, kk)*sphalf(kk, jj)
            end do
         end do
      end do

      do ii = 1, nn
         do jj = 1, nn
            do kk = 1, nn
               asclst(ii, jj) = asclst(ii, jj) + asclt(kk, jj)*sphalf(kk, ii)
               asost(ii, jj) = asost(ii, jj) + asot(kk, jj)*sphalf(kk, ii)
            end do
         end do
      end do

      ! find eigenvalues and eigenvectors of the A* matrices

      !      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

      call dsyev('V', 'U', nn, asclst, 2*mxprj, wsclst, work, 10*mxprj, info)

      if (info .ne. 0) then
         write (6, '(a,i4)') 'sr_so_r: A* matrix eigenvalue ERROR, info=', info
         stop
      end if

      call dsyev('V', 'U', nn, asost, 2*mxprj, wsost, work, 10*mxprj, info)

      if (info .ne. 0) then
         write (6, '(a,i4)') 'sr_so_r: A* matrix eigenvalue ERROR, info=', info
         stop
      end if

      ! take linear combinations to form orthonormal projectors

      vsr(:, :, l1) = 0.0d0
      vso(:, :, l1) = 0.0d0
      esr(:, l1) = 0.0d0
      eso(:, l1) = 0.0d0

      do ii = 1, nn
         esr(ii, l1) = wsclst(ii)
         eso(ii, l1) = wsost(ii)
         do jj = 1, nn
            vsr(:, ii, l1) = vsr(:, ii, l1) + asclst(jj, ii)*vkbst(:, jj)
            vso(:, ii, l1) = vso(:, ii, l1) + asost(jj, ii)*vkbst(:, jj)
         end do
      end do

      ! bubble-sort on coefficient magnitudes for scalar and then s-o
      ! (Yes, I know bubble-sort is the least-efficient sorting algorithm.)

      do ii = 1, 100
         sorted = .true.
         do jj = 2, nn
            if (abs(esr(jj - 1, l1)) < abs(esr(jj, l1))) then
               tt = esr(jj, l1)
               vkbt(:, 1) = vsr(:, jj, l1)
               esr(jj, l1) = esr(jj - 1, l1)
               vsr(:, jj, l1) = vsr(:, jj - 1, l1)
               esr(jj - 1, l1) = tt
               vsr(:, jj - 1, l1) = vkbt(:, 1)
               sorted = .false.
            end if
         end do
         if (sorted) exit
      end do

      do ii = 1, 100
         sorted = .true.
         do jj = 2, nn
            if (abs(eso(jj - 1, l1)) < abs(eso(jj, l1))) then
               tt = eso(jj, l1)
               vkbt(:, 1) = vso(:, jj, l1)
               eso(jj, l1) = eso(jj - 1, l1)
               vso(:, jj, l1) = vso(:, jj - 1, l1)
               eso(jj - 1, l1) = tt
               vso(:, jj - 1, l1) = vkbt(:, 1)
               sorted = .false.
            end if
         end do
         if (sorted) exit
      end do

      write (6, '(/a,i2)') &
      &         ' Orthonormal scalar projector coefficients, l = ', ll
      write (6, '(1p,6e12.4)') (esr(jj, l1), jj=1, nn)
      write (6, '(/a,i2)') &
      &         ' Orthonormal spin-orbit projector coefficients, l = ', ll
      write (6, '(1p,6e12.4)') (eso(jj, l1), jj=1, nn)

      ! Set sign of projectors (physically irrelevant) so that they are positive
      ! at their peak (needed for compaisons apparently)

      do jj = 1, nn
         apk = 0.0d0
         do ii = 1, mmax
            if (abs(vso(ii, jj, l1)) > apk) then
               apk = abs(vso(ii, jj, l1))
               ipk = ii
            end if
         end do
         if (vso(ipk, jj, l1) < 0.0d0) then
            vso(:, jj, l1) = -vso(:, jj, l1)
         end if
         apk = 0.0d0
         do ii = 1, mmax
            if (abs(vsr(ii, jj, l1)) > apk) then
               apk = abs(vsr(ii, jj, l1))
               ipk = ii
            end if
         end do
         if (vsr(ipk, jj, l1) < 0.0d0) then
            vsr(:, jj, l1) = -vsr(:, jj, l1)
         end if
      end do

   end do ! l1

   deallocate (vkbt, vkbst)
   return
end subroutine sr_so_r
end module vkb_m

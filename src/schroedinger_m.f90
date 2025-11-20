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
module schroedinger_m
   use adams_m, only: aeo, aio, aei, aii
   use vpinteg_m, only: vpinteg
   implicit none
   private
   public :: lschfb
   public :: lschfs
   public :: lschkb
   public :: lschpb
   public :: lschps
   public :: lschpsbar
   public :: lschpse
   public :: lschvkbb
   public :: lschvkbbe
   public :: lschvkbs
   public :: vkboutwf
contains
!> Finds bound states of an all-electron atomic potential using
!> Pauli-type  scalar-relativistic Schroedinger equation
subroutine lschfb(nn, ll, ierr, ee, rr, vv, uu, up, zz, mmax, mch, srel)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> mmax  size of log grid
   integer, intent(in) :: mmax
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> vv  local atomic potential
   real(dp), intent(in) :: vv(mmax)
   !> zz  atomic number
   real(dp), intent(in) :: zz
   !> nn  principal quantum number
   integer, intent(in) :: nn
   !> ll  angular-momentum quantum number
   integer, intent(in) :: ll
   !> srel .true. for scalar-relativistic, .false. for non-relativistic
   logical, intent(in) :: srel

   !Output variables
   !> uu  output radial wave function (*rr)
   real(dp) :: uu(mmax)
   !> up  d(uu)/dr
   real(dp) :: up(mmax)
   !> ee  bound-state energy, input guess and output calculated value
   real(dp) :: ee
   !> ierr  non-zero return if error
   integer :: ierr
   !> mch matching mesh point for inward-outward integrations
   integer :: mch

   !Local variables
   real(dp) :: als
   real(dp) :: de, emax, emin
   real(dp) :: eps, fss, tfss, gamma, ro, sc
   real(dp) :: sls, sn, cn, uout, upin, upout, xkap
   real(dp) :: amesh, al, xx
   integer :: ii, it, nint, node, nin

   real(dp), allocatable :: upp(:), cf(:), dv(:), fr(:), frp(:)
   allocate (upp(mmax), cf(mmax), dv(mmax), fr(mmax), frp(mmax))

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)

   ! convergence factor for solution of schroedinger eq.  if calculated
   ! correction to eigenvalue is smaller in magnitude than eps times
   ! the magnitude of the current guess, the current guess is not changed.
   eps = 1.0d-10
   ierr = 60

   ! check arguments
   if (ll > nn - 1) then
      write (6, '(/a,i4,a,i4)') 'lschfb: ERROR ll =', ll, ' > nn =', nn
      ierr = 1
      return
   end if

   ! relativistic - non-relativistic switch
   if (srel) then
      fss = (1.0d0/137.036d0)**2
   else
      fss = 1.0d-12
   end if

   if (ll == 0) gamma = dsqrt(1.0d0 - fss*zz**2)
   if (ll > 0) gamma = (ll*dsqrt(ll**2 - fss*zz**2) + &
   & (ll + 1)*dsqrt((ll + 1)**2 - fss*zz**2))/(2*ll + 1)

   sls = ll*(ll + 1)

   emax = vv(mmax) + 0.5d0*sls/rr(mmax)**2
   emin = emax
   do ii = 1, mmax
      if (ll == 0) then
         ! ad-hock step to eliminate probelms from absurd emin for ll=0
         emin = dmin1(emin, vv(ii) + 0.25d0/rr(ii)**2)
      else
         emin = dmin1(emin, vv(ii) + 0.5d0*sls/rr(ii)**2)
      end if
   end do
   if (ee == 0.0d0) ee = 0.5d0*(emax + emin)
   if (ee < emin) ee = 0.5d0*(emax + emin)
   if (ee > emax) ee = 0.5d0*(emax + emin)

   ! null arrays to remove leftover garbage
   uu(:) = 0.0d0
   up(:) = 0.0d0
   upp(:) = 0.0d0

   als = al**2

   ! return point for bound state convergence
   do nint = 1, 60

      ! coefficient array for u in differential eq.
      do ii = 1, mmax
         cf(ii) = als*sls + 2.0d0*als*(vv(ii) - ee)*rr(ii)**2
      end do

      ! calculate dv/dr for darwin correction
      dv(:) = 0.0d0

      dv(1) = (-50.d0*vv(1) + 96.d0*vv(2) - 72.d0*vv(3) + 32.d0*vv(4) &
      &         - 6.d0*vv(5))/(24.d0*al*rr(1))
      dv(2) = (-6.d0*vv(1) - 20.d0*vv(2) + 36.d0*vv(3) - 12.d0*vv(4) &
      &         + 2.d0*vv(5))/(24.d0*al*rr(2))

      do ii = 3, mmax - 2
         dv(ii) = (2.d0*vv(ii - 2) - 16.d0*vv(ii - 1) + 16.d0*vv(ii + 1) &
         &          - 2.d0*vv(ii + 2))/(24.d0*al*rr(ii))
      end do

      !  relativistic coefficient arrays for u (fr) and up (frp).
      do ii = 1, mmax
         tfss = fss
         fr(ii) = als*(rr(ii)**2)*(-tfss*(vv(ii) - ee)**2 + 0.5d0*tfss*dv(ii)/ &
         &     (rr(ii)*(1.0d0 + 0.5d0*tfss*(ee - vv(ii)))))
         frp(ii) = -al*rr(ii)*0.5d0*tfss*dv(ii)/(1.0d0 + 0.5d0*tfss*(ee - vv(ii)))
      end do

      ! find classical turning point for matching
      mch = 0
      do ii = mmax, 2, -1
         if (cf(ii - 1) <= 0.d0 .and. cf(ii) > 0.d0) then
            mch = ii
            exit
         end if
      end do

      if (mch == 0) then
         ierr = -1
         return
      end if

      ! start wavefunction with series

      do ii = 1, 4
         uu(ii) = rr(ii)**gamma
         up(ii) = al*gamma*rr(ii)**gamma
         upp(ii) = (al + frp(ii))*up(ii) + (cf(ii) + fr(ii))*uu(ii)
      end do

      ! outward integration using predictor once, corrector
      ! twice
      node = 0

      do ii = 4, mch - 1
         uu(ii + 1) = uu(ii) + aeo(up, ii)
         up(ii + 1) = up(ii) + aeo(upp, ii)
         do it = 1, 2
            upp(ii + 1) = (al + frp(ii + 1))*up(ii + 1) + (cf(ii + 1) + fr(ii + 1))*uu(ii + 1)
            up(ii + 1) = up(ii) + aio(upp, ii)
            uu(ii + 1) = uu(ii) + aio(up, ii)
         end do
         if (uu(ii + 1)*uu(ii) .le. 0.0d0) node = node + 1
      end do

      uout = uu(mch)
      upout = up(mch)

      if (node - nn + ll + 1 == 0) then

         ! start inward integration at 10*classical turning
         ! point with simple exponential

         nin = mch + 2.3d0/al
         if (nin + 4 > mmax) nin = mmax - 4
         xkap = dsqrt(sls/rr(nin)**2 + 2.0d0*(vv(nin) - ee))

         do ii = nin, nin + 4
            uu(ii) = exp(-xkap*(rr(ii) - rr(nin)))
            up(ii) = -rr(ii)*al*xkap*uu(ii)
            upp(ii) = (al + frp(ii))*up(ii) + (cf(ii) + fr(ii))*uu(ii)
         end do

         ! integrate inward

         do ii = nin, mch + 1, -1
            uu(ii - 1) = uu(ii) + aei(up, ii)
            up(ii - 1) = up(ii) + aei(upp, ii)
            do it = 1, 2
               upp(ii - 1) = (al + frp(ii - 1))*up(ii - 1) + (cf(ii - 1) + fr(ii - 1))*uu(ii - 1)
               up(ii - 1) = up(ii) + aii(upp, ii)
               uu(ii - 1) = uu(ii) + aii(up, ii)
            end do
         end do

         ! scale outside wf for continuity

         sc = uout/uu(mch)

         do ii = mch, nin
            up(ii) = sc*up(ii)
            uu(ii) = sc*uu(ii)
         end do

         upin = up(mch)

         ! perform normalization sum

         ro = rr(1)/dsqrt(amesh)
         sn = ro**(2.0d0*gamma + 1.0d0)/(2.0d0*gamma + 1.0d0)

         do ii = 1, nin - 3
            sn = sn + al*rr(ii)*uu(ii)**2
         end do

         sn = sn + al*(23.0d0*rr(nin - 2)*uu(nin - 2)**2 &
         &              + 28.0d0*rr(nin - 1)*uu(nin - 1)**2 &
         &              + 9.0d0*rr(nin)*uu(nin)**2)/24.0d0

         ! normalize u

         cn = 1.0d0/dsqrt(sn)
         uout = cn*uout
         upout = cn*upout
         upin = cn*upin

         do ii = 1, nin
            up(ii) = cn*up(ii)
            uu(ii) = cn*uu(ii)
         end do
         do ii = nin + 1, mmax
            uu(ii) = 0.0d0
         end do

         ! perturbation theory for energy shift

         de = 0.5d0*uout*(upout - upin)/(al*rr(mch))

         ! convergence test and possible exit

         if (dabs(de) < dmax1(dabs(ee), 0.2d0)*eps) then
            ierr = 0
            exit
         end if

         if (de > 0.0d0) then
            emin = ee
         else
            emax = ee
         end if
         ee = ee + de
         if (ee > emax .or. ee < emin) ee = 0.5d0*(emax + emin)

      else if (node - nn + ll + 1 < 0) then
         ! too few nodes
         emin = ee
         ee = 0.5d0*(emin + emax)

      else
         ! too many nodes
         emax = ee
         ee = 0.5d0*(emin + emax)
      end if

   end do

   !fix sign to be positive at rr->oo
   if (uu(mch) < 0.0d0) then
      uu(:) = -uu(:)
      up(:) = -up(:)
   end if

   deallocate (upp, cf, dv, fr, frp)
   return

end subroutine lschfb
!> integrates radial Pauli-type scalar-relativistic equation
!> on a logarithmic mesh
!> modified routine to be used in finding norm-conserving
!> pseudopotential
subroutine lschfs(nn, ll, ierr, ee, rr, vv, uu, up, zz, mmax, mch, srel)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> mmax  size of log grid
   integer, intent(in) :: mmax
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> vv  local atomic potential
   real(dp), intent(in) :: vv(mmax)
   !> zz  atomic number
   real(dp), intent(in) :: zz
   !> ll  angular-momentum quantum number
   integer, intent(in) :: ll
   !> mch matching mesh point for inward-outward integrations
   integer, intent(in) :: mch
   !> srel .true. for scalar-relativistic, .false. for non-relativistic
   logical, intent(in) :: srel

   !Output variables
   !> uu  output radial wave function (*rr)
   real(dp) :: uu(mmax)
   !> up  d(uu)/dr
   real(dp) :: up(mmax)
   !> ee  bound-state energy, input guess and output calculated value
   real(dp) :: ee
   !> ierr  non-zero return if error
   integer :: ierr
   !> nn  effective principal quantum number based on nodes inside mch (output)
   integer :: nn

   !Local variables
   real(dp) :: amesh, al
   real(dp) :: als, cn
   real(dp) :: fss, tfss, gamma, ro, xx
   real(dp) :: sls, sn, uout, upout
   integer :: ii, it, node

   real(dp), allocatable :: upp(:), cf(:), dv(:), fr(:), frp(:)

   allocate (upp(mmax), cf(mmax), dv(mmax), fr(mmax), frp(mmax))

   al = 0.01d0*log(rr(101)/rr(1))
   amesh = exp(al)

   ierr = 0

   ! relativistic - non-relativistic switch
   if (srel) then
      fss = (1.0d0/137.036d0)**2
   else
      fss = 1.0d-12
   end if

   if (ll == 0) gamma = sqrt(1.0d0 - fss*zz**2)
   if (ll > 0) gamma = (ll*sqrt(ll**2 - fss*zz**2) &
   & + (ll + 1)*sqrt((ll + 1)**2 - fss*zz**2))/(2*ll + 1)

   sls = ll*(ll + 1)

   ! null arrays to remove leftover garbage

   uu(:) = 0.0d0
   up(:) = 0.0d0
   upp(:) = 0.0d0

   node = 0

   als = al**2

   ! coefficient array for u in differential eq.
   do ii = 1, mmax
      cf(ii) = als*sls + 2.0d0*als*(vv(ii) - ee)*rr(ii)**2
   end do

   ! calculate dv/dr for darwin correction
   dv(:) = 0.0d0

   dv(1) = (-50.d0*vv(1) + 96.d0*vv(2) - 72.d0*vv(3) + 32.d0*vv(4) &
   &       - 6.d0*vv(5))/(24.d0*al*rr(1))
   dv(2) = (-6.d0*vv(1) - 20.d0*vv(2) + 36.d0*vv(3) - 12.d0*vv(4) &
   &       + 2.d0*vv(5))/(24.d0*al*rr(2))

   do ii = 3, mmax - 2
      dv(ii) = (2.d0*vv(ii - 2) - 16.d0*vv(ii - 1) + 16.d0*vv(ii + 1) &
      &         - 2.d0*vv(ii + 2))/(24.d0*al*rr(ii))
   end do

   !  relativistic coefficient arrays for u (fr) and up (frp).
   do ii = 1, mmax
      tfss = fss
      fr(ii) = als*(rr(ii)**2)*(-tfss*(vv(ii) - ee)**2 + 0.5d0*tfss*dv(ii)/ &
      &   (rr(ii)*(1.0d0 + 0.5d0*tfss*(ee - vv(ii)))))
      frp(ii) = -al*rr(ii)*0.5d0*tfss*dv(ii)/(1.0d0 + 0.5d0*tfss*(ee - vv(ii)))
   end do

   ! start wavefunction with series

   do ii = 1, 4
      uu(ii) = rr(ii)**gamma
      up(ii) = al*gamma*rr(ii)**gamma
      upp(ii) = (al + frp(ii))*up(ii) + (cf(ii) + fr(ii))*uu(ii)
   end do

   ! outward integration using predictor once, corrector
   ! twice

   do ii = 4, mch - 1
      uu(ii + 1) = uu(ii) + aeo(up, ii)
      up(ii + 1) = up(ii) + aeo(upp, ii)
      do it = 1, 2
         upp(ii + 1) = (al + frp(ii + 1))*up(ii + 1) + (cf(ii + 1) + fr(ii + 1))*uu(ii + 1)
         up(ii + 1) = up(ii) + aio(upp, ii)
         uu(ii + 1) = uu(ii) + aio(up, ii)
      end do
      if (uu(ii + 1)*uu(ii) .le. 0.0d0) node = node + 1
   end do

   uout = uu(mch)
   upout = up(mch)

   !perform normalization sum

   ro = rr(1)/dsqrt(amesh)
   sn = ro**(2.0d0*gamma + 1.0d0)/(2.0d0*gamma + 1.0d0)

   do ii = 1, mch - 3
      sn = sn + al*rr(ii)*uu(ii)**2
   end do

   sn = sn + al*(23.0d0*rr(mch - 2)*uu(mch - 2)**2 &
   &           + 28.0d0*rr(mch - 1)*uu(mch - 1)**2 &
   &          + 9.0d0*rr(mch)*uu(mch)**2)/24.0d0

   !normalize u

   cn = 1.0d0/dsqrt(sn)
   uout = cn*uout
   upout = cn*upout

   do ii = 1, mch
      up(ii) = cn*up(ii)
      uu(ii) = cn*uu(ii)
   end do
   do ii = mch + 1, mmax
      uu(ii) = 0.0d0
   end do

   !calculate effective principal quantum number as if this were a bound
   !state with a barrier at mch
   nn = node + ll + 1

   deallocate (upp, cf, dv, fr, frp)

   return
end subroutine lschfs
!> outward integration of the inhomogeneous radial Schroedinger equation
!> on a logarithmic mesh with local potential and one proector term
subroutine lschkb(nn, ll, ierr, ee, vkb, rr, vv, uu, up, mmax, mch)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> mmax  size of log grid
   integer, intent(in) :: mmax
   !> mch matching mesh point for inward-outward integrations
   integer, intent(in) :: mch
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> vv  local pseudopotential
   real(dp), intent(in) :: vv(mmax)
   !> vkb Vanderbilt-Kleinman-bylander projector
   real(dp), intent(in) :: vkb(mmax)
   !> nn  principal quantum number (not used)
   integer, intent(in) :: nn
   !> ll  angular-momentum quantum number
   integer, intent(in) :: ll

   !Output variables
   !> uu  output radial wave function (*rr)
   real(dp), intent(out) :: uu(mmax)
   !> up  d(uu)/dr
   real(dp), intent(out) :: up(mmax)
   !> ee  bound-state energy, input guess and output calculated value
   real(dp), intent(out) :: ee
   !> ierr  non-zero return if error
   integer, intent(out) :: ierr

   !Local variables
   real(dp) :: amesh, al
   real(dp) :: als, cn
   real(dp) :: sls, uout, upout
   real(dp) :: akb, ckb
   integer :: ii, it

   real(dp), allocatable :: upp(:), cf(:)

   allocate (upp(mmax), cf(mmax))

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)

   ierr = 0

   sls = ll*(ll + 1)

   ! null arrays to remove leftover garbage

   uu(:) = 0.0d0
   up(:) = 0.0d0
   upp(:) = 0.0d0

   als = al**2

   ! coefficient array for uu in differential eq.
   do ii = 1, mmax
      cf(ii) = als*sls + 2.0d0*als*(vv(ii) - ee)*rr(ii)**2
   end do

   ! start wavefunction with series based on projector

   ckb = 2.0d0*vkb(1)/rr(1)**(ll + 1)
   akb = ckb/(6.0d0 + 4.0d0*ll)
   do ii = 1, 4
      uu(ii) = akb*rr(ii)**(ll + 3)
      up(ii) = al*(ll + 3)*uu(ii)
      upp(ii) = als*(ll + 3)**2*uu(ii)
   end do

   ! outward integration using predictor once, corrector
   ! twice

   do ii = 4, mch - 1
      uu(ii + 1) = uu(ii) + aeo(up, ii)
      up(ii + 1) = up(ii) + aeo(upp, ii)
      do it = 1, 2
         upp(ii + 1) = al*up(ii + 1) + cf(ii + 1)*uu(ii + 1) &
         &            + 2.0d0*als*vkb(ii + 1)*rr(ii + 1)**2
         up(ii + 1) = up(ii) + aio(upp, ii)
         uu(ii + 1) = uu(ii) + aio(up, ii)
      end do
   end do

   uout = uu(mch)
   upout = up(mch)

   deallocate (upp, cf)

   return
end subroutine lschkb
!> Finds bound states of a semi-local pseudopotential
subroutine lschpb(nn, ll, ierr, ee, rr, vv, uu, up, mmax, mch)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> mmax  size of log grid
   integer, intent(in) :: mmax
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> vv  local psp
   real(dp), intent(in) :: vv(mmax)
   !> nn  principal quantum number
   integer, intent(in) :: nn
   !> ll  angular-momentum quantum number
   integer, intent(in) :: ll

   !Output variables
   !> uu  output radial wave function (*rr)
   real(dp), intent(out) :: uu(mmax)
   !> up  d(uu)/dr
   real(dp), intent(out) :: up(mmax)
   !> ee  bound-state energy, input guess and output calculated value
   real(dp), intent(out) :: ee
   !> ierr  non-zero return if error
   integer, intent(out) :: ierr
   !> mch matching mesh point for inward-outward integrations
   integer, intent(out) :: mch

   !Local Variables
   real(dp) :: de, emax, emin
   real(dp) :: eps, exp, ro, sc
   real(dp) :: sls, sn, cn, uout, upin, upout, xkap
   real(dp) :: amesh, al, als
   integer :: ii, it, nint, node, nin

   real(dp), allocatable :: upp(:), cf(:)
   allocate (upp(mmax), cf(mmax))

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)

   ! convergence factor for solution of schroedinger eq.  if calculated
   ! correction to eigenvalue is smaller in magnitude than eps times
   ! the magnitude of the current guess, the current guess is not changed.
   eps = 1.0d-10
   ierr = 60

   sls = ll*(ll + 1)

   emax = vv(mmax) + 0.5d0*sls/rr(mmax)**2
   emin = 0.0d0
   do ii = 1, mmax
      emin = dmin1(emin, vv(ii) + 0.5d0*sls/rr(ii)**2)
   end do
   emin = 4.0d0*emin
   if (ee > emax) ee = 1.25d0*emax
   if (ee < emin) ee = 0.75d0*emin
   if (ee > emax) ee = 0.5d0*(emax + emin)

   ! null arrays to remove leftover garbage
   uu(:) = 0.0d0
   up(:) = 0.0d0
   upp(:) = 0.0d0

   als = al**2

   ! return point for bound state convergence
   do nint = 1, 60

      ! coefficient array for u in differential eq.
      do ii = 1, mmax
         cf(ii) = als*sls + 2.0d0*als*(vv(ii) - ee)*rr(ii)**2
      end do

      ! find classical turning point for matching
      mch = 0
      do ii = mmax, 2, -1
         if (cf(ii - 1) <= 0.d0 .and. cf(ii) > 0.d0) then
            mch = ii
            exit
         end if
      end do
      if (mch == 0) then
         write (6, '(/a)') 'lschpb: ERROR no classical turning point'
         stop
      end if

      ! start wavefunction with series

      do ii = 1, 4
         uu(ii) = rr(ii)**(ll + 1)
         up(ii) = al*(ll + 1)*rr(ii)**(ll + 1)
         upp(ii) = al*up(ii) + cf(ii)*uu(ii)
      end do

      ! outward integration using predictor once, corrector
      ! twice
      node = 0

      do ii = 4, mch - 1
         uu(ii + 1) = uu(ii) + aeo(up, ii)
         up(ii + 1) = up(ii) + aeo(upp, ii)
         do it = 1, 2
            upp(ii + 1) = al*up(ii + 1) + cf(ii + 1)*uu(ii + 1)
            up(ii + 1) = up(ii) + aio(upp, ii)
            uu(ii + 1) = uu(ii) + aio(up, ii)
         end do
         if (uu(ii + 1)*uu(ii) .le. 0.0d0) node = node + 1
      end do

      uout = uu(mch)
      upout = up(mch)

      if (node - nn + ll + 1 == 0) then

         ! start inward integration at 10*classical turning
         ! point with simple exponential

         nin = mch + 2.3d0/al
         if (nin + 4 > mmax) nin = mmax - 4
         xkap = dsqrt(sls/rr(nin)**2 + 2.0d0*(vv(nin) - ee))

         do ii = nin, nin + 4
            uu(ii) = exp(-xkap*(rr(ii) - rr(nin)))
            up(ii) = -rr(ii)*al*xkap*uu(ii)
            upp(ii) = al*up(ii) + cf(ii)*uu(ii)
         end do

         ! integrate inward

         do ii = nin, mch + 1, -1
            uu(ii - 1) = uu(ii) + aei(up, ii)
            up(ii - 1) = up(ii) + aei(upp, ii)
            do it = 1, 2
               upp(ii - 1) = al*up(ii - 1) + cf(ii - 1)*uu(ii - 1)
               up(ii - 1) = up(ii) + aii(upp, ii)
               uu(ii - 1) = uu(ii) + aii(up, ii)
            end do
         end do

         ! scale outside wf for continuity

         sc = uout/uu(mch)

         do ii = mch, nin
            up(ii) = sc*up(ii)
            uu(ii) = sc*uu(ii)
         end do

         upin = up(mch)

         ! perform normalization sum

         ro = rr(1)/dsqrt(amesh)
         sn = ro**(2*ll + 3)/dfloat(2*ll + 3)

         do ii = 1, nin - 3
            sn = sn + al*rr(ii)*uu(ii)**2
         end do

         sn = sn + al*(23.0d0*rr(nin - 2)*uu(nin - 2)**2 &
         &              + 28.0d0*rr(nin - 1)*uu(nin - 1)**2 &
         &              + 9.0d0*rr(nin)*uu(nin)**2)/24.0d0

         ! normalize u

         cn = 1.0d0/dsqrt(sn)
         uout = cn*uout
         upout = cn*upout
         upin = cn*upin

         do ii = 1, nin
            up(ii) = cn*up(ii)
            uu(ii) = cn*uu(ii)
         end do
         do ii = nin + 1, mmax
            uu(ii) = 0.0d0
         end do

         ! perturbation theory for energy shift

         de = 0.5d0*uout*(upout - upin)/(al*rr(mch))

         ! convergence test and possible exit

         if (dabs(de) < dmax1(dabs(ee), 0.2d0)*eps) then
            ierr = 0
            exit
         end if

         if (de > 0.0d0) then
            emin = ee
         else
            emax = ee
         end if
         ee = ee + de
         if (ee > emax .or. ee < emin) ee = 0.5d0*(emax + emin)

      else if (node - nn + ll + 1 < 0) then
         ! too few nodes
         emin = ee
         ee = 0.5d0*(emin + emax)

      else
         ! too many nodes
         emax = ee
         ee = 0.5d0*(emin + emax)
      end if

   end do

   deallocate (upp, cf)
   return

end subroutine lschpb
!> outward integration of Srcroedinger equation for semi-local pseudopotential
!> on a logarithmic mesh
subroutine lschps(nn, ll, ierr, ee, rr, vv, uu, up, mmax, mch)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> mmax  size of log grid
   integer, intent(in) :: mmax
   !> mch matching mesh point for inward-outward integrations
   integer, intent(in) :: mch
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> vv  semi-local atomic pseudopotential
   real(dp), intent(in) :: vv(mmax)
   !> nn  principal quantum number (not used actually)
   integer, intent(in) :: nn
   !> ll  angular-momentum quantum number
   integer, intent(in) :: ll

   !Output variables
   !> uu  output radial wave function (*rr)
   real(dp), intent(out) :: uu(mmax)
   !> up  d(uu)/dr
   real(dp), intent(out) :: up(mmax)
   !> ee  bound-state energy, input guess and output calculated value
   real(dp), intent(out) :: ee
   !> ierr  non-zero return if error
   integer, intent(out) :: ierr

   !Local variables
   real(dp) :: amesh, al
   real(dp) :: als, cn
   real(dp) :: ro, xx
   real(dp) :: sls, sn, uout, upout
   integer :: ii, it

   real(dp), allocatable :: upp(:), cf(:)

   allocate (upp(mmax), cf(mmax))

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)

   ierr = 0

   sls = ll*(ll + 1)

   ! null arrays to remove leftover garbage

   uu(:) = 0.0d0
   up(:) = 0.0d0
   upp(:) = 0.0d0

   als = al**2

   ! coefficient array for u in differential eq.
   do ii = 1, mmax
      cf(ii) = als*sls + 2.0d0*als*(vv(ii) - ee)*rr(ii)**2
   end do

   ! start wavefunction with series

   do ii = 1, 4
      uu(ii) = rr(ii)**(ll + 1)
      up(ii) = al*(ll + 1)*rr(ii)**(ll + 1)
      upp(ii) = al*up(ii) + cf(ii)*uu(ii)
   end do

   ! outward integration using predictor once, corrector
   ! twice

   do ii = 4, mch - 1
      uu(ii + 1) = uu(ii) + aeo(up, ii)
      up(ii + 1) = up(ii) + aeo(upp, ii)
      do it = 1, 2
         upp(ii + 1) = al*up(ii + 1) + cf(ii + 1)*uu(ii + 1)
         up(ii + 1) = up(ii) + aio(upp, ii)
         uu(ii + 1) = uu(ii) + aio(up, ii)
      end do
   end do

   uout = uu(mch)
   upout = up(mch)

   !perform normalization sum

   ro = rr(1)/dsqrt(amesh)
   sn = ro**(2*ll + 3)/(2*ll + 3)

   do ii = 1, mch - 3
      sn = sn + al*rr(ii)*uu(ii)**2
   end do

   sn = sn + al*(23.0d0*rr(mch - 2)*uu(mch - 2)**2 &
   &           + 28.0d0*rr(mch - 1)*uu(mch - 1)**2 &
   &          + 9.0d0*rr(mch)*uu(mch)**2)/24.0d0

   !normalize u

   cn = 1.0d0/dsqrt(sn)
   uout = cn*uout
   upout = cn*upout

   do ii = 1, mch
      up(ii) = cn*up(ii)
      uu(ii) = cn*uu(ii)
   end do
   do ii = mch + 1, mmax
      uu(ii) = 0.0d0
   end do

   deallocate (upp, cf)

   return
end subroutine lschps
!> Finds bound states of a semi-local pseudopotential
subroutine lschpsbar(nn, ll, ierr, ee, emin, emax, rr, vv, uu, up, mmax, mbar, tht)

   implicit none
   integer, parameter :: dp = kind(1.0d0)
   real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp

   !Input variables
   !> mmax  size of log grid
   integer, intent(in) :: mmax
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> vv  local psp
   real(dp), intent(in) :: vv(mmax)
   !> emin  lower bound, potential minimum if ==0.0
   real(dp), intent(in out) :: emin
   !> emax  upper bound
   real(dp), intent(in out) :: emax
   !> tht  phast-shift angle for boundary condition, units of pi
   real(dp), intent(in) :: tht
   !> nn  principal quantum number
   integer, intent(in) :: nn
   !> ll  angular-momentum quantum number
   integer, intent(in) :: ll
   !> mbar  mesh point for infinite barrier
   integer, intent(in) :: mbar

   !Output variables
   !> uu  output radial wave function (*rr)
   real(dp), intent(out) :: uu(mmax)
   !> up  d(uu)/dr
   real(dp), intent(out) :: up(mmax)
   !> ee  bound-state energy, input guess and output calculated value (in/our)
   real(dp), intent(in out) :: ee
   !> ierr  non-zero return if error/ input for box boundary condition (0 or 1)
   integer, intent(out) :: ierr

   !Local Variables
   real(dp) :: de, uumax
   real(dp) :: eps, ro, sc, stht, ctht
   real(dp) :: sls, sn, cn, uout, upin, upout, ur, xkap, xkap2, xx
   real(dp) :: amesh, al, als
   integer :: ii, imin, it, mch, mchb, nint, node, nin
   logical :: maxset

   real(dp), allocatable :: upp(:), cf(:)
   allocate (upp(mmax), cf(mmax))

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)
   stht = sin(pi*tht)
   ctht = cos(pi*tht)

   nin = 0
   mch = 0
   mchb = 0
   ! convergence factor for solution of schroedinger eq.  if calculated
   ! correction to eigenvalue is smaller in magnitude than eps times
   ! the magnitude of the current guess, the current guess is not changed.
   eps = 1.0d-8
   ierr = 60

   sls = ll*(ll + 1)

   if (emin == 0.0d0) then
      emin = emax
      do ii = 1, mbar
         if (emin > vv(ii) + 0.5d0*sls/rr(ii)**2) then
            emin = vv(ii) + 0.5d0*sls/rr(ii)**2
            imin = ii
         end if
      end do
   end if
   ! avoid error in classical turning point search
   imin = max(imin, 2)

   if (ee < emin) ee = dmin1(emin + 5.0d0, 0.5d0*(emin + emax))
   node = 0
   de = 0.0d0
   mch = 0

   ! null arrays to remove leftover garbage
   uu(:) = 0.0d0
   up(:) = 0.0d0
   upp(:) = 0.0d0
   cf(:) = 0.0d0
   !maxset=.false.

   als = al**2
   ! return point for bound state convergence
   do nint = 1, 60

      uu(:) = 0.0d0; up(:) = 0.0d0; upp(:) = 0.0d0

      do ii = 1, mmax
         cf(ii) = als*sls + 2.0d0*als*(vv(ii) - ee)*rr(ii)**2
      end do

      ! find classical turning point for matching
      mch = mbar - 1
      if (cf(mbar) > 0.0d0) then
         do ii = mbar - 1, imin, -1
            if (cf(ii) <= 0.d0) then
               mch = ii
               exit
            end if
         end do
      end if

      de = 0.0d0
      ! start wavefunction with series

      do ii = 1, 4
         uu(ii) = rr(ii)**(ll + 1)
         up(ii) = al*(ll + 1)*rr(ii)**(ll + 1)
         upp(ii) = al*up(ii) + cf(ii)*uu(ii)
      end do

      ! outward integration using predictor once, corrector
      ! twice
      node = 0
      uumax = 0.0d0
      do ii = 4, mch - 1
         uu(ii + 1) = uu(ii) + aeo(up, ii)
         up(ii + 1) = up(ii) + aeo(upp, ii)
         do it = 1, 2
            upp(ii + 1) = al*up(ii + 1) + cf(ii + 1)*uu(ii + 1)
            up(ii + 1) = up(ii) + aio(upp, ii)
            uu(ii + 1) = uu(ii) + aio(up, ii)
         end do
         uumax = dmax1(uumax, abs(uu(ii + 1)))
      end do

      do ii = 4, mch - 1
         if (abs(uu(ii + 1)) > eps*uumax .or. abs(uu(ii)) > eps*uumax) then
            if (uu(ii + 1)*uu(ii) < 0.0d0) node = node + 1
         end if
      end do

      mchb = mch
      uout = uu(mchb)
      upout = up(mchb)

      de = 0.0d0
      if (node - nn + ll + 1 == 0) then

         ! start inward integration at 10*classical turning
         ! point with simple exponential

         nin = mch + 2.3d0/al
         if (nin < mbar - 4) then
            if (nin + 4 > mmax) nin = mmax - 4
            xkap = dsqrt(sls/rr(nin)**2 + 2.0d0*(vv(nin) - ee))

            do ii = nin, nin + 4
               uu(ii) = exp(-xkap*(rr(ii) - rr(nin)))
               up(ii) = -rr(ii)*al*xkap*uu(ii)
               upp(ii) = al*up(ii) + cf(ii)*uu(ii)
            end do

         else

            nin = mbar

            !   start inward integration at barrier
            !   point with sinh or sin

            xkap2 = sls/rr(mbar)**2 + 2.0d0*(vv(mbar) - ee)
            xkap = dsqrt(dabs(xkap2))

            do ii = mbar, mbar + 4
               xx = xkap*(rr(ii) - rr(mbar))

               if (xkap2 > 0) then
                  uu(ii) = rr(ii)*(ctht*sinh(-xx)/xkap + stht*cosh(-xx))

                  up(ii) = al*rr(ii)**2*(-ctht*cosh(-xx) - xkap*stht*sinh(-xx)) &
                  &                 + al*uu(ii)
               else
                  uu(ii) = rr(ii)*(ctht*sin(-xx)/xkap + stht*cos(-xx))

                  up(ii) = al*rr(ii)**2*(-ctht*cos(-xx) + xkap*stht*sin(-xx)) &
                  &                 + al*uu(ii)
               end if
               upp(ii) = al*up(ii) + cf(mbar)*uu(ii)
            end do

         end if !nin<mbar

         ! integrate inward

         do ii = nin, mchb + 1, -1
            uu(ii - 1) = uu(ii) + aei(up, ii)
            up(ii - 1) = up(ii) + aei(upp, ii)
            do it = 1, 2
               upp(ii - 1) = al*up(ii - 1) + cf(ii - 1)*uu(ii - 1)
               up(ii - 1) = up(ii) + aii(upp, ii)
               uu(ii - 1) = uu(ii) + aii(up, ii)
            end do
         end do

         do ii = nin + 1, mmax
            uu(ii) = 0.0d0
            up(ii) = 0.0d0
            upp(ii) = 0.0d0
         end do

         ! scale outside wf for continuity

         sc = uout/uu(mchb)

         do ii = mchb, nin
            up(ii) = sc*up(ii)
            uu(ii) = sc*uu(ii)
         end do

         upin = up(mchb)

         ! perform normalization sum

         ro = rr(1)/dsqrt(amesh)
         sn = ro**(2*ll + 3)/dfloat(2*ll + 3)
         do ii = 1, nin - 3
            sn = sn + al*rr(ii)*uu(ii)**2
         end do

         sn = sn + al*(23.0d0*rr(nin - 2)*uu(nin - 2)**2 &
         &              + 28.0d0*rr(nin - 1)*uu(nin - 1)**2 &
         &              + 9.0d0*rr(nin)*uu(nin)**2)/24.0d0

         ! normalize u

         cn = 1.0d0/dsqrt(sn)
         uout = cn*uout
         upout = cn*upout
         upin = cn*upin

         do ii = 1, nin
            up(ii) = cn*up(ii)
            uu(ii) = cn*uu(ii)
         end do
         do ii = nin + 1, mmax
            uu(ii) = 0.0d0
         end do

         ! perturbation theory for energy shift

         de = 0.5d0*uout*(upout - upin)/(al*rr(mchb))

         ! convergence test and possible exit

         if (dabs(de) < dmax1(dabs(ee), 0.2d0)*eps) then
            ierr = 0
            exit
         end if

         if (de > 0.0d0) then
            emin = ee
         else
            emax = ee
         end if
         ee = ee + de
         if (ee > emax .or. ee < emin) ee = 0.5d0*(emax + emin)
      else if (node - nn + ll + 1 < 0) then
         ! too few nodes
         emin = ee
         ee = 0.5d0*(emin + emax)
      else
         ! too many nodes
         emax = ee
         ee = 0.5d0*(emin + emax)
      end if

   end do !nint

   deallocate (upp, cf)
   return

end subroutine lschpsbar
!> integrates radial Schroedinger equation for pseudopotential
!> on a logarithmic mesh finding energy at which desired log derivative
!> uld is matched at point mch
subroutine lschpse(nn, ll, ierr, ee, uld, rr, vv, uu, up, mmax, mch)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> mmax  size of log grid
   integer, intent(in) :: mmax
   !> mch matching mesh point for inward-outward integrations
   integer, intent(in) :: mch
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> vv  semi-local pseudopotential
   real(dp), intent(in) :: vv(mmax)
   !> uld  log derivative to be matched
   real(dp), intent(in) :: uld
   !> nn  principal quantum number
   integer, intent(in) :: nn
   !> ll  angular-momentum quantum number
   integer, intent(in) :: ll

   !Output variables
   !> uu  output radial wave function (*rr)
   real(dp), intent(out) :: uu(mmax)
   !> up  d(uu)/dr
   real(dp), intent(out) :: up(mmax)
   !> ee  bound-state energy, input guess and output calculated value
   real(dp), intent(in out) :: ee
   !> ierr  non-zero return if error
   integer, intent(out) :: ierr

   !Local variables
   real(dp) :: als, cn
   real(dp) :: de, emax, emin
   real(dp) :: eps, ro, sc
   real(dp) :: amesh, al, xx
   real(dp) :: sls, sn, uout, upin, upout
   real(dp) :: xkap
   integer :: ii, it, nin, nint, node

   real(dp), allocatable :: upp(:), cf(:)
   allocate (upp(mmax), cf(mmax))

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)

   ! convergence factor for solution of schroedinger eq.  if calculated
   ! correction to eigenvalue is smaller in magnitude than eps times
   ! the magnitude of the current guess, the current guess is not changed.
   eps = 1.0d-10
   ierr = 60

   sls = ll*(ll + 1)

   emin = 0.0d0
   do ii = 1, mmax
      emin = dmin1(emin, rr(ii) + 0.5d0*sls/rr(ii)**2)
   end do
   emax = ee + 10.0d0
   emin = emin - 10.0d0

   ! null arrays to remove leftover garbage
   uu(:) = 0.0d0
   up(:) = 0.0d0
   upp(:) = 0.0d0

   als = al**2

   ! return point for bound state convergence
   do nint = 1, 60

      ! coefficient array for u in differential eq.
      do ii = 1, mmax
         cf(ii) = als*sls + 2.0d0*als*(vv(ii) - ee)*rr(ii)**2
      end do

      nin = mch

      ! start wavefunction with series

      do ii = 1, 4
         uu(ii) = rr(ii)**(ll + 1)
         up(ii) = al*(ll + 1)*rr(ii)**(ll + 1)
         upp(ii) = al*up(ii) + cf(ii)*uu(ii)
      end do

      ! outward integration using predictor once, corrector
      ! twice
      node = 0

      do ii = 4, mch - 1
         uu(ii + 1) = uu(ii) + aeo(up, ii)
         up(ii + 1) = up(ii) + aeo(upp, ii)
         do it = 1, 2
            upp(ii + 1) = al*up(ii + 1) + cf(ii + 1)*uu(ii + 1)
            up(ii + 1) = up(ii) + aio(upp, ii)
            uu(ii + 1) = uu(ii) + aio(up, ii)
         end do
         if (uu(ii + 1)*uu(ii) .le. 0.0d0) node = node + 1
      end do

      uout = uu(mch)
      upout = up(mch)

      if (node - nn + ll + 1 == 0) then

         upin = uld*uout

         ! perform normalization sum

         ro = rr(1)/dsqrt(amesh)
         sn = ro**(2*ll + 3)/dfloat(2*ll + 3)

         do ii = 1, nin - 3
            sn = sn + al*rr(ii)*uu(ii)**2
         end do

         sn = sn + al*(23.0d0*rr(nin - 2)*uu(nin - 2)**2 &
         &              + 28.0d0*rr(nin - 1)*uu(nin - 1)**2 &
         &              + 9.0d0*rr(nin)*uu(nin)**2)/24.0d0

         ! normalize u

         cn = 1.0d0/dsqrt(sn)
         uout = cn*uout
         upout = cn*upout
         upin = cn*upin

         do ii = 1, nin
            up(ii) = cn*up(ii)
            uu(ii) = cn*uu(ii)
         end do
         do ii = nin + 1, mmax
            uu(ii) = 0.0d0
         end do

         ! perturbation theory for energy shift

         de = 0.5d0*uout*(upout - upin)/(al*rr(mch))

         ! convergence test and possible exit

         if (dabs(de) < dmax1(dabs(ee), 0.2d0)*eps) then
            ierr = 0
            exit
         end if

         if (de > 0.0d0) then
            emin = ee
         else
            emax = ee
         end if
         ee = ee + de
         if (ee > emax .or. ee < emin) ee = 0.5d0*(emax + emin)

      else if (node - nn + ll + 1 < 0) then
         ! too few nodes
         emin = ee
         ee = 0.5d0*(emin + emax)

      else
         ! too many nodes
         emax = ee
         ee = 0.5d0*(emin + emax)
      end if

   end do

   deallocate (upp, cf)
   return

end subroutine lschpse
!> Finds bound states of a  pseudopotential with
!> Vanderbilt-Kleinman-Bylander non-local projectors
subroutine lschvkbb(nn, ll, nvkb, ierr, ee, emin, emax, &
&   rr, vloc, vkb, evkb, uu, up, mmax, mch)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input Variables
   !> nvkb  = number of VKB projectors to be used
   integer, intent(in) :: nvkb
   !> mmax  size of log grid
   integer, intent(in) :: mmax
   !> emin  externally generaated estimate of lower bound for ee
   real(dp), intent(in out) :: emin
   !> emax  externally generaated estimate of upper bound for ee
   real(dp), intent(in out) :: emax
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> vloc  local part of psp
   real(dp), intent(in) :: vloc(mmax)
   !> vkb  VKB projectors
   real(dp), intent(in) :: vkb(mmax, nvkb)
   !> evkb coefficients of BKB projectors
   real(dp), intent(in) :: evkb(nvkb)
   !> nn  principal quantum number
   integer, intent(in) :: nn
   !> ll  angular-momentum quantum number
   integer, intent(in) :: ll

   !Output variables
   !> uu  output radial wave function (*rr)
   real(dp), intent(out) :: uu(mmax)
   !> up  d(uu)/dr
   real(dp), intent(out) :: up(mmax)
   !> ee  bound-state energy, input guess and output calculated value
   real(dp), intent(in out) :: ee  !in/out, really - needs starting guess
   !> ierr  non-zero return if error
   integer, intent(out) :: ierr
   !> mch matching mesh point for inward-outward integrations
   integer, intent(out) :: mch

   !Local variables
   real(dp) :: cn
   real(dp) :: de
   real(dp) :: eps, ro, sc
   real(dp) :: sls, sn, uout, upin, upout
   real(dp) :: amesh, al, als
   real(dp) :: rc, xkap
   integer :: ii, it, jj, nin, nint, node

   real(dp), allocatable :: upp(:), cf(:)
   allocate (upp(mmax), cf(mmax))

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)
   node = 0

   ! convergence factor for solution of schroedinger eq.  if calculated
   ! correction to eigenvalue is smaller in magnitude than eps times
   ! the magnitude of the current guess, the current guess is not changed.
   !eps=1.0d-10
   eps = 1.0d-8
   ierr = 100

   sls = ll*(ll + 1)

   if (ee > emax) ee = 1.25d0*emax
   if (ee < emin) ee = 0.75d0*emin
   if (ee > emax) ee = 0.5d0*(emax + emin)

   ! null arrays to remove leftover garbage
   uu(:) = 0.0d0
   up(:) = 0.0d0
   upp(:) = 0.0d0

   als = al**2

   ! return point for bound state convergence
   do nint = 1, 100

      ! coefficient array for u in differential eq.
      do ii = 1, mmax
         cf(ii) = als*sls + 2.0d0*als*(vloc(ii) - ee)*rr(ii)**2
      end do

      ! find classical turning point for matching
      mch = 1
      do ii = mmax, 2, -1
         if (cf(ii - 1) <= 0.d0 .and. cf(ii) > 0.d0) then
            mch = ii
            exit
         end if
      end do

      if (mch == 1 .and. nvkb == 0) then
         write (6, '(/a)') 'lschvkbb: ERROR no classical turning point'
         write (6, '(a,i2,a,i6)') 'lschvkbb: ERROR nvkb=', nvkb, 'mch=', mch
         write (6, '(a,i2,a,f8.4,a)') 'lschvkbb: ERROR l=', ll, '  e=', ee
         !   stop
      end if

      if (nvkb > 0) then
         ! find cutoff radius for projectors
         rc = 0.0d0
         do ii = mmax, 1, -1
            if (abs(vkb(ii, 1)) > 0.0d0) then
               rc = rr(ii)
               exit
            end if
         end do

         if (mch == 1) rc = 1.25d0*rc

         ! adjust matching radius if necessary
         if (rr(mch) < rc) then
            do ii = mch, mmax
               if (rr(ii) > rc) then
                  mch = ii
                  exit
               end if
            end do
         end if
      end if
      ! outward integration
      call vkboutwf(ll, nvkb, ee, vkb, evkb, rr, vloc, uu, up, node, mmax, mch)

      uout = uu(mch)
      upout = up(mch)

      if (node - nn + ll + 1 == 0) then

         ! start inward integration at 10*classical turning
         ! point with simple exponential

         nin = mch + 2.3d0/al
         if (nin + 4 > mmax) nin = mmax - 4
         xkap = dsqrt(sls/rr(nin)**2 + 2.0d0*(vloc(nin) - ee))

         do ii = nin, nin + 4
            uu(ii) = exp(-xkap*(rr(ii) - rr(nin)))
            up(ii) = -rr(ii)*al*xkap*uu(ii)
            upp(ii) = al*up(ii) + cf(ii)*uu(ii)
         end do

         ! integrate inward

         do ii = nin, mch + 1, -1
            uu(ii - 1) = uu(ii) + aei(up, ii)
            up(ii - 1) = up(ii) + aei(upp, ii)
            do it = 1, 2
               upp(ii - 1) = al*up(ii - 1) + cf(ii - 1)*uu(ii - 1)
               up(ii - 1) = up(ii) + aii(upp, ii)
               uu(ii - 1) = uu(ii) + aii(up, ii)
            end do
         end do

         ! scale outside wf for continuity

         sc = uout/uu(mch)

         do ii = mch, nin
            up(ii) = sc*up(ii)
            uu(ii) = sc*uu(ii)
         end do

         upin = up(mch)

         ! perform normalization sum

         ro = rr(1)/dsqrt(amesh)
         sn = ro**(2*ll + 3)/dfloat(2*ll + 3)

         do ii = 1, nin - 3
            sn = sn + al*rr(ii)*uu(ii)**2
         end do

         sn = sn + al*(23.0d0*rr(nin - 2)*uu(nin - 2)**2 &
         &              + 28.0d0*rr(nin - 1)*uu(nin - 1)**2 &
         &              + 9.0d0*rr(nin)*uu(nin)**2)/24.0d0

         ! normalize u

         cn = 1.0d0/dsqrt(sn)
         uout = cn*uout
         upout = cn*upout
         upin = cn*upin

         do ii = 1, nin
            up(ii) = cn*up(ii)
            uu(ii) = cn*uu(ii)
         end do
         do ii = nin + 1, mmax
            uu(ii) = 0.0d0
         end do

         ! perturbation theory for energy shift

         de = 0.5d0*uout*(upout - upin)/(al*rr(mch))

         ! convergence test and possible exit

         if (dabs(de) < dmax1(dabs(ee), 0.2d0)*eps) then
            ierr = 0
            exit
         end if

         if (de > 0.0d0) then
            !      emin=ee
            emin = 0.5d0*(emin + ee)
         else
            !      emax=ee
            emax = 0.5d0*(emax + ee)
         end if
         ee = ee + de
         if (ee > emax .or. ee < emin) ee = 0.5d0*(emax + emin)

      else if (node - nn + ll + 1 < 0) then
         ! too few nodes
         !    emin=ee
         emin = 0.5d0*(emin + ee)
         ee = 0.5d0*(emin + emax)

      else
         ! too many nodes
         !    emax=ee
         emax = 0.5d0*(emax + ee)
         ee = 0.5d0*(emin + emax)
      end if

   end do
   !fix sign to be positive at rr->oo
   if (uu(mch) < 0.0d0) then
      uu(:) = -uu(:)
      up(:) = -up(:)
   end if

   deallocate (upp, cf)
   return

end subroutine lschvkbb
!> integrates radial schroedinger equation for pseudopotential with
!> Vanderbilt-Kleinman-Bylander non-local projectors finding energy at
!> which desired log derivative uld is matched at point mch
subroutine lschvkbbe(nn, ll, nvkb, ierr, ee, uld, emin, emax, &
&   rr, vloc, vkb, evkb, uu, up, mmax, mch)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input Variables
   !> nn  principal quantum number
   integer, intent(in) :: nn
   !> ll  angular-momentum quantum number
   integer, intent(in) :: ll
   !> nvkb  = number of VKB projectors to be used
   integer, intent(in) :: nvkb
   !> mmax  size of log grid
   integer, intent(in) :: mmax
   !> emin  externally generaated estimate of lower bound for ee
   real(dp), intent(in out) :: emin
   !> emax  externally generaated estimate of upper bound for ee
   real(dp), intent(in out) :: emax
   !> uld  bound-state energy, input guess and output calculated value
   real(dp), intent(in) :: uld
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> vloc  local part of psp
   real(dp), intent(in) :: vloc(mmax)
   !> vkb  VKB projectors
   real(dp), intent(in) :: vkb(mmax, nvkb)
   !> evkb coefficients of VKB projectors
   real(dp), intent(in) :: evkb(nvkb)

   !Output variables
   !> uu  output radial wave function (*rr)
   real(dp), intent(out) :: uu(mmax)
   !> up  d(uu)/dr
   real(dp), intent(out) :: up(mmax)
   !> ee  bound-state energy, input guess and output calculated value
   real(dp), intent(in out) :: ee  !in/out, really - needs starting guess
   !> ierr  non-zero return if error
   integer, intent(out) :: ierr
   !> mch matching mesh point for inward-outward integrations
   integer, intent(out) :: mch

   !Local variables
   real(dp) :: cn
   real(dp) :: de
   real(dp) :: eps, ro, sc
   real(dp) :: sls, sn, uout, upin, upout
   real(dp) :: amesh, al, als
   real(dp) :: rc, xkap
   integer :: ii, it, nin, nint, node

   real(dp), allocatable :: upp(:), cf(:)
   allocate (upp(mmax), cf(mmax))

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)

   ! convergence factor for solution of schroedinger eq.  if calculated
   ! correction to eigenvalue is smaller in magnitude than eps times
   ! the magnitude of the current guess, the current guess is not changed.
   eps = 1.0d-10
   ierr = 60

   sls = ll*(ll + 1)

   if (ee > emax) ee = 1.25d0*emax
   if (ee < emin) ee = 0.75d0*emin
   if (ee > emax) ee = 0.5d0*(emax + emin)

   ! null arrays to remove leftover garbage
   uu(:) = 0.0d0
   up(:) = 0.0d0
   upp(:) = 0.0d0

   als = al**2

   ! return point for bound state convergence
   do nint = 1, 60

      ! coefficient array for u in differential eq.
      do ii = 1, mmax
         cf(ii) = als*sls + 2.0d0*als*(vloc(ii) - ee)*rr(ii)**2
      end do

      nin = mch

      ! outward integration
      call vkboutwf(ll, nvkb, ee, vkb, evkb, rr, vloc, uu, up, node, mmax, mch)

      uout = uu(mch)
      upout = up(mch)

      if (node - nn + ll + 1 == 0) then

         upin = uld*uout

         ! perform normalization sum

         ro = rr(1)/dsqrt(amesh)
         sn = ro**(2*ll + 3)/dfloat(2*ll + 3)

         do ii = 1, nin - 3
            sn = sn + al*rr(ii)*uu(ii)**2
         end do

         sn = sn + al*(23.0d0*rr(nin - 2)*uu(nin - 2)**2 &
         &              + 28.0d0*rr(nin - 1)*uu(nin - 1)**2 &
         &              + 9.0d0*rr(nin)*uu(nin)**2)/24.0d0

         ! normalize u

         cn = 1.0d0/dsqrt(sn)
         uout = cn*uout
         upout = cn*upout
         upin = cn*upin

         do ii = 1, nin
            up(ii) = cn*up(ii)
            uu(ii) = cn*uu(ii)
         end do
         do ii = nin + 1, mmax
            uu(ii) = 0.0d0
         end do

         ! perturbation theory for energy shift

         de = 0.5d0*uout*(upout - upin)/(al*rr(mch))

         ! convergence test and possible exit

         if (dabs(de) < dmax1(dabs(ee), 0.2d0)*eps) then
            ierr = 0
            exit
         end if

         if (de > 0.0d0) then
            emin = ee
         else
            emax = ee
         end if
         ee = ee + de
         if (ee > emax .or. ee < emin) ee = 0.5d0*(emax + emin)

      else if (node - nn + ll + 1 < 0) then
         ! too few nodes
         emin = ee
         ee = 0.5d0*(emin + emax)

      else
         ! too many nodes
         emax = ee
         ee = 0.5d0*(emin + emax)
      end if
   end do

   deallocate (upp, cf)
   return

end subroutine lschvkbbe
!> Normalized scattering state for pseudopotential, fully non-local unless
!> ivkb=0, which should only happen if ll=lloc
subroutine lschvkbs(ll, ivkb, ee, rr, vloc, vkb, evkb, uu, up, mmax, mch)
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> mmax  size of log grid
   integer, intent(in) :: mmax
   !> mch  index of radius to which uu is computed
   integer :: mch
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> vloc  local pseudopotential
   real(dp) :: vloc(mmax)
   !> vkb  VKB projectors
   real(dp) :: vkb(mmax, 2)
   !> evkb  coefficients of VKB projectors
   real(dp) :: evkb(2)
   !> ee  scattering state energy
   real(dp) :: ee
   !> ivkb  = 0, 1 or 2 VKB proectors to be used
   integer :: ivkb
   !> ll  angular-momentum quantum number
   integer, intent(in) :: ll

   !Output variables
   !> uu  output radial wave function (*rr)
   real(dp) :: uu(mmax)
   !> up  d(uu)/dr
   real(dp) :: up(mmax)

   !Local variables
   real(dp) :: amesh, al
   real(dp) :: cn, ro, sn
   integer :: ii, node

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)

   ! null arrays to remove leftover garbage

   uu(:) = 0.0d0
   up(:) = 0.0d0

   call vkboutwf(ll, ivkb, ee, vkb, evkb, rr, vloc, uu, up, node, mmax, mch)

   !perform normalization sum

   ro = rr(1)/dsqrt(amesh)
   sn = ro**(2*ll + 3)/(2*ll + 3)

   do ii = 1, mch - 3
      sn = sn + al*rr(ii)*uu(ii)**2
   end do

   sn = sn + al*(23.0d0*rr(mch - 2)*uu(mch - 2)**2 &
   &           + 28.0d0*rr(mch - 1)*uu(mch - 1)**2 &
   &          + 9.0d0*rr(mch)*uu(mch)**2)/24.0d0

   !normalize u

   cn = 1.0d0/dsqrt(sn)

   do ii = 1, mch
      up(ii) = cn*up(ii)
      uu(ii) = cn*uu(ii)
   end do
   do ii = mch + 1, mmax
      uu(ii) = 0.0d0
   end do

   return
end subroutine lschvkbs
!> computes Vanderbilt / Kleinman-Bylander outward-integrated wave functions
subroutine vkboutwf(ll, nvkb, ep, vkb, evkb, rr, vloc, uu, up, node, mmax, mch)

   implicit none
   integer, parameter :: dp = kind(1.0d0)

   !Input variables
   !> mmax  dimension of log mesh
   integer, intent(in) :: mmax
   !> nvkb  switch for 1 or 2 projedtors
   integer ::  nvkb
   !> rr  log radial mesh
   real(dp), intent(in) :: rr(mmax)
   !> vloc  local pseudopotential
   real(dp) :: vloc(mmax)
   !> vkb  Vanderbilt-Kleinman-Bylander projectors for this l
   real(dp) :: vkb(mmax, nvkb)
   !> evkb  projector coefficients
   real(dp) :: evkb(nvkb)
   !> ep  energy at which wave function is to be calculated
   real(dp) :: ep
   !> ll  angular momentum
   integer, intent(in) :: ll
   !> mch  index of radius to which wave function is to be integrated
   integer :: mch

   !Output variables
   !> uu  wave function
   real(dp) :: uu(mmax)
   !> up  1st derivative of uu
   real(dp) :: up(mmax)
   !> node  count of number of nodes from 0 to rr(mch)
   integer node

   !Local variables
   real(dp) :: rc, rn
   real(dp), allocatable :: phi(:, :), phip(:, :), phi0(:), phi0p(:)
   real(dp), allocatable :: gg0(:), gg(:, :)
   integer, allocatable :: ipiv(:)

   integer ii, jj, kk, ierr, info

   uu(:) = 0.0d0
   up(:) = 0.0d0

   ! homogeneous solution
   call lschps(ll + 1, ll, ierr, ep, rr, vloc, uu, up, mmax, mch)

   rc = 0.0d0
   if (nvkb /= 0) then
      ! find cutoff radius for projectors
      do ii = mmax, 1, -1
         if (abs(vkb(ii, 1)) > 0.0d0) then
            rc = rr(ii)
            exit
         end if
      end do

      allocate (phi(mmax, nvkb), phip(mmax, nvkb))
      allocate (gg(nvkb, nvkb), gg0(nvkb))
      allocate (ipiv(nvkb))

      phi(:, :) = 0.0d0
      phip(:, :) = 0.0d0
      gg(:, :) = 0.0d0
      gg0(:) = 0.0d0

      ! inhomogeneous solutions
      do ii = 1, nvkb
         call lschkb(ll + 1, ll, ierr, ep, vkb(1, ii), rr, vloc, phi(1, ii), phip(1, ii), mmax, mch)
      end do

      ! projector matrix elements and coefficient matrix
      do jj = 1, nvkb
         call vpinteg(uu, vkb(1, jj), mch, 2*ll + 2, gg0(jj), rr)
         gg0(jj) = evkb(jj)*gg0(jj)
         do ii = 1, nvkb
            call vpinteg(phi(1, ii), vkb(1, jj), mch, 2*ll + 2, gg(jj, ii), rr)
            gg(jj, ii) = -evkb(jj)*gg(jj, ii)
         end do
         gg(jj, jj) = 1.0d0 + gg(jj, jj)
      end do

      ! solve linear equations for coefficients

      !    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      call dgesv(nvkb, 1, gg, nvkb, ipiv, gg0, nvkb, info)
      if (info /= 0) then
         write (6, '(/a,i4)') 'vkbout: dgesv ERROR, stopping info =', info
         stop
      end if

      ! output wave functions
      do jj = 1, nvkb
         uu(:) = uu(:) + gg0(jj)*phi(:, jj)
         up(:) = up(:) + gg0(jj)*phip(:, jj)
      end do

      deallocate (phi, phip)
      deallocate (gg, gg0)
      deallocate (ipiv)
   end if

   ! lower cutoff for node counting to avoid small-r noise
   rn = dmax1(0.10d0*rc, 0.05d0)

   node = 0
   do ii = 6, mch
      ! note historic evolution!
      !  if(rr(ii)>0.5d0 .and. uu(ii-1)*uu(ii)<0.0d0) then
      !  if(rr(ii)>0.25d0 .and. uu(ii-1)*uu(ii)<0.0d0) then
      !  if(rr(ii)>0.10d0 .and. uu(ii-1)*uu(ii)<0.0d0) then
      !  if(rr(ii)>0.07d0 .and. uu(ii-1)*uu(ii)<0.0d0) then
      if (rr(ii) > rn .and. uu(ii - 1)*uu(ii) < 0.0d0) then
         node = node + 1
      end if
   end do

   return
end subroutine vkboutwf
end module schroedinger_m

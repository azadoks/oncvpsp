!
! Copyright (c) 1989-2026 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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
! interpolates various arrays onto linear radial mesh and writes requested
! pseudopotential files

subroutine psfileout(lmax, lloc, rc, vkb, evkb, nproj, rr, vpuns, rho, rhomod, &
   irct, srel, &
   rhotae, rhoc, zz, zion, mmax, mxprj, iexc, icmod, nrl, drl, atsym, epstot, &
   na, la, ncon, nbas, nvcnf, nacnf, lacnf, nc, nv, lpopt, ncnf, &
   fa, rc0, ep, qcut, debl, facnf, dvloc0, fcfact, rcfact, &
   epsh1, epsh2, depsh, rlmax, psfile, uupsa, ea, &
   upffile, upfgrid, psp8file, psmlfile)
   use, intrinsic :: iso_fortran_env, only: stdout => output_unit
   use m_psmlout, only: psmlout
   implicit none
   integer, parameter :: dp=kind(1.0d0)

   ! Input variables
   !> Maximum angular momentum
   integer, intent(in) :: lmax
   !> l for local potential
   integer, intent(in) :: lloc
   !> Core radii
   real(dp), intent(in) :: rc(lmax + 1)
   !> VKB projectors
   real(dp), intent(in) :: vkb(mmax, mxprj, 4)
   !> Coefficients of VKB projectors
   real(dp), intent(in) :: evkb(mxprj, 4)
   !> Number of vkb projectors for each l
   integer, intent(in) :: nproj(lmax + 1)
   !> Log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> Unscreened semi-local pseudopotentials
   real(dp), intent(in) :: vpuns(mmax, 5)
   !> Valence pseudocharge
   real(dp), intent(in) :: rho(mmax)
   !> Model core charge
   real(dp), intent(in) :: rhomod(mmax, 5)
   !> All-electron valence charge
   real(dp), intent(in) :: rhotae(mmax)
   !> All-electron core charge
   real(dp), intent(in) :: rhoc(mmax)
   !> Atomic number
   real(dp), intent(in) :: zz
   !> Total valence charge
   real(dp), intent(in) :: zion
   !> Size of log radial grid
   integer, intent(in) :: mmax
   !> Maximum number of projectors
   integer, intent(in) :: mxprj
   !> Type of exchange-correlation
   integer, intent(in) :: iexc
   !> Model core charge type
   integer, intent(in) :: icmod
   !> Size of linear radial grid
   integer, intent(in) :: nrl
   !> Linear radial grid spacing
   real(dp), intent(in) :: drl
   !> Atomic symbol
   character(len=*), intent(in) :: atsym
   ! Additional input variables to be echoed
   integer, intent(in) :: na(30), la(30), ncon(6), nbas(6)
   integer, intent(in) :: nvcnf(5), nacnf(30, 5), lacnf(30, 5)
   integer, intent(in) :: nc, nv, lpopt, ncnf
   real(dp), intent(in) :: fa(30), rc0(6), ep(6), qcut(6), debl(6)
   real(dp), intent(in) :: facnf(30, 5), dvloc0, fcfact, rcfact
   real(dp), intent(in) :: epsh1, epsh2, depsh, rlmax
   character(len=*), intent(in) :: psfile
   !> UPF output file name
   character(len=*), intent(in) :: upffile
   !> UPF grid/mesh type
   character(len=*), intent(in) :: upfgrid
   !> PSP8 output file name
   character(len=*), intent(in) :: psp8file
   !> PSML output file name
   character(len=*), intent(in) :: psmlfile
   !> Pseudo-atomic valence wavefunctions
   real(dp), intent(in) :: uupsa(mmax, nv)
   !> Eigenvalues
   real(dp), intent(in) :: ea(30)
   !> Index where rhoc is matched
   integer, intent(in) :: irct
   !> Scalar-relativistic flag
   logical, intent(in) :: srel
   !> Pseudo-atom total energy
   real(dp), intent(in) :: epstot

   ! Local variables
   integer :: ii, jj, iprj, ll, l1
   integer :: unit, ios
   integer :: nrl_upf
   real(dp) :: al, uurcut, nrmsum
   real(dp), allocatable :: rrl(:)
   real(dp), allocatable :: drrl(:)
   real(dp), allocatable :: vkbl(:, :, :)
   real(dp), allocatable :: vpunsl(:, :)
   real(dp), allocatable :: rhotael(:)
   real(dp), allocatable :: rhocl(:)
   real(dp), allocatable :: rhol(:)
   real(dp), allocatable :: rhomodl(:, :)
   real(dp), allocatable :: uuaeal(:, :)
   real(dp), allocatable :: uupsal(:, :)

   if (trim(psfile) == 'psp8' .or. trim(psfile) == 'both') then
      allocate(rrl(nrl), drrl(nrl), vkbl(nrl, mxprj, 4), vpunsl(nrl, 5), &
         rhotael(nrl), rhocl(nrl), rhol(nrl), rhomodl(nrl, 5), &
         uuaeal(nrl, nc + nv), uupsal(nrl, nv))
      call linear_interpolate(nc, nv, lmax, lloc, la, mxprj, nproj, icmod, mmax, &
         nrl, drl, rr, vkb, vpuns, rhotae, rhoc, rho, rhomod, uuaea, uupsa, &
         rrl, drrl, vkbl, vpunsl, rhotael, rhocl, rhol, rhomodl, uuaeal, uupsal)
      if (trim(psp8file) /= "stdout") then
         open(newunit=unit, file=psp8file, status='new', action='write', iostat=ios)
      else
         unit = stdout
      end if
      call linout(lmax,lloc,rc,vkbl,evkb,nproj,rr,vpunsl,rhol,rhomodl, &
         rhotael,rhocl,zz,zion,mmax,mxprj,iexc,icmod,nrl,drl,atsym, &
         na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
         fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact, &
         epsh1,epsh2,depsh,rlmax,psfile)
      if (trim(psp8file) /= "stdout") then
         close(unit)
      end if
      deallocate(rrl, drrl, vkbl, vpunsl, rhotael, rhocl, rhol, rhomodl, uuaeal, uupsal)
   end if

   if (trim(psfile) == 'upf' .or. trim(psfile) == 'both') then
      if (trim(upfgrid) == 'log') then
         allocate(drrl(mmax))
         do ii = 1, mmax - 1
            drrl(ii) = ii * rr(1) * exp(0.01_dp * log(rr(101) / rr(1)) * (ii - 1))
         end do
         call upfout(unit,lmax,lloc,rc,vkb,evkb,nproj,rr,drrl,vpuns,rho,rhomod, &
            zz,zion,mmax,mxprj,iexc,icmod,drl,atsym,epstot, &
            na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
            fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact, &
            epsh1,epsh2,depsh,rlmax,psfile,uupsa,ea,upfgrid)
         deallocate(drrl)
      else
         al = log(rr(2) / rr(1))
         uurcut = 0.0_dp
         do ii = 1, nv
            nrmsum = 0.0_dp
            do jj = mmax - 1, 1, -1
               nrmsum = nrmsum + (uupsa(jj, ii)**2) * rr(jj) * al
               if (nrmsum > 1.0e-6_dp) then
                  exit  ! Found cutoff radius such that uu norm accurate to 10^-6
               end if
            end do
            if (rr(jj) > uurcut) uurcut = rr(jj)
         end do
         if (uurcut > drl * real(nrl - 1, dp)) then
            nrl_upf = 1 + int(uurcut / drl)
            if (mod(nrl_upf, 2) /= 0) nrl_upf = nrl_upf + 1
            write(6, '(a,i6,a,f10.5)') "Updating nrl = ", nrl_upf, " for uurcut = ", uurcut
         else
            nrl_upf = nrl
         end if
         allocate(rrl(nrl), drrl(nrl), vkbl(nrl, mxprj, 4), vpunsl(nrl, 5), &
            rhotael(nrl), rhocl(nrl), rhol(nrl), rhomodl(nrl, 5), &
            uuaeal(nrl, nc + nv), uupsal(nrl, nv))
         call linear_interpolate(nc, nv, lmax, lloc, la, mxprj, nproj, icmod, mmax, &
            nrl, drl, rr, vkb, vpuns, rhotae, rhoc, rho, rhomod, uuaeal, uupsal, &
            rrl, drl, vkbl, vpunsl, rhotael, rhocl, rhol, rhomodl, uuaeal, uupsal)
         if (trim(upffile) /= "stdout") then
            open(newunit=unit, file=upffile, status='new', action='write', iostat=ios)
         else
            unit = stdout
         end if
         call upfout(unit,lmax,lloc,rc,vkb,evkb,nproj,rr,drrl,vpuns,rho,rhomod, &
            zz,zion,mmax,mxprj,iexc,icmod,drl,atsym,epstot, &
            na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
            fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact, &
            epsh1,epsh2,depsh,rlmax,psfile,uupsa,ea,upfgrid)
         if (trim(upffile) /= "stdout") then
            close(unit)
         end if
         deallocate(rrl, drrl, vkbl, vpunsl, rhotael, rhocl, rhol, rhomodl, uuaeal, uupsal)
      end if
   end if

   if (trim(psfile) == 'psml' .or. trim(psfile) == 'both') then
      call psmlout(lmax,lloc,rc,vkb,evkb,mxprj,nproj,rr,vpuns,rho,rhomod, &
         irct, srel, &
         zz,zion,mmax,iexc,icmod,nrl,drl,atsym,epstot, &
         na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
         fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact, &
         epsh1,epsh2,depsh,rlmax,psfile)
   end if
end subroutine psfileout

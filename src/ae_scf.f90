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
! Convergence measure and energy terms of the all-electron SCF of sratom
! and relatom, which mix the screening (Hartree plus exchange-correlation)
! potential and add the nuclear -zz/r to it for the radial equations

 subroutine ae_residual(rho,vsi,vso,sf,rr,mmax,res)

!! root mean square over the electrons of the change of the screening
!! potential from input to output, the size of the first-order change of
!! the eigenvalues (the convergence measure of RadialAtomicSolvers.jl)
!!
!! @param rho  total charge density (4 pi r**2 rho integrates to sf)
!! @param vsi  input screening potential
!! @param vso  output screening potential
!! @param sf  number of electrons
!! @param rr  log radial mesh
!! @param mmax  dimension of log radial grid
!! @param res  residual (Ha)

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rho(mmax),vsi(mmax),vso(mmax),rr(mmax)
 real(dp), intent(in) :: sf

!Output variables
 real(dp), intent(out) :: res

!Local variables
 real(dp) :: al,sn
 integer :: ii

 al = 0.01d0 * log(rr(101) / rr(1))

 sn=0.0d0
 do ii=1,mmax
   sn=sn + rho(ii)*rr(ii)**3*(vso(ii)-vsi(ii))**2
 end do
 res=sqrt(al*sn/sf)

 return
 end subroutine ae_residual

 subroutine ae_energies(rho,vsi,vxc,zz,eeig,eeel,eexc,rr,mmax,etot,eterm)

!! total energy and its terms as RadialAtomicSolvers.jl and NIST define
!! them, the kinetic energy from the eigenvalue sum less the potential
!! energy in the input potential the orbitals are eigenstates of, so that
!! the total energy is stationary in that potential
!!
!! @param rho  total charge density of the orbitals of vsi
!! @param vsi  input screening potential
!! @param vxc  exchange-correlation potential of rho (vout)
!! @param zz  atomic number
!! @param eeig  sum of the occupied eigenvalues
!! @param eeel  electron-electron interaction of rho (vout)
!! @param eexc  exchange-correlation correction of rho (vout)
!! @param rr  log radial mesh
!! @param mmax  dimension of log radial grid
!! @param etot  total energy (Ha)
!! @param eterm  kinetic, electron-nucleus, Hartree and exchange-correlation
!!               energies (Ha)

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rho(mmax),vsi(mmax),vxc(mmax),rr(mmax)
 real(dp), intent(in) :: zz,eeig,eeel,eexc

!Output variables
 real(dp), intent(out) :: etot,eterm(4)

!Local variables
 real(dp) :: al,ekin,enuc,ecoul,exc,escr
 integer :: ii

 real(dp), allocatable :: ff(:)

 allocate(ff(mmax))

 al = 0.01d0 * log(rr(101) / rr(1))

! electron-nucleus energy; the integrand goes as rr**2 below rr(1)
 do ii=1,mmax
   ff(ii)=-zz*rho(ii)*rr(ii)**2
 end do
 enuc=radint(ff,2.0d0)

! potential energy in the input screening potential
 do ii=1,mmax
   ff(ii)=rho(ii)*vsi(ii)*rr(ii)**3
 end do
 escr=radint(ff,3.0d0)

! exchange-correlation energy from vout's correction to the eigenvalue sum
 do ii=1,mmax
   ff(ii)=rho(ii)*vxc(ii)*rr(ii)**3
 end do
 exc=eexc + radint(ff,3.0d0)

 ecoul=0.5d0*eeel
 ekin=eeig - enuc - escr

 etot=ekin + enuc + ecoul + exc
 eterm(1)=ekin
 eterm(2)=enuc
 eterm(3)=ecoul
 eterm(4)=exc

 deallocate(ff)
 return

 contains

 function radint(ff,pw)

!! integral of ff/rr over rr on the log mesh, by vout's rule; the part
!! below rr(1) is that of an integrand going as rr**(pw-1)

!Input variables
 real(dp), intent(in) :: ff(mmax)
 real(dp), intent(in) :: pw

!Output variables
 real(dp) :: radint

 radint=(9.0d0*ff(1) + 28.0d0*ff(2) + 23.0d0*ff(3))/24.0d0
 do ii=4,mmax
   radint=radint + ff(ii)
 end do
 radint=al*radint + ff(1)/pw

 return
 end function radint

 end subroutine ae_energies

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

 module aeatom_m

 implicit none
 private
 public :: aeatom

 contains

 subroutine aeatom(irel,c_light,nkap,na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho, &
&                  rr,vi,zz,mmax,iexc,etot,eterm,res,ierr,vext,warm,iout,ifail)

!! self-consistent all-electron atom: non-relativistic or scalar-relativistic
!! (sratom) or fully-relativistic (relatom)
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param nkap  number of kappa per level: 1, or 2 for irel=2
!! @param na  principal quantum number array, dimension ncv
!! @param la  angular-momenta
!! @param ea  eigenvalues, input guess (irel=2 or warm) and output;
!!            ea(:,1) for j=l+1/2, ea(:,2) for j=l-1/2 when irel=2
!! @param fa  occupancies
!! @param rpk  radius of outermost peak of wave function, indexed as ea
!! @param nc  number of core states
!! @param ncv  number of core+valence states
!! @param it  number of iterations (output)
!! @param rhoc  core charge density (output)
!! @param rho  total charge density (output)
!! @param rr  log radial mesh
!! @param vi  all-electron potential (output; input with warm)
!! @param zz  atomic number
!! @param mmax  size of log grid
!! @param iexc  exchange-correlation function to be used
!! @param etot  all-electron total energy (output)
!! @param eterm  kinetic, electron-nucleus, Hartree and exchange-correlation
!!               energies (output)
!! @param res  residual of the screening potential (output, ae_residual)
!! @param ierr  error flag
!! @param vext  external potential the orbitals also see (a confining wall,
!!              oncvpsp_confine.x), zero for the free atom; vi is without it
!! @param warm  start from vi and ea (the solution of a nearby potential,
!!              oncvpsp_confine.x) instead of the Thomas-Fermi potential
!!              and hydrogenic energies
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp, mxstate
 use relatom_m, only: relatom
 use sratom_m, only: sratom
 implicit none

!Input variables
 integer, intent(in) :: irel,nkap,mmax,iexc,nc,ncv
 real(dp), intent(in) :: c_light
 integer, intent(in) :: na(ncv),la(ncv)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: fa(ncv),rr(mmax),vext(mmax)

!Input/output variables
 integer, intent(inout) :: ierr
 real(dp), intent(inout) :: ea(mxstate,nkap),rpk(mxstate,nkap)

!Output variables
 integer, intent(out) :: it
 real(dp), intent(out) :: etot,eterm(4),res
 real(dp), intent(out) :: rho(mmax),rhoc(mmax)
 real(dp), intent(inout) :: vi(mmax)
 logical, intent(in) :: warm
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

 ifail=0
 if(irel==2) then
   call relatom(na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho,rr,vi,zz,mmax,iexc,c_light, &
&               etot,eterm,res,ierr,vext,warm,iout,ifail)
   if(ifail/=0) return
 else
   call sratom(na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho,rr,vi,zz,mmax,iexc,etot, &
&              eterm,res,ierr,irel,c_light,vext,warm,iout,ifail)
   if(ifail/=0) return
 end if
 return
 end subroutine aeatom

 end module aeatom_m

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
 subroutine aeatom(irel,nkap,na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho, &
&                  rr,vi,zz,mmax,iexc,etot,ierr)

!self-consistent all-electron atom: non-relativistic or scalar-relativistic
!(sratom) or fully-relativistic (relatom)

!irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!nkap  number of kappa per level: 1, or 2 for irel=2
!na  principal quantum number array, dimension ncv
!la  angular-momenta
!ea  eigenvalues, input guess (irel=2) and output; ea(:,1) for j=l+1/2,
!    ea(:,2) for j=l-1/2 when irel=2
!fa  occupancies
!rpk  radius of outermost peak of wave function, indexed as ea
!nc  number of core states
!ncv  number of core+valence states
!it  number of iterations (output)
!rhoc  core charge density (output)
!rho  total charge density (output)
!rr  log radial mesh
!vi  all-electron potential (output)
!zz  atomic number
!mmax  size of log grid
!iexc  exchange-correlation function to be used
!etot  all-electron total energy (output)
!ierr  error flag

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: irel,nkap,mmax,iexc,nc,ncv
 integer, intent(in) :: na(ncv),la(ncv)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: fa(ncv),rr(mmax)

!Input/output variables
 integer, intent(inout) :: ierr
 real(dp), intent(inout) :: ea(30,nkap),rpk(30,nkap)

!Output variables
 integer, intent(out) :: it
 real(dp), intent(out) :: etot
 real(dp), intent(out) :: rho(mmax),rhoc(mmax),vi(mmax)

 if(irel==2) then
   call relatom(na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho, &
&               rr,vi,zz,mmax,iexc,etot,ierr)
 else
   call sratom(na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho, &
&              rr,vi,zz,mmax,iexc,etot,ierr,irel)
 end if
 return
 end subroutine aeatom

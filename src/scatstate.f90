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
 subroutine scatstate(irel,nn,ll,kap,ierr,ee,rr,zz,vv,uu,up,mmax,mch)

!Integrates a scattering state of an all-electron atomic potential out to
!mch: with the non-relativistic or scalar-relativistic Schroedinger equation
!(lschfs), or with the Dirac equation (ldiracfs)

!irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!nn  effective principal quantum number based on nodes inside mch (output)
!ll  angular-momentum quantum number
!kap  =l, -(l+1) for j=l -/+ 1/2 (irel=2 only)
!ierr  non-zero return if error
!ee  energy of the state
!rr  log radial mesh
!zz  atomic number
!vv  local atomic potential
!uu  output radial wave function (*rr): uu(:,1) the Schroedinger function,
!    or the Dirac large and small components
!up  d(uu)/dr
!mmax  size of log grid
!mch  matching mesh point for inward-outward integrations

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: irel,ll,kap,mmax,mch
 real(dp), intent(in) :: rr(mmax),vv(mmax)
 real(dp), intent(in) :: zz,ee

!Input/output variables
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)
 integer, intent(inout) :: nn

!Output variables
 integer, intent(out) :: ierr

 if(irel==2) then
   call ldiracfs(nn,ll,kap,ierr,ee,rr,zz,vv,uu,up,mmax,mch)
 else
   call lschfs(nn,ll,ierr,ee,rr,vv,uu,up,zz,mmax,mch,irel)
 end if
 return
 end subroutine scatstate

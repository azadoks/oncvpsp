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
 subroutine boundstate(irel,nn,ll,kap,ierr,ee,rr,zz,vv,uu,up,mmax,mch)

!! Finds a bound state of an all-electron atomic potential: with the
!! non-relativistic or scalar-relativistic Schroedinger equation (lschfb), or
!! with the Dirac equation (ldiracfb)
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param nn  principal quantum number
!! @param ll  angular-momentum quantum number
!! @param kap  =l, -(l+1) for j=l -/+ 1/2 (irel=2 only)
!! @param ierr  non-zero return if error
!! @param ee  bound-state energy, input guess and output calculated value
!! @param rr  log radial mesh
!! @param zz  atomic number
!! @param vv  local atomic potential
!! @param uu  output radial wave function (*rr): uu(:,1) the Schroedinger
!!            function, or the Dirac large and small components
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of log grid
!! @param mch  matching mesh point for inward-outward integrations

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: irel,nn,ll,kap,mmax
 real(dp), intent(in) :: rr(mmax),vv(mmax)
 real(dp), intent(in) :: zz

!Input/output variables
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)
 real(dp), intent(inout) :: ee
 integer, intent(inout) :: mch

!Output variables
 integer, intent(out) :: ierr

 if(irel==2) then
   call ldiracfb(nn,ll,kap,ierr,ee,rr,zz,vv,uu,up,mmax,mch)
 else
   call lschfb(nn,ll,ierr,ee,rr,vv,uu,up,zz,mmax,mch,irel)
 end if
 return
 end subroutine boundstate

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

 module scatterstate_m

 implicit none
 private
 public :: scatterstate

 contains

 subroutine scatterstate(irel,c_light,nn,ll,kap,ee,rr,zz,vv,uu,up,mmax,adams, &
&                        mch)

!! Integrates a scattering state of an all-electron atomic potential out to
!! mch: with the non-relativistic or scalar-relativistic Schroedinger equation
!! (lschfs), or with the Dirac equation (ldiracfs)
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param nn  effective principal quantum number based on nodes inside mch
!!            (output)
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2 (irel=2 only)
!! @param ee  energy of the state
!! @param rr  log radial mesh (bohr)
!! @param zz  atomic number
!! @param vv  local atomic potential
!! @param uu  output radial wave function (*rr): uu(:,1) the Schroedinger
!!            function, or the Dirac large and small components
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param mch  matching mesh point for inward-outward integrations

 use constants_m, only: dp
 use ldiracfs_m, only: ldiracfs
 use lschfs_m, only: lschfs
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: irel,ll,kap,mmax,mch
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: rr(mmax),vv(mmax)
 real(dp), intent(in) :: zz,ee

!Input/output variables
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)
 integer, intent(inout) :: nn

 if(irel==2) then
   call ldiracfs(nn,ll,kap,c_light,ee,rr,zz,vv,uu,up,mmax,adams,mch)
 else
   call lschfs(nn,ll,ee,rr,vv,uu,up,zz,mmax,adams,mch,irel,c_light)
 end if
 return
 end subroutine scatterstate

 end module scatterstate_m

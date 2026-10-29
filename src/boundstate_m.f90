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

 module boundstate_m

 implicit none
 private
 public :: boundstate

 contains

 subroutine boundstate(irel,c_light,nn,ll,kap,ierr,ee,rr,zz,vv,uu,up,mmax, &
&      adams,mch, &
&      iout,ifail)

!! Finds a bound state of an all-electron atomic potential: with the
!! non-relativistic or scalar-relativistic Schroedinger equation (lschfb), or
!! with the Dirac equation (ldiracfb)
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param nn  principal quantum number
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2 (irel=2 only)
!! @param ierr  non-zero return if error
!! @param ee  bound-state energy, input guess and output calculated value
!! @param rr  log radial mesh (bohr)
!! @param zz  atomic number
!! @param vv  local atomic potential
!! @param uu  output radial wave function (*rr): uu(:,1) the Schroedinger
!!            function, or the Dirac large and small components
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param mch  matching mesh point for inward-outward integrations
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: dp
 use ldiracfb_m, only: ldiracfb
 use lschfb_m, only: lschfb
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: irel,nn,ll,kap,mmax
 type(adams_t), intent(in) :: adams
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: rr(mmax),vv(mmax)
 real(dp), intent(in) :: zz
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)
 real(dp), intent(inout) :: ee
 integer, intent(inout) :: mch

!Output variables
 integer, intent(out) :: ierr
 integer, intent(out) :: ifail

 ifail=0
 if(irel==2) then
   call ldiracfb(c_light,nn,ll,kap,ierr,ee,rr,zz,vv,uu,up,mmax,adams,mch,iout, &
&                ifail)
   if(ifail/=0) return
 else
   call lschfb(irel,c_light,nn,ll,ierr,ee,rr,zz,vv,uu,up,mmax,adams,mch,iout, &
&              ifail)
   if(ifail/=0) return
 end if
 return
 end subroutine boundstate

 end module boundstate_m

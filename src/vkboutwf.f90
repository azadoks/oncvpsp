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
 subroutine vkboutwf(ll,nvkb,ep,vkbpw,evkb,rr,vloc,vlocpw,uu,up,node,mmax, &
&                    mch,mnode)

!! computes Vanderbilt / Kleinman-Bylander outward-integrated wave functions:
!! the homogeneous solution of the local potential plus the inhomogeneous
!! ones of the projectors with the coefficients that make the projections
!! self-consistent, each segment between the kinks of the potential and
!! the projectors integrated on its own grid (outwf_m, doc/
!! smoothness_plan.md R-10)
!!
!! @param ll  angular momentum
!! @param nvkb  number of projectors (0 for the local potential alone)
!! @param ep  energy at which wave function is to be calculated
!! @param vkbpw  Vanderbilt-Kleinman-Bylander projectors for this l, as
!!               functions of pieces
!! @param evkb  projector coefficients
!! @param rr  log radial mesh
!! @param vloc  local pseudopotential on the mesh (the screened local
!!              potential, or that plus a smooth screening change)
!! @param vlocpw  the screened local potential's pieces
!! @param uu  wave function
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param node  count of number of nodes from 0 to rr(mnode)
!! @param mmax  dimension of log mesh
!! @param mch  index of radius to which wave function is to be integrated
!! @param mnode  last point of the node count, <= mch

 use constants_m, only: dp
 use pieces_m, only: piecewise_t
 use outwf_m, only: outwf_t, outwf_setup, outwf_solve, outwf_free
 implicit none

!Input variables
 real(dp), intent(in) :: rr(mmax),vloc(mmax),evkb(nvkb)
 real(dp), intent(in) :: ep
 integer, intent(in) :: nvkb,ll,mmax,mch,mnode
 type(piecewise_t), intent(in) :: vkbpw(nvkb),vlocpw

!Output variables
 real(dp), intent(out) :: uu(mmax),up(mmax)
 integer, intent(out) :: node

!Local variables
 type(outwf_t) :: ow

 call outwf_setup(ll,rr,mmax,vloc,vlocpw,nvkb,vkbpw,ow)
 call outwf_solve(ow,ep,evkb,mch,mnode,uu,up,node)
 call outwf_free(ow)
 return
 end subroutine vkboutwf

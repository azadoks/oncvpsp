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
!> The polynomial local potential (lloc=4) of vploc, at a continuous
!> matching radius (doc/smoothness_plan.md, R-7): the polynomial of lpopt
!> matches the all-electron potential and its first two or three
!> derivatives at rc, which lm_derivs gives at any rc (the Floater-Hormann
!> interpolant of fixed 11-point stencils, smooth in rc; for V = -2/r +
!> exp(-r) at 1.7 bohr on the default mesh its V''' is good to 1e-9, where
!> a least-squares polynomial of degree 8 on the nodes within 0.25 bohr,
!> weighted by a bump, was good to 2e-6: the truncation of a fit of fixed
!> width is the larger error for a potential with the derivatives of -Z/r).
!> The polynomial is evaluated analytically at any radius (its extension
!> past rc serves the integrals of pieces_m).

 module vploc_m
 use constants_m, only: dp
 use logmesh_m, only: lm_derivs
 implicit none
 private
 public :: vploc_t, vploc_fit, vploc_eval, vploc_mesh

 type vploc_t
   integer :: lpopt=0
   real(dp) :: rc=0.0d0,dvloc0=0.0d0
!> the all-electron potential and its first three derivatives at rc
   real(dp) :: vder(0:3)=0.0d0
!> the polynomial's coefficients: aco + bco r**p1 + cco r**p2 + dco r**p3
   real(dp) :: aco=0.0d0,bco=0.0d0,cco=0.0d0,dco=0.0d0
 end type vploc_t

 contains

 subroutine vploc_fit(rr,vv,mmax,rc,dvloc0,lpopt,cf)

!! the polynomial potential matched at rc
!!
!! @param rr  log radial mesh
!! @param vv  all-electron potential
!! @param mmax  size of the log radial mesh
!! @param rc  matching radius (bohr)
!! @param dvloc0  amplitude at rr==0 of the smooth "null space" function added
!! @param lpopt  1-5, the polynomial:
!!               1) match 2 derivatives, r**0,2,4
!!               2) match 2 derivatives, r**0,4,6
!!               3) match 3 derivatives, r**0,4,5,6
!!               4) match 3 derivatives, r**0,4,6,8
!!               5) match 3 derivatives, r**0,2,4,6
!! @param cf  the potential's coefficients

!Input variables
 integer, intent(in) :: mmax,lpopt
 real(dp), intent(in) :: rr(mmax),vv(mmax),rc,dvloc0

!Output variables
 type(vploc_t), intent(out) :: cf

!Local variables
 real(dp) :: vd(4),v0,dvv,d2vv,d3vv

 call lm_derivs(rr,vv,mmax,rc,3,vd)
 cf%vder(0:3)=vd(1:4)
 v0=cf%vder(0)
 dvv=cf%vder(1)
 d2vv=cf%vder(2)
 d3vv=cf%vder(3)

 cf%lpopt=lpopt
 cf%rc=rc
 cf%dvloc0=dvloc0
 cf%dco=0.0d0
 if(lpopt == 1) then
   cf%bco=(3*dvv-d2vv*rc)/(4*rc)
   cf%cco=( -dvv+d2vv*rc)/(8*rc**3)
   cf%aco=v0-cf%bco*rc**2-cf%cco*rc**4
 else if(lpopt == 2) then
   cf%bco=( 5*dvv-d2vv*rc)/(8*rc**3)
   cf%cco=(-3*dvv+d2vv*rc)/(12*rc**5)
   cf%aco=v0-cf%bco*rc**4-cf%cco*rc**6
 else if(lpopt == 3) then
   cf%bco=( 20*dvv-8*d2vv*rc+d3vv*rc**2)/(8*rc**3)
   cf%cco=(-15*dvv+7*d2vv*rc-d3vv*rc**2)/(5*rc**4)
   cf%dco=( 12*dvv-6*d2vv*rc+d3vv*rc**2)/(12*rc**5)
   cf%aco=v0-cf%bco*rc**4-cf%cco*rc**5-cf%dco*rc**6
 else if(lpopt == 4) then
   cf%bco=(35*dvv-11*d2vv*rc+d3vv*rc**2)/(32*rc**3)
   cf%cco=(-21*dvv+9*d2vv*rc-d3vv*rc**2)/(24*rc**5)
   cf%dco=(15*dvv-7*d2vv*rc+d3vv*rc**2)/(64*rc**7)
   cf%aco=v0-cf%bco*rc**4-cf%cco*rc**6-cf%dco*rc**8
 else if(lpopt == 5) then
   cf%bco=(15*dvv-7*d2vv*rc+d3vv*rc**2)/(16*rc)
   cf%cco=(-5*dvv+5*d2vv*rc-d3vv*rc**2)/(16*rc**3)
   cf%dco=( 3*dvv-3*d2vv*rc+d3vv*rc**2)/(48*rc**5)
   cf%aco=v0-cf%bco*rc**2-cf%cco*rc**4-cf%dco*rc**6
 else
   write(6,'(a,i0)') 'vploc_fit: ERROR lpopt ',lpopt
   error stop 1
 end if
 return
 end subroutine vploc_fit

 real(dp) function vploc_eval(cf,rq)

!! the polynomial potential at any radius: the matched polynomial plus the
!! supplementary function of dvloc0
!!
!! @param cf  the coefficients
!! @param rq  radius (bohr)

!Input variables
 type(vploc_t), intent(in) :: cf
 real(dp), intent(in) :: rq

!Local variables
 real(dp) :: xx

 select case(cf%lpopt)
 case(1)
   vploc_eval=cf%aco+cf%bco*rq**2+cf%cco*rq**4
 case(2)
   vploc_eval=cf%aco+cf%bco*rq**4+cf%cco*rq**6
 case(3)
   vploc_eval=cf%aco+cf%bco*rq**4+cf%cco*rq**5+cf%dco*rq**6
 case(4)
   vploc_eval=cf%aco+cf%bco*rq**4+cf%cco*rq**6+cf%dco*rq**8
 case(5)
   vploc_eval=cf%aco+cf%bco*rq**2+cf%cco*rq**4+cf%dco*rq**6
 case default
   vploc_eval=0.0d0
 end select
 xx=rq/cf%rc
 select case(cf%lpopt)
 case(1)
   vploc_eval=vploc_eval + cf%dvloc0*(1.0d0 - xx**2)**3
 case(2)
   vploc_eval=vploc_eval + cf%dvloc0*(1.0d0 - xx**4)**3
 case(3,4,5)
   vploc_eval=vploc_eval + cf%dvloc0*(1.0d0 - xx**4)**4
 case default
 end select
 return
 end function vploc_eval

 subroutine vploc_mesh(cf,rr,vv,mmax,vp)

!! the local potential on the mesh: the polynomial at the nodes up to rc,
!! the all-electron potential beyond
!!
!! @param cf  the coefficients
!! @param rr  log radial mesh
!! @param vv  all-electron potential
!! @param mmax  size of the log radial mesh
!! @param vp  the local potential

!Input variables
 integer, intent(in) :: mmax
 type(vploc_t), intent(in) :: cf
 real(dp), intent(in) :: rr(mmax),vv(mmax)

!Output variables
 real(dp), intent(out) :: vp(mmax)

!Local variables
 integer :: ii

 do ii=1,mmax
   if(rr(ii)<=cf%rc) then
     vp(ii)=vploc_eval(cf,rr(ii))
   else
     vp(ii)=vv(ii)
   end if
 end do
 return
 end subroutine vploc_mesh

 end module vploc_m

!
! Copyright (c) 2026 by Austin Zadoks
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

 module ae_delta_m

 implicit none
 private
 public :: ae_delta

 contains

 subroutine ae_delta(irel,c_light,ll,kap,ee,rr,vv,uu,up,mmax,i1,i2,delta)

!! the relativistic term of the all-electron radial equation at the mesh
!! points i1..i2: Delta = u'' - [l(l+1)/r**2 + 2(V-E)] u, the difference of
!! the exact second derivative of the large component from the Schroedinger
!! one (doc/smoothness_plan.md, D-4, R-11; PseudoPotentialGenerators.jl
!! relativistic_taper). u'' comes from the first-order system the solvers
!! integrate, with u and u' = up/(al r) at the point and V' there from
!! lm_derivs: no stencil on u, no division by u (vrel divided by u and
!! needed a threshold at its nodes)
!!
!! @param irel  0 non-relativistic (Delta = 0), 1 scalar-relativistic
!!              (the flux form of sr_flux_m), 2 fully-relativistic (Dirac,
!!              uu the large component)
!! @param c_light  speed of light (Hartree atomic units)
!! @param ll  angular-momentum quantum number
!! @param kap  =l, -(l+1) for j=l -/+ 1/2 (irel=2 only)
!! @param ee  energy of uu (Ha)
!! @param rr  log radial mesh
!! @param vv  the potential uu solves the equation in (Ha)
!! @param uu  the wave function (*rr)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of the log radial mesh
!! @param i1,i2  the mesh points wanted
!! @param delta  Delta at the points i1..i2, in delta(i1:i2)

 use constants_m, only: dp
 use logmesh_m, only: lm_derivs

!Input variables
 integer, intent(in) :: irel,ll,kap,mmax,i1,i2
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: ee
 real(dp), intent(in) :: rr(mmax),vv(mmax),uu(mmax),up(mmax)

!Output variables
 real(dp), intent(out) :: delta(mmax)

!Local variables
 real(dp) :: al,cc,fss,sls,vd(2),rq,u0,u1,u2,xm,xmp,qq,qp,aa,ap,ff,fp
 integer :: ii

 delta(i1:i2)=0.0d0
 if(irel==0) return

 al=log(rr(mmax)/rr(1))/dble(mmax-1)
 cc=c_light
 fss=1.0d0/cc**2
 sls=dble(ll*(ll+1))
 do ii=i1,i2
   rq=rr(ii)
   call lm_derivs(rr,vv,mmax,rq,1,vd)
   u0=uu(ii)
   u1=up(ii)/(al*rq)
   if(irel==1) then
! flux form: u'=u/r + M q, q'=(l(l+1)/(M r**2) + 2(V-E)) u - q/r,
! M = 1 + (E-V)/(2c**2)
     xm=1.0d0 + 0.5d0*fss*(ee - vd(1))
     xmp=-0.5d0*fss*vd(2)
     qq=(u1 - u0/rq)/xm
     qp=(sls/(xm*rq**2) + 2.0d0*(vd(1)-ee))*u0 - qq/rq
     u2=u1/rq - u0/rq**2 + xmp*qq + xm*qp
   else
! Dirac: g'=-kap g/r + A f, f'=kap f/r + (V-E) g/c, A=(2c**2 + E - V)/c
     aa=(2.0d0*cc**2 + ee - vd(1))/cc
     ap=-vd(2)/cc
     ff=(u1 + dble(kap)*u0/rq)/aa
     fp=dble(kap)*ff/rq + (vd(1)-ee)*u0/cc
     u2=dble(kap)*u0/rq**2 - dble(kap)*u1/rq + ap*ff + aa*fp
   end if
   delta(ii)=u2 - (sls/rq**2 + 2.0d0*(vd(1)-ee))*u0
 end do
 return
 end subroutine ae_delta

 end module ae_delta_m

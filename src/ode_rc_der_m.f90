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

 module ode_rc_der_m

 implicit none
 private
 public :: ode_rc_der

 contains

 subroutine ode_rc_der(irel,c_light,ll,kap,ee,rr,vv,uu,up,mmax,rc,uorder,iout, &
&                      ifail)

!! value and four derivatives at any radius rc of uu/r, uu an all-electron
!! wave function, from the radial equation it solves rather than from
!! stencils on uu (wf_rc_der): uu and uu'=up/(al r) at rc interpolated
!! by logmesh_m, then the Taylor coefficients of the first-order system
!! y'=A(r)y about rc, with A from the potential and its first three
!! derivatives at rc (see doc/smoothness_plan.md). Smooth in rc; every
!! derivative carries the solver's own error in uu, and the 4th the
!! round-off of V''' (eps/al**3) instead of that of a fourth derivative of
!! uu (eps/al**4): on hydrogen 1e-10 against 8e-7 for the stencils. uu and
!! up must agree near rc: if uu' from up differs from uu' by stencils by
!! more than 1e-6+100*(al*max(1,k rc))**4 relative, k the local wave
!! number at rc, the program stops. The solvers are 4th order in the phase
!! per mesh step, k al r: they agree to 1e-9 on the default mesh and 6e-6
!! at amesh=1.05 for valence states, and to 1e-5 for well states 5-18 Ha
!! up at rc 3-4 on amesh=1.004 (k al rc up to 0.08)
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic (the flux form
!!              of sr_flux_m), 2 fully-relativistic (Dirac, uu the large
!!              component)
!! @param c_light  speed of light (Hartree atomic units)
!! @param ll  angular momentum
!! @param kap  Dirac kappa, l or -(l+1) for j = l -/+ 1/2 (irel=2 only)
!! @param ee  energy of uu (Ha)
!! @param rr  log radial mesh (bohr)
!! @param vv  the potential uu solves the equation in (Ha)
!! @param uu  the wave function (*rr)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr), as the
!!            solvers return it
!! @param mmax  size of the log radial mesh
!! @param rc  core radius (bohr)
!! @param uorder  d^k(uu/r)/dr^k at rc in uorder(k+1), k=0..4
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp
 use logmesh_m, only: lm_interp, lm_derivs, lm_step

!Input variables
 integer, intent(in) :: irel,ll,kap,mmax
 real(dp), intent(in) :: c_light
 real(dp), intent(in) :: ee,rc
 real(dp), intent(in) :: rr(mmax),vv(mmax),uu(mmax),up(mmax)
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: uorder(5)
 integer, intent(out) :: ifail

!Local variables
 integer, parameter :: nord=4
 real(dp) :: aa(2,2,0:nord-1),yy(2,0:nord),ud(2),udst(2),vd(nord)
 real(dp) :: vk(0:nord-1),rinv(0:nord),rinv2(0:nord-1),xm(0:nord-1)
 real(dp) :: xminv(0:nord-1),cc,fss,sls,rk,fact,al,xk
 real(dp), allocatable :: dudr(:)
 integer :: kk,jj

 ifail=0
! Taylor coefficients about rc of V, 1/r, 1/r**2
 call lm_derivs(rr,vv,mmax,rc,nord-1,vd)
 fact=1.0d0
 do kk=0,nord-1
   if(kk>0) fact=fact*dble(kk)
   vk(kk)=vd(kk+1)/fact
 end do
 do kk=0,nord
   rinv(kk)=(-1.0d0)**kk/rc**(kk+1)
 end do
 do kk=0,nord-1
   rinv2(kk)=(-1.0d0)**kk*dble(kk+1)/rc**(kk+2)
 end do

! uu and uu' at rc, uu' from the solver's d(uu)/dn, checked against the
! stencils on uu so that inconsistent arrays from a solver are not used
 al=lm_step(rr,mmax)
 allocate(dudr(mmax))
 dudr(:)=up(:)/(al*rr(:))
 call lm_interp(rr,uu,mmax,rc,ud(1))
 call lm_interp(rr,dudr,mmax,rc,ud(2))
 deallocate(dudr)
 call lm_derivs(rr,uu,mmax,rc,1,udst)
 xk=sqrt(max(2.0d0*(ee-vk(0))-dble(ll*(ll+1))/rc**2,0.0d0))*rc
 if(abs(ud(2)-udst(2))>(1.0d-6+1.0d2*(al*max(1.0d0,xk))**4) &
&                      *(abs(udst(2))+abs(udst(1))/rc)) then
   write(iout,'(/a,i0,a,f12.6/a,2es12.4)') &
&        'ode_rc_der: ERROR uu and up disagree for l=',ll,' at rc=',rc, &
&        'uu'' from up and from stencils on uu:',ud(2),udst(2)
   ifail=1
   return
 end if

 aa(:,:,:)=0.0d0
 if(irel<=1) then
! flux form: y=(uu,qq), uu'=uu/r + M qq, qq'=(l(l+1)/(M r**2)
! + 2(V-E)) uu - qq/r, M = 1 + (E-V)/(2c**2) (1 for irel=0)
   fss=0.0d0
   if(irel==1) fss=1.0d0/c_light**2
   sls=dble(ll*(ll+1))
   xm(0)=1.0d0 + 0.5d0*fss*(ee - vk(0))
   do kk=1,nord-1
     xm(kk)=-0.5d0*fss*vk(kk)
   end do
   xminv(0)=1.0d0/xm(0)
   do kk=1,nord-1
     xminv(kk)=0.0d0
     do jj=1,kk
       xminv(kk)=xminv(kk) - xm(jj)*xminv(kk-jj)
     end do
     xminv(kk)=xminv(kk)*xminv(0)
   end do
   do kk=0,nord-1
     aa(1,1,kk)=rinv(kk)
     aa(1,2,kk)=xm(kk)
     do jj=0,kk
       aa(2,1,kk)=aa(2,1,kk) + sls*xminv(jj)*rinv2(kk-jj)
     end do
     aa(2,1,kk)=aa(2,1,kk) + 2.0d0*vk(kk)
     aa(2,2,kk)=-rinv(kk)
   end do
   aa(2,1,0)=aa(2,1,0) - 2.0d0*ee
   yy(1,0)=ud(1)
   yy(2,0)=(ud(2) - ud(1)/rc)/xm(0)
 else
! Dirac: y=(g,f), g'=-kap g/r + (2c**2 + E - V) f/c,
! f'=kap f/r + (V - E) g/c
   cc=c_light
   do kk=0,nord-1
     aa(1,1,kk)=-dble(kap)*rinv(kk)
     aa(1,2,kk)=-vk(kk)/cc
     aa(2,1,kk)=vk(kk)/cc
     aa(2,2,kk)=dble(kap)*rinv(kk)
   end do
   aa(1,2,0)=aa(1,2,0) + (2.0d0*cc**2 + ee)/cc
   aa(2,1,0)=aa(2,1,0) - ee/cc
   yy(1,0)=ud(1)
   yy(2,0)=(ud(2) + dble(kap)*ud(1)/rc)/aa(1,2,0)
 end if

! (k+1) y_k+1 = sum_j A_j y_k-j
 do kk=0,nord-1
   yy(:,kk+1)=0.0d0
   do jj=0,kk
     yy(:,kk+1)=yy(:,kk+1) + matmul(aa(:,:,jj),yy(:,kk-jj))
   end do
   yy(:,kk+1)=yy(:,kk+1)/dble(kk+1)
 end do

! uu/r = sum, and k! times its Taylor coefficients
 fact=1.0d0
 do kk=0,nord
   if(kk>0) fact=fact*dble(kk)
   rk=0.0d0
   do jj=0,kk
     rk=rk + yy(1,jj)*rinv(kk-jj)
   end do
   uorder(kk+1)=fact*rk
 end do
 return
 end subroutine ode_rc_der

 end module ode_rc_der_m

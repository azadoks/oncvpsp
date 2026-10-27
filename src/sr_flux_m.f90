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
!> The all-electron radial equation, non-relativistic or scalar-relativistic
!> (Koelling-Harmon), in its first-order flux form, as in
!> RadialAtomicSolvers.jl, for uu and qq = (uu' - uu/r)/M with the
!> relativistic mass M = 1 + (ee - vv)/(2c**2):
!> uu' = uu/r + M*qq,  qq' = (l(l+1)/(M*r**2) + 2*(vv - ee))*uu - qq/r.
!> It needs no dv/dr, whose finite differences of a GGA potential near the
!> nucleus made the scalar-relativistic SCF diverge on fine meshes. Shared
!> by the bound states (lschfb) and the scattering states (lschfs).

 module sr_flux_m
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: sr_gamma, sr_flux_coefs, sr_outward, sr_inward, sr_derivs

 contains

 function sr_gamma(ll,zz,fss)

!! the indicial exponent of the radial equation, uu ~ r**gamma at the
!! origin: that of the Koelling-Harmon equation, or l+1 for fss=0
!!
!! @param ll  angular-momentum quantum number
!! @param zz  atomic number
!! @param fss  1/c**2, or 0 for the non-relativistic equation

!Input variables
 integer, intent(in) :: ll
 real(dp), intent(in) :: zz,fss

!Output variables
 real(dp) :: sr_gamma

 if(fss==0.0d0) then
   sr_gamma=ll+1
 else
   sr_gamma=sqrt(ll*(ll+1)+1.0d0-fss*zz**2)
 end if

 return
 end function sr_gamma

 subroutine sr_flux_coefs(ll,ee,fss,al,rr,vv,mmax,xm,pm,qm)

!! the relativistic mass and the coefficients of the flux form on the
!! mesh index, d(uu)/dn = al*uu + pm*qq, d(qq)/dn = qm*uu - al*qq
!!
!! @param ll  angular-momentum quantum number
!! @param ee  energy
!! @param fss  1/c**2, or 0 for the non-relativistic equation
!! @param al  log of the mesh ratio
!! @param rr  log radial mesh
!! @param vv  local atomic potential
!! @param mmax  size of log grid
!! @param xm  relativistic mass M
!! @param pm  al*r*M
!! @param qm  al*r*(l(l+1)/(M*r**2) + 2*(vv - ee))

!Input variables
 integer, intent(in) :: ll,mmax
 real(dp), intent(in) :: ee,fss,al
 real(dp), intent(in) :: rr(mmax),vv(mmax)

!Output variables
 real(dp), intent(out) :: xm(mmax),pm(mmax),qm(mmax)

!Local variables
 real(dp) :: sls
 integer :: ii

 sls=ll*(ll+1)
 do ii=1,mmax
   xm(ii)=1.0d0+0.5d0*fss*(ee-vv(ii))
   pm(ii)=al*rr(ii)*xm(ii)
   qm(ii)=al*rr(ii)*(sls/(xm(ii)*rr(ii)**2) + 2.0d0*(vv(ii)-ee))
 end do

 return
 end subroutine sr_flux_coefs

 subroutine sr_outward(ll,ee,zz,fss,al,rr,vv,pm,qm,mmax,mlast, &
&                      uu,up,qq,qp,nnode)

!! outward integration of the flux form from the series at the origin to
!! mlast, counting the nodes; the coefficients pm, qm are those of
!! sr_flux_coefs at the energy ee
!!
!! @param ll  angular-momentum quantum number
!! @param ee  energy
!! @param zz  atomic number
!! @param fss  1/c**2, or 0 for the non-relativistic equation
!! @param al  log of the mesh ratio
!! @param rr  log radial mesh
!! @param vv  local atomic potential
!! @param pm  coefficient of qq in d(uu)/dn
!! @param qm  coefficient of uu in d(qq)/dn
!! @param mmax  size of log grid
!! @param mlast  last point of the integration
!! @param uu  radial wave function (*rr) on 1..mlast
!! @param up  d(uu)/dn on 1..mlast
!! @param qq  flux function on 1..mlast
!! @param qp  d(qq)/dn on 1..mlast
!! @param nnode  number of nodes of uu inside rr(mlast)

!Input variables
 integer, intent(in) :: ll,mmax,mlast
 real(dp), intent(in) :: ee,zz,fss,al
 real(dp), intent(in) :: rr(mmax),vv(mmax),pm(mmax),qm(mmax)

!Input/output variables
 real(dp), intent(inout) :: uu(mmax),up(mmax),qq(mmax),qp(mmax)

!Output variables
 integer, intent(out) :: nnode

!Local variables
 real(dp) :: gamma,aa,bb,ww,sls
 integer :: ii,it

 gamma=sr_gamma(ll,zz,fss)

! the model potential -zz/r + d of the start, d the screening at rr(1):
! M = (aa + bb*r)/r, and ww = d - ee
 sls=ll*(ll+1)
 ww=vv(1)+zz/rr(1)-ee
 aa=0.5d0*fss*zz
 bb=1.0d0-0.5d0*fss*ww

! start wavefunction with series

 call origin_start()
 do ii=1,4
   call sr_derivs(al,pm,qm,mmax,uu,up,qq,qp,ii)
 end do

! outward integration using predictor once, corrector
! twice
 nnode=0

 do ii=4,mlast-1
   uu(ii+1)=uu(ii)+aeo(up,ii)
   qq(ii+1)=qq(ii)+aeo(qp,ii)
   do it=1,2
     call sr_derivs(al,pm,qm,mmax,uu,up,qq,qp,ii+1)
     uu(ii+1)=uu(ii)+aio(up,ii)
     qq(ii+1)=qq(ii)+aio(qp,ii)
   end do
   call sr_derivs(al,pm,qm,mmax,uu,up,qq,qp,ii+1)
   if(uu(ii+1)*uu(ii) <= 0.0d0) nnode=nnode+1
 end do

! the first points again, inward from the outward solution: the series of
! the start is that of a model potential, and the kink where it meets the
! integration is amplified by the derivatives of a GGA density
 if(mlast>=16) call sr_inward(al,pm,qm,mmax,12,1,uu,up,qq,qp)

 return

 contains

 subroutine origin_start()

!! uu and qq at the first four points, from the series at the origin of
!! the flux form for the potential -zz/r + d (origin_series). With
!! M = (a + b*r)/r the series converges only inside r = a/b ~ zz/(2c**2),
!! for every l; points beyond an eighth of that (rr(1) of the default mesh
!! for H) are reached from there by RK4 in ln r on the same model equation

 real(dp) :: rs,r0,px,qx
 integer :: jj

 if(fss==0.0d0) then
   do jj=1,4
     call origin_series(rr(jj),uu(jj),qq(jj))
   end do
   return
 end if
! the radius is a/|b| for either sign of b, which a trial energy far below
! the bound states can make negative. b = 1 + (ee - d)/(2c**2) is at most 1
! for every energy below the screening at the origin (all bound and
! scattering states), so rs = a/8 there: which points the series reaches
! and the RK4 steps to the others do not depend on the energy, and the
! solution is smooth in it (doc/smoothness_plan.md, S-9); a/8 is inside a
! quarter of the radius for |b| <= 2
 rs=0.125d0*aa/max(1.0d0,abs(bb))
 r0=min(rs,rr(1))
 call origin_series(r0,px,qx)
 do jj=1,4
   if(rr(jj)<=rs) then
     call origin_series(rr(jj),px,qx)
   else
     call model_rk4(r0,rr(jj),px,qx)
   end if
   r0=rr(jj)
   uu(jj)=px
   qq(jj)=qx
 end do

 return
 end subroutine origin_start

 subroutine model_rk4(ra,rb,px,qx)

!! carries uu (px) and qq (qx) of the model equation of origin_series from
!! ra to rb by RK4 in x = ln r, in steps of at most 2e-3:
!! dP/dx = P + (a + b*r)*q,
!! dq/dx = (l(l+1)/(a + b*r) - 2*zz + 2*w*r)*P - q

 real(dp), intent(in) :: ra,rb
 real(dp), intent(inout) :: px,qx
 real(dp) :: hh,xa,p1,p2,p3,p4,q1,q2,q3,q4
 integer :: kk,nst

 nst=max(1,ceiling(log(rb/ra)/2.0d-3))
 hh=log(rb/ra)/nst
 xa=log(ra)
 do kk=1,nst
   call model_derivs(xa+(kk-1)*hh,px,qx,p1,q1)
   call model_derivs(xa+(kk-0.5d0)*hh,px+0.5d0*hh*p1,qx+0.5d0*hh*q1,p2,q2)
   call model_derivs(xa+(kk-0.5d0)*hh,px+0.5d0*hh*p2,qx+0.5d0*hh*q2,p3,q3)
   call model_derivs(xa+kk*hh,px+hh*p3,qx+hh*q3,p4,q4)
   px=px+hh*(p1+2.0d0*(p2+p3)+p4)/6.0d0
   qx=qx+hh*(q1+2.0d0*(q2+q3)+q4)/6.0d0
 end do

 return
 end subroutine model_rk4

 subroutine model_derivs(xx,pp,qv,dpx,dqx)

!! dP/dx (dpx) and dq/dx (dqx) of the model equation at x = xx

 real(dp), intent(in) :: xx,pp,qv
 real(dp), intent(out) :: dpx,dqx
 real(dp) :: rx,am

 rx=exp(xx)
 am=aa+bb*rx
 dpx=pp+am*qv
 dqx=(sls/am-2.0d0*zz+2.0d0*ww*rx)*pp-qv

 return
 end subroutine model_derivs

 subroutine origin_series(rx,px,qx)

!! uu (px) and qq (qx) at rx from their series at the origin for the
!! potential -zz/r + d, d the screening at rr(1), at the energy ee.
!! With M = (a + b*r)/r and w = d - ee, the flux form is
!! r*P' = P + (a + b*r)*q and
!! (a + b*r)*r*q' = (l(l+1) + (a + b*r)*(-2*zz + 2*w*r))*P - (a + b*r)*q.
!! Inside zz/(2c**2) the mass is set by the nuclear potential, so r**gamma
!! alone, or with a first-order term, is not the solution at rr(1) of the
!! default mesh for light atoms; its error in the logarithmic derivative
!! acts as a contact potential (1e-5 Ha for Na). The terms are summed as
!! p_k*r**k, q_k*r**k, whose coefficients alone overflow. The relativistic
!! series is used only inside an eighth of its radius a/b (origin_start)

 real(dp), intent(in) :: rx
 real(dp), intent(out) :: px,qx
 real(dp) :: pk,qk,pk1,qk1,pk2,rk
 integer :: kk

 if(fss==0.0d0) then
! P = sum p_k r**(l+1+k), q = P' - P/r
   pk2=0.0d0; pk1=1.0d0
   px=1.0d0; qx=ll
   do kk=1,80
     pk=(-2.0d0*zz*pk1*rx+2.0d0*ww*pk2*rx**2)/(kk*(2*ll+1+kk))
     px=px+pk
     qx=qx+(ll+kk)*pk
     if(abs(pk)<1.0d-17*abs(px)) exit
     pk2=pk1; pk1=pk
   end do
   if(kk>80) call series_error(rx)
   px=px*rx**(ll+1)
   qx=qx*rx**ll
   return
 end if
! P = sum p_k r**(gamma+k), q = sum q_k r**(gamma+k)
 pk2=0.0d0; pk1=1.0d0
 qk1=(gamma-1.0d0)/aa
 px=pk1; qx=qk1
 do kk=1,200
   rk=(-bb*(gamma+kk)*qk1 + (2.0d0*ww*aa-2.0d0*zz*bb)*pk1)*rx &
&    + 2.0d0*ww*bb*pk2*rx**2
   pk=((gamma+kk+1)*bb*qk1*rx + rk)/(kk*(2.0d0*gamma+kk))
   qk=((gamma+kk-1)*rk + (sls-2.0d0*zz*aa)*bb*qk1*rx) &
&    /(aa*kk*(2.0d0*gamma+kk))
   px=px+pk
   qx=qx+qk
   if(abs(pk)<1.0d-17*abs(px) .and. abs(qk)<1.0d-17*abs(qx)) exit
   pk2=pk1; pk1=pk; qk1=qk
 end do
 if(kk>200) call series_error(rx)
 px=px*rx**gamma
 qx=qx*rx**gamma

 return
 end subroutine origin_series

 subroutine series_error(rx)
!! the origin series at rx not converged: an error
 real(dp), intent(in) :: rx
 write(error_unit,'(/a,i0,a,es12.4,a,es12.4)') &
&     'sr_outward: ERROR origin series not converged, l=',ll, &
&     ', ee=',ee,', r=',rx
 error stop 1
 end subroutine series_error

 end subroutine sr_outward

 subroutine sr_inward(al,pm,qm,mmax,ifrom,ito,uu,up,qq,qp)

!! inward Adams integration of the flux form from ifrom down to ito, two
!! corrector passes; uu, qq and their derivatives are given on
!! ifrom..ifrom+3
!!
!! @param al  log of the mesh ratio
!! @param pm  coefficient of qq in d(uu)/dn
!! @param qm  coefficient of uu in d(qq)/dn
!! @param mmax  size of the arrays
!! @param ifrom  first point, given
!! @param ito  last point integrated
!! @param uu  radial wave function (*rr)
!! @param up  d(uu)/dn
!! @param qq  flux function
!! @param qp  d(qq)/dn

!Input variables
 integer, intent(in) :: mmax,ifrom,ito
 real(dp), intent(in) :: al,pm(mmax),qm(mmax)

!Input/output variables
 real(dp), intent(inout) :: uu(mmax),up(mmax),qq(mmax),qp(mmax)

!Local variables
 integer :: ii,it

 do ii=ifrom,ito+1,-1
   uu(ii-1)=uu(ii)+aei(up,ii)
   qq(ii-1)=qq(ii)+aei(qp,ii)
   do it=1,2
     call sr_derivs(al,pm,qm,mmax,uu,up,qq,qp,ii-1)
     uu(ii-1)=uu(ii)+aii(up,ii)
     qq(ii-1)=qq(ii)+aii(qp,ii)
   end do
   call sr_derivs(al,pm,qm,mmax,uu,up,qq,qp,ii-1)
 end do

 return
 end subroutine sr_inward

 subroutine sr_derivs(al,pm,qm,mmax,uu,up,qq,qp,jj)

!! derivatives on the mesh index of uu and qq at the point jj
!!
!! @param al  log of the mesh ratio
!! @param pm  coefficient of qq in d(uu)/dn
!! @param qm  coefficient of uu in d(qq)/dn
!! @param mmax  size of the arrays
!! @param uu  radial wave function (*rr)
!! @param up  d(uu)/dn
!! @param qq  flux function
!! @param qp  d(qq)/dn
!! @param jj  the point

!Input variables
 integer, intent(in) :: mmax,jj
 real(dp), intent(in) :: al,pm(mmax),qm(mmax),uu(mmax),qq(mmax)

!Input/output variables
 real(dp), intent(inout) :: up(mmax),qp(mmax)

 up(jj)=al*uu(jj)+pm(jj)*qq(jj)
 qp(jj)=qm(jj)*uu(jj)-al*qq(jj)

 return
 end subroutine sr_derivs

 include 'adams.inc'

 end module sr_flux_m

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
!> The outward integration of the radial Dirac equation for the large (g)
!> and small (f) components (times r),
!> g' = -kap*g/r + (2c**2 + ee - vv)*f/c,  f' = kap*f/r + (vv - ee)*g/c,
!> shared by the bound states (ldiracfb) and the scattering states
!> (ldiracfs), as sr_flux_m is for the scalar-relativistic equation: the
!> series at the origin for the potential -zz/r + d (d the screening at
!> rr(1)) at the first four points, Adams predictor-corrector with the
!> derivatives evaluated after each correction (PECECE), and the first
!> points integrated again inward from the outward solution
!> (doc/numerics_audit.md, C-4). The leading term r**gamma alone at the
!> first five points, as before, missed the solution's first-order terms
!> by about zz*rr(5); the scalar-relativistic solver once had a contact
!> error of 1e-5 Ha from the same start.

 module dirac_m
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: dirac_outward

 contains

 subroutine dirac_outward(ll,kap,ee,zz,al,rr,vv,mmax,mlast,gu,fu,gup,fup, &
&                         nnode)

!! outward integration from the series at the origin to mlast at the
!! energy ee, counting the nodes of the large component from the 5th point
!!
!! @param ll  angular-momentum quantum number
!! @param kap  =l, -(l+1) for j=l -/+ 1/2
!! @param ee  energy
!! @param zz  atomic number
!! @param al  log of the mesh ratio
!! @param rr  log radial mesh
!! @param vv  local atomic potential
!! @param mmax  size of log grid
!! @param mlast  last point of the integration
!! @param gu  large component (*rr) on 1..mlast
!! @param fu  small component (*rr) on 1..mlast
!! @param gup  d(gu)/dn on 1..mlast
!! @param fup  d(fu)/dn on 1..mlast
!! @param nnode  number of nodes of gu between the 5th point and mlast

 use relativity_m, only: c_light

!Input variables
 integer, intent(in) :: ll,kap,mmax,mlast
 real(dp), intent(in) :: ee,zz,al
 real(dp), intent(in) :: rr(mmax),vv(mmax)

!Input/output variables
 real(dp), intent(inout) :: gu(mmax),fu(mmax),gup(mmax),fup(mmax)

!Output variables
 integer, intent(out) :: nnode

!Local variables
 real(dp) :: cc,cci,gam,cof,dd
 integer :: ii,it

 cc=c_light
 cci=1.0d0/cc
 gam=sqrt(kap**2-(zz*cci)**2)
 cof=(kap+gam)*(cc/zz)
 dd=vv(1)+zz/rr(1)

 do ii=1,4
   call origin_series(rr(ii),gu(ii),fu(ii))
   call derivs_at(ii)
 end do

 nnode=0
 do ii=4,mlast-1
   fu(ii+1)=fu(ii)+aeo(fup,ii)
   gu(ii+1)=gu(ii)+aeo(gup,ii)
   do it=1,2
     call derivs_at(ii+1)
     fu(ii+1)=fu(ii)+aio(fup,ii)
     gu(ii+1)=gu(ii)+aio(gup,ii)
   end do
   call derivs_at(ii+1)
   if(ii>=5 .and. gu(ii+1)*gu(ii)<=0.0d0) nnode=nnode+1
 end do

! the first points again, inward from the outward solution: the series is
! that of a model potential (as in sr_flux_m)
 if(mlast>=16) then
   do ii=12,2,-1
     fu(ii-1)=fu(ii)+aei(fup,ii)
     gu(ii-1)=gu(ii)+aei(gup,ii)
     do it=1,2
       call derivs_at(ii-1)
       fu(ii-1)=fu(ii)+aii(fup,ii)
       gu(ii-1)=gu(ii)+aii(gup,ii)
     end do
     call derivs_at(ii-1)
   end do
 end if

 return

 contains

 subroutine origin_series(rx,gx,fx)

!! gu (gx) and fu (fx) at rx from their series at the origin for the
!! potential -zz/r + d: g = r**gamma sum a_k r**k, f = r**gamma sum b_k r**k,
!! a_0 = 1, b_0 = cof, and for k >= 1, with s = gamma + k,
!! (s + kap) a_k - (zz/c) b_k = (2c**2 + ee - d) b_(k-1)/c,
!! (zz/c) a_k + (s - kap) b_k = (d - ee) a_(k-1)/c,
!! whose determinant is k (2 gamma + k). The terms are summed as a_k r**k,
!! b_k r**k

 real(dp), intent(in) :: rx
 real(dp), intent(out) :: gx,fx
 real(dp) :: pk,qk,pk1,qk1,aa,bb,sk
 integer :: kk

 pk1=1.0d0
 qk1=cof
 gx=pk1
 fx=qk1
 do kk=1,100
   sk=gam+kk
   aa=(2.0d0*cc**2+ee-dd)*cci*qk1*rx
   bb=(dd-ee)*cci*pk1*rx
   pk=((sk-kap)*aa + zz*cci*bb)/(kk*(2.0d0*gam+kk))
   qk=((sk+kap)*bb - zz*cci*aa)/(kk*(2.0d0*gam+kk))
   gx=gx+pk
   fx=fx+qk
   if(abs(pk)<=1.0d-17*abs(gx) .and. abs(qk)<=1.0d-17*abs(fx)) exit
   if(kk==100) then
     write(error_unit,'(/a,i0,a,i0,a,es12.4,a,es12.4)') &
&          'dirac_outward: ERROR origin series not converged, l=',ll, &
&          ', kap=',kap,', ee=',ee,', r=',rx
     error stop 1
   end if
   pk1=pk
   qk1=qk
 end do
 gx=gx*rx**gam
 fx=fx*rx**gam

 return
 end subroutine origin_series

 subroutine derivs_at(jj)

!! derivatives on the mesh index of gu and fu at the point jj

 integer, intent(in) :: jj

 fup(jj)=al*rr(jj)*(kap*fu(jj)/rr(jj) + cci*(vv(jj)-ee)*gu(jj))
 gup(jj)=al*rr(jj)*(-kap*gu(jj)/rr(jj) &
&        + cci*(2.0d0*cc**2+ee-vv(jj))*fu(jj))

 return
 end subroutine derivs_at

 end subroutine dirac_outward

 include 'adams.inc'

 end module dirac_m

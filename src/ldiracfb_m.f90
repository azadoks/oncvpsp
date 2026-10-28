!
! Copyright (c) 1989-2012 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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

 module ldiracfb_m

 implicit none
 private
 public :: ldiracfb

 contains

 subroutine ldiracfb(nn,ll,kap,c_light,ierr,ee,rr,zz,vv,uu,up,mmax,mch,iout, &
&                    ifail)

!! Finds relativistic bound states of an al-electron potential
!!
!! @param nn  principal quantum number
!! @param ll  angular-momentum quantum number
!! @param kap  =l, -(l+1) for j=l -/+ 1/2
!! @param c_light  speed of light (Hartree atomic units)
!! @param ierr  non-zero return if error, 10 if there is no bound state
!!              and the box state (large component 0 at rr(mmax)) is
!!              returned
!! @param ee  bound-state energy, input guess and output calculated value
!! @param rr  log radial mesh
!! @param zz  atomic number
!! @param vv  local psp
!! @param uu(mmax,jj)  output radial wave functions (*rr) jj=1,2 for large,
!!                     small
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of log grid
!! @param mch  mesh point just beyond the classical turning point, where
!!             the inward and outward integrations are matched (at the
!!             continuous turning point, match_m)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp, eig_tol
 use boundsolve_m, only: bs_box_energy, bs_solve, bs_turning_point
 use match_m, only: match_k, match_point, match_weights, match_blend
 use dirac_m, only: dirac_outward
 use quad_m, only: log_norm
 use inward_start_m, only: inward_start
 use logmesh_m, only: lm_step

!Input variables
 integer, intent(in) :: nn,ll,kap,mmax
 real(dp), intent(in) :: rr(mmax),vv(mmax)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: c_light
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)
 real(dp), intent(inout) :: ee
 integer, intent(inout) :: mch

!Output variables
 integer, intent(out) :: ierr
 integer, intent(out) :: ifail

!Local variables

 real(dp), allocatable :: gu(:),fu(:),gup(:),fup(:),cf(:)

 real(dp) :: cc,cci,gam,cof
 real(dp) :: emax,emin
 real(dp) :: eps
 real(dp) :: sls,sn,cn
 real(dp) :: amesh,al,eecur
 integer :: ii

 ifail=0
 cc=c_light
 cci=1.0d0/cc

 al=lm_step(rr,mmax)
 amesh = exp(al)

! convergence: once the energy correction de is below eps*max(|ee|,0.2),
! passes go on at the corrected energy until de is at round-off, 8 ulps of
! max(|ee|,0.2), or stops halving (the round-off of de itself); the energy
! and function returned are those of the last pass, whose correction is
! not applied, so they are the fixed point to round-off and do not depend
! on the start guess or the path (doc/smoothness_plan.md, S-2 to S-4, and
! lschfb)
 eps=eig_tol
 ierr = 60

! check arguments
 if(ll>nn-1) then
   write(iout,'(/a,i4,a,i4)') 'ldiracfb: ERROR ll =',ll,' > nn =',nn
   ierr=1
   return
 end if
 if(kap/=ll .and. kap/=-(ll+1)) then
   write(iout,'(/a,i4,a,i4)') 'ldiracfb: ERROR kap =',kap,' ll =',ll
   ierr=2
   return
 end if
 if(zz<1.0d0) then
   write(iout,'(/a,f12.8)') 'ldiracfb: ERROR zz =',zz
   ierr=3
   return
 end if

 sls=ll*(ll+1)

! usual limits from non-relativistic Schroedinger eq. should be OK
 emax=vv(mmax)+0.5d0*sls/rr(mmax)**2
 emin=emax
 do ii=1,mmax
   emin=min(emin,vv(ii)+0.5d0*sls/rr(ii)**2)
 end do
 emin=max(emin,-zz**2/nn**2)
 if(ee>emax) ee=0.5d0*(emax+emin)
 if(ee<emin) ee=0.5d0*(emax+emin)
 if(ee>emax) ee=0.5d0*(emax+emin)

 allocate(gu(mmax),fu(mmax),gup(mmax),fup(mmax),cf(mmax))

! null arrays
 gu(:)=0.0d0; fu(:)=0.0d0; gup(:)=0.0d0; fup(:)=0.0d0; cf(:)=0.0d0

 gam=sqrt(kap**2-(zz*cci)**2)

 cof=(kap+gam)*(cc/zz)

! the energy iteration (boundsolve_m)
 call bs_solve(nn,ll,rr,vv,mmax,eps,setup,outward_nodes,match,boxstate, &
&              ee,emin,emax,mch,ierr,ifail)
 if(ifail/=0) return

! copy local arrays for output
 uu(:,1)=gu(:)
 up(:,1)=gup(:)
 uu(:,2)=fu(:)
 up(:,2)=fup(:)

!fix sign to be positive at rr->oo
 if(uu(mch,1)<0.0d0) then
   uu(:,:)=-uu(:,:)
   up(:,:)=-up(:,:)
 end if

 deallocate(gu,fu,gup,fup)
 deallocate(cf)
 return

 contains

 subroutine setup(et,mchs)

!! the coefficient for the turning point at et, and the classical turning
!! point mchs

 real(dp), intent(in) :: et
 integer, intent(out) :: mchs

 eecur=et
! coefficient array for u in Schroedinger eq.
 do ii=1,mmax
   cf(ii)=sls + 2.0d0*(vv(ii)-et)*rr(ii)**2
 end do
 call bs_turning_point(cf,mmax,mchs)

 return
 end subroutine setup

 integer function outward_nodes(mchs,mlast)

!! the outward solution up to mlast and the nodes of its large component
!! inside mchs

 integer, intent(in) :: mchs,mlast

 integer :: jj,node

 call outward(mlast,node)
 node=0
 do jj=5,mchs-1
   if(gu(jj+1)*gu(jj)<=0.0d0) node=node+1
 end do
 outward_nodes=node

 return
 end function outward_nodes

 subroutine match(et,mchs,nin,lo,hi,de,ifail)

!! the inward solution from nin at et, matched to the outward one across
!! lo..hi and normalized, and the energy correction de

 real(dp), intent(in) :: et
 integer, intent(in) :: mchs,nin,lo,hi
 real(dp), intent(out) :: de

 real(dp) :: sc,xkap,xmch
 real(dp) :: go(-match_k-1:match_k),fo(-match_k-1:match_k)
 real(dp) :: gpo(-match_k-1:match_k),fpo(-match_k-1:match_k)
 real(dp) :: ww(-match_k-1:match_k),ss(-match_k-1:match_k)
 integer :: jj,kk
 integer, intent(out) :: ifail

 ifail=0
 call match_point(rr,cf,mmax,mchs,xmch,iout,ifail)
 if(ifail/=0) return
 go(:)=gu(lo:hi)
 fo(:)=fu(lo:hi)
 gpo(:)=gup(lo:hi)
 fpo(:)=fup(lo:hi)

 xkap=sqrt(sls/rr(nin)**2 + 2.0d0*(vv(nin)-et))

 do jj=nin,nin+4
   gu(jj)=exp(-xkap*(rr(jj)-rr(nin)))
   fu(jj)=0.5d0*cci*(-xkap+kap/rr(jj))*gu(jj)

   fup(jj)= al*rr(jj)*(kap*fu(jj)/rr(jj) &
&           + cci*(vv(jj) - et)*gu(jj))
   gup(jj)= al*rr(jj)*(-kap*gu(jj)/rr(jj) &
&           + cci*(2.0d0*cc**2 + et - vv(jj))*fu(jj))
 end do

! integrate inward
 do jj=nin,lo+1,-1

   fu(jj-1)=fu(jj)+aei(fup,jj)
   gu(jj-1)=gu(jj)+aei(gup,jj)

   do kk=1,2
     fup(jj-1)= al*rr(jj-1)*(kap*fu(jj-1)/rr(jj-1) &
&             - cci*(et - vv(jj-1))*gu(jj-1))

     gup(jj-1)= al*rr(jj-1)*(-kap*gu(jj-1)/rr(jj-1) &
&             + cci*(2.0d0*cc**2 + et - vv(jj-1))*fu(jj-1))

     fu(jj-1)=fu(jj)+aii(fup,jj)
     gu(jj-1)=gu(jj)+aii(gup,jj)
   end do
   fup(jj-1)= al*rr(jj-1)*(kap*fu(jj-1)/rr(jj-1) &
&           - cci*(et - vv(jj-1))*gu(jj-1))
   gup(jj-1)= al*rr(jj-1)*(-kap*gu(jj-1)/rr(jj-1) &
&           + cci*(2.0d0*cc**2 + et - vv(jj-1))*fu(jj-1))
 end do

! perturbation theory for the energy shift based on the large component
! at each point of the window, g_o (g'_o - g'_i g_o/g_i)/2 for the
! normalized function (continuity of the small component should follow),
! averaged with the bump weights; the inward solution scaled to the
! outward one by the ratio of their weighted sums, and the two blended
! across the window (match_m)

 call match_weights(xmch,lo,hi,ww,ss)
 de=0.0d0
 do jj=lo,hi
   de=de+0.5d0*ww(jj-mchs)*go(jj-mchs) &
&     *(gpo(jj-mchs)-gup(jj)*go(jj-mchs)/gu(jj))/(al*rr(jj))
 end do
 sc=sum(ww*go)/sum(ww*gu(lo:hi))

 do jj=lo,nin
   gup(jj)=sc*gup(jj)
   gu(jj)=sc*gu(jj)
   fup(jj)=sc*fup(jj)
   fu(jj)=sc*fu(jj)
 end do
 call match_blend(lo,hi,mchs,ss,go,gu)
 call match_blend(lo,hi,mchs,ss,fo,fu)
 call match_blend(lo,hi,mchs,ss,gpo,gup)
 call match_blend(lo,hi,mchs,ss,fpo,fup)

! normalize u
 call normalize(nin)
 de=cn**2*de

 return
 end subroutine match

 subroutine outward(mlast,nnode)

!! outward integration from the series at the origin to mlast, at the
!! energy ee, counting the nodes of the large component (dirac_m)

 integer, intent(in) :: mlast
 integer, intent(out) :: nnode

 call dirac_outward(ll,kap,c_light,eecur,zz,al,rr,vv,mmax,mlast,gu,fu,gup,fup, &
&                   nnode)

 return
 end subroutine outward

 subroutine normalize(mlast)

!! normalize the large and small components on 1..mlast, zero beyond

 integer, intent(in) :: mlast

! perform normalization sum

 sn=log_norm(rr,gu(1:mlast)**2+fu(1:mlast)**2,mlast,ll,gam)

 cn=1.0d0/sqrt(sn)

 do ii=1,mlast
   gup(ii)=cn*gup(ii)
   gu(ii)=cn*gu(ii)
   fup(ii)=cn*fup(ii)
   fu(ii)=cn*fu(ii)
 end do
 do ii=mlast+1,mmax
   gu(ii)=0.0d0
   fu(ii)=0.0d0
   gup(ii)=0.0d0
   fup(ii)=0.0d0
 end do

 return
 end subroutine normalize

 subroutine boxstate(elow,etop,eb)

!! the box state above elow, widened above etop (bs_box_energy), and its
!! energy eb

 real(dp), intent(in) :: elow,etop
 real(dp), intent(out) :: eb

 integer :: nb

 call bs_box_energy(elow,etop,nn-ll-1,eps,box_nodes,eb)
 eecur=eb
 call outward(mmax,nb)
 call normalize(mmax)

 return
 end subroutine boxstate

 integer function box_nodes(ee)

!! the node count of the outward solution over the whole mesh at ee

 real(dp), intent(in) :: ee

 eecur=ee
 call outward(mmax,box_nodes)
 return
 end function box_nodes

 include 'adams.inc'

 end subroutine ldiracfb

 end module ldiracfb_m

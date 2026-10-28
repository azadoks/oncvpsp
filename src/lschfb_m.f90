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

 module lschfb_m

 implicit none
 private
 public :: lschfb

 contains

 subroutine lschfb(irel,c_light,nn,ll,ierr,ee,rr,zz,vv,uu,up,mmax,mch,iout, &
&                  ifail)

!! Finds bound states of an all-electron atomic potential using
!! Pauli-type  scalar-relativistic Schroedinger equation
!!
!! The Koelling-Harmon equation is integrated in its first-order flux form
!! (sr_flux_m)
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param nn  principal quantum number
!! @param ll  angular momentum
!! @param ierr  non-zero return if error, 10 if there is no bound state
!!              and the box state (u=0 at rr(mmax)) is returned
!! @param ee  bound-state energy, input guess and output calculated value
!! @param rr  log radial mesh (bohr)
!! @param zz  atomic number
!! @param vv  local atomic potential
!! @param uu  output radial wave function (*rr)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of the log radial mesh
!! @param mch  mesh point just beyond the classical turning point, where
!!             the inward and outward integrations are matched (at the
!!             continuous turning point, match_m)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp, eig_tol
 use boundsolve_m, only: bs_box_energy, bs_solve, bs_turning_point
 use sr_flux_m, only: sr_gamma, sr_flux_coefs, sr_outward, sr_inward, &
&                     sr_derivs
 use match_m, only: match_k, match_point, match_weights, match_blend
 use quad_m, only: quad_log_norm
 use inward_start_m, only: inward_start
 use logmesh_m, only: lm_step

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),vv(mmax)
 real(dp), intent(in) :: zz
 integer, intent(in) :: nn,ll
 integer, intent(in) :: irel
 real(dp), intent(in) :: c_light
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: uu(mmax),up(mmax)
 real(dp), intent(inout) :: ee
 integer, intent(inout) :: mch

!Output variables
 integer, intent(out) :: ierr
 integer, intent(out) :: ifail

!Local variables

 real(dp) :: als
 real(dp) :: emax,emin
 real(dp) :: eps,fss,gamma
 real(dp) :: sls,sn,cn
 real(dp) :: amesh,al,eecur
 integer :: ii

 real(dp), allocatable :: qq(:),qp(:),cf(:),pm(:),qm(:),xm(:),r2(:)
 ifail=0
 allocate(qq(mmax),qp(mmax),cf(mmax),pm(mmax),qm(mmax),xm(mmax),r2(mmax))

 al=lm_step(rr,mmax)
 amesh = exp(al)

! convergence: once the energy correction de is below eps*max(|ee|,0.2),
! passes go on at the corrected energy until de is at round-off, 8 ulps of
! max(|ee|,0.2), or stops halving (the round-off of de itself); the energy
! and function returned are those of the last pass, whose correction is
! not applied, so they are the fixed point to round-off and do not depend
! on the start guess or the path (doc/smoothness_plan.md, S-2 to S-4).
! Stopping at the second pass below eps left that pass's correction,
! 1e-15 to 1e-14 relative, as noise in a well state's energy as a function
! of its well, which its well width and its tail carried into the
! projectors at 1e-9. 1e-10 left 3.6e-7 Ha in the 1s of U, which the
! eigenvalue sum carries into the total energy, and which depended on the
! SCF history
 eps=eig_tol
 ierr = 60

! check arguments
 if(ll>nn-1) then
   write(iout,'(/a,i4,a,i4)') 'lschfb: ERROR ll =',ll,' > nn =',nn
   ierr=1
   return
 end if

! relativistic - non-relativistic switch
 if(irel/=0) then
   fss=(1.0d0/c_light)**2
 else
   fss=0.0d0
 end if

! the indicial exponent, uu ~ r**gamma at the origin
 gamma=sr_gamma(ll,zz,fss)

 sls=ll*(ll+1)

 emax=vv(mmax)+0.5d0*sls/rr(mmax)**2
 emin=emax
 do ii=1,mmax
   if(ll==0) then
! ad-hock step to eliminate probelms from absurd emin for ll=0
     emin=min(emin,vv(ii)+0.25d0/rr(ii)**2)
   else
     emin=min(emin,vv(ii)+0.5d0*sls/rr(ii)**2)
   end if
 end do
 if(ee==0.0d0) ee=0.5d0*(emax+emin)
 if(ee<emin) ee=0.5d0*(emax+emin)
 if(ee>emax) ee=0.5d0*(emax+emin)

! null arrays to remove leftover garbage
 uu(:)=0.0d0
 up(:)=0.0d0
 qq(:)=0.0d0
 qp(:)=0.0d0

 als=al**2

! squared radii, which do not depend on the energy
 do ii=1,mmax
   r2(ii)=rr(ii)**2
 end do

! the energy iteration (boundsolve_m)
 call bs_solve(nn,ll,rr,vv,mmax,eps,setup,outward_nodes,match,boxstate, &
&              ee,emin,emax,mch,ierr,ifail)
 if(ifail/=0) return

!fix sign to be positive at rr->oo
 if(uu(mch)<0.0d0) then
   uu(:)=-uu(:)
   up(:)=-up(:)
 end if

 deallocate(qq,qp,cf,pm,qm,xm,r2)
 return

 contains

 subroutine coefs(et)

!! coefficient arrays of the radial equation at the energy et

 real(dp), intent(in) :: et

! coefficient array for u in the second-order equation, whose sign
! change gives the classical turning point
 do ii=1,mmax
   cf(ii)=als*sls + 2.0d0*als*(vv(ii)-et)*r2(ii)
 end do

 call sr_flux_coefs(ll,et,fss,al,rr,vv,mmax,xm,pm,qm)
 eecur=et

 return
 end subroutine coefs

 subroutine setup(et,mchs)

!! the coefficients at et and the classical turning point mchs

 real(dp), intent(in) :: et
 integer, intent(out) :: mchs

 call coefs(et)
 call bs_turning_point(cf,mmax,mchs)

 return
 end subroutine setup

 integer function outward_nodes(mchs,mlast)

!! the outward solution up to mlast and its nodes inside mchs

 integer, intent(in) :: mchs,mlast

 integer :: jj,node

 call outward(mlast,node)
 node=0
 do jj=4,mchs-1
   if(uu(jj+1)*uu(jj)<=0.0d0) node=node+1
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
 real(dp) :: uo(-match_k-1:match_k),qo(-match_k-1:match_k)
 real(dp) :: upo(-match_k-1:match_k),qpo(-match_k-1:match_k)
 real(dp) :: ww(-match_k-1:match_k),ss(-match_k-1:match_k)
 integer :: jj
 integer, intent(out) :: ifail

 ifail=0
 call match_point(rr,cf,mmax,mchs,xmch,iout,ifail)
 if(ifail/=0) return
 uo(:)=uu(lo:hi)
 qo(:)=qq(lo:hi)
 upo(:)=up(lo:hi)
 qpo(:)=qp(lo:hi)

 xkap=sqrt(sls/rr(nin)**2 + 2.0d0*(vv(nin)-et))

 do jj=nin,nin+4
   uu(jj)=exp(-xkap*(rr(jj)-rr(nin)))
   qq(jj)=-(xkap+1.0d0/rr(jj))*uu(jj)/xm(jj)
   call sr_derivs(al,pm,qm,mmax,uu,up,qq,qp,jj)
 end do

! integrate inward
 call sr_inward(al,pm,qm,mmax,nin,lo,uu,up,qq,qp)

! perturbation theory for the energy shift at each point of the window,
! u_o (q_o - q_i u_o/u_i)/2 for the normalized function, averaged with
! the bump weights; the inward solution scaled to the outward one by the
! ratio of their weighted sums, and the two blended across the window

 call match_weights(xmch,lo,hi,ww,ss)
 de=0.0d0
 do jj=lo,hi
   de=de+0.5d0*ww(jj-mchs)*uo(jj-mchs) &
&     *(qo(jj-mchs)-qq(jj)*uo(jj-mchs)/uu(jj))
 end do
 sc=sum(ww*uo)/sum(ww*uu(lo:hi))

 do jj=lo,nin
   up(jj)=sc*up(jj)
   uu(jj)=sc*uu(jj)
   qq(jj)=sc*qq(jj)
   qp(jj)=sc*qp(jj)
 end do
 call match_blend(lo,hi,mchs,ss,uo,uu)
 call match_blend(lo,hi,mchs,ss,qo,qq)
 call match_blend(lo,hi,mchs,ss,upo,up)
 call match_blend(lo,hi,mchs,ss,qpo,qp)

! normalize u
 call normalize(nin)
 de=cn**2*de

 return
 end subroutine match

 subroutine outward(mlast,nnode)

!! outward integration from the series at the origin to mlast, at the
!! energy of the last coefs, counting the nodes

 integer, intent(in) :: mlast
 integer, intent(out) :: nnode

 call sr_outward(ll,eecur,zz,fss,al,rr,vv,pm,qm,mmax,mlast, &
&                uu,up,qq,qp,nnode)

 return
 end subroutine outward

 subroutine normalize(mlast)

!! normalize u on 1..mlast, zero beyond

 integer, intent(in) :: mlast

! perform normalization sum

 sn=quad_log_norm(rr,uu(1:mlast)**2,mlast,ll,gamma)

 cn=1.0d0/sqrt(sn)

 do ii=1,mlast
   up(ii)=cn*up(ii)
   uu(ii)=cn*uu(ii)
   qq(ii)=cn*qq(ii)
   qp(ii)=cn*qp(ii)
 end do
 do ii=mlast+1,mmax
   uu(ii)=0.0d0
   up(ii)=0.0d0
   qq(ii)=0.0d0
   qp(ii)=0.0d0
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
 call coefs(eb)
 call outward(mmax,nb)
 call normalize(mmax)

 return
 end subroutine boxstate

 integer function box_nodes(ee)

!! the node count of the outward solution over the whole mesh at ee

 real(dp), intent(in) :: ee

 call coefs(ee)
 call outward(mmax,box_nodes)
 return
 end function box_nodes

 end subroutine lschfb

 end module lschfb_m

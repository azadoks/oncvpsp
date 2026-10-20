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
 subroutine ldiracfb(nn,ll,kap,ierr,ee,rr,zz,vv,uu,up,mmax,mch)

!! Finds relativistic bound states of an al-electron potential
!!
!! @param nn  principal quantum number
!! @param ll  angular-momentum quantum number
!! @param kap  =l, -(l+1) for j=l -/+ 1/2
!! @param ierr  non-zero return if error, 10 if there is no bound state
!!              and the box state (large component 0 at rr(mmax)) is
!!              returned
!! @param ee  bound-state energy, input guess and output calculated value
!! @param rr  log radial mesh
!! @param zz  atomic number
!! @param vv  local psp
!! @param uu(mmax,jj)  output radial wave functions (*rr) jj=1,2 for large,
!!                     small
!! @param up  d(uu)/dr
!! @param mmax  size of log grid
!! @param mch  matching mesh point for inward-outward integrations

 use constants_m, only: dp
 use relativity_m, only: c_light
 implicit none

!Input variables
 real(dp), intent(in) :: rr(mmax),vv(mmax)
 real(dp), intent(in) :: zz
 integer, intent(in) :: nn,ll,kap,mmax

!Input/output variables
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)
 real(dp), intent(inout) :: ee
 integer, intent(inout) :: mch

!Output variables
 integer, intent(out) :: ierr

!Local variables

 real(dp), allocatable :: gu(:),fu(:),gup(:),fup(:),cf(:)

 real(dp) :: cc,cci,gam,cof
 real(dp) :: de,emax,emin,etop,elow
 real(dp) :: eps,ro,sc
 real(dp) :: sls,sn,cn,uout,upin,upout,xkap
 real(dp) :: amesh,al,als
 integer :: ii,kk,nint,node,nin,nsmall

 cc=c_light
 cci=1.0d0/cc

 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

! convergence: two passes in a row whose energy correction de is smaller
! in magnitude than eps*max(|ee|,0.2), the second at the energy with the
! first's correction applied, so that the energy and function returned are
! those of a pass at the corrected energy, whose own correction is far
! below eps: the result does not depend on the start guess or the path
! (doc/smoothness_plan.md, S-2 to S-4). A first correction already at
! round-off (a guess from a converged SCF) needs no second pass
 eps=1.0d-12
 nsmall=-1
 ierr = 60

! check arguments
 if(ll>nn-1) then
   write(6,'(/a,i4,a,i4)') 'ldiracfb: ERROR ll =',ll,' > nn =',nn
   ierr=1
   return
 end if
 if(kap/=ll .and. kap/=-(ll+1)) then
   write(6,'(/a,i4,a,i4)') 'ldiracfb: ERROR kap =',kap,' ll =',ll
   ierr=2
   return
 end if
 if(zz<1.0d0) then
   write(6,'(/a,f12.8)') 'ldiracfb: ERROR zz =',zz
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

! top of the bracket: no bound state lies above it on the mesh;
! highest energy found to have too few nodes
 etop=emax
 elow=etop

! return point for bound state convergence
 do nint=1,60

! a bracket closed without a converged energy shift: there is no bound
! state, which an occupied state can lack in an early, overshooting SCF
! potential; the box state stands in for it
   if(emax-emin<max(abs(ee),0.2d0)*eps) then
     call boxstate()
     ierr=10
     exit
   end if

! coefficient array for u in Schroedinger eq.
   do ii=1,mmax
     cf(ii)=sls + 2.0d0*(vv(ii)-ee)*rr(ii)**2
   end do

! find classical turning point for matching
   mch=0
   do ii=mmax,2,-1
     if(cf(ii-1)<=0.d0 .and. cf(ii)>0.d0) then
       mch=ii
       exit
     end if
   end do

   if(mch==0 .and. cf(mmax)>0.0d0) then
! classically forbidden everywhere: too low
     emin=ee
     elow=ee
     ee=0.5d0*(emin+emax)
     cycle
   else if(mch==0) then
! classically allowed at rr(mmax): too high
     emax=ee
     ee=0.5d0*(emin+emax)
     cycle
   end if

   call outward(mch,node)

   uout=gu(mch)
   upout=gup(mch)

! start inward integration with simple exponential where the
! wave function has decayed to machine precision
   call inward_start(ll,ee,rr,vv,mch,mmax,nin)

   if(node-nn+ll+1==0 .and. nin<=mch+4) then
! no room for the inward integration: too high
     emax=ee
     ee=0.5d0*(emin+emax)

   else if(node-nn+ll+1==0) then

     xkap=sqrt(sls/rr(nin)**2 + 2.0d0*(vv(nin)-ee))

     do ii=nin,nin+4
       gu(ii)=exp(-xkap*(rr(ii)-rr(nin)))
       fu(ii)=0.5d0*cci*(-xkap+kap/rr(ii))*gu(ii)

       fup(ii)= al*rr(ii)*(kap*fu(ii)/rr(ii) &
&               + cci*(vv(ii) - ee)*gu(ii))
       gup(ii)= al*rr(ii)*(-kap*gu(ii)/rr(ii) &
&               + cci*(2.0d0*cc**2 + ee - vv(ii))*fu(ii))
     end do

! integrate inward
     do ii=nin,mch+1,-1

       fu(ii-1)=fu(ii)+aei(fup,ii)
       gu(ii-1)=gu(ii)+aei(gup,ii)

       do kk=1,2
         fup(ii-1)= al*rr(ii-1)*(kap*fu(ii-1)/rr(ii-1) &
&                 - cci*(ee - vv(ii-1))*gu(ii-1))

         gup(ii-1)= al*rr(ii-1)*(-kap*gu(ii-1)/rr(ii-1) &
&                 + cci*(2.0d0*cc**2 + ee - vv(ii-1))*fu(ii-1))

         fu(ii-1)=fu(ii)+aii(fup,ii)
         gu(ii-1)=gu(ii)+aii(gup,ii)
       end do
     end do

! scale outside wf for continuity

     sc=uout/gu(mch)

     do ii=mch,nin
       gup(ii)=sc*gup(ii)
       gu(ii)=sc*gu(ii)
       fup(ii)=sc*fup(ii)
       fu(ii)=sc*fu(ii)
     end do

     upin=gup(mch)

! normalize u
     call normalize(nin)
     uout=cn*uout
     upout=cn*upout
     upin=cn*upin

! perturbation theory for energy shift based on large component
! continuity of small component should follow

     de=0.5d0*uout*(upout-upin)/(al*rr(mch))

! convergence test and possible exit

     if(abs(de)<max(abs(ee),0.2d0)*eps) then
       if(nsmall==nint-1 .or. &
&         abs(de)<=8.0d0*epsilon(1.0d0)*max(abs(ee),0.2d0)) then
         ierr = 0
         exit
       end if
       nsmall=nint
     end if

     if(de>0.0d0) then
       emin=ee
     else
       emax=ee
     end if
     ee=ee+de
     if(ee>emax .or. ee<emin) ee=0.5d0*(emax+emin)

   else if(node-nn+ll+1<0) then
! too few nodes
     emin=ee
     elow=ee
     ee=0.5d0*(emin+emax)

   else
! too many nodes
     emax=ee
     ee=0.5d0*(emin+emax)
   end if

 end do

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

 subroutine outward(mlast,nnode)

!! outward integration from the series at the origin to mlast, at the
!! energy ee, counting the nodes of the large component

 integer, intent(in) :: mlast
 integer, intent(out) :: nnode

! start wavefunctions with series
 do ii=1,5
   gu(ii)=rr(ii)**gam
   fu(ii)=cof*gu(ii)

   fup(ii)= al*rr(ii)*(kap*fu(ii)/rr(ii) &
&           + cci*(vv(ii) - ee)*gu(ii))

   gup(ii)= al*rr(ii)*(-kap*gu(ii)/rr(ii) &
&           + cci*(2.0d0*cc**2 + ee - vv(ii))*fu(ii))
 end do

! outward integration using predictor once, corrector
! twice
 nnode=0
 do ii=5,mlast-1

   fu(ii+1)=fu(ii)+aeo(fup,ii)
   gu(ii+1)=gu(ii)+aeo(gup,ii)

   do kk=1,2
     fup(ii+1)= al*rr(ii+1)*(kap*fu(ii+1)/rr(ii+1) &
&             - cci*(ee - vv(ii+1))*gu(ii+1))

     gup(ii+1)= al*rr(ii+1)*(-kap*gu(ii+1)/rr(ii+1) &
&             + cci*(2.0d0*cc**2 + ee - vv(ii+1))*fu(ii+1))

     fu(ii+1)=fu(ii)+aio(fup,ii)
     gu(ii+1)=gu(ii)+aio(gup,ii)
   end do
   if(gu(ii+1)*gu(ii) <= 0.0d0) nnode=nnode+1
 end do

 return
 end subroutine outward

 subroutine normalize(mlast)

!! normalize the large and small components on 1..mlast, zero beyond

 integer, intent(in) :: mlast

! perform normalization sum

 ro=rr(1)/sqrt(amesh)
 sn=((1.0d0+cof**2)*ro**(2.0d0*gam+1.0d0))/(2.0d0*gam+1.0d0)

 do ii=1,mlast-3
   sn=sn+al*rr(ii)*(gu(ii)**2 + fu(ii)**2)
 end do

 sn=sn + al*(23.0d0*rr(mlast-2)*(gu(mlast-2)**2 + fu(mlast-2)**2) &
&            + 28.0d0*rr(mlast-1)*(gu(mlast-1)**2 + fu(mlast-1)**2) &
&            +  9.0d0*rr(mlast  )*(gu(mlast  )**2 + fu(mlast  )**2))/24.0d0

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

 subroutine boxstate()

!! the state with nn-ll-1 nodes and a large component zero at rr(mmax),
!! above elow (the box state of RadialAtomicSolvers.jl): the outward
!! solution over the whole mesh gains a node each time it crosses zero
!! at rr(mmax), so the energy is bisected on its node count

 real(dp) :: eblo,ebhi
 integer :: nb,jj

 eblo=elow
 ebhi=etop+1.0d0
 do jj=1,60
   ee=ebhi
   call outward(mmax,nb)
   if(nb>nn-ll-1) exit
   eblo=ebhi
   ebhi=etop+2.0d0*(ebhi-etop)
 end do

 do jj=1,200
   ee=0.5d0*(eblo+ebhi)
   call outward(mmax,nb)
   if(nb>nn-ll-1) then
     ebhi=ee
   else
     eblo=ee
   end if
   if(ebhi-eblo<max(abs(ee),0.2d0)*eps) exit
 end do

 ee=eblo
 call outward(mmax,nb)
 call normalize(mmax)
 mch=mmax

 return
 end subroutine boxstate

 include 'adams.inc'

 end subroutine ldiracfb

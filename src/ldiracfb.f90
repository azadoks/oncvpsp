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
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of log grid
!! @param mch  mesh point just beyond the classical turning point, where
!!             the inward and outward integrations are matched (at the
!!             continuous turning point, match_m)

 use constants_m, only: dp
 use relativity_m, only: c_light
 use match_m, only: match_k, match_point, match_weights
 use dirac_m, only: dirac_outward
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
 real(dp) :: de,deprev,emax,emin,etop,elow
 real(dp) :: eps,ro,sc
 real(dp) :: sls,sn,cn,xkap,xmch
 real(dp) :: amesh,al,als
 real(dp) :: go(-match_k-1:match_k),fo(-match_k-1:match_k)
 real(dp) :: gpo(-match_k-1:match_k),fpo(-match_k-1:match_k)
 real(dp) :: ww(-match_k-1:match_k),ss(-match_k-1:match_k)
 integer :: ii,kk,nint,node,nin,lo,hi
 logical :: polish

 cc=c_light
 cci=1.0d0/cc

 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

! convergence: once the energy correction de is below eps*max(|ee|,0.2),
! passes go on at the corrected energy until de is at round-off, 8 ulps of
! max(|ee|,0.2), or stops halving (the round-off of de itself); the energy
! and function returned are those of the last pass, whose correction is
! not applied, so they are the fixed point to round-off and do not depend
! on the start guess or the path (doc/smoothness_plan.md, S-2 to S-4, and
! lschfb)
 eps=1.0d-12
 deprev=huge(1.0d0)
 polish=.false.
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
! potential; the box state stands in for it. Once a shift is below eps
! the passes polish a converging energy, whose alternating shifts close
! the bracket too
   if(.not. polish .and. emax-emin<max(abs(ee),0.2d0)*eps) then
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

! the outward solution, carried through the matching window
! mch-1-match_k..mch+match_k (match_m), and the nodes of its large
! component inside the turning point
   lo=mch-1-match_k
   hi=mch+match_k
   call outward(min(hi,mmax),node)
   node=0
   do ii=5,mch-1
     if(gu(ii+1)*gu(ii)<=0.0d0) node=node+1
   end do

! start inward integration with simple exponential where the
! wave function has decayed to machine precision
   call inward_start(ll,ee,rr,vv,mch,mmax,nin)

   if(node-nn+ll+1==0 .and. nin<=hi+4) then
! no room for the inward integration: too high
     emax=ee
     ee=0.5d0*(emin+emax)

   else if(node-nn+ll+1==0) then

     call match_point(rr,cf,mmax,mch,xmch)
     go(:)=gu(lo:hi)
     fo(:)=fu(lo:hi)
     gpo(:)=gup(lo:hi)
     fpo(:)=fup(lo:hi)

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
     do ii=nin,lo+1,-1

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
       fup(ii-1)= al*rr(ii-1)*(kap*fu(ii-1)/rr(ii-1) &
&               - cci*(ee - vv(ii-1))*gu(ii-1))
       gup(ii-1)= al*rr(ii-1)*(-kap*gu(ii-1)/rr(ii-1) &
&               + cci*(2.0d0*cc**2 + ee - vv(ii-1))*fu(ii-1))
     end do

! perturbation theory for the energy shift based on the large component
! at each point of the window, g_o (g'_o - g'_i g_o/g_i)/2 for the
! normalized function (continuity of the small component should follow),
! averaged with the bump weights; the inward solution scaled to the
! outward one by the ratio of their weighted sums, and the two blended
! across the window (match_m)

     call match_weights(xmch,lo,hi,ww,ss)
     de=0.0d0
     do ii=lo,hi
       de=de+0.5d0*ww(ii-mch)*go(ii-mch) &
&         *(gpo(ii-mch)-gup(ii)*go(ii-mch)/gu(ii))/(al*rr(ii))
     end do
     sc=sum(ww*go)/sum(ww*gu(lo:hi))

     do ii=lo,nin
       gup(ii)=sc*gup(ii)
       gu(ii)=sc*gu(ii)
       fup(ii)=sc*fup(ii)
       fu(ii)=sc*fu(ii)
     end do
     do ii=lo,hi
       gu(ii)=(1.0d0-ss(ii-mch))*go(ii-mch)+ss(ii-mch)*gu(ii)
       fu(ii)=(1.0d0-ss(ii-mch))*fo(ii-mch)+ss(ii-mch)*fu(ii)
       gup(ii)=(1.0d0-ss(ii-mch))*gpo(ii-mch)+ss(ii-mch)*gup(ii)
       fup(ii)=(1.0d0-ss(ii-mch))*fpo(ii-mch)+ss(ii-mch)*fup(ii)
     end do

! normalize u
     call normalize(nin)
     de=cn**2*de

! convergence test and possible exit

     if(abs(de)<max(abs(ee),0.2d0)*eps) then
       if(abs(de)<=8.0d0*epsilon(1.0d0)*max(abs(ee),0.2d0) .or. &
&         abs(de)>=0.5d0*abs(deprev)) then
         ierr = 0
         exit
       end if
       polish=.true.
     end if
     deprev=de

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
!! energy ee, counting the nodes of the large component (dirac_m)

 integer, intent(in) :: mlast
 integer, intent(out) :: nnode

 call dirac_outward(ll,kap,ee,zz,al,rr,vv,mmax,mlast,gu,fu,gup,fup,nnode)

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

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
 subroutine lschvkbb(nn,ll,nvkb,ierr,ee,emin,emax, &
&   rr,vloc,vlocpw,vkbpw,evkb,uu,up,mmax,mch)

!! Finds bound states of a  pseudopotential with
!! Vanderbilt-Kleinman-Bylander non-local projectors, to round-off
!! (lschvkbb_tol with no tolerance; the arguments are its)

 use constants_m, only: dp
 use pieces_m, only: piecewise_t
 implicit none
 external :: lschvkbb_tol

 real(dp), intent(in) :: rr(mmax),vloc(mmax),evkb(nvkb)
 integer, intent(in) :: nn,ll,nvkb,mmax
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(nvkb)
 real(dp), intent(inout) :: emin,emax,ee
 real(dp), intent(out) :: uu(mmax),up(mmax)
 integer, intent(out) :: ierr,mch

 call lschvkbb_tol(nn,ll,nvkb,ierr,ee,emin,emax, &
&                  rr,vloc,vlocpw,vkbpw,evkb,uu,up,mmax,mch,0.0d0)
 return
 end subroutine lschvkbb

 subroutine lschvkbb_tol(nn,ll,nvkb,ierr,ee,emin,emax, &
&   rr,vloc,vlocpw,vkbpw,evkb,uu,up,mmax,mch,etol)

!! Finds bound states of a  pseudopotential with
!! Vanderbilt-Kleinman-Bylander non-local projectors
!!
!! @param nn  principal quantum number
!! @param ll  angular-momentum quantum number
!! @param nvkb  = number of VKB projectors to be used
!! @param ierr  non-zero return if error
!! @param ee  bound-state energy, input guess and output calculated value
!! @param emin  externally generaated estimate of lower bound for ee
!! @param emax  externally generaated estimate of upper bound for ee
!! @param rr  log radial mesh
!! @param vloc  local part of psp on the mesh (the screened local
!!              potential, or that plus a smooth screening change)
!! @param vlocpw  the screened local potential's pieces
!! @param vkbpw  VKB projectors as functions of pieces
!! @param evkb  coefficients of BKB projectors
!! @param uu  output radial wave function (*rr)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of log grid
!! @param mch  mesh point just beyond the matching point of the inward and
!!             outward integrations: the continuous classical turning
!!             point (match_m), or, if that is less, a point match_k+2
!!             steps beyond the last kink of the potential and the
!!             projectors (a smooth maximum of the two), so that the inward
!!             integration crosses no kink
!! @param etol  0 for the energy to round-off (below); otherwise the
!!              relative energy correction at which to stop at once, for
!!              the passes of an SCF whose potential is still moving

 use constants_m, only: dp
 use match_m, only: match_k, match_point, match_weights, match_smax
 use pieces_m, only: piecewise_t
 use outwf_m, only: outwf_t, outwf_setup, outwf_solve, outwf_free
 implicit none

!Input variables
 real(dp), intent(in) :: rr(mmax),vloc(mmax),evkb(nvkb),etol
 integer, intent(in) :: nn,ll,nvkb,mmax
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(nvkb)

!Input/output variables
 real(dp), intent(inout) :: emin,emax
 real(dp), intent(inout) :: ee  !in/out, really - needs starting guess

!Output variables
 real(dp), intent(out) :: uu(mmax),up(mmax)
 integer, intent(out) :: ierr,mch


!Local variables
 real(dp) :: cn
 real(dp) :: de,deprev
 real(dp) :: eps,ro,sc
 real(dp) :: sls,sn
 real(dp) :: amesh,al,als
 real(dp) :: xkap,xmch
 real(dp) :: uo(-match_k-1:match_k),upo(-match_k-1:match_k)
 real(dp) :: ww(-match_k-1:match_k),ss(-match_k-1:match_k)
 integer :: ii,it,nin,nint,neig,lo,hi
 logical :: room
 type(outwf_t) :: ow

 real(dp), allocatable :: upp(:),cf(:)
 allocate(upp(mmax),cf(mmax))

 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)
 neig=0
 nin=mmax-4

! convergence: once the energy correction de is below eps*max(|ee|,0.2),
! passes go on at the corrected energy until de is at round-off, 8 ulps of
! max(|ee|,0.2), or stops halving (the round-off of de itself); the energy
! and function returned are those of the last pass, whose correction is
! not applied, so they are the fixed point to round-off and do not depend
! on the start guess or the path (doc/smoothness_plan.md, S-2 to S-4, and
! lschfb)
 eps=1.0d-12
 deprev=huge(1.0d0)
 ierr = 100

 sls=ll*(ll+1)


 if(ee>emax) ee=1.25d0*emax
 if(ee<emin) ee=0.75d0*emin
 if(ee>emax) ee=0.5d0*(emax+emin)

! null arrays to remove leftover garbage
 uu(:)=0.0d0
 up(:)=0.0d0
 upp(:)=0.0d0

 als=al**2

! the grids of the outward integration (outwf_m), once for this potential
 call outwf_setup(ll,rr,mmax,vloc,vlocpw,nvkb,vkbpw,ow)

! return point for bound state convergence
 do nint=1,100

! coefficient array for u in differential eq.
   do ii=1,mmax
     cf(ii)=als*sls + 2.0d0*als*(vloc(ii)-ee)*rr(ii)**2
   end do

! find classical turning point for matching
   mch=1
   do ii=mmax,2,-1
     if(cf(ii-1)<=0.d0 .and. cf(ii)>0.d0) then
       mch=ii
       exit
     end if
   end do

   if(mch==1 .and. cf(mmax)<=0.0d0) then
! classically allowed at rr(mmax): too high
     emax=0.5d0*(emax+ee)
     ee=0.5d0*(emin+emax)
     cycle
   else if(mch==1 .and. nvkb==0) then
! classically forbidden everywhere with no projectors: too low
     emin=0.5d0*(emin+ee)
     ee=0.5d0*(emin+emax)
     cycle
   end if

! the matching point: the continuous turning point, kept clear of the
! kinks by a smooth maximum so that the inward integration in the local
! potential alone is valid over the whole window
   xmch=0.0d0
   if(mch>1) call match_point(rr,cf,mmax,mch,xmch)
   if(ow%tlast>0.0d0) xmch=match_smax(xmch,ow%tlast+dble(match_k+2))
   mch=ceiling(xmch)
   lo=mch-1-match_k
   hi=mch+match_k
   room=hi+4<=mmax-4

! outward integration through the matching window, and the number of
! eigenvalues below ee on (0, rr(mch)] (outwf_index; mch is beyond the last
! kink, so beyond the projectors): a node count would also count nodes
! that enter at r = 0 with the energy
   if(room) then
     call outwf_solve(ow,ee,evkb,hi,mch,uu,up,neig)
   else
     call outwf_solve(ow,ee,evkb,mmax,mch,uu,up,neig)
   end if

! start inward integration with simple exponential where the
! wave function has decayed to machine precision
   if(neig-nn+ll+1==0 .and. room) then
     call inward_start(ll,ee,rr,vloc,mch,mmax,nin)
   end if

   if(neig-nn+ll+1==0 .and. (.not. room .or. nin<=hi+4)) then
! no room for the inward integration: too high
     emax=0.5d0*(emax+ee)
     ee=0.5d0*(emin+emax)

   else if(neig-nn+ll+1==0) then

     uo(:)=uu(lo:hi)
     upo(:)=up(lo:hi)

     xkap=sqrt(sls/rr(nin)**2 + 2.0d0*(vloc(nin)-ee))

     do ii=nin,nin+4
       uu(ii)=exp(-xkap*(rr(ii)-rr(nin)))
       up(ii)=-rr(ii)*al*xkap*uu(ii)
       upp(ii)=al*up(ii)+cf(ii)*uu(ii)
     end do

! integrate inward

     do ii=nin,lo+1,-1
       uu(ii-1)=uu(ii)+aei(up,ii)
       up(ii-1)=up(ii)+aei(upp,ii)
       do it=1,2
         upp(ii-1)=al*up(ii-1)+cf(ii-1)*uu(ii-1)
         up(ii-1)=up(ii)+aii(upp,ii)
         uu(ii-1)=uu(ii)+aii(up,ii)
       end do
     end do

! perturbation theory for the energy shift at each point of the window,
! u_o (u'_o - u'_i u_o/u_i)/2 for the normalized function, averaged with
! the bump weights; the inward solution scaled to the outward one by the
! ratio of their weighted sums, and the two blended across the window
! (match_m)

     call match_weights(xmch,lo,hi,ww,ss)
     de=0.0d0
     do ii=lo,hi
       de=de+0.5d0*ww(ii-mch)*uo(ii-mch) &
&         *(upo(ii-mch)-up(ii)*uo(ii-mch)/uu(ii))/(al*rr(ii))
     end do
     sc=sum(ww*uo)/sum(ww*uu(lo:hi))

     do ii=lo,nin
       up(ii)=sc*up(ii)
       uu(ii)=sc*uu(ii)
     end do
     do ii=lo,hi
       uu(ii)=(1.0d0-ss(ii-mch))*uo(ii-mch)+ss(ii-mch)*uu(ii)
       up(ii)=(1.0d0-ss(ii-mch))*upo(ii-mch)+ss(ii-mch)*up(ii)
     end do

! perform normalization sum

     ro=rr(1)/sqrt(amesh)
     sn=ro**(2*ll+3)/real(2*ll+3,dp)

     do ii=1,nin-3
       sn=sn+al*rr(ii)*uu(ii)**2
     end do

     sn=sn + al*(23.0d0*rr(nin-2)*uu(nin-2)**2 &
&              + 28.0d0*rr(nin-1)*uu(nin-1)**2 &
&              +  9.0d0*rr(nin  )*uu(nin  )**2)/24.0d0

! normalize u

     cn=1.0d0/sqrt(sn)

     do ii=1,nin
       up(ii)=cn*up(ii)
       uu(ii)=cn*uu(ii)
     end do
     do ii=nin+1,mmax
       uu(ii)=0.0d0
       up(ii)=0.0d0
     end do
     de=cn**2*de

! convergence test and possible exit

     if(abs(de)<max(abs(ee),0.2d0)*etol) then
       ierr = 0
       exit
     end if
     if(abs(de)<max(abs(ee),0.2d0)*eps) then
       if(abs(de)<=8.0d0*epsilon(1.0d0)*max(abs(ee),0.2d0) .or. &
&         abs(de)>=0.5d0*abs(deprev)) then
         ierr = 0
         exit
       end if
     end if
     deprev=de

     if(de>0.0d0) then
       emin=0.5d0*(emin+ee)
     else
       emax=0.5d0*(emax+ee)
     end if
     ee=ee+de
     if(ee>emax .or. ee<emin) ee=0.5d0*(emax+emin)

   else if(neig-nn+ll+1<0) then
! too few eigenvalues below ee
     emin=0.5d0*(emin+ee)
     ee=0.5d0*(emin+emax)

   else
! too many eigenvalues below ee
     emax=0.5d0*(emax+ee)
     ee=0.5d0*(emin+emax)
   end if

 end do
!fix sign to be positive at rr->oo
 if(uu(mch)<0.0d0) then
   uu(:)=-uu(:)
   up(:)=-up(:)
 end if

 call outwf_free(ow)
 deallocate(upp,cf)
 return

 contains

 include 'adams.inc'

 end subroutine lschvkbb_tol

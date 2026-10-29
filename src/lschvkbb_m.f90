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

 module lschvkbb_m

 implicit none
 private
 public :: lschvkbb,lschvkbb_tol

 contains

 subroutine lschvkbb(nn,ll,nvkb,ierr,ee,emin,emax, &
&   rr,vloc,vlocpw,vkbpw,evkb,uu,up,mmax,adams,mch,iout,ifail)

!! Finds bound states of a  pseudopotential with
!! Vanderbilt-Kleinman-Bylander non-local projectors, to round-off
!! (lschvkbb_tol with no tolerance; the arguments are its)
!!
!! @param nn  principal quantum number
!! @param ll  angular momentum
!! @param nvkb  number of projectors used
!! @param ierr  0, or non-zero after a failure; 100 not converged, 200 a
!!              state the node count does not expect (a ghost) below emin
!! @param ee  the bound-state energy (Ha): the first guess, and the
!!            eigenvalue on return
!! @param emin  an estimate of a lower bound of ee (Ha)
!! @param emax  an estimate of an upper bound of ee (Ha)
!! @param rr  log radial mesh (bohr)
!! @param vloc  the local potential on the mesh (Ha)
!! @param vlocpw  the screened local potential as pieces
!! @param vkbpw  the projectors as pieces
!! @param evkb  coefficients of the projectors (Ha)
!! @param uu  the radial wave function (times r)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param mch  the matching point (lschvkbb_tol)
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: dp
 use pieces_m, only: piecewise_t
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: nn,ll,nvkb,mmax
 real(dp), intent(in) :: rr(mmax),vloc(mmax),evkb(nvkb)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(nvkb)
 type(adams_t), intent(in) :: adams
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: emin,emax,ee

!Output variables
 real(dp), intent(out) :: uu(mmax),up(mmax)
 integer, intent(out) :: ierr,mch
 integer, intent(out) :: ifail

 ifail=0
 call lschvkbb_tol(nn,ll,nvkb,ierr,ee,emin,emax, &
&                  rr,vloc,vlocpw,vkbpw,evkb,uu,up,mmax,adams,mch,0.0d0,iout, &
&                  ifail)
 if(ifail/=0) return
 return
 end subroutine lschvkbb

 subroutine lschvkbb_tol(nn,ll,nvkb,ierr,ee,emin,emax, &
&   rr,vloc,vlocpw,vkbpw,evkb,uu,up,mmax,adams,mch,etol,iout,ifail)

!! Finds bound states of a  pseudopotential with
!! Vanderbilt-Kleinman-Bylander non-local projectors
!!
!! @param nn  principal quantum number
!! @param ll  angular momentum
!! @param nvkb  = number of VKB projectors to be used
!! @param ierr  non-zero return if error: 100 not converged; 200 every
!!              energy tried had more than nn-ll-1 states of l below it,
!!              so a state the count does not expect (a ghost) lies below
!!              emin
!! @param ee  bound-state energy, input guess and output calculated value
!! @param emin  externally generaated estimate of lower bound for ee
!! @param emax  externally generaated estimate of upper bound for ee
!! @param rr  log radial mesh (bohr)
!! @param vloc  local part of psp on the mesh (the screened local
!!              potential, or that plus a smooth screening change)
!! @param vlocpw  the screened local potential as pieces
!! @param vkbpw  the projectors as pieces
!! @param evkb  coefficients of the projectors
!! @param uu  output radial wave function (*rr)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param mch  mesh point just beyond the matching point of the inward and
!!             outward integrations: the continuous classical turning
!!             point (match_m), or, if that is less, a point match_k+2
!!             steps beyond the last kink of the potential and the
!!             projectors (a smooth maximum of the two), so that the inward
!!             integration crosses no kink
!! @param etol  0 for the energy to round-off (below); otherwise the
!!              relative energy correction at which to stop at once, for
!!              the passes of an SCF whose potential is still moving
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: dp, eig_tol, eig_floor
 use boundsolve_m, only: bs_converged, bs_correct, bs_bisect
 use match_m, only: match_k, match_point, match_weights, match_blend, &
&                    match_smax
 use pieces_m, only: piecewise_t
 use adams_m, only: adams_t
 use outwf_m, only: outwf_t, outwf_setup, outwf_solve, outwf_free
 use quad_m, only: quad_log_norm
 use inward2_m, only: inward2
 use inward_start_m, only: inward_start
 use logmesh_m, only: lm_step

!Input variables
 integer, intent(in) :: nn,ll,nvkb,mmax
 real(dp), intent(in) :: rr(mmax),vloc(mmax),evkb(nvkb),etol
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(nvkb)
 type(adams_t), intent(in) :: adams
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: emin,emax
 real(dp), intent(inout) :: ee  !in/out, really - needs starting guess

!Output variables
 real(dp), intent(out) :: uu(mmax),up(mmax)
 integer, intent(out) :: ierr,mch
 integer, intent(out) :: ifail


!Local variables
 real(dp) :: cn
 real(dp) :: de,deprev
 real(dp) :: eps,sc
 real(dp) :: sls,sn
 real(dp) :: amesh,al,als
 real(dp) :: xkap,xmch
 real(dp) :: uo(-match_k-1:match_k),upo(-match_k-1:match_k)
 real(dp) :: ww(-match_k-1:match_k),ss(-match_k-1:match_k)
 integer :: ii,nin,nint,neig,lo,hi
 logical :: room,toomany
 type(outwf_t) :: ow

 real(dp), allocatable :: upp(:),cf(:)
 ifail=0
 allocate(upp(mmax),cf(mmax))

 al=lm_step(rr,mmax)
 amesh = exp(al)
 neig=0
 nin=mmax-adams%order

! convergence: once the energy correction de is below eps*max(|ee|,0.2),
! passes go on at the corrected energy until de is at round-off, 8 ulps of
! max(|ee|,0.2), or stops halving (the round-off of de itself); the energy
! and function returned are those of the last pass, whose correction is
! not applied, so they are the fixed point to round-off and do not depend
! on the start guess or the path (doc/smoothness_plan.md, S-2 to S-4, and
! lschfb)
 eps=eig_tol
 deprev=huge(1.0d0)
 ierr = 100
 toomany=.true.

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
 call outwf_setup(ll,rr,mmax,adams,vloc,vlocpw,nvkb,vkbpw,ow,iout,ifail)
 if(ifail/=0) return

! return point for bound state convergence
 do nint=1,100

! coefficient array for u in differential eq.
   do ii=1,mmax
     cf(ii)=als*sls + 2.0d0*als*(vloc(ii)-ee)*rr(ii)**2
   end do

! find classical turning point for matching
   mch=1
   do ii=mmax,2,-1
     if(cf(ii-1)<=0.0d0 .and. cf(ii)>0.0d0) then
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
     toomany=.false.
     emin=0.5d0*(emin+ee)
     ee=0.5d0*(emin+emax)
     cycle
   end if

! the matching point: the continuous turning point, kept clear of the
! kinks by a smooth maximum so that the inward integration in the local
! potential alone is valid over the whole window
   xmch=0.0d0
   if(mch>1) then
     call match_point(rr,cf,mmax,mch,xmch,iout,ifail)
     if(ifail/=0) return
   end if
   if(ow%tlast>0.0d0) xmch=match_smax(xmch,ow%tlast+dble(match_k+2))
   mch=ceiling(xmch)
   lo=mch-1-match_k
   hi=mch+match_k
   room=hi+4<=mmax-adams%order

! outward integration through the matching window, and the number of
! eigenvalues below ee on (0, rr(mch)] (outwf_index; mch is beyond the last
! kink, so beyond the projectors): a node count would also count nodes
! that enter at r = 0 with the energy
   if(room) then
     call outwf_solve(ow,ee,evkb,hi,mch,uu,up,neig,iout=iout,ifail=ifail)
     if(ifail/=0) return
   else
     call outwf_solve(ow,ee,evkb,mmax,mch,uu,up,neig,iout=iout,ifail=ifail)
     if(ifail/=0) return
   end if

! start inward integration with simple exponential where the
! wave function has decayed to machine precision
   if(neig-nn+ll+1<=0) toomany=.false.
   if(neig-nn+ll+1==0 .and. room) then
     call inward_start(ll,ee,rr,vloc,mch,mmax,adams,nin)
   end if

   if(neig-nn+ll+1==0 .and. (.not. room .or. nin<=hi+4)) then
! no room for the inward integration: too high
     emax=0.5d0*(emax+ee)
     ee=0.5d0*(emin+emax)

   else if(neig-nn+ll+1==0) then

     uo(:)=uu(lo:hi)
     upo(:)=up(lo:hi)

     xkap=sqrt(sls/rr(nin)**2 + 2.0d0*(vloc(nin)-ee))

     do ii=nin,nin+adams%order
       uu(ii)=exp(-xkap*(rr(ii)-rr(nin)))
       up(ii)=-rr(ii)*al*xkap*uu(ii)
       upp(ii)=al*up(ii)+cf(ii)*uu(ii)
     end do

! integrate inward
     call inward2(al,cf,mmax,adams,nin,lo,uu,up,upp)

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
     call match_blend(lo,hi,mch,ss,uo,uu)
     call match_blend(lo,hi,mch,ss,upo,up)

! perform normalization sum

     sn=quad_log_norm(rr,uu(1:nin)**2,nin,ll,dble(ll+1))

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

     if(abs(de)<max(abs(ee),eig_floor)*etol) then
       ierr = 0
       exit
     end if
     if(bs_converged(de,deprev,ee,eps)) then
       ierr = 0
       exit
     end if
     deprev=de

     call bs_correct(de,ee,emin,emax,.true.)

   else if(neig-nn+ll+1<0) then
! too few eigenvalues below ee
     call bs_bisect(ee,emin,emax,.true.,.true.)

   else
! too many eigenvalues below ee
     call bs_bisect(ee,emin,emax,.false.,.true.)
   end if

 end do
 if(ierr/=0 .and. toomany) ierr=200
!fix sign to be positive at rr->oo
 if(uu(mch)<0.0d0) then
   uu(:)=-uu(:)
   up(:)=-up(:)
 end if

 call outwf_free(ow)
 deallocate(upp,cf)
 return
 end subroutine lschvkbb_tol

 end module lschvkbb_m

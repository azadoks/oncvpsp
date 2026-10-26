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
 subroutine lschpsbar(nn,ll,ierr,ee,emin,emax,rr,vv,uu,up,mmax,mbar,tht)

!! Finds bound states of a semi-local pseudopotential
!!
!! @param nn  principal quantum number
!! @param ll  angular-momentum quantum number
!! @param ierr  non-zero return if error
!! @param ee  bound-state energy, input guess and output calculated value
!!            (in/our)
!! @param emin  lower bound, potential minimum if ==0.0
!! @param emax  upper bound
!! @param rr  log radial mesh
!! @param vv  local psp
!! @param uu  output radial wave function (*rr)
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr)
!! @param mmax  size of log grid
!! @param mbar  mesh point for infinite barrier
!! @param tht  phast-shift angle for boundary condition, units of pi

 use constants_m, only: dp, pi
 use boundsolve_m, only: bs_converged, bs_correct, bs_bisect
 use quad_m, only: log_norm
 use inward2_m, only: inward2
 implicit none


!Input variables
 real(dp), intent(in) :: rr(mmax),vv(mmax)
 real(dp), intent(in) :: tht
 integer, intent(in) :: nn,ll,mmax,mbar

!Input/output variables
 real(dp), intent(inout) :: emin,emax
 real(dp), intent(inout) :: ee

!Output variables
 real(dp), intent(out) :: uu(mmax),up(mmax)
 integer, intent(out) :: ierr

!Local variables

 real(dp) :: de,deprev,uumax
 real(dp) :: eps,epsnode,sc,stht,ctht
 real(dp) :: sls,sn,cn,uout,upin,upout,xkap,xkap2,xx
 real(dp) :: amesh,al,als
 integer :: ii,imin,it,mch,mchb,nint,node,nin

 real(dp), allocatable :: upp(:),cf(:)
 allocate(upp(mmax),cf(mmax))

 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)
 stht=sin(pi*tht)
 ctht=cos(pi*tht)

 nin=0
 mch=0
 mchb=0
! convergence: once the energy correction de is below eps*max(|ee|,0.2),
! passes go on at the corrected energy until de is at round-off, 8 ulps of
! max(|ee|,0.2), or stops halving (the round-off of de itself); the energy
! and function returned are those of the last pass, whose correction is
! not applied, so they are the fixed point to round-off and do not depend
! on the start guess or the path (doc/smoothness_plan.md, S-2 to S-4, and
! lschfb)
 eps=1.0d-12
 deprev=huge(1.0d0)
! amplitude below which sign changes of uu are not counted as nodes
 epsnode=1.0d-8
 ierr = 60

 sls=ll*(ll+1)

 imin=2
 if(emin==0.0d0) then
   emin=emax
   do ii=1,mbar
     if(emin>vv(ii)+0.5d0*sls/rr(ii)**2) then
       emin=vv(ii)+0.5d0*sls/rr(ii)**2
       imin=ii
     end if
   end do
 end if
! avoid error in classical turning point search
 imin=max(imin,2)

 if(ee<emin) ee=min(emin+5.0d0, 0.5d0*(emin+emax))
 node=0
 de=0.0d0
 mch=0

! null arrays to remove leftover garbage
 uu(:)=0.0d0
 up(:)=0.0d0
 upp(:)=0.0d0
 cf(:)=0.0d0

 als=al**2
! return point for bound state convergence
 do nint=1,60

   uu(:)=0.0d0 ; up(:)=0.0d0 ; upp(:)=0.0d0

   do ii=1,mmax
     cf(ii)=als*sls + 2.0d0*als*(vv(ii)-ee)*rr(ii)**2
   end do

! find classical turning point for matching
   mch=mbar-1
   if(cf(mbar)>0.0d0) then
     do ii=mbar-1,imin,-1
       if(cf(ii)<=0.d0) then
         mch=ii
         exit
       end if
     end do
   end if

   de=0.0d0
! start wavefunction with series

   do ii=1,4
     uu(ii)=rr(ii)**(ll+1)
     up(ii)=al*(ll+1)*rr(ii)**(ll+1)
     upp(ii)=al*up(ii)+cf(ii)*uu(ii)
   end do

! outward integration using predictor once, corrector
! twice
   node=0
   uumax=0.0d0
   do ii=4,mch-1
     uu(ii+1)=uu(ii)+aeo(up,ii)
     up(ii+1)=up(ii)+aeo(upp,ii)
     do it=1,2
       upp(ii+1)=al*up(ii+1)+cf(ii+1)*uu(ii+1)
       up(ii+1)=up(ii)+aio(upp,ii)
       uu(ii+1)=uu(ii)+aio(up,ii)
     end do
     uumax=max(uumax,abs(uu(ii+1)))
   end do

   do ii=4,mch-1
     if(abs(uu(ii+1))>epsnode*uumax .or. abs(uu(ii))>epsnode*uumax) then
       if(uu(ii+1)*uu(ii) < 0.0d0) node=node+1
     end if
   end do

   mchb=mch
   uout=uu(mchb)
   upout=up(mchb)

   de=0.0d0
   if(node-nn+ll+1==0) then

! start inward integration with simple exponential where the
! wave function has decayed to machine precision, if that is inside the
! barrier

     call inward_start(ll,ee,rr,vv,mch,mmax,nin)
     if(nin<mbar-4) then
       xkap=sqrt(sls/rr(nin)**2 + 2.0d0*(vv(nin)-ee))

       do ii=nin,nin+4
         uu(ii)=exp(-xkap*(rr(ii)-rr(nin)))
         up(ii)=-rr(ii)*al*xkap*uu(ii)
         upp(ii)=al*up(ii)+cf(ii)*uu(ii)
       end do

     else

       nin=mbar

!   start inward integration at barrier
!   point with sinh or sin

       xkap2=sls/rr(mbar)**2 + 2.0d0*(vv(mbar)-ee)
       xkap=sqrt(abs(xkap2))


       do ii=mbar,mbar+4
         xx=xkap*(rr(ii)-rr(mbar))

         if(xkap2>0) then
           uu(ii)=rr(ii)*(ctht*sinh(-xx)/xkap + stht*cosh(-xx))

           up(ii)=al*rr(ii)**2*(-ctht*cosh(-xx) - xkap*stht*sinh(-xx)) &
&                 +al*uu(ii)
         else
           uu(ii)=rr(ii)*(ctht*sin(-xx)/xkap + stht*cos(-xx))

           up(ii)=al*rr(ii)**2*(-ctht*cos(-xx) + xkap*stht*sin(-xx)) &
&                 +al*uu(ii)
         end if
         upp(ii)=al*up(ii)+cf(ii)*uu(ii)
       end do

     end if !nin<mbar

! integrate inward
     call inward2(al,cf,mmax,nin,mchb,uu,up,upp)

     do ii=nin+1,mmax
       uu(ii)=0.0d0
       up(ii)=0.0d0
       upp(ii)=0.0d0
     end do

! scale outside wf for continuity

     sc=uout/uu(mchb)

     do ii=mchb,nin
       up(ii)=sc*up(ii)
       uu(ii)=sc*uu(ii)
     end do

     upin=up(mchb)

! perform normalization sum

     sn=log_norm(rr,uu(1:nin)**2,nin,ll,dble(ll+1))

! normalize u

     cn=1.0d0/sqrt(sn)
     uout=cn*uout
     upout=cn*upout
     upin=cn*upin

     do ii=1,nin
       up(ii)=cn*up(ii)
       uu(ii)=cn*uu(ii)
     end do
     do ii=nin+1,mmax
       uu(ii)=0.0d0
     end do

! perturbation theory for energy shift

     de=0.5d0*uout*(upout-upin)/(al*rr(mchb))

! convergence test and possible exit


     if(bs_converged(de,deprev,ee,eps)) then
       ierr = 0
       exit
     end if
     deprev=de

     call bs_correct(de,ee,emin,emax,.false.)
   else if(node-nn+ll+1<0) then
! too few nodes
     call bs_bisect(ee,emin,emax,.true.,.false.)
   else
! too many nodes
     call bs_bisect(ee,emin,emax,.false.,.false.)
   end if

 end do !nint

 deallocate(upp,cf)
 return

 contains

 include 'adams.inc'

 end subroutine lschpsbar

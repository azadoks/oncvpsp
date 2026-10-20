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
 subroutine lschfb(nn,ll,ierr,ee,rr,vv,uu,up,zz,mmax,mch,irel)

!! Finds bound states of an all-electron atomic potential using
!! Pauli-type  scalar-relativistic Schroedinger equation
!!
!! @param nn  principal quantum number
!! @param ll  angular-momentum quantum number
!! @param ierr  non-zero return if error, 10 if there is no bound state
!!              and the box state (u=0 at rr(mmax)) is returned
!! @param ee  bound-state energy, input guess and output calculated value
!! @param rr  log radial mesh
!! @param vv  local atomic potential
!! @param uu  output radial wave function (*rr)
!! @param up  d(uu)/dr
!! @param zz  atomic number
!! @param mmax  size of log grid
!! @param mch  matching mesh point for inward-outward integrations
!! @param irel  0 non-relativistic, 1 scalar-relativistic

 use constants_m, only: dp
 use relativity_m, only: c_light
 implicit none

!Input variables
 real(dp), intent(in) :: rr(mmax),vv(mmax)
 real(dp), intent(in) :: zz
 integer, intent(in) :: nn,ll
 integer, intent(in) :: mmax
 integer, intent(in) :: irel

!Input/output variables
 real(dp), intent(inout) :: uu(mmax),up(mmax)
 real(dp), intent(inout) :: ee
 integer, intent(inout) :: mch

!Output variables
 integer, intent(out) :: ierr

!Local variables

 real(dp) :: als
 real(dp) :: de,emax,emin,etop,elow
 real(dp) :: eps,fss,tfss,gamma,ro,sc
 real(dp) :: sls,sn,cn,uout,upin,upout,xkap
 real(dp) :: amesh,al,xx
 integer :: ii,it,nint,node,nin

 real(dp), allocatable :: upp(:),cf(:),dv(:),fr(:),frp(:),r2(:)
 allocate(upp(mmax),cf(mmax),dv(mmax),fr(mmax),frp(mmax),r2(mmax))

 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

! convergence factor for solution of schroedinger eq.  if calculated
! correction to eigenvalue is smaller in magnitude than eps times
! the magnitude of the current guess, the current guess is not changed.
 eps=1.0d-10
 ierr = 60

! check arguments
 if(ll>nn-1) then
   write(6,'(/a,i4,a,i4)') 'lschfb: ERROR ll =',ll,' > nn =',nn
   ierr=1
   return
 end if

! relativistic - non-relativistic switch
 if(irel/=0) then
   fss=(1.0d0/c_light)**2
 else
   fss=1.0d-12
 end if

 if(ll==0) gamma=sqrt(1.0d0-fss*zz**2)
 if(ll>0) gamma=(ll*sqrt(ll**2-fss*zz**2) + &
& (ll+1)*sqrt((ll+1)**2-fss*zz**2))/(2*ll+1)

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
 upp(:)=0.0d0

 als=al**2

! squared radii and dv/dr for the darwin correction, which do not depend
! on the energy
 do ii=1,mmax
   r2(ii)=rr(ii)**2
 end do

 dv(:)=0.0d0

 dv(1)=(-50.d0*vv(1)+96.d0*vv(2)-72.d0*vv(3)+32.d0*vv(4) &
&       -6.d0*vv(5))/(24.d0*al*rr(1))
 dv(2)=(-6.d0*vv(1)-20.d0*vv(2)+36.d0*vv(3)-12.d0*vv(4) &
&       +2.d0*vv(5))/(24.d0*al*rr(2))

 do ii=3,mmax-2
   dv(ii)=(2.d0*vv(ii-2)-16.d0*vv(ii-1)+16.d0*vv(ii+1) &
&        -2.d0*vv(ii+2))/(24.d0*al*rr(ii))
 end do

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

   call coefs(ee)

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

   uout=uu(mch)
   upout=up(mch)

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
       uu(ii)=exp(-xkap*(rr(ii)-rr(nin)))
       up(ii)=-rr(ii)*al*xkap*uu(ii)
       upp(ii)=(al+frp(ii))*up(ii)+(cf(ii)+fr(ii))*uu(ii)
     end do

! integrate inward

     do ii=nin,mch+1,-1
       uu(ii-1)=uu(ii)+aei(up,ii)
       up(ii-1)=up(ii)+aei(upp,ii)
       do it=1,2
         upp(ii-1)=(al+frp(ii-1))*up(ii-1)+(cf(ii-1)+fr(ii-1))*uu(ii-1)
         up(ii-1)=up(ii)+aii(upp,ii)
         uu(ii-1)=uu(ii)+aii(up,ii)
       end do
     end do

! scale outside wf for continuity

     sc=uout/uu(mch)

     do ii=mch,nin
       up(ii)=sc*up(ii)
       uu(ii)=sc*uu(ii)
     end do

     upin=up(mch)

! normalize u
     call normalize(nin)
     uout=cn*uout
     upout=cn*upout
     upin=cn*upin

! perturbation theory for energy shift

     de=0.5d0*uout*(upout-upin)/(al*rr(mch))

! convergence test and possible exit

     if(abs(de)<max(abs(ee),0.2d0)*eps) then
       ierr = 0
       exit
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

!fix sign to be positive at rr->oo
 if(uu(mch)<0.0d0) then
   uu(:)=-uu(:)
   up(:)=-up(:)
 end if

 deallocate(upp,cf,dv,fr,frp,r2)
 return

 contains

 subroutine coefs(et)

!! coefficient arrays of the radial equation at the energy et

 real(dp), intent(in) :: et

! coefficient array for u in differential eq.
 do ii=1,mmax
   cf(ii)=als*sls + 2.0d0*als*(vv(ii)-et)*r2(ii)
 end do

!  relativistic coefficient arrays for u (fr) and up (frp).
 do ii=1,mmax
   tfss=fss
   fr(ii)=als*r2(ii)*(-tfss*(vv(ii)-et)**2 + 0.5d0*tfss*dv(ii)/ &
&   (rr(ii)*(1.0d0+0.5d0*tfss*(et-vv(ii)))))
   frp(ii)=-al*rr(ii)*0.5d0*tfss*dv(ii)/(1.0d0+0.5d0*tfss*(et-vv(ii)))
 end do

 return
 end subroutine coefs

 subroutine outward(mlast,nnode)

!! outward integration from the series at the origin to mlast, at the
!! energy of the last coefs, counting the nodes

 integer, intent(in) :: mlast
 integer, intent(out) :: nnode

! start wavefunction with series

 do ii=1,4
   uu(ii)=rr(ii)**gamma
   up(ii)=al*gamma*rr(ii)**gamma
   upp(ii)=(al+frp(ii))*up(ii)+(cf(ii)+fr(ii))*uu(ii)
 end do

! outward integration using predictor once, corrector
! twice
 nnode=0

 do ii=4,mlast-1
   uu(ii+1)=uu(ii)+aeo(up,ii)
   up(ii+1)=up(ii)+aeo(upp,ii)
   do it=1,2
     upp(ii+1)=(al+frp(ii+1))*up(ii+1)+(cf(ii+1)+fr(ii+1))*uu(ii+1)
     up(ii+1)=up(ii)+aio(upp,ii)
     uu(ii+1)=uu(ii)+aio(up,ii)
   end do
   if(uu(ii+1)*uu(ii) <= 0.0d0) nnode=nnode+1
 end do

 return
 end subroutine outward

 subroutine normalize(mlast)

!! normalize u on 1..mlast, zero beyond

 integer, intent(in) :: mlast

! perform normalization sum

 ro=rr(1)/sqrt(amesh)
 sn=ro**(2.0d0*gamma+1.0d0)/(2.0d0*gamma+1.0d0)

 do ii=1,mlast-3
   sn=sn+al*rr(ii)*uu(ii)**2
 end do

 sn=sn + al*(23.0d0*rr(mlast-2)*uu(mlast-2)**2 &
&            + 28.0d0*rr(mlast-1)*uu(mlast-1)**2 &
&            +  9.0d0*rr(mlast  )*uu(mlast  )**2)/24.0d0

 cn=1.0d0/sqrt(sn)

 do ii=1,mlast
   up(ii)=cn*up(ii)
   uu(ii)=cn*uu(ii)
 end do
 do ii=mlast+1,mmax
   uu(ii)=0.0d0
 end do

 return
 end subroutine normalize

 subroutine boxstate()

!! the state with nn-ll-1 nodes and u=0 at rr(mmax), above elow (the
!! box state of RadialAtomicSolvers.jl): the outward solution over the
!! whole mesh gains a node each time u(mmax) crosses zero, so the
!! energy is bisected on its node count

 real(dp) :: eblo,ebhi
 integer :: nb,kk

 eblo=elow
 ebhi=etop+1.0d0
 do kk=1,60
   call coefs(ebhi)
   call outward(mmax,nb)
   if(nb>nn-ll-1) exit
   eblo=ebhi
   ebhi=etop+2.0d0*(ebhi-etop)
 end do

 do kk=1,200
   ee=0.5d0*(eblo+ebhi)
   call coefs(ee)
   call outward(mmax,nb)
   if(nb>nn-ll-1) then
     ebhi=ee
   else
     eblo=ee
   end if
   if(ebhi-eblo<max(abs(ee),0.2d0)*eps) exit
 end do

 ee=eblo
 call coefs(ee)
 call outward(mmax,nb)
 call normalize(mmax)
 mch=mmax

 return
 end subroutine boxstate

 include 'adams.inc'

 end subroutine lschfb

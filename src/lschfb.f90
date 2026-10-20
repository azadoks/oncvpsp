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
!! The Koelling-Harmon equation is integrated in its first-order flux form,
!! as in RadialAtomicSolvers.jl, for uu and qq = (uu' - uu/r)/M with the
!! relativistic mass M = 1 + (ee - vv)/(2c**2):
!! uu' = uu/r + M*qq,  qq' = (l(l+1)/(M*r**2) + 2*(vv - ee))*uu - qq/r.
!! It needs no dv/dr, whose finite differences of a GGA potential near the
!! nucleus made the scalar-relativistic SCF diverge on fine meshes
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
 real(dp) :: eps,fss,gamma,ro,sc
 real(dp) :: sls,sn,cn,uout,upin,upout,xkap
 real(dp) :: amesh,al,eecur
 integer :: ii,it,nint,node,nin

 real(dp), allocatable :: qq(:),qp(:),cf(:),pm(:),qm(:),xm(:),r2(:)
 allocate(qq(mmax),qp(mmax),cf(mmax),pm(mmax),qm(mmax),xm(mmax),r2(mmax))

 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

! convergence factor for solution of schroedinger eq.  if calculated
! correction to eigenvalue is smaller in magnitude than eps times
! the magnitude of the current guess, the current guess is not changed.
! 1e-10 left 3.6e-7 Ha in the 1s of U, which the eigenvalue sum carries
! into the total energy, and which depended on the SCF history
 eps=1.0d-12
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
   fss=0.0d0
 end if

! the indicial exponent of the Koelling-Harmon equation, uu ~ r**gamma
 gamma=sqrt(ll*(ll+1)+1.0d0-fss*zz**2)
! the non-relativistic equation (fss=0): u ~ r**(l+1) at the origin
 if(irel==0) gamma=ll+1

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
   upout=qq(mch)

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
       qq(ii)=-(xkap+1.0d0/rr(ii))*uu(ii)/xm(ii)
       call derivs_at(ii)
     end do

! integrate inward

     do ii=nin,mch+1,-1
       uu(ii-1)=uu(ii)+aei(up,ii)
       qq(ii-1)=qq(ii)+aei(qp,ii)
       do it=1,2
         call derivs_at(ii-1)
         uu(ii-1)=uu(ii)+aii(up,ii)
         qq(ii-1)=qq(ii)+aii(qp,ii)
       end do
       call derivs_at(ii-1)
     end do

! scale outside wf for continuity

     sc=uout/uu(mch)

     do ii=mch,nin
       up(ii)=sc*up(ii)
       uu(ii)=sc*uu(ii)
       qq(ii)=sc*qq(ii)
       qp(ii)=sc*qp(ii)
     end do

     upin=qq(mch)

! normalize u
     call normalize(nin)
     uout=cn*uout
     upout=cn*upout
     upin=cn*upin

! perturbation theory for energy shift

     de=0.5d0*uout*(upout-upin)

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

 deallocate(qq,qp,cf,pm,qm,xm,r2)
 return

 contains

 subroutine coefs(et)

!! coefficient arrays of the radial equation at the energy et

 real(dp), intent(in) :: et

 eecur=et

! coefficient array for u in the second-order equation, whose sign
! change gives the classical turning point
 do ii=1,mmax
   cf(ii)=als*sls + 2.0d0*als*(vv(ii)-et)*r2(ii)
 end do

! relativistic mass and the coefficients of the flux form,
! d(uu)/dn = al*uu + pm*qq, d(qq)/dn = qm*uu - al*qq
 do ii=1,mmax
   xm(ii)=1.0d0+0.5d0*fss*(et-vv(ii))
   pm(ii)=al*rr(ii)*xm(ii)
   qm(ii)=al*rr(ii)*(sls/(xm(ii)*r2(ii)) + 2.0d0*(vv(ii)-et))
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
   call origin_series(rr(ii),uu(ii),qq(ii))
   call derivs_at(ii)
 end do

! outward integration using predictor once, corrector
! twice
 nnode=0

 do ii=4,mlast-1
   uu(ii+1)=uu(ii)+aeo(up,ii)
   qq(ii+1)=qq(ii)+aeo(qp,ii)
   do it=1,2
     call derivs_at(ii+1)
     uu(ii+1)=uu(ii)+aio(up,ii)
     qq(ii+1)=qq(ii)+aio(qp,ii)
   end do
   call derivs_at(ii+1)
   if(uu(ii+1)*uu(ii) <= 0.0d0) nnode=nnode+1
 end do

! the first points again, inward from the outward solution: the series of
! the start is that of a model potential, and the kink where it meets the
! integration is amplified by the derivatives of a GGA density
 if(mlast>=16) then
   do ii=12,2,-1
     uu(ii-1)=uu(ii)+aei(up,ii)
     qq(ii-1)=qq(ii)+aei(qp,ii)
     do it=1,2
       call derivs_at(ii-1)
       uu(ii-1)=uu(ii)+aii(up,ii)
       qq(ii-1)=qq(ii)+aii(qp,ii)
     end do
     call derivs_at(ii-1)
   end do
 end if

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
   qq(ii)=cn*qq(ii)
   qp(ii)=cn*qp(ii)
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

 subroutine origin_series(rx,px,qx)

!! uu (px) and qq (qx) at rx from their series at the origin for the
!! potential -zz/r + d, d the screening at rr(1), at the energy of the last
!! coefs. With M = (a + b*r)/r and w = d - ee, the flux form is
!! r*P' = P + (a + b*r)*q and
!! (a + b*r)*r*q' = (l(l+1) + (a + b*r)*(-2*zz + 2*w*r))*P - (a + b*r)*q.
!! Inside zz/(2c**2) the mass is set by the nuclear potential, so r**gamma
!! alone, or with a first-order term, is not the solution at rr(1) of the
!! default mesh for light atoms; its error in the logarithmic derivative
!! acts as a contact potential (1e-5 Ha for Na). The terms are summed as
!! p_k*r**k, q_k*r**k, whose coefficients alone overflow. For l=0 the
!! series is entire; for l>0 it converges inside r = a/b ~ zz/(2c**2), and
!! beyond half of that the leading term is used, which for uu ~ r**(l+1)
!! leaves no contact error

 real(dp), intent(in) :: rx
 real(dp), intent(out) :: px,qx
 real(dp) :: aa,bb,ww,dd,pk,qk,pk1,qk1,pk2,rk
 integer :: kk

 dd=vv(1)+zz/rr(1)
 ww=dd-eecur
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
   px=px*rx**(ll+1)
   qx=qx*rx**ll
   return
 end if
 aa=0.5d0*fss*zz
 bb=1.0d0+0.5d0*fss*(eecur-dd)
 if(ll>0 .and. bb*rx>0.5d0*aa) then
   px=rx**gamma
   qx=(gamma-1.0d0)*rx**gamma/(aa+bb*rx)
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
 px=px*rx**gamma
 qx=qx*rx**gamma

 return
 end subroutine origin_series

 subroutine derivs_at(jj)

!! derivatives on the mesh index of uu and qq at the point jj

 integer, intent(in) :: jj

 up(jj)=al*uu(jj)+pm(jj)*qq(jj)
 qp(jj)=qm(jj)*uu(jj)-al*qq(jj)

 return
 end subroutine derivs_at

 include 'adams.inc'

 end subroutine lschfb

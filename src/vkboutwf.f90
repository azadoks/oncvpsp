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
 subroutine vkboutwf(ll,nvkb,ep,vkb,evkb,rr,vloc,uu,up,node,mmax,mch)

! computes Vanderbilt / Kleinman-Bylander outward-integrated wave functions

!ll  angular momentum
!nvkb  switch for 1 or 2 projedtors
!ep  energy at which wave function is to be calculated
!vkb  Vanderbilt-Kleinman-Bylander projectors for this l
!evkb  projector coefficients
!rr  log radial mesh
!vloc  local pseudopotential
!uu  wave function
!up  1st derivative of uu
!node  count of number of nodes from 0 to rr(mch)
!mmax  dimension of log mesh
!mch  index of radius to which wave function is to be integrated

 implicit none
 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp) :: rr(mmax),vloc(mmax),vkb(mmax,nvkb),evkb(nvkb)
 real(dp) :: ep
 integer nvkb,ll,mmax,mch

!Output variables
 real(dp) :: uu(mmax),up(mmax)
 integer node

!Local variables
real(dp) :: rc,rn
 real(dp), allocatable ::  phi(:,:),phip(:,:),phi0(:),phi0p(:)
 real(dp), allocatable ::  gg0(:),gg(:,:)
 integer, allocatable :: ipiv(:)

 integer ii,jj,kk,info

 uu(:)=0.0d0
 up(:)=0.0d0

! homogeneous solution, and the inhomogeneous solutions below, integrated
! together by outwf_all
 if(nvkb==0) call outwf_all(nvkb,uu,up,uu,up)

 rc=0.0d0
 if(nvkb/=0) then
! find cutoff radius for projectors
  do ii=mmax,1,-1
    if(abs(vkb(ii,1))>0.0d0) then
      rc=rr(ii)
      exit
    end if
  end do

  allocate(phi(mmax,nvkb),phip(mmax,nvkb))
  allocate(gg(nvkb,nvkb),gg0(nvkb))
  allocate(ipiv(nvkb))

  phi(:,:)=0.0d0
  phip(:,:)=0.0d0
  gg(:,:)=0.0d0
  gg0(:)=0.0d0

! inhomogeneous solutions
  call outwf_all(nvkb,uu,up,phi,phip)


! projector matrix elements and coefficient matrix
  do jj=1,nvkb
   call vpinteg(uu,vkb(1,jj),mch,2*ll+2,gg0(jj),rr)
   gg0(jj)=evkb(jj)*gg0(jj)
   do ii=1,nvkb
    call vpinteg(phi(1,ii),vkb(1,jj),mch,2*ll+2,gg(jj,ii),rr)
    gg(jj,ii)=-evkb(jj)*gg(jj,ii)
   end do
   gg(jj,jj)=1.0d0+gg(jj,jj)
  end do

! solve linear equations for coefficients

!    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )

   call dgesv(nvkb, 1, gg, nvkb, ipiv, gg0, nvkb, info)
   if(info/=0) then
    write(6,'(/a,i4)') 'vkbout: dgesv ERROR, stopping info =',info
    stop
   end if

! output wave functions
   do jj=1,nvkb
    uu(:)=uu(:)+gg0(jj)*phi(:,jj)
    up(:)=up(:)+gg0(jj)*phip(:,jj)
   end do

  deallocate(phi,phip)
  deallocate(gg,gg0)
  deallocate(ipiv)
 end if

! lower cutoff for node counting to avoid small-r noise
 rn=dmax1(0.10d0*rc, 0.05d0)

 node=0
 do ii=6,mch
! note historic evolution!
!  if(rr(ii)>0.5d0 .and. uu(ii-1)*uu(ii)<0.0d0) then
!  if(rr(ii)>0.25d0 .and. uu(ii-1)*uu(ii)<0.0d0) then
!  if(rr(ii)>0.10d0 .and. uu(ii-1)*uu(ii)<0.0d0) then
!  if(rr(ii)>0.07d0 .and. uu(ii-1)*uu(ii)<0.0d0) then
   if(rr(ii)>rn .and. uu(ii-1)*uu(ii)<0.0d0) then
    node=node+1
   end if
 end do

 return

 contains

 subroutine outwf_all(nkb,uh,uph,uk,upk)

! lschps (uh, uph) and lschkb for each of the nkb projectors (uk, upk),
! integrated in one loop so the independent predictor-corrector
! recurrences can overlap; each solution is computed exactly as those
! routines compute it

 integer :: nkb
 real(dp) :: uh(mmax),uph(mmax),uk(mmax,*),upk(mmax,*)

 real(dp) :: al,als,amesh,cn,ro,sls,sn,akb,ckb
 integer :: ii,it,kk,nsol,nn
 real(dp), allocatable :: cf(:),src(:,:),u(:,:),upr(:,:),upp(:,:)

 nsol=nkb+1
!the start-up points 1..4 are set even when mch<4
 nn=max(mch,4)
 allocate(cf(nn),src(nn,nsol),u(nn,nsol),upr(nn,nsol),upp(nn,nsol))

 al = 0.01d0 * dlog(rr(101) / rr(1))
 amesh = dexp(al)
 sls=ll*(ll+1)
 als=al**2

 do ii=1,nn
   cf(ii)=als*sls + 2.0d0*als*(vloc(ii)-ep)*rr(ii)**2
 end do

! solution 1 homogeneous (lschps), kk+1 inhomogeneous for projector kk
! (lschkb)
 src(:,1)=0.0d0
 do ii=1,4
   u(ii,1)=rr(ii)**(ll+1)
   upr(ii,1)=al*(ll+1)*rr(ii)**(ll+1)
   upp(ii,1)=al*upr(ii,1)+cf(ii)*u(ii,1)
 end do
 do kk=1,nkb
   do ii=1,nn
     src(ii,kk+1)=2.0d0*als*vkb(ii,kk)*rr(ii)**2
   end do
   ckb = 2.0d0*vkb(1,kk)/rr(1)**(ll+1)
   akb = ckb/(6.0d0+4.0d0*ll)
   do ii=1,4
     u(ii,kk+1)=akb*rr(ii)**(ll+3)
     upr(ii,kk+1)= al*(ll+3)*u(ii,kk+1)
     upp(ii,kk+1)= als*(ll+3)**2*u(ii,kk+1)
   end do
 end do

 do ii=4,mch-1
   do kk=1,nsol
     u(ii+1,kk)=u(ii,kk)+aeo(upr(1,kk),ii)
     upr(ii+1,kk)=upr(ii,kk)+aeo(upp(1,kk),ii)
   end do
   do it=1,2
     upp(ii+1,1)=al*upr(ii+1,1)+cf(ii+1)*u(ii+1,1)
     do kk=2,nsol
       upp(ii+1,kk)=al*upr(ii+1,kk)+cf(ii+1)*u(ii+1,kk) &
&                  + src(ii+1,kk)
     end do
     do kk=1,nsol
       upr(ii+1,kk)=upr(ii,kk)+aio(upp(1,kk),ii)
       u(ii+1,kk)=u(ii,kk)+aio(upr(1,kk),ii)
     end do
   end do
 end do

! normalize the homogeneous solution as lschps does
 ro=rr(1)/dsqrt(amesh)
 sn=ro**(2*ll+3)/(2*ll+3)
 do ii=1,mch-3
   sn=sn+al*rr(ii)*u(ii,1)**2
 end do
 sn =sn + al*(23.0d0*rr(mch-2)*u(mch-2,1)**2 &
&           + 28.0d0*rr(mch-1)*u(mch-1,1)**2 &
&          +  9.0d0*rr(mch  )*u(mch  ,1)**2)/24.0d0
 cn=1.0d0/dsqrt(sn)

 uh(:)=0.0d0
 uph(:)=0.0d0
 do ii=1,mch
   uph(ii)=cn*upr(ii,1)
   uh(ii)=cn*u(ii,1)
 end do
 do kk=1,nkb
   uk(:,kk)=0.0d0
   upk(:,kk)=0.0d0
   uk(1:mch,kk)=u(1:mch,kk+1)
   upk(1:mch,kk)=upr(1:mch,kk+1)
 end do

 deallocate(cf,src,u,upr,upp)
 end subroutine outwf_all

 include 'adams.inc'

 end subroutine vkboutwf

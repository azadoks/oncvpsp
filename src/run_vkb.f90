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
 subroutine run_vkb(lmax,lloc,lpopt,dvloc0,irc,nproj,rr,mmax,mxprj,pswf, &
&                   vfull,vp,evkb,vkb,nlim,vr,irel,nkap)

!! computes Vanderbilt / Kleinman-Bylander non-local potentials
!!
!! @param lmax  maximum angular momentum
!! @param lloc  l for local potential (lloc==4 => use linear combination)
!! @param lpopt  choice of polynomial for lloc==4
!! @param dvloc0  amplitude at rr==0 to be smoothely added for lloc==4
!! @param irc  core radii indices
!! @param nproj  number of projectors for each l
!! @param rr  log radial grid
!! @param mmax  size of radial grid
!! @param mxprj  dimension of number of projectors
!! @param pswf  pseudo wave functions
!! @param vfull  all-electron potential
!! @param vp  semi-local pseudopotentials for first projectors
!!            (vp(:,5) is local potential if lloc=4)
!! @param vkb  semi-local "potentials"*pswf (input); VKB projectors (output)
!! @param evkb  coefficients of VKB projectors (output)
!! @param nlim  index of maximum rc including that of vlocal (output)
!! @param vr  effective scalar-relativistic "potential" calculated in vrel
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param nkap  number of kappa per l: 1, or 2 for irel=2
!!              (for pswf, vp, vkb, evkb and vr the last index is 1 for j=l+1/2,
!!              2 for j=l-1/2 when nkap=2)

 use constants_m, only: dp, mxl, mxlin
 implicit none

!Input variables
 integer, intent(in) :: lmax,lloc,lpopt,mmax,mxprj,irel,nkap
 integer, intent(in) :: irc(mxlin+1),nproj(mxlin+1)
 real(dp), intent(in) :: rr(mmax),pswf(mmax,mxprj,mxl+1,nkap),vfull(mmax)
 real(dp), intent(in) :: vr(mmax,mxprj,mxlin+1,nkap)
 real(dp), intent(in) :: dvloc0

!Input/output variables
 real(dp), intent(inout) :: vp(mmax,mxl+2,nkap)
 real(dp), intent(inout) :: vkb(mmax,mxprj,mxl+1,nkap)

!Output variables
 real(dp), intent(out) :: evkb(mxprj,mxl+1,nkap)
 integer, intent(out) :: nlim

!Local variables
 integer :: lwork
 integer :: ii,ipk,jj,kk,ikap,ll,l1,info,np,kap,mkap
 real(dp) :: apk,sn,tt,qq12,xx,ff
 real(dp) :: bb(mxprj,mxprj),bbev(mxprj),bbi(mxprj,mxprj)
 real(dp) :: sovl(mxprj,mxprj),sovlev(mxprj),smhalf(mxprj,mxprj)
 real(dp) :: sphalf(mxprj,mxprj)
 real(dp) :: bbist(mxprj,mxprj),bbistev(mxprj),bbit(mxprj,mxprj)
 real(dp), allocatable :: vloc(:),vkbt(:,:),vkbst(:,:),work(:)
 logical :: sorted

 evkb(:,:,:)=0.0d0

 lwork=3*mxprj
 allocate(vloc(mmax),vkbt(mmax,mxprj),vkbst(mmax,mxprj),work(lwork))


! calculate local potential
! for nkap=2 use "scalar-relativistic" weighted average if real ll/=0 is
! specified, and reset semi-local potential to this average (ie., there is
! no SO for this l)

 if(lloc==4) then
   call vploc(rr,vfull,vp(1,1,1),dvloc0,irc(5),mmax,lpopt)
   vloc(:)=vp(:,5,1)
   if(nkap==2) vp(:,5,2)=vp(:,5,1)
 else if(nkap==1) then
   vloc(:)=vp(:,lloc+1,1)
 else
   l1=lloc+1
   vloc(:)=(dble(l1)*vp(:,l1,1) + dble(l1-1)*vp(:,l1,2))/dble(2*l1-1)
 end if

! a relativistic semi-local potential is vae + vr at rc and vae beyond: the
! taper below takes out vr as it does from the projectors, so that vloc,
! and with it every projector, is continuous at rc (it jumped by vr(rc),
! about 1e-4 Ha, which the radial solvers integrated across at O(al))
 if(lloc<4) then
   l1=lloc+1
   do ii=1,irc(l1)
     xx=20.0d0*((rr(ii)/rr(irc(l1)))-1.0d0)
     if(xx>-1.0d0) then
       ff=(1.0d0-xx**2)**2
       if(nkap==1) then
         tt=vr(ii,1,l1,1)
       else
         tt=(dble(l1)*vr(ii,1,l1,1) + dble(l1-1)*vr(ii,1,l1,2))/dble(2*l1-1)
       end if
       vloc(ii)=vloc(ii)-ff*tt
     end if
   end do
   vp(:,l1,1)=vloc(:)
   if(nkap==2) vp(:,l1,2)=vloc(:)
 end if


! Vanderbilt / Kleinman-Bylander projector construction
 nlim=irc(lloc+1)
 do l1=1,lmax+1
   if(l1==lloc+1) cycle
   ll=l1-1
   if(nkap==1 .or. l1==1) then
     mkap=1
   else
     mkap=2
   end if
! loop on J = ll +/- 1/2
   do ikap=1,mkap
     if(ikap==1) kap=-(ll+1)
     if(ikap==2) kap=  ll

     nlim=max(nlim,irc(l1))
     do jj=1,nproj(l1)
!incorporates scalar-relativistic "potential" in last 5% of [0:rc] to smoothly
!force projectors to zero.
       do ii=1,irc(l1)
         xx=20.0d0*((rr(ii)/rr(irc(l1)))-1.0d0)
         if(xx>-1.0d0) then
           ff=(1.0d0-xx**2)**2
         else
           ff=0.0d0
         end if

         vkb(ii,jj,l1,ikap)=vkb(ii,jj,l1,ikap) &
&                           -(vloc(ii)+ff*vr(ii,jj,l1,ikap))*pswf(ii,jj,l1,ikap)
       end do
       do ii=irc(l1)+1,mmax
         vkb(ii,jj,l1,ikap)=0.0d0
       end do
     end do
   end do !ikap
 end do !l1

! Vanderbilt B-matrix construction

 do l1=1,lmax+1
   if(l1==lloc+1) cycle
   if(nproj(l1)==0) cycle
   ll=l1-1
   if(nkap==1 .or. l1==1) then
     mkap=1
   else
     mkap=2
   end if
! loop on J = ll +/- 1/2
   do ikap=1,mkap
     if(ikap==1) kap=-(ll+1)
     if(ikap==2) kap=  ll
     np=nproj(l1)
     bb(:,:)=0.0d0
     call vpinteg(pswf(1,1,l1,ikap),vkb(1,1,l1,ikap),irc(l1),2*l1,bb(1,1),rr)
     evkb(1,l1,ikap)=1.0d0/bb(1,1)
! this is Kleinman-Bylander result if nproj==1

     if(nproj(l1)>=2) then

       do jj=1,np
         do ii=1,np
           call vpinteg(pswf(1,ii,l1,ikap),vkb(1,jj,l1,ikap),irc(l1),2*l1, &
&                       bb(ii,jj),rr)
         end do
       end do

! symmetrize exactly
       if(irel==2) then
         write(6,'(/a,i3,a,i3)') 'B matrix Hermiticity error, ll=',l1-1, &
&              '  kap=',kap
       else
         write(6,'(/a,i4)') 'B matrix Hermiticity error, ll=',l1-1
       end if
       do jj=2,np
         do ii=1,jj-1
           if(irel==2) then
             write(6,'(2i4,1p,d14.4)') ii,jj,bb(ii,jj)-bb(jj,ii)
           else
             write(6,'(2i4,1p,e14.4)') ii,jj,bb(ii,jj)-bb(jj,ii)
           end if
           tt=0.5d0*(bb(ii,jj)+bb(jj,ii))
           bb(ii,jj)=tt
           bb(jj,ii)=tt
         end do
       end do

! find eigenvalues and eigenvectors of the B matrix

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

       call dsyev( 'V', 'U', np, bb, mxprj, bbev, work, lwork, info )
       if(info /= 0) then
         write(6,'(a,i4)') 'run_vkb: B matrix eigenvalue ERROR, info=',info
         error stop 1
       end if

! take linear combinations to form diagonal projectors

       do jj=1,np
         vkbt(:,jj)=vkb(:,jj,l1,ikap)
       end do

       do jj=1,np
         vkb(:,jj,l1,ikap)=0.0d0
         evkb(jj,l1,ikap)=1.0d0/bbev(jj)
         do ii=1,np
           vkb(:,jj,l1,ikap)=vkb(:,jj,l1,ikap)+bb(ii,jj)*vkbt(:,ii)
         end do
       end do

     end if !nproj>=2

! normalize projectors
     do jj=1,np
       call vpinteg(vkb(1,jj,l1,ikap),vkb(1,jj,l1,ikap),irc(l1),2*l1,sn,rr)
       if(vkb(irc(l1)-10,jj,l1,ikap)>=0.0d0) then
         vkb(:,jj,l1,ikap)=vkb(:,jj,l1,ikap)/sqrt(sn)
       else
         vkb(:,jj,l1,ikap)=-vkb(:,jj,l1,ikap)/sqrt(sn)
       end if
       evkb(jj,l1,ikap)=sn*evkb(jj,l1,ikap)
     end do


! calculate overlap
     if(nproj(l1)>=2) then


       bbi(:,:)=0.0d0
       do jj=1,np
         sovl(jj,jj)=1.0d0
         bbi(jj,jj)=evkb(jj,l1,ikap)
         do ii=jj+1,np
           call vpinteg(vkb(1,ii,l1,ikap),vkb(1,jj,l1,ikap),irc(l1),2*l1,sn,rr)
           sovl(ii,jj)=sn
           sovl(jj,ii)=sn
         end do
       end do
       call vpinteg(vkb(1,1,l1,ikap),vkb(1,2,l1,ikap),irc(l1),2*l1,sn,rr)

! find eigenvalues and eigenvectors of the overlap matrix

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

       call dsyev( 'V', 'U', np, sovl, mxprj, sovlev, work, lwork, info )
       if(info /= 0) then
         write(6,'(a,i4)') 'run_vkb: sovl matrix eigenvalue ERROR, info=',info
         error stop 1
       end if

! construct S^(-1/2) AND s^(1/2)

       do jj=1,np
         tt=sqrt(sovlev(jj))
         do ii=1,np
           sphalf(ii,jj)=tt*sovl(ii,jj)
           smhalf(ii,jj)=sovl(ii,jj)/tt
         end do
       end do

! construct B^(-1)* = S^(1/2)^T B^(-1) S^(1/2)

       bbit(:,:)=0.0d0
       bbist(:,:)=0.0d0

       do ii=1,np
         do jj=1,np
           do kk=1,np
             bbit(ii,jj)=bbit(ii,jj)+bbi(ii,kk)*sphalf(kk,jj)
           end do
         end do
       end do

       do ii=1,np
         do jj=1,np
           do kk=1,np
             bbist(ii,jj)=bbist(ii,jj)+bbit(kk,jj)*sphalf(kk,ii)
           end do
         end do
       end do

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

       call dsyev( 'V', 'U', np , bbist, mxprj, bbistev, work, lwork, info )
       if(info /= 0) then
         write(6,'(a,i4)') &
&              'run_vkb: B^(-1)* matrix eigenvalue ERROR, info=',info
         error stop 1
       end if


! take linear combinations to form orthonormal basis functions

       vkbt(:,:)=vkb(:,:,l1,ikap)
       vkbst(:,:)=0.0d0

       do jj=1,np
         do ii=1,np
           vkbst(:,jj)=vkbst(:,jj) + smhalf(ii,jj)*vkbt(:,ii)
         end do
       end do


! take linear combinations to form orthonormal projectors

       vkb(:,:,L1,ikap)=0.0d0

       do ii=1,np
         evkb(ii,l1,ikap)=bbistev(ii)
         do jj=1,np
           vkb(:,ii,l1,ikap)=vkb(:,ii,l1,ikap) + bbist(jj,ii)*vkbst(:,jj)
         end do
       end do

! re-order if necessary so first projector has largest magnitude coefficient
! which is consistent with relativistic sr_so_r.f90 output

       do ii=1,100
         sorted=.true.
         do jj=2,np

           if(abs(evkb(jj,l1,ikap))>abs(evkb(jj-1,l1,ikap))) then
             tt=evkb(jj-1,l1,ikap)
             vkbt(:,1)=vkb(:,jj-1,l1,ikap)
             evkb(jj-1,l1,ikap)=evkb(jj,l1,ikap)
             vkb(:,jj-1,l1,ikap)=vkb(:,jj,l1,ikap)
             evkb(jj,l1,ikap)=tt
             vkb(:,jj,l1,ikap)=vkbt(:,1)
             sorted=.false.
           end if
         end do
         if(sorted) exit
       end do

       write(6,'(/a,1p,5e12.4)') '  Orthonormal projector coefficients',&
&          (evkb(jj,l1,ikap),jj=1,np)

! Set sign of projectors (physically irrelevant) so that they are positive
! at their peak (needed for compaisons apparently)

       do jj=1,np
         apk=0.0d0
         ipk=1
         do ii=1,mmax
           if(abs(vkb(ii,jj,l1,ikap))>apk) then
             apk=abs(vkb(ii,jj,l1,ikap))
             ipk=ii
           end if
         end do
         if(vkb(ipk,jj,l1,ikap)<0.0d0) then
           vkb(:,jj,l1,ikap)=-vkb(:,jj,l1,ikap)
         end if
       end do

     end if !nproj(l1)>=2

     if(nproj(l1)<mxprj) then
       do jj=nproj(l1)+1, mxprj
         evkb(jj,l1,ikap)=0.0d0
         vkb(:,jj,l1,ikap)=0.0d0
       end do
     end if
   end do !ikap
 end do !l1

 deallocate(vloc,vkbt,vkbst,work)

 return
 end subroutine run_vkb

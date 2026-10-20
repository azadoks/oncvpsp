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
 subroutine run_plot(lmax,npa,epa,lloc,irc, &
&                    vkb,evkb,nproj,rr,vfull,vp,vpuns,zz,mmax,mxprj,drl,nrl, &
&                    rho,rhoc,rhomod,irel,nkap,cvgplt)

!! write output for plotting pseudopotentials, model core charge, and
!! all-electron and pseudo wave functions
!!
!! @param lmax  maximum angular momentum
!! @param npa  principal quantum number for corresponding all-electron state
!! @param epa  bound-state or scattering state reference energies for vkb
!!             potentials (last index 1 for j=l+1/2, 2 for j=l-1/2, as in vkb,
!!             evkb, vp, cvgplt)
!! @param lloc  l for local potential
!! @param irc  indices of core radii
!! @param vkb  VKB projectors
!! @param evkb  coefficients of VKB projectors
!! @param nproj  number of vkb projectors for each l
!! @param rr  log radial grid
!! @param vfull  all-electron potential
!! @param vp  semi-local pseudopotentials
!! @param vpuns  unscreened vp
!! @param zz  atomic number
!! @param mmax  size of radial grid
!! @param mxprj  dimension of number of projectors
!! @param drl  spacing of linear radial mesh
!! @param nrl  number of points in linear radial mesh
!! @param rho  valence pseudocharge
!! @param rhoc  all-electron core charge
!! @param rhomod  model core charge
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param nkap  number of kappa per l: 1, or 2 for irel=2
!! @param cvgplt  energy per electron error vs. cutoff (convergence profiles)

 use constants_m, only: dp, mxl, mxlin, nrhomod, nerrlv
 implicit none

!Input variables
 integer, intent(in) :: lmax,lloc,mmax,mxprj,nrl,irel,nkap
 integer, intent(in) :: npa(mxprj,mxlin+1),irc(mxlin+1),nproj(mxlin+1)
 real(dp), intent(in) :: zz,drl
 real(dp), intent(in) :: rr(mmax),vp(mmax,mxl+2,nkap),vpuns(mmax,mxl+2)
 real(dp), intent(in) :: vfull(mmax)
 real(dp), intent(in) :: vkb(mmax,mxprj,mxl+1,nkap)
 real(dp), intent(in) :: rho(mmax),rhoc(mmax),rhomod(mmax,nrhomod)
 real(dp), intent(in) :: epa(mxprj,mxlin+1,nkap),evkb(mxprj,mxl+1,nkap)
 real(dp), intent(in) :: cvgplt(2,nerrlv,mxprj,mxl+1,nkap)

!Local variables
 integer :: ll,l1,ii,jj,ierr,mch,n1,n2,n3,n4,nn
 integer :: iprj,npr,lout,iout
 integer :: ikap,kap,mkap
 real(dp) :: al,cnorm,emax,emin,etest,rmx,sgnae,sgnps
 real(dp) :: r0,dr
 real(dp), allocatable :: u2(:),up(:),ur(:,:),urp(:,:)

 allocate(u2(mmax),up(mmax),ur(mmax,2),urp(mmax,2))

 al = 0.01d0 * log(rr(101)/rr(1))

 write(6,'(/a)') 'DATA FOR PLOTTING'

 n1 = int(log(drl/rr(1))/al+1.0d0)
 n2 = int(log(real(nrl,dp)*drl/rr(1))/al+1.0d0)
 n3=0
 do l1=1,lmax+1
   n3=max(n3,irc(l1)-1)
 end do
 n4=min(n2,int(log((rr(n3)+1.0d0)/rr(1))/al))
 n3=int(log(1.1d0*rr(n3)/rr(1))/al+1.0d0)

!the fully-relativistic code writes every point of the potentials
 dr=min(drl,rr(n4)/200)
 write(6,'(/2a/)') ' radii, charge,',' pseudopotentials (ll=0, 1, lmax)'
 r0=rr(n1)-dr
 do ii=n1,n4
   if(rr(ii)>r0 .or. irel==2) then
     write(6,'(a,6(f12.7,1x))') '!p',rr(ii),rho(ii),(vpuns(ii,l1),l1=1,lmax+1)
     r0=rr(ii)+dr
   end if
 end do

 if(lloc == 4) then
   r0=rr(n1)-dr
   do ii=n1,n4
     if(rr(ii)>r0) then
       write(6,'(a, 2(f12.7,1x))') ' !L',rr(ii),vpuns(ii,lloc + 1)
       r0=rr(ii)+dr
     end if
   end do
 end if

! output for analysis of core charge and models
 rmx=0
 do ii=n1,n2
   rmx=max(rmx,5.0d0*rho(ii),rhomod(ii,1))
 end do

 write(6,'(//2a/)') ' radii, charge,',' core charge, model core charge'
 dr=min(drl,rr(n2)/200.0d0)
 r0=rr(n1)-dr
 do ii=n1,n2
   if(rr(ii)>r0) then
     if(rhoc(ii)<rmx)then
       write(6,'(a,8(f12.7,1x))') '!r',rr(ii),rho(ii),rhoc(ii), &
&       rhomod(ii,1)
     else if(rr(ii)>0.01d0) then
       write(6,'(a,8(f12.7,1x))') '!r',rr(ii),rho(ii),rmx, &
&      rhomod(ii,1)
     end if
     r0=rr(ii)+dr
   end if
 end do

! loop for wave function output, the all-electron functions (their large
! components when irel=2) normalized as Schroedinger functions; for
! irel=2, the j=l+1/2 functions and projectors are labelled with -iprj
! and -l

 do l1 = 1, lmax + 1
   ll = l1 - 1
   npr=nproj(l1)
   if(l1==lloc+1) npr=0
   if(nkap==1 .or. ll==0) then
     mkap=1
   else
     mkap=2
   end if
! loop on J = ll +/- 1/2 when nkap=2
   do ikap=1,mkap
     if(ikap==1) then
       kap=-(ll+1)
     else
       kap=  ll
     end if

     do iprj=1,nproj(l1)
       if(epa(iprj,l1,ikap)<0.0d0) then
         if(irel==2) then
           write(6,'(//a,i2,a,i2,a,i2,a,a/)') 'n=',npa(iprj,l1),',  l=', &
&                ll,', kap=',kap, ', all-electron wave function', &
&                ', pseudo w-f'
         else
           write(6,'(//a,i2,a,i2,a,a/)') 'n=',npa(iprj,l1),',  l=',ll, &
&                ', all-electron wave function', ', pseudo w-f'
         end if

         etest=epa(iprj,l1,ikap)
         call boundstate(irel,npa(iprj,l1),ll,kap,ierr,etest, &
&                        rr,zz,vfull,ur,urp,mmax,mch)
         if(ierr /= 0 .and. irel==2) then
           write(6,'(a,i2,a,i2,a,i2,a,i4)') &
&                'run_plot_r 409: ldiracfb convergence ERROR, n=',&
&                 npa(iprj,l1),'  l=',ll, 'kap-',kap,'  ierr=',ierr
           error stop 1
         else if(ierr /= 0) then
           write(6,'(/a,3i4)') &
&                'run_plot: lschfb convergence ERROR n,l,ierr=', &
&                npa(iprj,l1),ll,ierr
           error stop 1
         end if
         if(irel==2) call renorm_r(ur,rr,ll,kap,zz,mmax,cnorm)

         emax=0.9d0*etest
         emin=1.1d0*etest

         call lschvkbb(ll+iprj,ll,npr,ierr,etest,emin,emax, &
&                      rr,vp(1,lloc+1,ikap),vkb(1,1,l1,ikap), &
&                      evkb(1,l1,ikap),u2,up,mmax,mch)
         sgnae=1.0d0
         sgnps=1.0d0
       else
         if(irel==2) then
           write(6,'(//a,i2,a,i2,a,i2,a,a/)') 'scattering, iprj=',iprj, &
&               ',  l=',ll,', kap=',kap,', all-electron wave function', &
&               ', pseudo w-f'
         else
           write(6,'(//a,i2,a,i2,a,a/)') 'scattering, iprj=',iprj, &
&               ',  l=',ll,', all-electron wave function', ', pseudo w-f'
         end if
         call scatterstate(irel,nn,ll,kap,ierr,epa(iprj,l1,ikap), &
&                       rr,zz,vfull,ur,urp,mmax,n2)
         if(irel==2) then
           if(ierr/=0) then
             write(6,'(/2a,i0,a,i0,a,i0)') &
&                  'run_plot_r: ERROR ldiracfs failed', &
&                  ', ll=',ll,', kap=',kap,', ierr=',ierr
             error stop 1
           end if
           call renorm_r(ur,rr,ll,kap,zz,mmax,cnorm)
         end if

         call lschvkbs(ll,npr,epa(iprj,l1,ikap),rr,vp(1,lloc+1,ikap), &
&                      vkb(1,1,l1,ikap),evkb(1,l1,ikap),u2,up,mmax,n2)
         sgnae=sign(1.0d0,ur(n2,1))
         sgnps=sign(1.0d0,u2(n2))
       end if

!the fully-relativistic code starts a step further out
       dr=min(drl,rr(n2)/200.0d0)
       if(irel==2) then
         r0=rr(n1)+dr
       else
         r0=rr(n1)-dr
       end if
       iout=iprj
       if(nkap==2 .and. ikap==1) iout=-iprj
       do ii = n1, n2
         if(rr(ii)>r0) then
           write(6,'(a,i5,i1,3(f12.6,1x))') '&',iout,ll,rr(ii), &
&                sgnae*ur(ii,1),sgnps*u2(ii)
           r0=rr(ii)+dr
         end if
       end do
     end do !iproj

! orthonormal projector plots

     write(6,'(/)')
     dr=min(drl,rr(n3)/200.0d0)
     r0=rr(n1)-dr
     lout=ll
     if(nkap==2 .and. ikap==1) lout=-ll
     do ii=n1,n3
       if(rr(ii)>r0) then
         write(6,'(a,i6,6(f12.6,1x))') '!J',lout,rr(ii), &
&              (vkb(ii,jj,l1,ikap),jj=1,nproj(l1))
         r0=rr(ii)+dr
       end if
     end do
   end do !ikap
 end do !l1
!
! convergence profile plots

 if(irel==2) then
! the average over kappa for l>0
   write(6,'(/a)') 'convergence profiles, (ll=0,lmax), kappa average'
   l1=1
   ll=l1-1
   do jj=1,nerrlv
     if(cvgplt(1,jj,1,l1,1)/=0.d0) then
       write(6,'(a,i6,3f12.6)') '!C',ll,cvgplt(1,jj,1,l1,1),cvgplt(2,jj,1,l1,1)
     end if
   end do
   do l1=2,lmax+1
     ll=l1-1
     do jj=1,nerrlv
       if(cvgplt(1,jj,1,l1,1)/=0.d0 .and. cvgplt(1,jj,1,l1,nkap)/=0.0d0) then
         write(6,'(a,i6,3f12.6)') '!C',ll,cvgplt(1,jj,1,l1,1), &
&         0.5d0*(cvgplt(2,jj,1,l1,1)+cvgplt(2,jj,1,l1,nkap))
       end if
     end do
   end do
 else
   write(6,'(/a)') 'convergence profiles, (ll=0,lmax)'
   write(6,*) 'lmax',lmax
   do l1=1,lmax+1
     ll=l1-1
     do jj=1,nerrlv
       if(cvgplt(1,jj,1,l1,1)/=0.d0) then
         write(6,'(a,i6,3(f12.6,1x))') '!C',ll,cvgplt(1,jj,1,l1,1), &
&              cvgplt(2,jj,1,l1,1)
       end if
     end do
   end do
 end if

 deallocate(u2,up,ur,urp)
 return
 end subroutine run_plot

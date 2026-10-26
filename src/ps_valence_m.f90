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
!> The valence of oncvpsp_main's pseudopotential: the all-electron and
!> pseudo valence states of the reference configuration, their densities
!> and eigenvalue sum, and the unscreening of the semi-local potentials.
 module ps_valence_m
 use constants_m, only: dp, mxprj, mxstate, mxl, mxlin, mxkap
 use kappa_m, only: kap_count, kap_of, j_share
 use pieces_m, only: piecewise_t
 implicit none
 private
 public :: valence_densities, unscreen_potentials

 contains

 subroutine valence_densities(irel,nkap,nc,nv,na,la,fa,ea,it,lloc,nproj, &
&                             rr,zz,vfull,vp,vlocpw,vkbpw,evkb,mmax,uu,up, &
&                             uupsa,eeig,zval,rho,rhops,rhoae,rhotae)
!! the valence states, all-electron (in vfull) and pseudo (with the VKB
!! projectors, for the most consistent unscreening), state by state: their
!! densities rhoae and rhops (of the share of each J when fully
!! relativistic), the pseudo-waves uupsa for upfout, the all-electron
!! valence density rhotae, the pseudo one rho, the valence charge zval and
!! the sum of the pseudo eigenvalues eeig
 integer, intent(in) :: irel,nkap,nc,nv,na(mxstate),la(mxstate),it,lloc
 integer, intent(in) :: nproj(mxlin+1),mmax
 real(dp), intent(in) :: fa(mxstate),ea(mxstate,mxkap),rr(mmax),zz
 real(dp), intent(in) :: vfull(mmax),vp(mmax,mxl+2,nkap)
 real(dp), intent(in) :: evkb(mxprj,mxl+1,nkap)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1,nkap)
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)
 real(dp), intent(out) :: uupsa(mmax,nkap,nv),eeig,zval,rho(mmax)
 real(dp), intent(out) :: rhops(mmax,nv),rhoae(mmax,nv),rhotae(mmax)
 integer :: nodes(mxl+1),kk,ll,l1,mkap,ikap,kap,ierr,mch
 real(dp) :: fj,et,emin,emax

! null charge and eigenvalue accumulators
 uupsa(:,:,:)=0.0d0
 eeig=0.0d0
 zval=0.0d0
 rhops(:,:)=0.0d0
 rho(:)=0.0d0
 nodes(:)=0
 rhoae(:,:)=0.0d0
 rhotae(:)=0.0d0
 do kk=1,nv
   ll=la(nc+kk)
   l1=ll+1

   mkap=kap_count(nkap,ll)
   do ikap=1,mkap
     kap=kap_of(ll,ikap)
     fj=1.0d0
     if(nkap==2) fj=j_share(ll,ikap,1.0d0)

     et=ea(nc+kk,ikap)

     call boundstate(irel,na(nc+kk),ll,kap,ierr,et, &
&                    rr,zz,vfull,uu,up,mmax,mch)
     if(ierr /= 0) then
       if(irel==2) then
         write(6,'(a,i2,a,i2,a,i2,a,i4)') &
&              'oncvpsp_r 541: ldiracfb convergence ERROR, n=',&
&               na(nc+kk),'  l=',ll, '  kap-',kap,'  ierr=',ierr
       else
         write(6,'(/a,3i4)') &
&              'oncvpsp-461: lschfb convergence ERROR n,l,iter=', &
&              na(nc+kk),la(nc+kk),it
       end if
       error stop 1
     end if

     if(irel==2) then
       rhoae(:,kk)=rhoae(:,kk)+fj*((uu(:,1)/rr(:))**2 + (uu(:,2)/rr(:))**2)
       emax=0.90d0*et
       emin=1.10d0*et
     else
       rhoae(:,kk)=(uu(:,1)/rr(:))**2
       emax=0.75d0*et
       emin=1.25d0*et
     end if

     call lschvkbb(ll+nodes(l1)+1,ll,nproj(l1),ierr,et,emin,emax, &
&                  rr,vp(1,lloc+1,ikap),vlocpw,vkbpw(1:nproj(l1),l1,ikap), &
&                  evkb(1,l1,ikap),uu,up,mmax,mch)
     if(ierr /= 0) then
       if(irel==2) then
         write(6,'(a,i2,a,i2)') 'oncvpsp_r: lschvkbb convergence ERROR, n=',&
&         ll+1,'  l=',ll
       else
         write(6,'(a,3i4)') 'oncvpsp: lschvkbb ERROR',ll+nodes(l1)+1,ll,ierr
         if(ierr==200) write(6,'(a,i2,a,f12.6,a/a)') &
&             'a state of l=',ll,' lies below',emin, &
&             ' Ha that the node count does not expect:', &
&             'a bound ghost of the projectors'
         flush(6)
       end if
       error stop 1
     end if

! save valence pseudo wave functions for upfout
     uupsa(:,ikap,kk)=uu(:,1)

     if(irel==2) then
       rhops(:,kk)=rhops(:,kk)+fj*(uu(:,1)/rr(:))**2
     else
       rhops(:,kk)=(uu(:,1)/rr(:))**2
     end if
     eeig=eeig+fa(nc+kk)*fj*et

   end do !ikap
   rhotae(:)=rhotae(:)+fa(nc+kk)*rhoae(:,kk)
   rho(:)=rho(:)+fa(nc+kk)*rhops(:,kk)
   zval=zval+fa(nc+kk)
   nodes(l1)=nodes(l1)+1
 end do !kk
 end subroutine valence_densities

 subroutine unscreen_potentials(irel,lmax,lloc,vp,vo,zion,rr,mmax,nkap,vpuns)
!! the semi-local potentials unscreened of the valence's screening vo,
!! the scalar-relativistic average when fully relativistic; and -zion/r
!! far out, where the all-electron charge outranges rho: switched over
!! the outer half of the mesh (O-5), not cut where rho underflows to 0
 use cinf_taper_m, only: cinf_taper
 integer, intent(in) :: irel,lmax,lloc,mmax,nkap
 real(dp), intent(in) :: vp(mmax,mxl+2,nkap),vo(mmax),zion,rr(mmax)
 real(dp), intent(inout) :: vpuns(mmax,mxl+2)
 real(dp) :: sw
 integer :: ii,ll,l1

 do l1=1,max(lmax+1,lloc+1)
   ll=l1-1
   if(irel/=2 .or. ll==0 .or. ll==4) then
     vpuns(:,l1)=vp(:,l1,1)-vo(:)
   else
     vpuns(:,l1)=(dble(ll+1)*vp(:,l1,1)+ dble(ll)*vp(:,l1,2))/dble(2*ll+1) &
&                -vo(:)
   end if
 end do

 do ii=1,mmax
   sw=cinf_taper(rr(ii),0.5d0*rr(mmax),rr(mmax))
   do l1=1,max(lmax+1,lloc+1)
     vpuns(ii,l1)=sw*vpuns(ii,l1)-(1.0d0-sw)*zion/rr(ii)
   end do
 end do
 end subroutine unscreen_potentials

 end module ps_valence_m

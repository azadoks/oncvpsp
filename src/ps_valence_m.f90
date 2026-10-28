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
 implicit none
 private
 public :: valence_densities, unscreen_potentials, coulomb_switch

 contains

 subroutine valence_densities(irel,c_light,nkap,nc,nv,na,la,fa,ea,it,lloc, &
&                             nproj,rr,zz,vfull,vp,vlocpw,vkbpw,evkb,mmax,uu, &
&                             up,uupsa,eeig,zval,rho,rhops,rhoae,rhotae,iout, &
&                             ifail)
!! the valence states, all-electron (in vfull) and pseudo (with the VKB
!! projectors, for the most consistent unscreening), state by state: their
!! densities rhoae and rhops (of the share of each J when fully
!! relativistic), the pseudo-waves uupsa for upfout, the all-electron
!! valence density rhotae, the pseudo one rho, the valence charge zval and
!! the sum of the pseudo eigenvalues eeig
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param na  principal quantum numbers
!! @param la  angular momenta
!! @param fa  occupancies
!! @param ea  the all-electron eigenvalues (Ha) of each state and kappa,
!!            the first guesses of both solvers
!! @param it  the iterations of the reference atom, for the error message
!! @param lloc  l of the local potential, 4 for the polynomial
!! @param nproj  number of projectors for each l
!! @param rr  log radial mesh (bohr)
!! @param zz  atomic number
!! @param vfull  all-electron potential (Ha)
!! @param vp  the screened semi-local pseudopotentials (Ha),
!!            vp(:,lloc+1,:) the local one
!! @param vlocpw  the screened local potential as pieces
!! @param vkbpw  the projectors as pieces
!! @param evkb  coefficients of the projectors (Ha)
!! @param mmax  size of the log radial mesh
!! @param uu  work space for the radial wave functions (times r); the last
!!            pseudo wave function on return
!! @param up  d(uu)/dn, its derivative in the mesh index (al*r*d(uu)/dr);
!!            work space as uu
!! @param uupsa  the pseudo wave function (times r) of each kappa and
!!               valence state
!! @param eeig  the sum of the pseudo eigenvalues (Ha), weighted by the
!!              occupancies
!! @param zval  the valence charge, the sum of the valence occupancies
!! @param rho  the pseudo valence charge density
!! @param rhops  the pseudo density (u/r)**2 of each valence state, of
!!               unit occupation
!! @param rhoae  the all-electron density of each valence state, as rhops
!! @param rhotae  all-electron valence charge density
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written
 use boundstate_m, only: boundstate
 use lschvkbb_m, only: lschvkbb
 use constants_m, only: dp, mxprj, mxstate, mxl, mxlin, mxkap
 use kappa_m, only: kap_count, kap_of, kap_j_share
 use pieces_m, only: piecewise_t

!Input variables
 integer, intent(in) :: irel,nkap,nc,nv,na(mxstate),la(mxstate),it,lloc
 real(dp), intent(in) :: c_light
 integer, intent(in) :: nproj(mxlin+1),mmax
 real(dp), intent(in) :: fa(mxstate),ea(mxstate,nkap),rr(mmax),zz
 real(dp), intent(in) :: vfull(mmax),vp(mmax,mxl+2,nkap)
 real(dp), intent(in) :: evkb(mxprj,mxl+1,nkap)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1,nkap)
 integer, intent(in) :: iout

!Input/output variables
 real(dp), intent(inout) :: uu(mmax,2),up(mmax,2)

!Output variables
 real(dp), intent(out) :: uupsa(mmax,nkap,nv),eeig,zval,rho(mmax)
 real(dp), intent(out) :: rhops(mmax,nv),rhoae(mmax,nv),rhotae(mmax)
 integer, intent(out) :: ifail

!Local variables
 integer :: nodes(mxl+1),kk,ll,l1,mkap,ikap,kap,ierr,mch
 real(dp) :: fj,et,emin,emax

 ifail=0
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
     if(nkap==2) fj=kap_j_share(ll,ikap,1.0d0)

     et=ea(nc+kk,ikap)

     call boundstate(irel,c_light,na(nc+kk),ll,kap,ierr,et, &
&                    rr,zz,vfull,uu,up,mmax,mch,iout,ifail)
     if(ifail/=0) return
     if(ierr /= 0) then
       if(irel==2) then
         write(iout,'(a,i2,a,i2,a,i2,a,i4)') &
&              'oncvpsp_r 541: ldiracfb convergence ERROR, n=',&
&               na(nc+kk),'  l=',ll, '  kap-',kap,'  ierr=',ierr
       else
         write(iout,'(/a,3i4)') &
&              'oncvpsp-461: lschfb convergence ERROR n,l,iter=', &
&              na(nc+kk),la(nc+kk),it
       end if
       ifail=1
       return
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
&                  evkb(1,l1,ikap),uu,up,mmax,mch,iout,ifail)
     if(ifail/=0) return
     if(ierr /= 0) then
       if(irel==2) then
         write(iout,'(a,i2,a,i2)') 'oncvpsp_r: lschvkbb convergence ERROR, n=',&
&         ll+1,'  l=',ll
       else
         write(iout,'(a,3i4)') 'oncvpsp: lschvkbb ERROR',ll+nodes(l1)+1,ll,ierr
         if(ierr==200) write(iout,'(a,i2,a,f12.6,a/a)') &
&             'a state of l=',ll,' lies below',emin, &
&             ' Ha that the node count does not expect:', &
&             'a bound ghost of the projectors'
         flush(6)
       end if
       ifail=1
       return
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
!!
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param lmax  maximum angular momentum
!! @param lloc  l of the local potential, 4 for the polynomial
!! @param vp  the screened semi-local pseudopotentials (Ha)
!! @param vo  the screening (Ha): the Hartree and exchange-correlation
!!            potential of the pseudo valence density and model core
!! @param zion  ion charge
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param vpuns  the unscreened semi-local pseudopotentials (Ha), set for
!!               l up to max(lmax,lloc)
 use constants_m, only: dp, mxl

!Input variables
 integer, intent(in) :: irel,lmax,lloc,mmax,nkap
 real(dp), intent(in) :: vp(mmax,mxl+2,nkap),vo(mmax),zion,rr(mmax)

!Input/output variables
 real(dp), intent(inout) :: vpuns(mmax,mxl+2)

!Local variables
 integer :: ll,l1

 do l1=1,max(lmax+1,lloc+1)
   ll=l1-1
   if(irel/=2 .or. ll==0 .or. ll==4) then
     vpuns(:,l1)=vp(:,l1,1)-vo(:)
   else
     vpuns(:,l1)=(dble(ll+1)*vp(:,l1,1)+ dble(ll)*vp(:,l1,2))/dble(2*ll+1) &
&                -vo(:)
   end if
 end do

 do l1=1,max(lmax+1,lloc+1)
   call coulomb_switch(rr,mmax,zion,vpuns(:,l1))
 end do
 end subroutine unscreen_potentials

 subroutine coulomb_switch(rr,mmax,zion,vv)
!! vv switched to -zion/r over the outer half of the mesh (O-5)
!!
!! @param rr  log radial mesh (bohr)
!! @param mmax  size of the log radial mesh
!! @param zion  ion charge
!! @param vv  a potential (Ha), switched in place
 use cinf_taper_m, only: cinf_taper
 use constants_m, only: dp

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: rr(mmax),zion

!Input/output variables
 real(dp), intent(inout) :: vv(mmax)

!Local variables
 real(dp) :: sw
 integer :: ii

 do ii=1,mmax
   sw=cinf_taper(rr(ii),0.5d0*rr(mmax),rr(mmax))
   vv(ii)=sw*vv(ii)-(1.0d0-sw)*zion/rr(ii)
 end do
 end subroutine coulomb_switch

 end module ps_valence_m

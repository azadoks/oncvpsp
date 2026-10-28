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

 module run_plot_m

 implicit none
 private
 public :: run_plot

 contains

 subroutine run_plot(lmax,npa,epa,wellprj,lloc,rc,vkb,vlocpw,vkbpw,evkb, &
&                    nproj,rr,vfull,vp,vpuns,zz,mmax,drl,nrl,rho,rhoc, &
&                    rhomod,irel,c_light,nkap,cvgplt,iout,ifail)

!! write output for plotting pseudopotentials, model core charge, and
!! all-electron and pseudo wave functions. The points are the mesh points
!! of one thinning, by the spacing dr=min(drl,rr(n2)/200) from rr(n1),
!! which no core radius moves, so each plotted value is a smooth function
!! of the radii; the plots of the potentials and projectors end at radii
!! of the continuous rc (doc/smoothness_plan.md, X-4)
!!
!! @param lmax  maximum angular momentum
!! @param npa  principal quantum number for corresponding all-electron state
!! @param epa  bound-state or scattering state reference energies for vkb
!!             potentials (last index 1 for j=l+1/2, 2 for j=l-1/2, as in vkb,
!!             evkb, vp, cvgplt)
!! @param wellprj  whether each projector is a well state (its energy then
!!                 of either sign), else a bound state
!! @param lloc  l of the local potential, 4 for the polynomial
!! @param rc  core radii (bohr)
!! @param vkb  VKB projectors on the mesh
!! @param vlocpw  the screened local potential as pieces
!! @param vkbpw  the projectors as pieces
!! @param evkb  coefficients of the projectors
!! @param nproj  number of projectors for each l
!! @param rr  log radial mesh (bohr)
!! @param vfull  all-electron potential (Ha)
!! @param vp  semi-local pseudopotentials
!! @param vpuns  unscreened vp
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param drl  spacing of linear radial mesh
!! @param nrl  size of the output mesh
!! @param rho  valence pseudocharge
!! @param rhoc  core charge density
!! @param rhomod  model core charge
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param cvgplt  energy per electron error vs. cutoff (convergence profiles)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp, mxl, mxlin, nrhomod, nerrlv, nerrdg, mxprj
 use kappa_m, only: kap_count, kap_of
 use pieces_m, only: piecewise_t
 use boundstate_m, only: boundstate
 use lschvkbb_m, only: lschvkbb
 use lschvkbs_m, only: lschvkbs
 use renorm_r_m, only: renorm_r
 use scatterstate_m, only: scatterstate
 use logmesh_m, only: lm_index

!Input variables
 integer, intent(in) :: lmax,lloc,mmax,nrl,irel,nkap
 real(dp), intent(in) :: c_light
 integer, intent(in) :: npa(mxprj,mxlin+1),nproj(mxlin+1)
 real(dp), intent(in) :: zz,drl,rc(mxlin+1)
 real(dp), intent(in) :: rr(mmax),vp(mmax,mxl+2,nkap),vpuns(mmax,mxl+2)
 real(dp), intent(in) :: vfull(mmax)
 real(dp), intent(in) :: vkb(mmax,mxprj,mxl+1,nkap)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1,nkap)
 real(dp), intent(in) :: rho(mmax),rhoc(mmax),rhomod(mmax,nrhomod)
 real(dp), intent(in) :: epa(mxprj,mxlin+1,nkap),evkb(mxprj,mxl+1,nkap)
 logical, intent(in) :: wellprj(mxprj,mxlin+1)
 real(dp), intent(in) :: cvgplt(2,nerrdg,mxprj,mxl+1,nkap)
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: ll,l1,n1,n2,n3,n4,nn
 integer :: ikap,mkap
 real(dp) :: dr,rcm
 real(dp), allocatable :: u2(:),up(:),ur(:,:),urp(:,:)
!the points that are written of a plot, each beyond the last by dr
 logical, allocatable :: keep(:)

 ifail=0
 allocate(u2(mmax),up(mmax),ur(mmax,2),urp(mmax,2),keep(mmax))

 write(iout,'(/a)') 'DATA FOR PLOTTING'

 n1=lm_index(rr,mmax,drl)
 n2=lm_index(rr,mmax,real(nrl,dp)*drl)
 rcm=maxval(rc(1:lmax+1))
 n4=min(n2,lm_index(rr,mmax,rcm+1.0d0)-1)
 n3=lm_index(rr,mmax,1.1d0*rcm)
 dr=min(drl,rr(n2)/200.0d0)

 call write_potentials()
 call write_charges()

! loop for wave function output, the all-electron functions (their large
! components when irel=2) normalized as Schroedinger functions; for
! irel=2, the j=l+1/2 functions and projectors are labelled with -iprj
! and -l

 do l1 = 1, lmax + 1
   ll = l1 - 1
   mkap=kap_count(nkap,ll)
! loop on J = ll +/- 1/2 when nkap=2
   do ikap=1,mkap
     call write_waves(l1,ikap,ifail)
     if(ifail/=0) return
     call write_projectors(l1,ikap)
   end do !ikap
 end do !l1

 call write_profiles()

 deallocate(u2,up,ur,urp,keep)
 return

 contains

 subroutine thin_points(i1,i2,dr,r0)
!! keep(i1:i2): the points of a plot from i1 to i2, each beyond r0, the
!! last kept point plus dr (r0 at first the start)
!!
!! @param i1  first point
!! @param i2  last point
!! @param dr  least spacing of the points kept
!! @param r0  start
 integer, intent(in) :: i1,i2
 real(dp), intent(in) :: dr
 real(dp), intent(inout) :: r0
 integer :: ii
 do ii=i1,i2
   keep(ii)=rr(ii)>r0
   if(keep(ii)) r0=rr(ii)+dr
 end do
 end subroutine thin_points

 subroutine write_potentials()
!! the !p lines of the charge and the potentials, and !L of the local
!! potential for lloc=4; the fully-relativistic code writes every point
!! of the potentials
 integer :: ii
 real(dp) :: r0
 write(iout,'(/2a/)') ' radii, charge,',' pseudopotentials (ll=0, 1, lmax)'
 r0=rr(n1)-dr
 call thin_points(n1,n4,dr,r0)
 do ii=n1,n4
   if(keep(ii) .or. irel==2) then
     write(iout,'(a,6(f12.7,1x))') '!p',rr(ii),rho(ii), &
&          (vpuns(ii,l1),l1=1,lmax+1)
   end if
 end do

 if(lloc == 4) then
   r0=rr(n1)-dr
   call thin_points(n1,n4,dr,r0)
   do ii=n1,n4
     if(keep(ii)) then
       write(iout,'(a, 2(f12.7,1x))') ' !L',rr(ii),vpuns(ii,lloc + 1)
     end if
   end do
 end if
 end subroutine write_potentials

 subroutine write_charges()
!! the !r lines of the charge, the core charge (cut at rmx) and the model
!! core charge, for the analysis of the model
 integer :: ii
 real(dp) :: r0,rmx
 rmx=0
 do ii=n1,n2
   rmx=max(rmx,5.0d0*rho(ii),rhomod(ii,1))
 end do

 write(iout,'(//2a/)') ' radii, charge,',' core charge, model core charge'
 r0=rr(n1)-dr
 call thin_points(n1,n2,dr,r0)
 do ii=n1,n2
   if(keep(ii)) then
     if(rhoc(ii)<rmx)then
       write(iout,'(a,8(f12.7,1x))') '!r',rr(ii),rho(ii),rhoc(ii), &
&       rhomod(ii,1)
     else if(rr(ii)>0.01d0) then
       write(iout,'(a,8(f12.7,1x))') '!r',rr(ii),rho(ii),rmx, &
&      rhomod(ii,1)
     end if
   end if
 end do
 end subroutine write_charges

 subroutine write_waves(l1,ikap,ifail)
!! the & lines of the all-electron and pseudo wave functions of each
!! projector of l1-1 and ikap: bound states as such, the others
!! scattering states at their energies, both of the sign of the
!! all-electron one at n2
 integer, intent(in) :: l1,ikap
 integer :: ii,iprj,npr,ierr,mch,kap,ilabel
 integer, intent(out) :: ifail
 real(dp) :: r0,cnorm,emax,emin,etest,sgnae,sgnps
 ifail=0
 kap=kap_of(ll,ikap)
 npr=nproj(l1)
 if(l1==lloc+1) npr=0

 do iprj=1,nproj(l1)
! a bound state is plotted as one; so is the j = l-1/2 partner of a well
! state below 0 (the eigenvalue of the j = l+1/2 well, D-7c), as before;
! a well state's own energy can be below 0 too (ep_as_given)
   if(epa(iprj,l1,ikap)<0.0d0 .and. &
&     .not. (wellprj(iprj,l1) .and. ikap==1)) then
     if(irel==2) then
       write(iout,'(//a,i2,a,i2,a,i2,a,a/)') 'n=',npa(iprj,l1),',  l=', &
&            ll,', kap=',kap, ', all-electron wave function', &
&            ', pseudo w-f'
     else
       write(iout,'(//a,i2,a,i2,a,a/)') 'n=',npa(iprj,l1),',  l=',ll, &
&            ', all-electron wave function', ', pseudo w-f'
     end if

     etest=epa(iprj,l1,ikap)
     call boundstate(irel,c_light,npa(iprj,l1),ll,kap,ierr,etest, &
&                    rr,zz,vfull,ur,urp,mmax,mch,iout,ifail)
     if(ifail/=0) return
     if(ierr /= 0 .and. irel==2) then
       write(iout,'(a,i2,a,i2,a,i2,a,i4)') &
&            'run_plot_r 409: ldiracfb convergence ERROR, n=',&
&             npa(iprj,l1),'  l=',ll, 'kap-',kap,'  ierr=',ierr
       ifail=1
       return
     else if(ierr /= 0) then
       write(iout,'(/a,3i4)') &
&            'run_plot: lschfb convergence ERROR n,l,ierr=', &
&            npa(iprj,l1),ll,ierr
       ifail=1
       return
     end if
     if(irel==2) then
       call renorm_r(ur,rr,ll,kap,c_light,zz,mmax,cnorm,iout=iout,ifail=ifail)
       if(ifail/=0) return
     end if

     emax=0.9d0*etest
     emin=1.1d0*etest

     call lschvkbb(ll+iprj,ll,npr,ierr,etest,emin,emax, &
&                  rr,vp(1,lloc+1,ikap),vlocpw,vkbpw(1:npr,l1,ikap), &
&                  evkb(1,l1,ikap),u2,up,mmax,mch,iout,ifail)
     if(ifail/=0) return
     if(ierr/=0) then
       write(iout,'(/a,i0,a,i0,a,i0)') &
&            'write_waves: WARNING lschvkbb did not converge, iprj=',iprj, &
&            ', ll=',ll,', ierr=',ierr
     end if
     sgnae=1.0d0
     sgnps=1.0d0
   else
     if(irel==2) then
       write(iout,'(//a,i2,a,i2,a,i2,a,a/)') 'scattering, iprj=',iprj, &
&           ',  l=',ll,', kap=',kap,', all-electron wave function', &
&           ', pseudo w-f'
     else
       write(iout,'(//a,i2,a,i2,a,a/)') 'scattering, iprj=',iprj, &
&           ',  l=',ll,', all-electron wave function', ', pseudo w-f'
     end if
     call scatterstate(irel,c_light,nn,ll,kap,epa(iprj,l1,ikap), &
&                   rr,zz,vfull,ur,urp,mmax,n2)
     if(irel==2) then
       call renorm_r(ur,rr,ll,kap,c_light,zz,mmax,cnorm,iout=iout,ifail=ifail)
       if(ifail/=0) return
     end if

     call lschvkbs(ll,npr,epa(iprj,l1,ikap),rr,vp(1,lloc+1,ikap), &
&                  vlocpw,vkbpw(1:npr,l1,ikap),evkb(1,l1,ikap),u2,up, &
&                  mmax,n2,iout,ifail)
     if(ifail/=0) return
     sgnae=sign(1.0d0,ur(n2,1))
     sgnps=sign(1.0d0,u2(n2))
   end if

!the fully-relativistic code starts a step further out
   if(irel==2) then
     r0=rr(n1)+dr
   else
     r0=rr(n1)-dr
   end if
   ilabel=iprj
   if(nkap==2 .and. ikap==1) ilabel=-iprj
   call thin_points(n1,n2,dr,r0)
   do ii = n1, n2
     if(keep(ii)) then
       write(iout,'(a,i5,i1,3(f12.6,1x))') '&',ilabel,ll,rr(ii), &
&            sgnae*ur(ii,1),sgnps*u2(ii)
     end if
   end do
 end do !iproj
 end subroutine write_waves

 subroutine write_projectors(l1,ikap)
!! the !J lines of the orthonormal projectors of l1-1 and ikap
 integer, intent(in) :: l1,ikap
 integer :: ii,jj,lout
 real(dp) :: r0
 write(iout,'(/)')
 r0=rr(n1)-dr
 lout=ll
 if(nkap==2 .and. ikap==1) lout=-ll
 call thin_points(n1,n3,dr,r0)
 do ii=n1,n3
   if(keep(ii)) then
     write(iout,'(a,i6,6(f12.6,1x))') '!J',lout,rr(ii), &
&          (vkb(ii,jj,l1,ikap),jj=1,nproj(l1))
   end if
 end do
 end subroutine write_projectors

 subroutine write_profiles()
!! the !C lines of the convergence profiles, the average over kappa for
!! l>0 when fully relativistic
 integer :: jj
 if(irel==2) then
   write(iout,'(/a)') 'convergence profiles, (ll=0,lmax), kappa average'
   l1=1
   ll=l1-1
   do jj=1,nerrlv
     if(cvgplt(1,jj,1,l1,1)/=0.0d0) then
       write(iout,'(a,i6,3f12.6)') '!C',ll,cvgplt(1,jj,1,l1,1), &
&            cvgplt(2,jj,1,l1,1)
     end if
   end do
   do l1=2,lmax+1
     ll=l1-1
     do jj=1,nerrlv
       if(cvgplt(1,jj,1,l1,1)/=0.0d0 .and. cvgplt(1,jj,1,l1,nkap)/=0.0d0) then
         write(iout,'(a,i6,3f12.6)') '!C',ll,cvgplt(1,jj,1,l1,1), &
&         0.5d0*(cvgplt(2,jj,1,l1,1)+cvgplt(2,jj,1,l1,nkap))
       end if
     end do
   end do
 else
   write(iout,'(/a)') 'convergence profiles, (ll=0,lmax)'
   write(iout,*) 'lmax',lmax
   do l1=1,lmax+1
     ll=l1-1
     do jj=1,nerrlv
       if(cvgplt(1,jj,1,l1,1)/=0.0d0) then
         write(iout,'(a,i6,3(f12.6,1x))') '!C',ll,cvgplt(1,jj,1,l1,1), &
&              cvgplt(2,jj,1,l1,1)
       end if
     end do
   end do
 end if
 end subroutine write_profiles

 end subroutine run_plot

 end module run_plot_m

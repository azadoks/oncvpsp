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
 subroutine run_diag(lmax,npa,epa,lloc,irc, &
&                    vkb,evkb,nproj,rr,vfull,vp,zz,mmax,mxprj,irel,nkap)

!! diagnostics for semi-local and Vanderbilt-Kleinman-bylander pseudopotentials
!! checks bound-state energies, norms, and slopes, and pseudo-bound state
!! quantities for positive-energy projectors by matching all-electron
!! log derivatives at rc.
!! Should test as essentially exact for non-relativitic calculations
!! Error for relativistic is B matrix Hermiticity error
!!
!! @param lmax  maximum angular momentum
!! @param npa  principal quantum number for corresponding all-electron state
!! @param epa  bound-state or scattering state reference energies for vkb
!!             potentials (last index 1 for j=l+1/2, 2 for j=l-1/2, as in vkb,
!!             evkb and vp)
!! @param lloc  l for local potential
!! @param irc  core radii indices
!! @param vkb  VKB projectors
!! @param evkb  coefficients of VKB projectors
!! @param nproj  number of vkb projectors for each l
!! @param rr  log radial grid
!! @param vfull  all-electron potential
!! @param vp  semi-local pseudopotentials (vp(:,5) is local potential if linear
!!            comb.)
!! @param zz  atomic number
!! @param mmax  size of radial grid
!! @param mxprj  dimension of number of projectors
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param nkap  number of kappa per l: 1, or 2 for irel=2

 use constants_m, only: dp, mxl, mxlin
 implicit none

!Input variables
 integer, intent(in) :: lmax,lloc,mmax,mxprj,irel,nkap
 integer, intent(in) :: npa(mxprj,mxlin+1),irc(mxlin+1),nproj(mxlin+1)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: rr(mmax),vp(mmax,mxl+2,nkap),vfull(mmax)
 real(dp), intent(in) :: vkb(mmax,mxprj,mxl+1,nkap)
 real(dp), intent(in) :: epa(mxprj,mxlin+1,nkap),evkb(mxprj,mxl+1,nkap)

!Local variables
 integer :: ll,l1,ikap,kap,mkap,ierr,mch,mchf
 integer :: iprj,nnae,nnp,npr
 real(dp) :: al,emax,emin,etest,umch,upmch,uldf,gam,gpr,cnorm
 real(dp), allocatable :: uu(:),up(:),ur(:,:),urp(:,:)

 allocate(uu(mmax),up(mmax),ur(mmax,2),urp(mmax,2))

 al = 0.01d0 * log(rr(101)/rr(1))

! loop for diagnostic output using Vanderbilt Kleinman-Bylander projectors
!
 if(irel==2) then
   write(6,'(/a/a)') &
&   'Diagnostic tests using Vanderbilt-Kleinman-Bylander pseudopotentials',&
&   '  relativistic with J = L +/- 1/2 projectors'
   write(6,'(/2a)')'   l  kap   rcore       rmatch      e in        ', &
&     'delta e    norm test   slope test'
 else
   write(6,'(/a/a)') &
&   'Diagnostic tests using Vanderbilt-Kleinman-Bylander pseudopotentials',&
&   '  1 or more projectors used as specified by nproj input data'
   write(6,'(/2a)')'   l    rcore       rmatch      e in        ', &
&     'delta e   norm test  slope test'
 end if
!
 do l1 = 1, lmax + 1
   write(6,'(a)') ''
   ll = l1 - 1
   nnp=l1
   do iprj=1,max(1,nproj(l1))
     if(nkap==1 .or. l1==1) then
       mkap=1
     else
       mkap=2
     end if
!  loop on J = ll +/- 1/2 when nkap=2
     do ikap=1,mkap
       if(ikap==1) then
         kap=-(ll+1)
       else
         kap=  ll
       end if

       npr=nproj(l1)

!   find cutoff radius for projectors
       mchf=max(irc(l1),irc(lloc+1))+5

!   all-electron state, its large component when irel=2, normalized as a
!   Schroedinger function
       etest=epa(iprj,l1,ikap)
       if(epa(iprj,l1,1)<0.0d0) then
         call boundstate(irel,npa(iprj,l1),ll,kap,ierr,etest, &
&                        rr,zz,vfull,ur,urp,mmax,mch)
         if(ierr /= 0 .and. irel==2) then
           write(6,'(/a,4i4)') 'run_diag_r: ldiracfb convergence ERROR &
&         n,l,,kap,ierr=', &
&          npa(iprj,l1),ll,kap,ierr
           error stop 1
         else if(ierr /= 0) then
           write(6,'(/a,3i4)') &
&                'run_diag: lschfb convergence ERROR n,l,ierr=', &
&                npa(iprj,l1),ll,ierr
           error stop 1
         end if
       else
         call scatstate(irel,nnae,ll,kap,ierr,etest, &
&                       rr,zz,vfull,ur,urp,mmax,mchf)
       end if
       cnorm=1.0d0
       if(irel==2) call renorm_r(ur,rr,ll,kap,zz,mmax,cnorm)
       umch=ur(mchf,1)
       upmch=cnorm*urp(mchf,1)
       uldf=upmch/umch

       if(l1==lloc+1) npr=0

       etest=epa(iprj,l1,ikap)
       emax=etest+0.1d0*abs(etest)+0.05d0
       emin=etest-0.1d0*abs(etest)-0.05d0

       if(epa(iprj,l1,1)<0.0d0) then

         emax=min(emax,0.0d0)

         call lschvkbb(ll+iprj,ll,npr,ierr,etest,emin,emax, &
&                      rr,vp(1,lloc+1,ikap),vkb(1,1,l1,ikap), &
&                      evkb(1,l1,ikap),uu,up,mmax,mch)
         if(ierr/=0 .and. irel==2) then
           write(6,'(a,4i4,1p,2e16.8)') 'run_diag_r: lschvkbb ERROR', &
&                iprj,ll,kap,ierr,epa(iprj,l1,ikap),etest
         else if(ierr/=0) then
           write(6,'(a,3i4,1p,2e16.8)') 'run_diag: lschvkbb ERROR', &
&                iprj,ll,ierr,epa(iprj,l1,ikap),etest
         end if

         nnp=nnp+1
       else
!     calculate effective pseudo wave function principal quantum number
!     from all-electron node count

         nnp=nnae-npa(1,l1)+ll+1
         call lschvkbbe(nnp,ll,npr,ierr,etest,uldf,emin,emax, &
&                       rr,vp(1,lloc+1,ikap),vkb(1,1,l1,ikap), &
&                       evkb(1,l1,ikap),uu,up,mmax,mchf)

         if(ierr/=0 .and. irel==2) then
           write(6,'(a,4i4,1p,2e16.8)') 'run_diag_r: lschvkbbe ERROR', &
&                iprj,ll,kap,ierr,epa(iprj,l1,ikap),etest
         else if(ierr/=0) then
           write(6,'(a,4i4,1p,2e16.8)') 'run_diag: lschvkbbe ERROR', &
&                iprj,nnp,ll,ierr,epa(iprj,l1,ikap),etest
         end if
       end if !epa<0 (bound or not)

       gam=abs(umch/uu(mchf))
       gpr=abs(upmch/up(mchf))
!
       if(irel==2) then
         write(6,'(2i4,6f12.7)') ll,kap,rr(irc(l1)),rr(mchf), &
&              epa(iprj,l1,ikap),etest-epa(iprj,l1,ikap),gam,gpr
       else
         write(6,'(i4,6f12.7)') ll,rr(irc(l1)),rr(mchf),epa(iprj,l1,ikap), &
&              etest-epa(iprj,l1,ikap),gam,gpr
       end if

     end do  !ikap
   end do  !iprj
 end do !l1

 deallocate(uu,up,ur,urp)
 return
 end subroutine run_diag

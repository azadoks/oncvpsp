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
 subroutine run_diag_sr_so_r(lmax,npa,epa,wellprj,lloc,irc, &
&                    vlocpw,vsrpw,esr,vsopw,eso,nproj,rr,vfull,vp,zz,mmax, &
&                    mxprj)

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
!!             potentials
!! @param wellprj  whether each projector is a well state (its energy then
!!                 of either sign), else a bound state
!! @param lloc  l for local potential
!! @param irc  core radii indices
!! @param vlocpw  the screened local potential's pieces
!! @param vsrpw  normalized scalar-relativistic projectors, as functions of
!!               pieces
!! @param esr  energy coefficients of vsr
!! @param vsopw  normalized spin-orbit projectors, as functions of pieces
!! @param eso  energy coefficients of vso (multiplied by <L dot S> here)
!! @param nproj  number of vkb projectors for each l
!! @param rr  log radial grid
!! @param vfull  all-electron potential
!! @param vp  semi-local pseudopotentials (vp(:,5) is local potential if linear
!!            comb.)
!! @param zz  atomic number
!! @param mmax  size of radial grid
!! @param mxprj  dimension of number of projectors

 use constants_m, only: dp, mxl, mxlin, mxkap
 use pieces_m, only: piecewise_t
 use outwf_m, only: outwf_index_at
 implicit none

!Input variables
 integer, intent(in) :: lmax,lloc,mmax,mxprj
 integer, intent(in) :: npa(mxprj,mxlin+1),irc(mxlin+1),nproj(mxlin+1)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: rr(mmax),vp(mmax,mxl+2,mxkap),vfull(mmax)
 real(dp), intent(in) :: esr(2*mxprj,mxl+1)
 real(dp), intent(in) :: eso(2*mxprj,mxl+1)
 real(dp), intent(in) :: epa(mxprj,mxlin+1,mxkap)
 logical, intent(in) :: wellprj(mxprj,mxlin+1)
 type(piecewise_t), intent(in) :: vlocpw,vsrpw(2*mxprj,mxl+1)
 type(piecewise_t), intent(in) :: vsopw(2*mxprj,mxl+1)

!Local variables
 integer :: ll,l1,ikap,kap,mkap,ii,jj,ierr,mch,mchf,nvkbt
 integer :: iprj,nnae,nnp,npr
 real(dp) :: al,emax,emin,etest,sls,umch,upmch,uldf,gam,gpr,cnorm,fso
 real(dp), allocatable :: uu(:),up(:),ur(:,:),urp(:,:)
 real(dp), allocatable :: evkbt(:)
 type(piecewise_t), allocatable :: vkbtpw(:)

 allocate(uu(mmax),up(mmax),ur(mmax,2),urp(mmax,2))
 allocate(evkbt(4*mxprj),vkbtpw(4*mxprj))


 al = 0.01d0 * log(rr(101)/rr(1))

! loop for diagnostic output using Vanderbilt Kleinman-Bylander projectors
!
 write(6,'(/a/a)') &
& 'Diagnostic tests using Vanderbilt-Kleinman-Bylander pseudopotentials',&
& '  relativistic with scalar and spin-orbit non-local projectors'
 write(6,'(/2a)')'   l  kap   rcore       rmatch      e in        ', &
&   'delta e    norm test   slope test'
!
 do l1 = 1, lmax + 1
   write(6,'(a)') ''
   ll = l1 - 1
   nnp=l1
   do iprj=1,max(1,nproj(l1))
     if(l1==1) then
       mkap=1
     else
       mkap=2
     end if
!  loop on J = ll +/- 1/2
     do ikap=1,mkap
       if(ikap==1) then
         kap=-(ll+1)
       else
         kap=  ll
       end if

       npr=nproj(l1)

! construct multi-projector potential
       evkbt(:)=0.0d0
       if(ikap==1) fso= 0.5d0*ll
       if(ikap==2) fso=-0.5d0*(ll+1)
       if(ll==0) then
         nvkbt=npr
         do ii=1,npr
           vkbtpw(ii)=vsrpw(ii,l1)
           evkbt(ii)=esr(ii,l1)
         end do
       else
         nvkbt=4*npr
         do ii=1,2*npr
           vkbtpw(ii)=vsrpw(ii,l1)
           evkbt(ii)=esr(ii,l1)
           vkbtpw(2*npr+ii)=vsopw(ii,l1)
           evkbt(2*npr+ii)=fso*eso(ii,l1)
         end do
       end if


!   find cutoff radius for projectors
       mchf=max(irc(l1),irc(lloc+1))+5

       etest=epa(iprj,l1,ikap)
       if(.not. wellprj(iprj,l1)) then
         call ldiracfb(npa(iprj,l1),ll,kap,ierr,etest, &
&                        rr,zz,vfull,ur,urp,mmax,mch)
       else
         call ldiracfs(nnae,ll,kap,ierr,etest, &
&                     rr,zz,vfull,ur,urp,mmax,mchf)
       end if
       call renorm_r(ur,rr,ll,kap,zz,mmax,cnorm)
       umch=ur(mchf,1)
       upmch=cnorm*urp(mchf,1)
       uldf=upmch/umch


       etest=epa(iprj,l1,ikap)
       emax=etest+0.1d0*abs(etest)+0.05d0
       emin=etest-0.1d0*abs(etest)-0.05d0

       if(.not. wellprj(iprj,l1)) then
         sls=(l1-1)*l1
         emax=min(emax,0.0d0)

         call lschvkbb(ll+iprj,ll,nvkbt,ierr,etest,emin,emax, &
&                      rr,vp(1,lloc+1,ikap),vlocpw,vkbtpw(1:nvkbt),evkbt, &
&                      uu,up,mmax,mch)
         if(ierr/=0) then
           write(6,'(a,4i4,1p,2e16.8)') 'run_diag_sr_so_r: lschvkbb ERROR', &
&                iprj,ll,kap,ierr,epa(iprj,l1,ikap),etest
         end if

         nnp=nnp+1
       else
!     the root with as many eigenvalues below it as the pseudo-atom has
!     at epa, where the log derivatives match by construction: the root
!     nearest epa. The all-electron node count is no target, as the
!     pseudo-wave function may have a node that entered at r = 0 without
!     an eigenvalue passing (outwf_index)
         nnp=ll+1+outwf_index_at(ll,rr,mmax,vp(1,lloc+1,ikap),vlocpw,nvkbt, &
&            vkbtpw(1:nvkbt),evkbt(1:nvkbt),epa(iprj,l1,ikap),mchf)
         call lschvkbbe(nnp,ll,nvkbt,ierr,etest,uldf,emin,emax, &
&                       rr,vp(1,lloc+1,ikap),vlocpw,vkbtpw(1:nvkbt),evkbt, &
&                       uu,up,mmax,mchf)

         if(ierr/=0) then
           write(6,'(a,4i4,1p,2e16.8)') 'run_diag_r: lschvkbbe ERROR', &
&                iprj,ll,kap,ierr,epa(iprj,l1,ikap),etest
         end if
       end if !epa<0 (bound or not)

       gam=abs(umch/uu(mchf))
       gpr=abs(upmch/up(mchf))
!
       write(6,'(2i4,6f12.7)') ll,kap,rr(irc(l1)),rr(mchf),epa(iprj,l1,ikap), &
&            etest-epa(iprj,l1,ikap),gam,gpr

     end do  !ikap
   end do  !iprj
 end do !l1

 deallocate(uu,up)
 return
 end subroutine run_diag_sr_so_r

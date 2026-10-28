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

 module run_diag_m

 implicit none
 private
 public :: run_diag

 contains

 subroutine run_diag(lmax,npa,epa,wellprj,lloc,rc,vlocpw,vkbpw,evkb,nproj,rr, &
&                    vfull,vp,zz,mmax,irel,c_light,nkap,iso,vsrpw,esr, &
&                    vsopw,eso,iout,ifail)

!! diagnostics for semi-local and Vanderbilt-Kleinman-bylander pseudopotentials
!! checks bound-state energies, norms, and slopes, and pseudo-bound state
!! quantities for positive-energy projectors by matching all-electron
!! log derivatives at the radius rm = max(rc(l), rc(lloc)) diag_match, the
!! functions interpolated there and a scattering state normalized to it,
!! so every number is a smooth function of the radii.
!! Should test as essentially exact for non-relativitic calculations
!! Error for relativistic is B matrix Hermiticity error
!! With iso=1 (fully relativistic) the projectors of each J are instead
!! the scalar-relativistic ones and the spin-orbit ones times <L dot S>
!! (sr_so_r)
!!
!! @param lmax  maximum angular momentum
!! @param npa  principal quantum number for corresponding all-electron state
!! @param epa  bound-state or scattering state reference energies for vkb
!!             potentials (last index 1 for j=l+1/2, 2 for j=l-1/2, as in vkb,
!!             evkb and vp)
!! @param wellprj  whether each projector is a well state (its energy then
!!                 of either sign), else a bound state
!! @param lloc  l for local potential
!! @param rc  core radii (bohr)
!! @param vlocpw  the screened local potential's pieces
!! @param vkbpw  VKB projectors as functions of pieces
!! @param evkb  coefficients of VKB projectors
!! @param nproj  number of vkb projectors for each l
!! @param rr  log radial grid
!! @param vfull  all-electron potential
!! @param vp  semi-local pseudopotentials (vp(:,5) is local potential if linear
!!            comb.)
!! @param zz  atomic number
!! @param mmax  size of radial grid
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param nkap  number of kappa per l: 1, or 2 for irel=2
!! @param iso  0 the projectors vkbpw, evkb; 1 those of vsrpw..eso
!! @param vsrpw  normalized scalar-relativistic projectors, as functions of
!!               pieces (iso=1)
!! @param esr  energy coefficients of vsr (iso=1)
!! @param vsopw  normalized spin-orbit projectors, as functions of pieces
!!               (iso=1)
!! @param eso  energy coefficients of vso (multiplied by <L dot S> here;
!!             iso=1)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp, mxl, mxlin, diag_match, mxprj
 use logmesh_m, only: lm_interp
 use phsft_m, only: phsft_mch
 use quad_m, only: log_norm_to
 use sr_flux_m, only: sr_gamma
 use kappa_m, only: kap_count, kap_of
 use pieces_m, only: piecewise_t
 use outwf_m, only: outwf_index_at
 use boundstate_m, only: boundstate
 use lschvkbb_m, only: lschvkbb
 use lschvkbbe_m, only: lschvkbbe
 use renorm_r_m, only: renorm_r
 use scatterstate_m, only: scatterstate

!Input variables
 integer, intent(in) :: lmax,lloc,mmax,irel,nkap,iso
 real(dp), intent(in) :: c_light
 integer, intent(in) :: npa(mxprj,mxlin+1),nproj(mxlin+1)
 real(dp), intent(in) :: zz,rc(mxlin+1)
 real(dp), intent(in) :: rr(mmax),vp(mmax,mxl+2,nkap),vfull(mmax)
 real(dp), intent(in) :: epa(mxprj,mxlin+1,nkap),evkb(mxprj,mxl+1,nkap)
 logical, intent(in) :: wellprj(mxprj,mxlin+1)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1,nkap)
 real(dp), intent(in) :: esr(2*mxprj,mxl+1),eso(2*mxprj,mxl+1)
 type(piecewise_t), intent(in) :: vsrpw(2*mxprj,mxl+1),vsopw(2*mxprj,mxl+1)
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: ll,l1,ikap,kap,mkap,ierr,mch,mchf,nvkbt
 integer :: iprj,nnae,nnp,npr
 real(dp) :: emax,emin,etest,umch,upmch,uldf,gam,gpr,cnorm
 real(dp) :: rm,fss,uq,upq
 real(dp), allocatable :: uu(:),up(:),ur(:,:),urp(:,:)

 real(dp), allocatable :: evkbt(:)
 type(piecewise_t), allocatable :: vkbtpw(:)
 character(len=16) :: name

 ifail=0
 allocate(uu(mmax),up(mmax),ur(mmax,2),urp(mmax,2))
 allocate(evkbt(4*mxprj),vkbtpw(4*mxprj))

 fss=0.0d0
 if(irel/=0) fss=(1.0d0/c_light)**2

! loop for diagnostic output using Vanderbilt Kleinman-Bylander projectors
!
 if(iso==1) then
   name='run_diag_sr_so_r'
   write(iout,'(/a/a)') &
& 'Diagnostic tests using Vanderbilt-Kleinman-Bylander pseudopotentials',&
& '  relativistic with scalar and spin-orbit non-local projectors'
   write(iout,'(/2a)')'   l  kap   rcore       rmatch      e in        ', &
&   'delta e    norm test   slope test'
 else if(irel==2) then
   name='run_diag_r'
   write(iout,'(/a/a)') &
&   'Diagnostic tests using Vanderbilt-Kleinman-Bylander pseudopotentials',&
&   '  relativistic with J = L +/- 1/2 projectors'
   write(iout,'(/2a)')'   l  kap   rcore       rmatch      e in        ', &
&     'delta e    norm test   slope test'
 else
   name='run_diag'
   write(iout,'(/a/a)') &
&   'Diagnostic tests using Vanderbilt-Kleinman-Bylander pseudopotentials',&
&   '  1 or more projectors used as specified by nproj input data'
   write(iout,'(/2a)')'   l    rcore       rmatch      e in        ', &
&     'delta e   norm test  slope test'
 end if
!
 do l1 = 1, lmax + 1
   write(iout,'(a)') ''
   ll = l1 - 1
   nnp=l1
   do iprj=1,max(1,nproj(l1))
     mkap=kap_count(nkap,l1-1)
!  loop on J = ll +/- 1/2 when nkap=2
     do ikap=1,mkap
       kap=kap_of(ll,ikap)

       npr=nproj(l1)

!   the matching radius, and the last mesh point its interpolants read
       rm=max(rc(l1),rc(lloc+1))*diag_match
       mchf=phsft_mch(rr,mmax,rm)
       if(mchf==0) then
         write(iout,'(/a,f12.6)') trim(name)// &
&              ': ERROR matching radius beyond the mesh, rm=',rm
         ifail=1
         return
       end if

!   all-electron state, its large component when irel=2, normalized as a
!   Schroedinger function
       etest=epa(iprj,l1,ikap)
       if(.not. wellprj(iprj,l1)) then
         call boundstate(irel,c_light,npa(iprj,l1),ll,kap,ierr,etest, &
&                        rr,zz,vfull,ur,urp,mmax,mch,iout,ifail)
         if(ifail/=0) return
         if(ierr /= 0 .and. irel==2) then
           write(iout,'(/a,4i4)') 'run_diag_r: ldiracfb convergence ERROR &
&         n,l,,kap,ierr=', &
&          npa(iprj,l1),ll,kap,ierr
           ifail=1
           return
         else if(ierr /= 0) then
           write(iout,'(/a,3i4)') &
&                'run_diag: lschfb convergence ERROR n,l,ierr=', &
&                npa(iprj,l1),ll,ierr
           ifail=1
           return
         end if
       else
         call scatterstate(irel,c_light,nnae,ll,kap,etest, &
&                       rr,zz,vfull,ur,urp,mmax,mchf)
       end if
       cnorm=1.0d0
       if(irel==2 .and. wellprj(iprj,l1)) then
         call renorm_r(ur,rr,ll,kap,c_light,zz,mmax,cnorm,iout,ifail,rnorm=rm)
         if(ifail/=0) return
       else if(irel==2) then
         call renorm_r(ur,rr,ll,kap,c_light,zz,mmax,cnorm,iout,ifail)
         if(ifail/=0) return
       else if(wellprj(iprj,l1)) then
!   a scattering state normalized to rm, as the pseudo-wave function
         cnorm=1.0d0/sqrt(log_norm_to(rr,ur(1:mchf,1)**2,mchf,ll, &
&              sr_gamma(ll,zz,fss),rm))
         ur(1:mchf,1)=cnorm*ur(1:mchf,1)
       end if
       call lm_interp(rr(1:mchf),ur(1:mchf,1),mchf,rm,umch)
       call lm_interp(rr(1:mchf),urp(1:mchf,1),mchf,rm,upmch)
       upmch=cnorm*upmch
       uldf=upmch/umch

       if(l1==lloc+1) npr=0
       call assemble_projectors(ll,l1,ikap,npr,nvkbt)

       etest=epa(iprj,l1,ikap)
       emax=etest+0.1d0*abs(etest)+0.05d0
       emin=etest-0.1d0*abs(etest)-0.05d0

       if(.not. wellprj(iprj,l1)) then

         emax=min(emax,0.0d0)

         call lschvkbb(ll+iprj,ll,nvkbt,ierr,etest,emin,emax, &
&                      rr,vp(1,lloc+1,ikap),vlocpw,vkbtpw(1:nvkbt),evkbt, &
&                      uu,up,mmax,mch,iout,ifail)
         if(ifail/=0) return
         if(ierr/=0 .and. irel==2) then
           write(iout,'(a,4i4,1p,2e16.8)') trim(name)//': lschvkbb ERROR', &
&                iprj,ll,kap,ierr,epa(iprj,l1,ikap),etest
         else if(ierr/=0) then
           write(iout,'(a,3i4,1p,2e16.8)') 'run_diag: lschvkbb ERROR', &
&                iprj,ll,ierr,epa(iprj,l1,ikap),etest
         end if

         nnp=nnp+1
       else
!     the root with as many eigenvalues below it as the pseudo-atom has
!     at epa, where the log derivatives match by construction: the root
!     nearest epa. The all-electron node count is no target, as the
!     pseudo-wave function may have a node that entered at r = 0 without
!     an eigenvalue passing (outwf_index)

         nnp=ll+1+outwf_index_at(ll,rr,mmax,vp(1,lloc+1,ikap),vlocpw,nvkbt, &
&            vkbtpw(1:nvkbt),evkbt(1:nvkbt),epa(iprj,l1,ikap),mchf,iout,ifail)
         if(ifail/=0) return
         call lschvkbbe(nnp,ll,nvkbt,ierr,etest,uldf,emin,emax, &
&                       rr,vp(1,lloc+1,ikap),vlocpw,vkbtpw(1:nvkbt),evkbt, &
&                       uu,up,mmax,mchf,rm,iout,ifail)
         if(ifail/=0) return

         if(ierr/=0 .and. irel==2) then
           write(iout,'(a,4i4,1p,2e16.8)') trim(name)//': lschvkbbe ERROR', &
&                iprj,ll,kap,ierr,epa(iprj,l1,ikap),etest
         else if(ierr/=0) then
           write(iout,'(a,4i4,1p,2e16.8)') 'run_diag: lschvkbbe ERROR', &
&                iprj,nnp,ll,ierr,epa(iprj,l1,ikap),etest
         end if
       end if !epa<0 (bound or not)

       call lm_interp(rr(1:mchf),uu(1:mchf),mchf,rm,uq)
       call lm_interp(rr(1:mchf),up(1:mchf),mchf,rm,upq)
       gam=abs(umch/uq)
       gpr=abs(upmch/upq)
!
       if(irel==2) then
         write(iout,'(2i4,6f12.7)') ll,kap,rc(l1),rm, &
&              epa(iprj,l1,ikap),etest-epa(iprj,l1,ikap),gam,gpr
       else
         write(iout,'(i4,6f12.7)') ll,rc(l1),rm,epa(iprj,l1,ikap), &
&              etest-epa(iprj,l1,ikap),gam,gpr
       end if

     end do  !ikap
   end do  !iprj
 end do !l1

 deallocate(uu,up,ur,urp,evkbt,vkbtpw)
 return

 contains

 subroutine assemble_projectors(ll,l1,ikap,npr,nvkbt)
!! the nvkbt projectors vkbtpw, evkbt of l and ikap: the npr of vkbpw,
!! evkb, or with iso=1 the scalar-relativistic ones and the spin-orbit
!! ones times <L dot S>
!!
!! @param ll  angular momentum
!! @param l1  ll+1
!! @param ikap  1 for j=l+1/2, 2 for j=l-1/2
!! @param npr  number of projectors of l
!! @param nvkbt  number of projectors of the channel
 integer, intent(in) :: ll,l1,ikap,npr
 integer, intent(out) :: nvkbt
 real(dp) :: fso
 integer :: ii
 evkbt(:)=0.0d0
 if(iso==0) then
   nvkbt=npr
   do ii=1,npr
     vkbtpw(ii)=vkbpw(ii,l1,ikap)
     evkbt(ii)=evkb(ii,l1,ikap)
   end do
   return
 end if
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
 end subroutine assemble_projectors

 end subroutine run_diag

 end module run_diag_m

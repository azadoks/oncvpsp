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
! self-consistent all-electron and pseudopotential atomic calculations
! compared for reference and tests atomic configurations


 subroutine run_config_ps(jj,cae,nc,rhov,rhomod,rr, &
&                 rcmax,mmax,mxprj,iexc,etot,epstot,nproj,vpuns, &
&                 lloc,vkb,evkb)

!jj  index of current configufation
!cae  all-electron results for this configuration, from run_config_ae
!nc  number of core states
!rhov valence pseudo-charge of reference configuration
!rhomod  model core charge
!rr  log radial mesh
!rcmax  maximum core radius for psp
!mmax  size of log grid
!mxprj  dimension of number of projectors
!iexc  exchange-correlation function to be used
!etot  reference configuration total energy
!epstot  pseudoatom total energy
!nproj  number of VKB projectors to use for each l
!vpuns  unscreened semi-local pseudopotentials (plus different vloc if lloc==4)
!lloc  index-1 of local potential
!vkb   Vanderbilt-Kleinman-Bylander projectors
!evkb VKB projector coefficients

 use config_ae_m, only: config_ae_t, config_ae_warning
 use constants_m, only: dp
 implicit none

!Input variables

 type(config_ae_t), intent(in) :: cae
 integer, intent(in) :: jj,mmax,mxprj,iexc,nc,lloc
 integer, intent(in) :: nproj(5)
 real(dp), intent(in) :: etot,epstot,rcmax
 real(dp), intent(in) :: rhov(mmax),rr(mmax)
 real(dp), intent(in) :: vpuns(mmax,5),vkb(mmax,mxprj,4),evkb(mxprj,4),rhomod(mmax,5)

!Output variables - printing only

!Local variables
 integer :: ii,it,kk,l1,ierr,nvt
 integer :: nat(30),lat(30),natp(30),latp(30),nav(4)
 real(dp) :: eaetst,etsttot
 real(dp) :: eat(30,3),fat(30,3),eatp(30),fatp(30,3)

 real(dp),allocatable :: rho(:),rhocps(:)

! atom tests comparing reference state and excited configuration for
! all-electron and pseudo atoms.  Total excitation energies and excited-
! -configuration eigenvalues are compared

! For multi-projector potentials, the bound-state solver lschvkbb has poor
! stability, and requires a trial energy sufficiently close to the correct
! result for a given local potential. This appears to be due to the fact
! that intermediate stages towards the solution of the radial Schroedinger
! equation with a non-local potential may not increase their node count
! monotonically with increasing energy (the same phenonenon which gives
! rise to ghost states).

! To aviod these problems, the pseudoatom calculation is started with the
! screened local potential and valence eigenvalues of the reference
! configuration.  This is adiabatically changed in two steps.  In the first,
! the occupation numbers are changed by 2% within the potential iteration
! loop for each of the first 50 steps to those of a maximally-ionized
! intermediate configuration.

! The second step is started with the local potential of this maximally-
! ionized psuedo-atom, and eigenvalues from an all-electron
! calculation of this configuration including unoccupied states that
! are to be adiabatically filled to produce the final configuration.
! The same 2% strategy is used filling pseudo-atom states in this step.

! The all-electron atoms for the maximally-ionized and excited
! configurations are computed by run_config_ae.

 if(cae%stat/=0) then
   call config_ae_warning(jj,cae)
   return
 end if

 nvt=cae%nvt
 nat(:)=cae%nat(:)
 lat(:)=cae%lat(:)
 fat(:,:)=cae%fat(:,:)
 eat(:,:)=cae%eat(:,:)
 eaetst=cae%eaetst

 allocate(rho(mmax),rhocps(mmax))

!first pseudoatom run from reference to maximally-ionized configuration
 do l1=1,4
   nav(l1)=l1-1
 end do

 do kk=1,nvt
   latp(kk)=lat(nc+kk)
   l1=latp(kk)+1
   nav(l1)=nav(l1)+1
   natp(kk)=nav(l1)
   eatp(kk)=eat(nc+kk,1)
   fatp(kk,1)=fat(nc+kk,1)
   fatp(kk,2)=fat(nc+kk,2)
 end do

 rho(:)=rhov(:)

 rhocps(:)=rhomod(:,1)

 call psatom(natp,latp,eatp,fatp,nvt,it,rhocps,rho, &
&           rr,rcmax,mmax,mxprj,iexc,etsttot,nproj,vpuns,lloc, &
&           vkb,evkb,ierr)

 if(ierr/=0) then
   write(6,'(a,a/a,i2)') 'run_config: WARNING for fully non-local PS atom,', &
&        ' stg. 1', &
&        ' WARNING no output for configuration',jj
   deallocate(rho,rhocps)
   return
 end if

!second pseudoatom run from maximally-ionized to excited configuration

 do kk=1,nvt
   eatp(kk)=eat(nc+kk,2)
   fatp(kk,1)=fat(nc+kk,2)
   fatp(kk,2)=fat(nc+kk,3)
 end do
 do ii=1,nvt
 end do

 call psatom(natp,latp,eatp,fatp,nvt,it,rhocps,rho, &
&           rr,rcmax,mmax,mxprj,iexc,etsttot,nproj,vpuns,lloc, &
&           vkb,evkb,ierr)

 if(ierr/=0) then
   write(6,'(a,a/a,i2)') 'run_config: WARNING for fully non-local PS atom,', &
&        ' stg. 2', &
&        ' WARNING no output for configuration',jj
   deallocate(rho,rhocps)
   return
 end if



 write(6,'(/a)') '   n   l     f        eae           eps        diff'
 do ii=1,nc
   write(6,'(2i4,f8.4,f14.8)') nat(ii),lat(ii),fat(ii,3),eat(ii,3)
 end do
 do ii=1,nvt
   if(fat(ii+nc,3)==0.0d0) cycle
   write(6,'(2i4,f8.4,2f14.8,1p,d12.2)') nat(ii+nc),lat(ii+nc),fat(ii+nc,3), &
&   eat(ii+nc,3),eatp(ii),eatp(ii)-eat(ii+nc,3)
 end do

 write(6,'(/a)') '    Total energies and differences'
 write(6,'(a,1p,d16.8,a,d16.8,a,d10.2)') '      AE_ref=',etot, &
& '  AE_tst=',eaetst,'  dif=',eaetst-etot
 write(6,'(a,1p,d16.8,a,d16.8,a,d10.2)') '      PS_ref=',epstot, &
& '  PS_tst=',etsttot,'  dif=',etsttot-epstot
 write(6,'(a,1p,d10.2)') '      PSP excitation error=', &
& eaetst-etot-etsttot+epstot

 deallocate(rho,rhocps)

 return
 end subroutine run_config_ps

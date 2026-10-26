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
! self-consistent pseudoatom calculation

 subroutine psatom(na,la,ea,fat,nv,it,rhoc,rho, &
&           rr,rcmax,mmax,mxprj,iexc,etot,nproj,vpuns,lloc,vlocpw,vkbpw, &
&           evkb,ierr,nkap,vext)

!! @param na  principal quantum number array, dimension nv
!! @param la  angular-momenta
!! @param ea  eigenvalues (input starting guess, output), ea(:,1) for j=l+1/2,
!!            ea(:,2) for j=l-1/2 when nkap=2
!! @param fat  occupancies, fat(:,1) reference and fat(:,2) new configuration
!! @param nv  number of valence states
!! @param it  number of iterations (output)
!! @param rhoc  core charge (model core charge of the pseudopotential)
!! @param rho  valence charge, input all-electron starting guess, output
!!             self-consistent pseudo valence charge
!! @param rr  log radial mesh
!! @param rcmax  maximum core radius for psp
!! @param mmax  size of log grid
!! @param mxprj  dimension of number of projectors
!! @param iexc  exchange-correlation function to be used
!! @param etot  pseudoatom total energy (output)
!! @param nproj  number of VKB projectora to use for each l
!! @param vpuns  unscreened semi-local pseudopotentials (plus differenv vloc if
!!               lloc==4)
!! @param lloc  index-1 of local potential
!! @param vlocpw  the screened local potential's pieces
!! @param vkbpw  Vanderbilt-Kleinman-Bylander projectors as functions of
!!               pieces, last index 1 for j=l+1/2, 2 for j=l-1/2 when nkap=2
!! @param evkb  VKB projector coefficients, last index as for vkbpw
!! @param ierr  non-zero return if error
!! @param nkap  number of kappa per l: 1, or 2 for the fully-relativistic code
!! @param vext  external potential the orbitals also see (a confining wall,
!!              oncvpsp_confine.x, smooth, so outwf_m takes it as part of
!!              the screening), zero for the free atom; nonzero, it is
!!              ramped in with the occupations, and the eigenvalues, which
!!              it can raise above 0, are bracketed around their last
!!              values in place of 10 % of them

 use constants_m, only: dp, mxstate, mxl, mxlin
 use pieces_m, only: piecewise_t
 implicit none

!Input variables

 integer, intent(in) :: mmax,mxprj,iexc,nv,lloc,nkap
 integer, intent(in) :: na(nv),la(nv),nproj(mxlin+1)
 real(dp), intent(in) :: rcmax
 real(dp), intent(in) :: fat(mxstate,2),rr(mmax)
 real(dp), intent(in) :: vpuns(mmax,mxl+2)
 real(dp), intent(in) :: evkb(mxprj,mxl+1,nkap)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1,nkap)
 real(dp), intent(in) :: rhoc(mmax),vext(mmax)

!Input/output variables
 real(dp), intent(inout) :: ea(mxstate,nkap)
 real(dp), intent(inout) :: rho(mmax)

!Output variables
 integer, intent(out) :: ierr
 integer, intent(out) :: it
 real(dp), intent(out) :: etot

!Local variables
 integer :: okb
 real(dp) :: vi(mmax)
 integer :: nin,mch
 real(dp) :: amesh,al
 real(dp) :: dr,eeel,eexc,et,emin,emax,rl,rl1,sd,sn,sls,eeig
 real(dp) :: thl,vn,zval,fj,res,wt,resbest
 real(dp) :: fa(mxstate)
 integer :: ii,jj,l1,itarget,nprj,nflat
 integer :: ikap,kap,ll,mkap
 logical :: convg,walled
 real(dp) :: wext

 real(dp), allocatable :: uu(:),up(:)
 real(dp), allocatable :: vo(:),vi1(:),vo1(:),vxc(:),vp(:)
 real(dp), allocatable :: vtot(:)


! blend parameter for Anderson iterative potential mixing
 real(dp), parameter ::  bl=0.5d0
! convergence threshold for the residual of the screening potential and
! for the change of each eigenvalue (Ha): near round-off, so that the
! result does not depend on the iteration count (doc/smoothness_plan.md,
! S-5)
 real(dp), parameter ::  tolscf=1.0d-12
! the relative energy correction at which the bound states stop while the
! occupations still move to the target: those iterations only lead to the
! target's start, and from it on the states are solved to round-off
 real(dp), parameter ::  tolramp=1.0d-4
! the iterations over which a wall is ramped in: in 50, those of the
! occupations, 50 of the 68 iterations of Bi in 0.1 Ha/bohr**2; 10 hold
! every wall of the test inputs that 50 did, up to 3 Ha/bohr**2 from the
! free atom (oncvpsp_confine.x starts each wall from the one before)
 integer, parameter :: nwramp=10
 real(dp) :: etol

 ierr=0

 allocate(uu(mmax),up(mmax))
 allocate(vo(mmax),vi1(mmax),vo1(mmax),vxc(mmax),vp(mmax))
 allocate(vtot(mmax))

! this seems necessary in sratom, so I might as well do it  here too
! don't ask why!
 uu(:)=0.d0; up(:)=0.d0
 vo(:)=0.d0; vi1(:)=0.d0; vo1(:)=0.d0; vxc(:)=0.d0; vp(:)=0.d0
 vtot(:)=0.d0
 dr=0.d0; eeel=0.d0; eexc=0.d0; et=0.d0; emin=0.d0; emax=0.d0
 rl=0.d0; rl1=0.d0; sd=0.d0; sn=0.d0; sls=0.d0; eeig=0.d0
 thl=0.d0; vn=0.d0; zval=0.d0
 nin=0; mch=0

 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

 nprj=0
 do l1=1,4
   nprj=max(nprj,nproj(l1))
 end do

! total valence charge, initial config
 zval=0.0d0
 do ii=1,nv
   fa(ii)=fat(ii,1)
   zval=zval+fa(ii)
 end do

! first iteration at the target occupations (none yet)
 itarget=huge(1)
 walled=any(vext(:)/=0.0d0)
 wext=1.0d0
 if(walled) wext=1.0d0/dble(nwramp)

! screening potential for pseudocharge
! input rho is assumed to be valence rho from all-electron calculation
! which shouldn''t be a bad approximation for a starting screened
! potential

 call vout(1,rho,rhoc,vi,vxc,zval,eeel,eexc,rr,mmax,iexc)


! big self  self-consietency loop
 do it=1,100

! convergence is only to be considered at the target occupations
   convg=all(fa(1:nv)==fat(1:nv,2)) .and. wext==1.0d0
   if(convg) itarget=min(itarget,it)
   etol=tolramp
   if(convg) etol=0.0d0

! solve for bound states in turn

   eeig=0.0d0

   vtot(:)=vpuns(:,lloc+1)+vi(:)
   if(walled) vtot(:)=vtot(:)+wext*vext(:)

   rho(:)=0.0d0

   do ii=1,nv

     if(fa(ii)==0.0d0) then
       cycle
     end if
     ierr = 0
     l1=la(ii)+1

     ll=la(ii)
     if(nkap==1 .or. ll==0) then
       mkap=1
     else
       mkap=2
     end if
! loop on J = ll +/- 1/2 when nkap=2, the occupation shared out as 2j+1
     do ikap=1,mkap
       if(nkap==1) then
         kap=0
         fj=fa(ii)
       else if(ikap==1) then
         kap=-(ll+1)
         fj=dble(ll+1)*fa(ii)/dble(2*ll+1)
       else
         kap=  ll
         fj=dble(ll)*fa(ii)/dble(2*ll+1)
       end if
       et=ea(ii,ikap)

       emin=1.10d0*ea(ii,ikap)
       emax=0.90d0*ea(ii,ikap)
       if(walled) then
         emin=ea(ii,ikap)-max(0.1d0*abs(ea(ii,ikap)),0.5d0)
         emax=ea(ii,ikap)+max(0.1d0*abs(ea(ii,ikap)),0.5d0)
       end if
       call lschvkbb_tol(na(ii),la(ii),nproj(l1),ierr,et,emin,emax, &
&                       rr,vtot,vlocpw,vkbpw(1:nproj(l1),l1,ikap), &
&                       evkb(1,l1,ikap),uu,up,mmax,mch,etol)

       if(ierr /= 0) then
         if(nkap==2) then
           write(6,'(/a,3i4)') &
&                'psatom_r: WARNING lschvkbb convergence error n,l,kap,iter=', &
&                 na(ii),la(ii),kap,it
         else
           write(6,'(/a,3i4)') &
&                'psatom: WARNING lschvkbb convergence error n,l,iter=', &
&                na(ii),la(ii),it
         end if
         write(6,'(a,1p,4e14.6)') 'ea(ii),et,emin,emax',ea(ii,ikap),et, &
&              emin,emax
         exit
       end if

! each eigenvalue settled to tolscf (see the convergence test below)
       if(abs(ea(ii,ikap)-et)>tolscf) convg=.false.
       ea(ii,ikap)=et

! accumulate charge and eigenvalues
       eeig = eeig + fj*ea(ii,ikap)
       rho(:)=rho(:) + fj*(uu(:)/rr(:))**2
     end do !ikap
     if(ierr/=0) exit

   end do !ii

   if(ierr/=0) exit

! total valence charge
   zval=0.0d0
   do ii=1,nv
     zval=zval+fa(ii)
   end do

! output potential
   call vout(1,rho,rhoc,vo,vxc,zval,eeel,eexc,rr,mmax,iexc)

! converged when no eigenvalue has moved by more than tolscf and the
! screening potential is self-consistent to tolscf, or its residual has
! stopped falling at its round-off floor (ae_converged): 1.4-3e-12 for
! Li with a near-bound s ghost, where 1e-12 was never reached
   call ae_residual(rho,vi,vo,zval,rr,mmax,res)
   call ae_converged(it,res,tolscf,resbest,nflat,convg)

! generate next iteration using d. g. anderson''s
! method, once the previous iteration was at the target occupations too
   thl=0.0d0
   if(it>itarget) then
     sn=0.0d0
     sd=0.0d0
     do ii=1,mmax
       rl=vo(ii)-vi(ii)
       rl1=vo1(ii)-vi1(ii)
       dr=rl-rl1
       sn=sn + rl*dr*rr(ii)**2
       sd=sd + dr*dr*rr(ii)**2
     end do
     thl=sn/sd
   end if

   do ii=1,mmax
     vn=(1.0d0-bl)*((1.0d0-thl)*vi(ii) + thl*vi1(ii)) &
&     + bl*((1.0d0-thl)*vo(ii) + thl*vo1(ii))
     vi1(ii)=vi(ii)
     vo1(ii)=vo(ii)
     vi(ii)=vn
   end do

   if(convg) exit

   if(it==100 .and. .not. convg) then
     if(nkap==2) then
       write(6,'(/a)') 'psatom_r: WARNING potential failed to converge'
     else
       write(6,'(/a)') 'psatom: potential failed to converge'
     end if
     ierr=101
   end if

! switch from reference to new configuration in 50 steps, exactly the
! new one from the 50th on (and throughout if they are the same)
! (assigned, since fat1 + 1*(fat2 - fat1) need not round to fat2); a wall
! in nwramp
   if(walled) wext=min(1.0d0,dble(it+1)/dble(nwramp))
   if(it>=50) then
     fa(1:nv)=fat(1:nv,2)
   else
     wt=dble(it)/50.0d0
     do ii=1,nv
       fa(ii)=fat(ii,1) + wt*(fat(ii,2)-fat(ii,1))
     end do
   end if


 end do !it

! total energy output

! output potential for e-e interactions

 call vout(1,rho,rhoc,vo,vxc,zval,eeel,eexc, &
&          rr,mmax,iexc)

 etot =  eeig + eexc - 0.5d0*eeel

 deallocate(uu,up)
 deallocate(vo,vi1,vo1,vxc,vp)
 deallocate(vtot)
 return

 end subroutine psatom

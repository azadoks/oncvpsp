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
 subroutine relatom(na,la,ea,fa,rpk,nc,ncv,it,rhoc,rho, &
&           rr,vi,zz,mmax,iexc,etot,eterm,res,ierr)

!! self-consistent fully-relativistic all-electron atom
!! calculation using log mesh
!! partially-occupied orbitals are weighted proportionally to 2j+1
!!
!! @param na  principal quantum number array, dimension ncv
!! @param la  angular-momenta
!! @param ea  eigenvalues (output), ea(:,1) for j=l+1/2, ea(:,2) for j=l-1/2
!! @param fa  occupancies
!! @param rpk  radius of outermost peak of wave function, indexed as ea
!! @param nc  number of core states
!! @param ncv  number of core+valence states
!! @param it  number of iterations (output)
!! @param rhoc  core charge (output)
!! @param rho  total core+valence charge (output)
!! @param rr  log radial mesh
!! @param vi  all-electron potential (output), the one the output orbitals
!!            are eigenstates of
!! @param zz  atomic number
!! @param mmax  size of log grid
!! @param iexc  exchange-correlation function to be used
!! @param etot  all-electron total energy (output)
!! @param eterm  kinetic, electron-nucleus, Hartree and exchange-correlation
!!               energies (output)
!! @param res  residual of the screening potential (output, ae_residual)
!! @param ierr  error flag

 use constants_m, only: dp, mxstate, mxkap
 implicit none

!Input variables

 integer, intent(in) :: mmax,iexc,nc,ncv
 integer, intent(in) :: na(ncv),la(ncv)
 real(dp), intent(in) :: zz
 real(dp), intent(in) :: fa(ncv),rr(mmax)

!Input/output variables
 integer, intent(inout) :: ierr
 real(dp), intent(inout) :: ea(mxstate,mxkap),rpk(mxstate,mxkap)

!Output variables
 integer, intent(out) :: it
 real(dp), intent(out) :: etot,eterm(4),res
 real(dp), intent(out) :: rho(mmax),rhoc(mmax),vi(mmax)

!Local variables
 real(dp) :: tfapot
 integer :: nin,mch
 real(dp) :: amesh,al
 real(dp) :: dr,eeel,eexc,et,rl,rl1,sd,sf,sn,eeig
 real(dp) :: thl,vn,zion,fj
 integer :: ii,jj,ll,ikap,kap,mkap
 logical :: convg

 real(dp), allocatable :: uu(:,:),up(:,:)
 real(dp), allocatable :: vsi(:),vso(:),vsi1(:),vso1(:),vxc(:)

! blend parameter for Anderson iterative potential mixing
 real(dp), parameter ::  bl=0.5d0
! convergence threshold for the residual of the screening potential (Ha)
 real(dp), parameter ::  tolres=1.0d-8

 allocate(uu(mmax,2),up(mmax,2))
 allocate(vsi(mmax),vso(mmax),vsi1(mmax),vso1(mmax),vxc(mmax))

! why all this is necessary is unclear, but it seems to be
 uu(:,:)=0.d0; up(:,:)=0.d0; vso(:)=0.d0; vsi1(:)=0.d0; vso1(:)=0.d0
 vxc(:)=0.d0
 dr=0.d0; eeel=0.d0; eexc=0.d0; et=0.d0; rl=0.d0; rl1=0.d0
 sd=0.d0; sf=0.d0; sn=0.d0; eeig=0.d0; thl=0.d0; vn=0.d0; zion=0.d0
 nin=0; mch=0

 al = 0.01d0 * log(rr(101) / rr(1))
 amesh = exp(al)

 do ii=1,mmax
   vi(ii)=tfapot(rr(ii),zz)
 end do

! the screening (Hartree plus exchange-correlation) potential is mixed,
! and the nuclear potential added to it for the radial equation
 do ii=1,mmax
   vsi(ii)=vi(ii) + zz/rr(ii)
 end do

! starting approximation for energies
 sf=0.0d0
 do ii=1,ncv
   sf=sf+fa(ii)
   zion=zz+1.0d0-sf
   do ikap=1,2
     ea(ii,ikap)=-0.5d0*(zion/na(ii))**2
     if(ea(ii,ikap)>vi(mmax)) ea(ii,ikap)=2.0d0*vi(mmax)
   end do
   if(la(ii)==0) ea(ii,2)=0.0d0
 end do

! big self  self-consietency loop

 do it=1,100
   convg=.true.

   rhoc(:) = 0.0d0
   rho(:)=0.0d0

! solve for bound states in turn
   eeig=0.0d0
   do ii=1,ncv
     ll = la(ii)
     if(ll==0) then
       mkap=1
     else
       mkap=2
     end if
! loop on J = ll +/- 1/2
     do ikap=1,mkap
       if(ikap==1) then
         kap=-(ll+1)
         fj=dble(ll+1)*fa(ii)/dble(2*ll+1)
       else
         kap=  ll
         fj=dble(ll)*fa(ii)/dble(2*ll+1)
       end if
       et=ea(ii,ikap)
       ierr = 0

       call ldiracfb(na(ii),ll,kap,ierr,et, &
&                  rr,zz,vi,uu,up,mmax,mch)
! a state with no bound solution in this potential gets its box state;
! the SCF goes on but cannot converge with it
       if(ierr == 10) then
         convg=.false.
         ierr=0
       end if

       if(ierr /= 0) then
         write(6,'(/2a,5i4)') 'relatom: ldiracfb convergence ERROR', &
&             ' n,l,kap,iter,ierr=',na(ii),ll,kap,it,ierr
         error stop 1
       end if

! overall convergence criterion based on eps within lschfb
       if(ea(ii,ikap)/= et) convg=.false.
       ea(ii,ikap)=et

! accumulate charge and eigenvalues
       eeig = eeig + fj*ea(ii,ikap)
       rho(:)=rho(:) + fj*((uu(:,1)/rr(:))**2 + (uu(:,2)/rr(:))**2)
       if(ii<=nc) then
         rhoc(:)=rhoc(:) + fj*((uu(:,1)/rr(:))**2 + (uu(:,2)/rr(:))**2)
       end if


! find outermost peak of wavefunction
       do jj=mch-1,1,-1
         if(up(jj,1)*up(jj+1,1)<0.0d0) then
           rpk(ii,ikap)=rr(jj)
           exit
         end if
       end do
     end do
   end do

   if(ierr/=0) then
     exit
   end if


! output screening potential, and the energies of the orbitals of the
! input one
   call vout(0,rho,rhoc,vso,vxc,sf,eeel,eexc, &
&            rr,mmax,iexc)
   call ae_energies(rho,vsi,vxc,zz,eeig,eeel,eexc,rr,mmax,etot,eterm)

! converged when no eigenvalue has moved and the screening potential is
! self-consistent; the eigenvalues alone can stand still from a start far
! from self-consistency
   call ae_residual(rho,vsi,vso,sf,rr,mmax,res)
   if(it<2 .or. res>tolres) convg=.false.
   if(convg) exit

! generate next iteration using d. g. anderson''s
! method
   thl=0.0d0
   if(it>1) then
     sn=0.0d0
     sd=0.0d0
     do ii=1,mmax
       rl=vso(ii)-vsi(ii)
       rl1=vso1(ii)-vsi1(ii)
       dr=rl-rl1
       sn=sn + rl*dr*rr(ii)**2
       sd=sd + dr*dr*rr(ii)**2
     end do
     thl=sn/sd
   end if

   do ii=1,mmax
     vn=(1.0d0-bl)*((1.0d0-thl)*vsi(ii) + thl*vsi1(ii)) &
&     + bl*((1.0d0-thl)*vso(ii) + thl*vso1(ii))
     vsi1(ii)=vsi(ii)
     vso1(ii)=vso(ii)
     vsi(ii)=vn
     vi(ii)=vn - zz/rr(ii)
   end do

   if(it==100 .and. .not. convg) then
     write(6,'(/a)') 'relatom: WARNING failed to converge'
   end if

 end do !it

 if(.not. convg .and. ierr==0) then
   ierr=100
 end if

 deallocate(uu,up)
 deallocate(vsi,vso,vsi1,vso1,vxc)
 return

 end subroutine relatom

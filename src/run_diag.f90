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

!> diagnostics for semi-local and Vanderbilt-Kleinman-bylander pseudopotentials
!> checks bound-state energies, norms, and slopes, and pseudo-bound state
!> quantities for positive-energy projectors by matching all-electron
!> log derivatives at rc.
!> Should test as essentially exact for non-relativitic calculations
!> Error for relativistic is B matrix Hermiticity error
subroutine run_diag(lmax,npa,epa,lloc,irc, &
&                    vkb,evkb,nproj,rr,vfull,vp,zz,mmax,mxprj,srel)
   implicit none
   integer, parameter :: dp=kind(1.0d0)

!Input variables
   !> maximum angular momentum
   integer, intent(in) :: lmax
   !> l for local potential
   integer, intent(in) :: lloc
   !> size of radial grid
   integer, intent(in) :: mmax
   !> dimension of number of projectors
   integer, intent(in) :: mxprj
   !> principal quantum number for corresponding all-electron state
   integer, intent(in) :: npa(mxprj,6)
   !> core radii indices
   integer, intent(in) :: irc(6)
   !> number of vkb projectors for each l
   integer, intent(in) :: nproj(6)
   !> atomic number
   real(dp), intent(in) :: zz
   !> log radial grid
   real(dp), intent(in) :: rr(mmax)
   !> semi-local pseudopotentials (vp(:,5) is local potential if linear comb.)
   real(dp), intent(in) :: vp(mmax,5)
   !> all-electron potential
   real(dp), intent(in) :: vfull(mmax)
   !> VKB projectors
   real(dp), intent(in) :: vkb(mmax,mxprj,4)
   !> bound-state or scattering state reference energies for vkb potentials
   real(dp), intent(in) :: epa(mxprj,6)
   !> coefficients of VKB projectors
   real(dp), intent(in) :: evkb(mxprj,4)
   !> .true. for scalar-relativistic, .false. for non-relativistic
   logical, intent(in) :: srel

!Output variables - printing only

!Local variables
   integer :: ll,l1,ierr,mch,mchf
   integer :: iprj,nnae,nnp,npr
   real(dp) :: al,emax,emin,etest,sls,umch,upmch,uldf,gam,gpr
   real(dp), allocatable :: uu(:),up(:)

   allocate(uu(mmax),up(mmax))

   al = 0.01d0 * dlog(rr(101)/rr(1))

! loop for diagnostic output using Vanderbilt Kleinman-Bylander projectors
!
   write(6,'(/a/a)') &
   & 'Diagnostic tests using Vanderbilt-Kleinman-Bylander pseudopotentials',&
   & '  1 or more projectors used as specified by nproj input data'
   write(6,'(/2a)')'   l    rcore       rmatch      e in        ', &
   &   'delta e   norm test  slope test'
!
   do l1 = 1, lmax + 1
      write(6,'(a)') ''
      ll = l1 - 1
      nnp=l1
      do iprj=1,max(1,nproj(l1))
         npr=nproj(l1)

!   find cutoff radius for projectors
         mchf=max(irc(l1),irc(lloc+1))+5

         etest=epa(iprj,l1)
         if(etest<0.0d0) then
            call lschfb(npa(iprj,l1),ll,ierr,etest, &
            &                rr,vfull,uu,up,zz,mmax,mch,srel)
            if(ierr /= 0) then
               write(6,'(/a,3i4)') 'run_diag: lschfb convergence ERROR n,l,ierr=', &
               &        npa(iprj,l1),ll,ierr
               stop
            end if

         else
            call lschfs(nnae,ll,ierr,etest, &
            &               rr,vfull,uu,up,zz,mmax,mchf,srel)
         end if

         umch=uu(mchf)
         upmch=up(mchf)
         uldf=upmch/umch
         etest=epa(iprj,l1)

         if(l1==lloc+1) npr=0

         etest=epa(iprj,l1)
         emax=etest+0.1d0*dabs(etest)+0.05d0
         emin=etest-0.1d0*dabs(etest)-0.05d0


         if(epa(iprj,l1)<0.0d0) then
            sls=(l1-1)*l1

            emax=dmin1(emax,0.0d0)

            call lschvkbb(ll+iprj,ll,npr,ierr,etest,emin,emax, &
            &                   rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1), &
            &                   uu,up,mmax,mch)
            if(ierr/=0) then
               write(6,'(a,3i4,1p,2e16.8)') 'run_diag: lschvkbb ERROR', &
               &             iprj,ll,ierr,epa(iprj,l1),etest
            end if

            nnp=nnp+1
         else
!     calculate effective pseudo wave function principal quantum number
!     from all-electron node count
            nnp=nnae-npa(1,l1)+ll+1
            call lschvkbbe(nnp,ll,npr,ierr,etest,uldf,emin,emax, &
            &                   rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1), &
            &                   uu,up,mmax,mchf)

            if(ierr/=0) then
               write(6,'(a,4i4,1p,2e16.8)') 'run_diag: lschvkbbe ERROR', &
               &             iprj,nnp,ll,ierr,epa(iprj,l1),etest

            end if
         end if  !epa<0 (bound or not)

         gam=dabs(umch/uu(mchf))
         gpr=dabs(upmch/up(mchf))
!
         write(6,'(i4,6f12.7)') ll,rr(irc(l1)),rr(mchf),epa(iprj,l1), &
         &         etest-epa(iprj,l1),gam,gpr

      end do  !iprj
   end do  !l1

   deallocate(uu,up)
   return
end subroutine run_diag

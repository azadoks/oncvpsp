!
! Copyright (c) 2026 by Austin Zadoks
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
! The machine-readable block of the diagnostics, for programs that read the
! output (an optimizer taking finite differences of them): the values that
! the tables above print to 6 or 2 decimals, to 17 significant digits.

 subroutine diagout(lmax,nproj,mxprj,nkap,cvgplt,npsh,epsh2,depsh,rphs, &
&                   pshae,pshps,ngh,ghn,ghe,ghp,wcut,errw,ecutw,wcurve)

!! writes the block "Begin DIAGNOSTICS" ... "End DIAGNOSTICS", one value
!! per line, each line tagged with its quantity:
!!
!! ecut l ikap iproj err ecut: the cutoff (Ha) at which the residual
!! kinetic energy per electron of the optimized wave of projector iproj
!! falls to err (Ha), for the levels err that it reaches (run_optimize)
!!
!! ld_r l ikap r: the radius (bohr) of the log derivatives of l, ikap
!!
!! ld l ikap e ae ps: atan(r u'/u) at that radius and energy e (Ha), of
!! the all-electron and the pseudo wave (run_phsft)
!!
!! ghost l ikap n e_ref e_basis: a bound state of the nonlocal
!! Schroedinger equation (Ha; 0 for the vacuum level when l has none) and
!! the eigenvalue of the barrier-basis diagonalization matched to it
!! (run_ghosts); a lower eigenvalue is a GHOST(-), reported below
!! min(e_ref - 1e-4, e_ref*(1 + 1e-4)), so e_basis - e_ref is the margin
!!
!! ghost_r l ikap e r: of the basis eigenstates more than 0.05 Ha above
!! the reference, the one of least <r>/rc: its energy and <r>/rc, a
!! GHOST(+) below 1
!!
!! ecut_w err ecut: the cutoff (Ha) at which the weighted mean of the
!! residual kinetic energies per electron falls to err (Ha), the cutoff of
!! the model core's budget with budget_cutoff = "weighted"
!! (modcore_budget_table); only then written, with
!!
!! ecut_w_weight l ikap iproj w: the weight of each pseudo-wave in the mean
!! (those with w > 0)
!!
!! ikap is 1 for j=l+1/2 (and for nkap=1), 2 for j=l-1/2
!!
!! @param lmax  maximum angular momentum
!! @param nproj  number of projectors for each l
!! @param mxprj  dimension of the number of projectors
!! @param nkap  number of kappa per l: 1, or 2 for irel=2
!! @param cvgplt  the crossings, (Ecut, err) per level, projector, l, ikap
!! @param npsh  number of energies of the log derivatives
!! @param epsh2  the first (highest) energy
!! @param depsh  the energy step (downwards)
!! @param rphs  the radius of the log derivatives of each l
!! @param pshae  the all-electron log derivatives
!! @param pshps  the pseudopotential log derivatives
!! @param ngh  rows of the bound-state ghost test per l and kap
!! @param ghn  their principal quantum numbers
!! @param ghe  their energies, e_ref and e_basis
!! @param ghp  the least-localized-ghost margin per l and kap, e and <r>/rc
!!             (0 if no state qualifies)
!! @param wcut  whether to write the weighted cutoff
!! @param errw  its error per electron (Ha)
!! @param ecutw  the weighted cutoff (Ha)
!! @param wcurve  the weights of the pseudo-waves

 use constants_m, only: dp, mxl, mxlin, nerrlv, mxstate, mxkap
 implicit none

!Input variables
 integer, intent(in) :: lmax,mxprj,nkap,npsh
 integer, intent(in) :: nproj(mxlin+1)
 real(dp), intent(in) :: cvgplt(2,nerrlv,mxprj,mxl+1,nkap)
 real(dp), intent(in) :: epsh2,depsh
 real(dp), intent(in) :: rphs(4),pshae(npsh,4,nkap),pshps(npsh,4,nkap)
 integer, intent(in) :: ngh(mxl+1,mxkap),ghn(mxstate,mxl+1,mxkap)
 real(dp), intent(in) :: ghe(2,mxstate,mxl+1,mxkap),ghp(2,mxl+1,mxkap)
 logical, intent(in) :: wcut
 real(dp), intent(in) :: errw,ecutw,wcurve(mxprj,mxl+1,nkap)

!Local variables
 integer :: ii,jj,l1,ikap,mkap,iprj

 write(6,'(/a)') 'Begin DIAGNOSTICS'
 write(6,'(a)') '# ecut l ikap iproj err(Ha) ecut(Ha)'
 do l1=1,lmax+1
   mkap=nkap
   if(l1==1) mkap=1
   do ikap=1,mkap
     do iprj=1,nproj(l1)
       do jj=1,nerrlv
         if(cvgplt(1,jj,iprj,l1,ikap)/=0.0d0) then
           write(6,'(a,3i4,1p,2e24.16)') 'ecut',l1-1,ikap,iprj, &
&                cvgplt(2,jj,iprj,l1,ikap),cvgplt(1,jj,iprj,l1,ikap)
         end if
       end do
     end do
   end do
 end do

 if(wcut) then
   write(6,'(a)') '# ecut_w err(Ha) ecut(Ha)'
   write(6,'(a,1p,2e24.16)') 'ecut_w',errw,ecutw
   write(6,'(a)') '# ecut_w_weight l ikap iproj w'
   do l1=1,lmax+1
     mkap=nkap
     if(l1==1) mkap=1
     do ikap=1,mkap
       do iprj=1,mxprj
         if(wcurve(iprj,l1,ikap)>0.0d0) then
           write(6,'(a,3i4,1p,e24.16)') 'ecut_w_weight',l1-1,ikap,iprj, &
&                wcurve(iprj,l1,ikap)
         end if
       end do
     end do
   end do
 end if

 write(6,'(a)') '# ld_r l ikap r(bohr)'
 write(6,'(a)') '# ld l ikap e(Ha) atan(r*u''/u)_ae atan(r*u''/u)_ps'
 do l1=1,4
   mkap=nkap
   if(l1==1) mkap=1
   do ikap=1,mkap
     write(6,'(a,2i4,1p,e24.16)') 'ld_r',l1-1,ikap,rphs(l1)
     do ii=1,npsh
       write(6,'(a,2i4,1p,3e24.16)') 'ld',l1-1,ikap, &
&            epsh2-depsh*real(ii-1,dp),pshae(ii,l1,ikap),pshps(ii,l1,ikap)
     end do
   end do
 end do
 write(6,'(a)') '# ghost l ikap n e_ref(Ha) e_basis(Ha)'
 write(6,'(a)') '# ghost_r l ikap e(Ha) <r>/rc'
 do l1=1,lmax+1
   mkap=nkap
   if(l1==1) mkap=1
   do ikap=1,mkap
     do jj=1,ngh(l1,ikap)
       write(6,'(a,3i4,1p,2e24.16)') 'ghost',l1-1,ikap,ghn(jj,l1,ikap), &
&            ghe(:,jj,l1,ikap)
     end do
     if(ghp(2,l1,ikap)>0.0d0) then
       write(6,'(a,2i4,1p,2e24.16)') 'ghost_r',l1-1,ikap,ghp(:,l1,ikap)
     end if
   end do
 end do
 write(6,'(a)') 'End DIAGNOSTICS'

 return
 end subroutine diagout

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
 subroutine run_phsft(lmax,lloc,nproj,epa,npsh,epsh2,depsh,rxpsh,vlocpw, &
&                      vkbpw,evkb,rr,vfull,vp,zz,mmax,mxprj,irc,irel,nkap, &
&                      rphs,pshae,pshps)

!! computes log derivatives, actually atan(r * ((d psi(r)/dr)/psi(r)))
!! at rr(irphs) comparing all-electron with Vanderbilt-Kleinman-Bylander
!! results for 1 and 2 projectors, or the semi-local pseudpotential
!! when that is the local potential for some l
!! the computed quantity is reminiscent of a scattering phase shift, but isn't
!! for irel=2, the data for j=l+1/2 and j=l-1/2 follow each other, the
!! former labelled with -l
!!
!! @param lmax  maximum angular momentum
!! @param lloc  l for local potential
!! @param nproj  number of VKB projectors for each l
!! @param epa  bound-state or scattering state reference energies for vkb
!!             potentials (last index 1 for j=l+1/2, 2 for j=l-1/2, as in vkb,
!!             evkb and vp)
!! @param npsh  number of energies, from epsh2 down
!! @param epsh2  high energy limit for "phase shift" calculation
!! @param depsh  energy increment
!! @param rxpsh  radius (bohr) at which log derivatives are computed; if
!!               negative, rr(irc(l)+2) is used, and rr(irc(lloc)) for l>lmax
!! @param vlocpw  the screened local potential's pieces
!! @param vkbpw  VKB projectors as functions of pieces
!! @param evkb  coefficients of VKB projectors
!! @param rr  log radial grid
!! @param vfull  all-electron potential
!! @param vp  semi-local pseudopotentials (vp(:,5) is local potential if linear
!!            comb.)
!! @param zz  atomic number
!! @param mmax  size of radial grid
!! @param mxprj  dimension of number of projectors
!! @param irc  core radii indices
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param nkap  number of kappa per l: 1, or 2 for irel=2
!! @param rphs  the radius of the log derivatives of each l (output)
!! @param pshae  the all-electron log derivatives, per energy, l and kappa
!!               (output; zero for l=0, ikap=2)
!! @param pshps  the pseudopotential log derivatives, likewise (output)

 use constants_m, only: dp, mxl, mxlin
 use pieces_m, only: piecewise_t
 implicit none

!Input variables
 integer, intent(in) :: lmax,lloc,mmax,mxprj,irel,nkap,npsh
 integer, intent(in) :: nproj(mxlin+1),irc(mxlin+1)
 real(dp), intent(in) :: epsh2,depsh,zz,rxpsh
 real(dp), intent(in) :: rr(mmax),vp(mmax,mxl+2,nkap),epa(mxprj,mxlin+1,nkap)
 real(dp), intent(in) :: vfull(mmax)
 real(dp), intent(in) :: evkb(mxprj,mxl+1,nkap)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1,nkap)

!Output variables
 real(dp), intent(out) :: rphs(4),pshae(npsh,4,nkap),pshps(npsh,4,nkap)

!Local variables
 integer :: ii,irphs,ll,l1
 integer :: ikap,kap,mkap,lout
 real(dp) :: epsh
 character(len=11) :: rname

 pshae(:,:,:)=0.0d0
 pshps(:,:,:)=0.0d0

!the routine named in the messages
 if(irel==2) then
   rname='run_phsft_r'
 else
   rname='run_phsft'
 end if

! loop for phase shift calculation -- full, then local or Kleinman-
! Bylander / Vanderbilt

 do l1 = 1, 4

   ll = l1 - 1
   if(ll <= lmax) then
     if(rxpsh < 0.0) then
       irphs = irc(l1) + 2
     else
       irphs = 0
       do ii = 1,mmax
         if(rr(ii) >= rxpsh) then
           irphs = ii
           exit
         end if
       end do
       if(irphs == 0) then
         write(6, '(/2a,f12.6,a,f12.6)') &
&           trim(rname)//': ERROR rxpsh beyond radial mesh', &
&           ', rxpsh=', rxpsh, ', rr(mmax)=', rr(mmax)
         error stop 1
       end if
! check_data keeps rxpsh beyond every rc, so only a direct call gets here
       if(irphs < irc(l1)) then
         write(6, '(/2a,i0,2(a,f12.6))') &
&           trim(rname)//': WARNING rxpsh inside the core radius', &
&           ', l=', ll, ', rxpsh=', rxpsh, ', rc=', rr(irc(l1))
       end if
     end if
   else
     if(rxpsh < 0.0) then
       irphs = irc(lloc + 1)
     else
       irphs = 0
       do ii = 1,mmax
         if(rr(ii) >= rxpsh) then
           irphs = ii
           exit
         end if
       end do
       if(irphs == 0) then
         write(6, '(/2a,f12.6,a,f12.6)') &
&           trim(rname)//': ERROR rxpsh beyond radial mesh', &
&           ', rxpsh=', rxpsh, ', rr(mmax)=', rr(mmax)
         error stop 1
       end if
! check_data keeps rxpsh beyond every rc, so only a direct call gets here
       if(irphs < irc(lloc+1)) then
         write(6, '(/2a,i0,2(a,f12.6))') &
&           trim(rname)//': WARNING rxpsh inside ' // &
&           'the local potential core radius', &
&           ', l=', ll, ', rxpsh=', rxpsh, ', rc=', rr(irc(lloc+1))
       end if
     end if
   end if

   rphs(l1)=rr(irphs)
   if(nkap==1 .or. l1==1) then
     mkap=1
   else
     mkap=2
   end if
! loop on J = ll +/- 1/2 when irel=2
   do ikap=1,mkap
     if(ikap==1) kap=-(ll+1)
     if(ikap==2) kap=  ll
     call fphsft(ll,kap,epsh2,depsh,pshae(1,l1,ikap),rr,vfull,zz,mmax, &
&                irphs,npsh,irel)
     if(ll == lloc) then
       call  vkbphsft(ll,0,epsh2,depsh,epa(1,l1,ikap),pshae(1,l1,ikap), &
&                          pshps(1,l1,ikap), &
&                          rr,vp(1,lloc+1,ikap),vlocpw,vkbpw(1:0,l1,ikap), &
&                          evkb(1,l1,ikap),mmax,irphs,npsh)
     else
       call  vkbphsft(ll,nproj(l1),epsh2,depsh,epa(1,l1,ikap), &
&                          pshae(1,l1,ikap),pshps(1,l1,ikap), &
&                          rr,vp(1,lloc+1,ikap),vlocpw, &
&                          vkbpw(1:nproj(l1),l1,ikap), &
&                          evkb(1,l1,ikap),mmax,irphs,npsh)
     end if

     lout=ll
     if(nkap==2 .and. ikap==1) lout=-ll
     write(6,'(/a,i2)') 'log derivativve data for plotting, l=',ll
     write(6,'(a,f6.2)') 'atan(r * ((d psi(r)/dr)/psi(r))), r=',rr(irphs)
     write(6,'(a/)') 'l, energy, all-electron, pseudopotential'
     do ii = 1, npsh
       epsh = epsh2 - depsh * real(ii - 1,dp)
       write(6,'(a, i6, 3f12.6)') '! ',lout, epsh, pshae(ii,l1,ikap), &
&            pshps(ii,l1,ikap)
     end do
   end do !ikap
 end do !l1
 return
 end subroutine run_phsft

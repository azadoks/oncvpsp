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

 module run_phsft_m

 implicit none
 private
 public :: run_phsft

 contains

 subroutine run_phsft(lmax,lloc,nproj,epa,npsh,epsh2,depsh,rxpsh,vlocpw,vkbpw, &
&                      evkb,rr,vfull,vp,zz,mmax,rc,irel,c_light,nkap, &
&                      rphs,pshae,pshps,iout,ifail)

!! computes log derivatives, actually atan(r * ((d psi(r)/dr)/psi(r)))
!! at the radius rphs comparing all-electron with Vanderbilt-Kleinman-Bylander
!! results for 1 and 2 projectors, or the semi-local pseudpotential
!! when that is the local potential for some l
!! the computed quantity is reminiscent of a scattering phase shift, but isn't
!! for irel=2, the data for j=l+1/2 and j=l-1/2 follow each other, the
!! former labelled with -l
!! the phases are continuous in the energy, the radius and the potentials
!! (phsft_m, phsft_pair); the pseudopotential's are shifted by the multiple
!! of pi that makes them equal the all-electron ones at the reference
!! energy epa
!!
!! @param lmax  maximum angular momentum
!! @param lloc  l of the local potential, 4 for the polynomial
!! @param nproj  number of projectors for each l
!! @param epa  bound-state or scattering state reference energies for vkb
!!             potentials (last index 1 for j=l+1/2, 2 for j=l-1/2, as in vkb,
!!             evkb and vp)
!! @param npsh  number of energies, from epsh2 down
!! @param epsh2  upper limit of the log derivative analysis energy range (Ha)
!! @param depsh  energy step of the log derivative analysis (Ha)
!! @param rxpsh  radius (bohr) at which log derivatives are computed; if
!!               negative, rc(l)*amesh**2 (the second mesh point beyond rc(l)
!!               when that is a mesh point), and rc(lloc) for l>lmax
!! @param vlocpw  the screened local potential as pieces
!! @param vkbpw  the projectors as pieces
!! @param evkb  coefficients of the projectors
!! @param rr  log radial mesh (bohr)
!! @param vfull  all-electron potential (Ha)
!! @param vp  semi-local pseudopotentials (vp(:,5) is local potential if linear
!!            comb.)
!! @param zz  atomic number
!! @param mmax  size of the log radial mesh
!! @param rc  core radii (bohr)
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param c_light  speed of light (Hartree atomic units)
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param rphs  the radius of the log derivatives of each l (output)
!! @param pshae  the all-electron log derivatives, per energy, l and kappa
!!               (output; zero for l=0, ikap=2)
!! @param pshps  the pseudopotential log derivatives, likewise (output)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp, mxl, mxlin, mxprj
 use kappa_m, only: kap_count, kap_of
 use pieces_m, only: piecewise_t
 use phsft_m, only: phsft_mch, phsft_pair

!Input variables
 integer, intent(in) :: lmax,lloc,mmax,irel,nkap,npsh
 real(dp), intent(in) :: c_light
 integer, intent(in) :: nproj(mxlin+1)
 real(dp), intent(in) :: epsh2,depsh,zz,rxpsh
 real(dp), intent(in) :: rr(mmax),vp(mmax,mxl+2,nkap),epa(mxprj,mxlin+1,nkap)
 real(dp), intent(in) :: vfull(mmax),rc(mxlin+1)
 real(dp), intent(in) :: evkb(mxprj,mxl+1,nkap)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1,nkap)
 integer, intent(in) :: iout

!Output variables
 real(dp), intent(out) :: rphs(4),pshae(npsh,4,nkap),pshps(npsh,4,nkap)
 integer, intent(out) :: ifail

!Local variables
 integer :: ii,ll,l1,lc,nvkb
 integer :: ikap,kap,mkap,lout
 real(dp) :: epsh,rq,amesh
 character(len=11) :: rname

 ifail = 0
 pshae(:,:,:)=0.0d0
 pshps(:,:,:)=0.0d0
 amesh=rr(2)/rr(1)

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
! the core radius of l, or of the local potential beyond lmax
   if(ll <= lmax) then
     lc = l1
   else
     lc = lloc + 1
   end if
   if(rxpsh < 0.0d0) then
     if(ll <= lmax) then
       rq = rc(lc)*amesh**2
     else
       rq = rc(lc)
     end if
   else
     rq = rxpsh
   end if
   if(phsft_mch(rr,mmax,rq) == 0) then
     write(iout, '(/2a,f12.6,a,f12.6)') &
&       trim(rname)//': ERROR rxpsh beyond radial mesh', &
&       ', rxpsh=', rq, ', rr(mmax)=', rr(mmax)
     ifail = 1
     return
   end if
! check_data keeps rxpsh beyond every rc, so only a direct call gets here
   if(rxpsh >= 0.0d0 .and. rq < rc(lc)) then
     if(ll <= lmax) then
       write(iout, '(/2a,i0,2(a,f12.6))') &
&         trim(rname)//': WARNING rxpsh inside the core radius', &
&         ', l=', ll, ', rxpsh=', rq, ', rc=', rc(lc)
     else
       write(iout, '(/2a,i0,2(a,f12.6))') &
&         trim(rname)//': WARNING rxpsh inside ' // &
&         'the local potential core radius', &
&         ', l=', ll, ', rxpsh=', rq, ', rc=', rc(lc)
     end if
   end if

   rphs(l1)=rq
   mkap=kap_count(nkap,l1-1)
! loop on J = ll +/- 1/2 when irel=2
   do ikap=1,mkap
     kap=kap_of(ll,ikap)
     if(ll == lloc) then
       nvkb = 0
     else
       nvkb = nproj(l1)
     end if
     call phsft_pair(irel,c_light,ll,kap,nvkb,rr,vfull,zz,vp(1,lloc+1,ikap), &
&                    vlocpw,vkbpw(1:nvkb,l1,ikap),evkb(1:nvkb,l1,ikap),mmax,1, &
&                    [rq],npsh,epsh2,depsh,epa(1,l1,ikap),pshae(:,l1,ikap), &
&                    pshps(:,l1,ikap),iout,ifail)
     if(ifail /= 0) return

     lout=ll
     if(nkap==2 .and. ikap==1) lout=-ll
     write(iout,'(/a,i2)') 'log derivativve data for plotting, l=',ll
     write(iout,'(a,f6.2)') 'atan(r * ((d psi(r)/dr)/psi(r))), r=',rq
     write(iout,'(a/)') 'l, energy, all-electron, pseudopotential'
     do ii = 1, npsh
       epsh = epsh2 - depsh * real(ii - 1,dp)
       write(iout,'(a, i6, 3f12.6)') '! ',lout, epsh, pshae(ii,l1,ikap), &
&            pshps(ii,l1,ikap)
     end do
   end do !ikap
 end do !l1
 return
 end subroutine run_phsft

 end module run_phsft_m

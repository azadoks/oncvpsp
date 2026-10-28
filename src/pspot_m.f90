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

 module pspot_m

 implicit none
 private
 public :: pspot

 contains

 subroutine pspot(ipr,ll,rr,irc,rc,mmax,al,nbas,qroot,eig,uu,pswfopt_sb, &
&                 psx,chix, &
&           psopt,vae,vpsp,vkb,ekin_num,iout,ifail)

!! calculates optimized pseudopotential from coefficients of optimized
!! pseudo wave function and all-electron wave function
!!
!! @param ipr  1 for first projector, 2 for second which usually has a node
!! @param ll  angular momentum
!! @param rr  log radial mesh
!! @param irc  last mesh point at or inside rc
!! @param rc  core radius (bohr)
!! @param mmax  last point of rr
!! @param al  log of mesh spacing factor
!! @param nbas  number of sbfs
!! @param qroot  q's for set of sbf
!! @param eig  eigenvalue of all-electron wave function
!! @param uu  rr*all-electron wave function
!! @param pswfopt_sb  sbf coefficients for optimized pseudo wave function
!! @param vae  all-electron potential
!! @param psopt  rr * pseudo wave function on log grid
!! @param vpsp  pseudopotential
!! @param vkb  (epsilon-T)*phi component of chi projector
!! @param ekin_num  total radial kinetic energy: the Bessel sum's to rc and
!!                  the all-electron wave function's beyond, each by the
!!                  rule of logmesh_m
!! @param psx  the Bessel sum rr * pseudo wave function on the nodes up to
!!             pw_reach+16 past irc, 0 beyond: the inside piece's smooth
!!             extension past rc (pieces_m)
!! @param chix  (eig - T) psx on the same nodes: the projector's seed
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp
 use pieces_m, only: pw_reach
 use logmesh_m, only: lm_derivs, lm_integral
 use sbf_rc_der_m, only: sbf_rc_der
 implicit none

!Input variables
 integer, intent(in) :: ipr,ll,irc,mmax,nbas
 real(dp), intent(in) :: rr(mmax),qroot(nbas),uu(mmax),pswfopt_sb(nbas)
 real(dp), intent(in) :: vae(mmax)
 real(dp), intent(in) :: al,eig,rc

!Input/output variables
 real(dp), intent(inout) :: vkb(mmax)

!Output variables
 real(dp), intent(out) :: psopt(mmax)
 real(dp), intent(out) :: vpsp(mmax)
 real(dp), intent(out) :: ekin_num
 real(dp), intent(out) :: psx(mmax),chix(mmax)
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
 real(dp), allocatable :: work(:,:)
 real(dp) :: sbfder(5)
 real(dp) :: tt,ddt,ske,uder(2)
 integer :: ii,jj,ll1,nx

 ifail=0
 allocate(work(mmax,5))
 work(:,:)=0.0d0


 ll1=ll+1
 do ii=1,mmax
   work(ii,1)=uu(ii)
 end do

!7-point numerical first derivative of the all-electron wave function
!(of uu here, not uu/rr as for sbf matching), for its kinetic energy
!beyond rc
 do ii=7,mmax-6
   work(ii,2)=(-work(ii-3,1)+ 9.d0*work(ii-2,1)&
&   -45.d0*work(ii-1,1)+45.d0*work(ii+1,1)&
&   -9.d0*work(ii+2,1)+work(ii+3,1))&
&   /(60.d0*al*rr(ii))
 end do

!the inside piece and its seed on the nodes past rc too, as smooth
!extensions for the integrals of pieces_m
 nx=min(mmax,irc+pw_reach+16)
 psx(:)=0.0d0
 chix(:)=0.0d0
 do ii=1,nx
   tt=0.0d0
   ddt=0.0d0
   do jj=1,nbas
     call sbf_rc_der(ll,qroot(jj),rr(ii),sbfder)
     tt=tt+pswfopt_sb(jj)*rr(ii)*sbfder(1)
     ddt=ddt+pswfopt_sb(jj)*qroot(jj)**2*rr(ii)*sbfder(1)
   end do
   psx(ii)=tt
   chix(ii)=eig*tt-0.5d0*ddt
 end do

!fill [0, rc] with rr * analytic pswf and its kinetic operator: for
!u = r j_l(q r), l(l+1) u/r**2 - u'' = q**2 u (Bessel's equation), exactly;
!the two terms of the left side go as r**(l-1) and cancel to r**(l+1),
!which lost eps l(l+1)/r**2 near the origin (doc/numerics_audit.md, C-11)
 do ii=1,irc
   tt=0.0d0
   ddt=0.0d0
   do jj=1,nbas
     call sbf_rc_der(ll,qroot(jj),rr(ii),sbfder)
     tt=tt+pswfopt_sb(jj)*rr(ii)*sbfder(1)
     ddt=ddt+pswfopt_sb(jj)*qroot(jj)**2*rr(ii)*sbfder(1)
   end do
   work(ii,1)=tt
   work(ii,4)=0.5d0*ddt
 end do
 do ii=1,mmax
   psopt(ii)=work(ii,1)
 end do

!pseudopotential and node test
 do ii=1,irc
   if(abs(work(ii,1))>0.0d0) then
     if(ipr<=1 .and. work(ii,1)*work(ii+1,1)<0.0d0) then
       write(iout,'(a)') ' ERROR pspot:  first pseudo wave function has node, &
&         program will stop'
       write(iout,'(a)') ' ERROR pspot: try changing psp parameters for this l'
       ifail=1
       return
     end if
     vpsp(ii)=-work(ii,4)/work(ii,1) + eig
     vkb(ii)=-work(ii,4) + eig*work(ii,1)
   else
     vpsp(ii)=0.0d0
   end if
 end do

 do ii=irc+1,mmax
   vpsp(ii)=vae(ii)
   vkb(ii)=vae(ii)*work(ii,1)
 end do

!radial kinetic energy, split at rc itself (a sum over the mesh that
!switched functions at irc, inside the cell of rc, was good to 3e-7 Ha in
!60_Nd, O(h**3)): to rc, that of the Bessel sum, psx T psx = psx (eig psx -
!chix), which goes as r**(2l+2) at the origin and is smooth past rc on
!the nodes of psx; beyond, that of the all-electron function u, by parts,
!u(rc) u'(rc)/2 + int_rc^inf (u'**2 + l(l+1) u**2/r**2)/2 dr, which needs
!only its first derivative (as tail_kinetic does for the free tail). The
!boundary term from psx, whose u'(rc) is the matched one of ode_rc_der,
!as in the function E_resid transforms: that of the mesh differs from it
!by the solver's error, O(h**4), 2e-8 Ha in 60_Nd
 work(:,5)=psx(:)*(eig*psx(:)-chix(:))
 call lm_integral(rr,work(:,5),mmax,2*ll+2,0.0d0,rc,ekin_num)
 do ii=1,mmax
   work(ii,5)=0.5d0*(work(ii,2)**2+dble(ll*ll1)*(uu(ii)/rr(ii))**2)
 end do
 call lm_integral(rr,work(:,5),mmax,2*ll,rc,rr(mmax),ske)
 call lm_derivs(rr,psx,mmax,rc,1,uder)
 ekin_num=ekin_num+ske+0.5d0*uder(1)*uder(2)

 deallocate(work)
 return
 end subroutine pspot

 end module pspot_m

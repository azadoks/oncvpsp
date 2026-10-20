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
! Creates monotonic polynomial model core charge matching all-electron
! core charge and 4 derivatives at "crossover" radius.
! Polynomial is 8th-order with no linear term.

! Performs analysis and based on "hardness" criterion described in
! Teter, Phys. Rev. B 48, 5031 (1993) , Appendix, as

 subroutine modcore(icmod,rhops,rhotps,rhoc,rhoae,rhotae,rhomod, &
&                   fcfact,rcfact,irps,mmax,rr,nc,nv,la,zion,iexc)

!! @param icmod  (not used)
!! @param rhops  state-by-state pseudocharge density
!! @param rhotps  total pseudocharge density
!! @param rhoc  core-charge density
!! @param rhoae  state-by-state all-electron valence charge density
!! @param rhotae  total all-electron valence charge density
!! @param rhomod  model core density and 4 derivatives
!! @param fcfact  the model is matched where rhoc equals fcfact times the total
!!                pseudocharge
!! @param rcfact  (not used)
!! @param irps  rr index of maximum rc
!! @param mmax  dimension of log grid
!! @param rr  log radial grid
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param la  angular-momenta
!! @param zion  ion charge
!! @param iexc  exchange-correlation function to be used

 use constants_m, only: dp, mxstate
 implicit none

!Input variables
 integer, intent(in) :: icmod,nv,nc,iexc,irps,mmax
 integer, intent(in) :: la(mxstate)
 real(dp), intent(in) :: rhoae(mmax,nv),rhops(mmax,nv),rhotae(mmax)
 real(dp), intent(in) :: rhotps(mmax),rhoc(mmax),rr(mmax)
 real(dp), intent(in) :: zion,fcfact,rcfact

!Output variables
 real(dp), intent(out) :: rhomod(mmax,5)

!Local variables
!convergence criterion
 real(dp), parameter :: eps=1.0d-7
 real(dp) :: a0,al,et,yy,gg,a0min,a0max,dermax,psum,sf,eeel,eexc
 real(dp) :: d2mdiff,rmatch,rhocmatch
 real(dp) :: emin,emax,sls,ss,hh,rgg,xx
 real(dp) :: aco(5),polym(5,5),work(5,5),constm(5,5),xpow(9),fmatch(5)
 real(dp), allocatable :: vxcae(:),vxcpsp(:),vo(:),d2excae(:,:),d2excps(:,:)
 real(dp), allocatable :: dvxcae(:,:),dvxcps(:,:),vxct(:)
 integer :: ii,ierr,ircc,irmod,iter,jj,kk,ll,l1,mch
 integer :: ipvt(5),nodes(4)


 allocate(vxcae(mmax),vxcpsp(mmax),vo(mmax))
 allocate(dvxcae(mmax,nv),dvxcps(mmax,nv),vxct(mmax))
 allocate(d2excae(nv,nv),d2excps(nv,nv))

 d2excae(:,:)=0.0d0
 d2excps(:,:)=0.0d0
 ircc = 0

 do ii = mmax,1,-1
   if(rhoc(ii) > fcfact*rhotps(ii)) then
     ircc=ii
     exit
   end if
 end do

 if(ircc == 0) then
   write(6,'(/a)') 'rhomod: ERROR ircc (core-valence charge crossover) &
&        not found'
   error stop 1
 end if

! the 7-point derivative stencils below reach down to ircc-22
 if(ircc<23) then
   write(6,'(/2a,i0)') 'modcore: ERROR core-valence crossover too close', &
&        ' to the origin, ircc=',ircc
   error stop 1
 end if

!set limit for d2Exc calculation to radius beyond which rhomod=rhoc
 irmod=max(irps,ircc)

 write(6,'(/a/a)') 'Model core correction analysis',&
&                  '  based on d2Exc/dn_idn_j'


! get derivatives of all-electron xc energy

 call der2exc(rhotae,rhoc,rhoae,rr,d2excae,d2excps,d2mdiff, &
&                   zion,iexc,nc,nv,la,irmod,mmax)

 write(6,'(/a/)') 'd2excae - all-electron derivatives'
 do kk=1,nv
   write(6,'(1p,4d16.6)') (d2excae(kk,jj),jj=1,nv)
 end do


! set model charge to zero
 rhomod(:,:)=0.0d0

! compute d2excps with no core correction
 rhomod(:,:)=0.0d0

 call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
&                   zion,iexc,nc,nv,la,irmod,mmax)

 write(6,'(/a/)') 'd2excps - pseudofunction derivatives with no core correction'
 do kk=1,nv
   write(6,'(1p,4d16.6)') (d2excps(kk,jj),jj=1,nv)
 end do
 write(6,'(/a,1p,e16.6)') 'rms 2nd derivative error',d2mdiff

 al = 0.01d0 * log(rr(101) / rr(1))

! core charge density for Louie-Froyen-Cohen correction

 rhomod(:,1)= rhoc(:)

 xx=rr(ircc)
 fmatch(1)=rhoc(ircc)

! core charge derivatives
 do jj=2,5

!7-point numerical first derivatives applied successively
   do ii=ircc-25+3*jj,mmax-3
     rhomod(ii,jj)=(-rhomod(ii-3,jj-1)+ 9.d0*rhomod(ii-2,jj-1)&
&     -45.d0*rhomod(ii-1,jj-1)+45.d0*rhomod(ii+1,jj-1)&
&     -9.d0*rhomod(ii+2,jj-1)+rhomod(ii+3,jj-1))&
&     /(60.d0*al*rr(ii))
   end do
   fmatch(jj)=rhomod(ircc,jj)
 end do

! constants for polynomial matrix
 do jj=1,5
   constm(1,jj)=1.0d0
 end do
 do jj=1,5
   constm(2,jj)=jj+2
 end do
 do ii=3,5
   do jj=1,5
     constm(ii,jj)=constm(ii-1,jj)*(constm(2,jj)-ii+2)
   end do
 end do

! powers of xx=rr(ircc)

 xpow(1)=0.0d0
 xpow(2)=1.0d0
 do ii=3,9
   xpow(ii)=xx*xpow(ii-1)
 end do

! polynomial matrix
 do jj=1,5
   do ii=1,5
     polym(ii,jj)=constm(ii,jj)*xpow(jj-ii+5)
   end do
 end do


! begin iteration for monotonic polynomial match

 a0=fmatch(1)
 a0min=a0
 a0max=0.0d0

 do iter=1,50

   work(:,:)=polym(:,:)

! fill target value vector
   aco(1)=fmatch(1)-a0
   do jj=2,5
     aco(jj)=fmatch(jj)
   end do

! solve linear equations for coefficients

   call dgesv(5,1,work,5,ipvt,aco,5,kk)
   if(kk/= 0) then
     if(kk>0) write(6,'(a,i4)') &
&      'modcore:ERROR stop - singular polym matrix',kk
     if(kk<0) write(6,'(a,i4)') &
&      'modcore:ERROR stop - dgesv input error',kk
     error stop 1
   end if

! find maximum derivative
! note that xpow end up properly re-written for continued iteation
   dermax=-1.0d20
   do kk=1,ircc
     do ii=3,9
       xpow(ii)=rr(kk)*xpow(ii-1)
     end do
     ii=2
     psum=0.0d0
     do jj=1,5
       psum=psum+constm(ii,jj)*xpow(jj-ii+5)*aco(jj)
     end do
     dermax=max(dermax,psum)
   end do

! test maximum derivative and adjust a0
! interval-halving search after achieving monotonicity to get
! barely monotonic model
   if(dermax > 0.0d0) then
     a0min=a0
   else
     a0max=a0
   end if

! success when maximum derivative bracketed just above zero
   if(dermax>0.0d0 .and. dermax<eps*abs(fmatch(2))) exit

   if(a0max==0.0d0) then
     a0=2.5d0*a0
   else
     a0=0.5d0*(a0max+a0min)
   end if

 end do !iterr
 if(iter>50) write(6,'(/a/)') 'WARNING - modcore not conveged'


! fill in model and derivatives with fitted polynomial
 do kk=1,ircc-1
   do ii=3,9
     xpow(ii)=rr(kk)*xpow(ii-1)
   end do
   do ii=1,5
     if(ii == 1) then
       psum=a0
     else
       psum=0.0d0
     end if
     do jj=1,5
       psum=psum+constm(ii,jj)*xpow(jj-ii+5)*aco(jj)
     end do
     rhomod(kk,ii)=psum
   end do
 end do


!test model
 call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
&                   zion,iexc,nc,nv,la,irmod,mmax)

 write(6,'(/a/)') 'Polynomial model core charge'
 write(6,'(/a/)') 'd2excps - pseudofunction derivatives with core correction'

 do kk=1,nv
   write(6,'(1p,4d16.6)') (d2excps(kk,jj),jj=1,nv)
 end do
 write(6,'(/a,1p,e16.6)') 'rms 2nd derivative error',d2mdiff

 deallocate(vxcae,vxcpsp,vo)
 deallocate(dvxcae,dvxcps)
 deallocate(d2excae,d2excps)

 return
 end subroutine modcore

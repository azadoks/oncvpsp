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
! Creates a Teter-function model core charge, amplitude fcfact times the
! core charge and scale rcfact times the radius at the core-valence
! crossover.  The prefactors are input (icmod=3) or optimized to reduce the
! Exc 2nd-derivative errors (icmod=4).
! The crossover is the root of ln(rhoc) - ln(rhotps), and the Teter function
! and its derivatives are analytic everywhere, tapered to 0 by a C-infinity
! switch between x = 3 and 4 instead of cut at 3 (modcore_m,
! doc/smoothness_plan.md M-1, M-5, M-7).  For icmod=4 the coarse grid and
! Nelder-Mead choose the basin and Newton finds its stationary point (M-6).

! Performs analysis and based on "hardness" criterion described in
! Teter, Phys. Rev. B 48, 5031 (1993) , Appendix, as

 subroutine modcore3(icmod,rhops,rhotps,rhoc,rhoae,rhotae,rhomod, &
&                   fcfact,rcfact,fcgrid,rcgrid,mmax,rr,nc,nv,la,zion,iexc)

!! @param icmod  3 for specified fcfact and rcfact, 4 for coefficient
!!               optimization
!! @param rhops  state-by-state pseudocharge density
!! @param rhotps  total pseudocharge density
!! @param rhoc  core-charge density
!! @param rhoae  state-by-state all-electron valence charge density
!! @param rhotae  total all-electron valence charge density
!! @param rhomod  model core density and 4 derivatives
!! @param fcfact  prefactor for model amplitude (multiplies crossover value),
!!                set to 1 for icmod=4
!! @param rcfact  prefactor for model scale (multiplies crossover radius)
!! @param fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!! @param rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!! @param mmax  dimension of log grid
!! @param rr  log radial grid
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param la  angular-momenta
!! @param zion  ion charge
!! @param iexc  exchange-correlation function to be used

 use constants_m, only: dp, mxstate, nrhomod
 use logmesh_m, only: lm_derivs
 use modcore_m, only: core_crossover, teter_core_derivs
 implicit none

!Input variables
 integer, intent(in) :: icmod,nv,nc,iexc,mmax
 integer, intent(in) :: la(mxstate)
 real(dp), intent(in) :: rhoae(mmax,nv),rhops(mmax,nv),rhotae(mmax)
 real(dp), intent(in) :: rhotps(mmax),rhoc(mmax),rr(mmax)
 real(dp), intent(in) :: zion,rcfact
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)

!Input/output variables
 real(dp), intent(inout) :: fcfact

!Output variables
 real(dp), intent(out) :: rhomod(mmax,nrhomod)

!Local variables
 real(dp) :: d2mdiff,rmatch,rhocmatch
 real(dp) :: dd(1),gd(0:4)
 real(dp), allocatable :: d2excae(:,:),d2excps(:,:)
 integer :: ii,ierr,jj,kk

!2-dimensional Nelder-Mead variables
 real(dp), parameter :: alpha_nm=1.0d0
 real(dp), parameter :: gamma_nm=2.0d0
 real(dp), parameter :: rho_nm=-0.5d0
 real(dp), parameter :: sigma_nm=0.5d0
 real(dp) :: xx(2,3),ff(3),xt(2),ft,x0(2),xr(2),fr,xe(2),fe,xc(2),fc
 real(dp) :: d2ref,d2min
 real(dp), allocatable :: fta(:)
 integer :: nfc,nrc
 character(len=24) :: rowfmt


 allocate(d2excae(nv,nv),d2excps(nv,nv))

 d2excae(:,:)=0.0d0
 d2excps(:,:)=0.0d0

 write(6,'(/a/a)') 'Model core correction analysis',&
&                  '  based on d2Exc/dn_idn_j'

 call der2exc(rhotae,rhoc,rhoae,rr,d2excae,d2excps,d2mdiff, &
&                   zion,iexc,nc,nv,la,rr(mmax),mmax)

 write(6,'(/a/)') 'd2excae - all-electron derivatives'
 do kk=1,nv
   write(6,'(1p,4d16.6)') (d2excae(kk,jj),jj=1,nv)
 end do


! compute d2excps with no core correction
 rhomod(:,:)=0.0d0

 call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
&                   zion,iexc,nc,nv,la,rr(mmax),mmax)

 write(6,'(/a/)') 'd2excps - pseudofunction derivatives with no core correction'
 do kk=1,nv
   write(6,'(1p,4d16.6)') (d2excps(kk,jj),jj=1,nv)
 end do
 write(6,'(/a,1p,e16.6)') 'rms 2nd derivative error',d2mdiff

 d2ref=d2mdiff

! find valence pseudocharge - core charge crossover
 call core_crossover(rr,rhoc,rhotps,mmax,1.0d0,rmatch,ierr)
 if(ierr/=0) then
   write(6,'(/a)') 'modcore3: ERROR ircc (core-valence charge crossover) &
&        not found'
   error stop 1
 end if
 call lm_derivs(rr,rhoc,mmax,rmatch,0,dd)
 rhocmatch=dd(1)
 write(6,'(/a,2f10.4)') 'rmatch, rhocmatch',rmatch,rhocmatch

!option for optimization
 if(icmod==4) then

   fcfact=1.0d0

!Coarse-grid search for minimum

   write(6,'(/a)') 'Coarse scan for minimum error'
   write(6,'(a)') '  matrix elements: rms 2nd-derivative errors (mHa)'
   write(6,'(a)') '  column index : amplitude prefactor to rhocmatch'
   write(6,'(a)') '  row index : scale prefactor to rmatch'

   nfc=nint((fcgrid(2)-fcgrid(1))/fcgrid(3))+1
   nrc=nint((rcgrid(2)-rcgrid(1))/rcgrid(3))+1
   allocate(fta(nfc))

!row labels need more than one decimal if the rcfact grid is finer than 0.1
   if(abs(10.0d0*rcgrid(1)-nint(10.0d0*rcgrid(1)))<1.0d-8 .and. &
&     abs(10.0d0*rcgrid(3)-nint(10.0d0*rcgrid(3)))<1.0d-8) then
     rowfmt='(f5.1,f9.3,*(f7.3))'
   else
     rowfmt='(f6.3,f8.3,*(f7.3))'
   end if

   d2min=10.0d0
   write(6,'(/7x,*(f7.3))') (fcgrid(1)+fcgrid(3)*(jj-1), jj=1,nfc)
   write(6,'(a)') ''
   do  kk=1,nrc
     xt(2)=(rcgrid(1)+rcgrid(3)*(kk-1))*rmatch
     do jj=1,nfc
       xt(1)=(fcgrid(1)+fcgrid(3)*(jj-1))*rhocmatch

       call trial(xt,d2mdiff)
       fta(jj)=1.0d3*d2mdiff

       if(d2mdiff<d2min) then
         xx(:,1)=xt(:)
         d2min=d2mdiff
       end if

     end do
     write(6,rowfmt) rcgrid(1)+rcgrid(3)*(kk-1),(fta(jj),jj=1,nfc)
   end do

!Initial Nelder-Mead simplex from coarse-search minimum, scaled to the
!grid steps

   xx(1,2)=xx(1,1)+0.5d0*fcgrid(3)*rhocmatch
   xx(2,2)=xx(2,1)
   xx(1,3)=xx(1,1)
   xx(2,3)=xx(2,1)+0.5d0*rcgrid(3)*rmatch
   xx(1,1)=xx(1,1)-0.25d0*fcgrid(3)*rhocmatch
   xx(2,1)=xx(2,1)-0.25d0*rcgrid(3)*rmatch

!Fill function values for initial simplex

   do kk=1,3
     call trial(xx(:,kk),ff(kk))

!  write(6,'(i4,a,2f10.4,1p,e14.4)') kk,'  xt',xt(1),xt(2),ff(kk)
   end do

!Nelder-Mead iteration loop
   write(6,'(/a)') 'Nelder-Mead iteration'
   do kk=1,101

!(1) Order (dumb bubble sort)
     do jj=1,3
       if(ff(3)<ff(2)) then
         ft=ff(2) ; xt(:)=xx(:,2)
         ff(2)=ff(3) ; xx(:,2)=xx(:,3)
         ff(3)=ft ; xx(:,3)=xt(:)
       end if
       if(ff(2)<ff(1)) then
         ft=ff(1) ; xt(:)=xx(:,1)
         ff(1)=ff(2) ; xx(:,1)=xx(:,2)
         ff(2)=ft ; xx(:,2)=xt(:)
       end if
     end do

!stopping criterion
     if(ff(3)-ff(1)<1.0d-4*d2ref) then
       write(6,'(a,i4,a)') ' converged in',kk-1,' steps'
       exit
     else if(kk==101) then
       write(6,'(a,i4,a)') ' WARNING: not fully converged in 100 steps'
       exit
     end if

!(2) Centroid
     x0(:)=0.5d0*(xx(:,1)+xx(:,2))

!(3) Reflection
     xr(:)=x0(:)+alpha_nm*(x0(:)-xx(:,3))

     call trial(xr,fr)

     if(ff(1)<= fr .and. fr<ff(2)) then
       ff(3)=fr ; xx(:,3)=xr(:)
       cycle !kk loop
     end if

!(4) Expansion
     if(fr<ff(1)) then
       xe(:)=x0(:)+gamma_nm*(x0(:)-xx(:,3))

       call trial(xe,fe)

       if(fe<fr) then
         ff(3)=fe ; xx(:,3)=xe(:)
         cycle !kk
       else
         ff(3)=fr ; xx(:,3)=xr(:)
         cycle !kk
       end if
     end if

!(5) Contraction
     xc(:)=x0(:)+rho_nm*(x0(:)-xx(:,3))

     call trial(xc,fc)
!  write(6,'(i4,a,2f10.4,1p,e14.4)') kk,'  xc',xc(1),xc(2),fc
     if(fc<ff(3)) then
       ff(3)=fc ; xx(:,3)=xc(:)
       cycle !kk
     end if

!(6) Reduction
     do jj=2,3
       xx(:,jj)=xx(:,1)+sigma_nm*(xx(:,jj)-xx(:,1))

       call trial(xx(:,jj),ff(jj))
!  write(6,'(i4,a,2f10.4,1p,e14.4)') kk,' xrd',xx(1,jj),xx(2,jj),ff(jj)
     end do !jj

   end do !kk

!Newton on the gradient from the best vertex, so the result is the
!stationary point itself rather than wherever the simplex stopped
!(doc/smoothness_plan.md M-6)
   call newton_refine(xx(:,1))
   write(6,'(a)') 'amplitude prefactor, scale prefactor'
   write(6,'(2f10.4)') xx(1,1)/rhocmatch,xx(2,1)/rmatch

   write(6,'(/a)') 'Optimized Teter model core charge'

!option for specifying prefactors as input
 else if(icmod==3) then
   xx(1,1)=fcfact*rhocmatch
   xx(2,1)=rcfact*rmatch

   write(6,'(/a/)') &
&        'Teter function model core charge with specified parameters'
 end if


! the Teter function and its derivatives, analytic everywhere (tapered)
 do ii=1,mmax
   call teter_core_derivs(rr(ii)/xx(2,1),4,gd)
   do kk=0,4
     rhomod(ii,kk+1)=xx(1,1)*gd(kk)/xx(2,1)**kk
   end do
 end do

 call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,d2mdiff, &
&                   zion,iexc,nc,nv,la,rr(mmax),mmax)

 write(6,'(/a/)') 'd2excps - pseudofunction derivatives with core correction'
 do kk=1,nv
   write(6,'(1p,4d16.6)') (d2excps(kk,jj),jj=1,nv)
 end do
 write(6,'(/a,1p,e16.6)') 'rms 2nd derivative error',d2mdiff

 deallocate(d2excae,d2excps)

 return

 contains

 subroutine trial(xm,fm)

!! the rms 2nd-derivative error of a trial Teter model core
!!
!! @param xm  its amplitude and scale (bohr)
!! @param fm  the rms 2nd-derivative error

!Input variables
 real(dp), intent(in) :: xm(2)

!Output variables
 real(dp), intent(out) :: fm

!Local variables
 integer :: im

 do im=1,mmax
   call teter_core_derivs(rr(im)/xm(2),0,gd)
   rhomod(im,1)=xm(1)*gd(0)
 end do
 call der2exc(rhotps,rhomod(1,1),rhops,rr,d2excps,d2excae,fm, &
&             zion,iexc,nc,nv,la,rr(mmax),mmax)
 end subroutine trial

 subroutine newton_refine(xm)

!! Newton's method on the gradient of the rms error in (ln amplitude,
!! ln scale), held inside a trust radius.  The valley of the error is
!! narrow and curved (Hessian condition ~600 for 40_Zr), and der2exc is
!! noisy at ~1e-17 on f ~ 1e-4, so the gradient is a central difference
!! along each eigenvector of the Hessian, with a step hk ~ 1/sqrt(lambda_k):
!! hj across the valley, longer along it.  The Hessian, for the directions
!! only, takes the step hb, as its noise turns them and the turn moves the
!! root through the valley's bend.  Axis differences at hj stop at ~1e-11
!! in ln(amplitude); these at ~5e-13, the noise over the soft curvature.
!! The steps and directions are smooth in the point, so the root is a
!! smooth function of the inputs.
!!
!! @param xm  amplitude and scale (bohr), in: start, out: stationary point

!Input/output variables
 real(dp), intent(inout) :: xm(2)

!Local variables
 real(dp), parameter :: hj=1.0d-3,hb=5.0d-3,hcap=5.0d-2,trmax=0.1d0
 real(dp), parameter :: tol=3.0d-12
 real(dp) :: uu(2),ut(2),fs(-1:1,-1:1),hm(2,2),bb(2,2),lam(2),hk(2)
 real(dp) :: f0,fp,fm,gk,ck,du(2),dnrm,tr,dt
 integer :: it,ia,ib,kk

 uu(:)=log(xm(:))
 do it=1,30
!Hessian at the step hb, for the directions and step lengths
   do ib=-1,1
     do ia=-1,1
       ut(1)=uu(1)+ia*hb
       ut(2)=uu(2)+ib*hb
       call trial(exp(ut),fs(ia,ib))
     end do
   end do
   f0=fs(0,0)
   hm(1,1)=(fs(1,0)-2.0d0*f0+fs(-1,0))/hb**2
   hm(2,2)=(fs(0,1)-2.0d0*f0+fs(0,-1))/hb**2
   hm(1,2)=(fs(1,1)-fs(1,-1)-fs(-1,1)+fs(-1,-1))/(4.0d0*hb**2)
   hm(2,1)=hm(1,2)
   tr=0.5d0*(hm(1,1)+hm(2,2))
   dt=sqrt((0.5d0*(hm(1,1)-hm(2,2)))**2+hm(1,2)**2)
   lam(1)=tr+dt
   lam(2)=tr-dt
   if(abs(hm(1,1)-lam(2))>=abs(hm(2,2)-lam(2))) then
     bb(:,2)=[hm(1,2),lam(2)-hm(1,1)]
   else
     bb(:,2)=[lam(2)-hm(2,2),hm(1,2)]
   end if
   if(dt==0.0d0) bb(:,2)=[0.0d0,1.0d0]
   bb(:,2)=bb(:,2)/sqrt(bb(1,2)**2+bb(2,2)**2)
   bb(:,1)=[bb(2,2),-bb(1,2)]

!steps of equal rise in f, at most hcap
   do kk=1,2
     if(lam(kk)>0.0d0 .and. lam(1)>0.0d0) then
       hk(kk)=min(hcap,hj*sqrt(lam(1)/lam(kk)))
     else
       hk(kk)=hj
     end if
   end do

!Newton step along each direction
   du(:)=0.0d0
   do kk=1,2
     ut(:)=uu(:)+hk(kk)*bb(:,kk)
     call trial(exp(ut),fp)
     ut(:)=uu(:)-hk(kk)*bb(:,kk)
     call trial(exp(ut),fm)
     gk=(fp-fm)/(2.0d0*hk(kk))
     ck=(fp-2.0d0*f0+fm)/hk(kk)**2
     if(ck>0.0d0) then
       du(:)=du(:)-(gk/ck)*bb(:,kk)
     else
!not convex along this direction: descend to the trust radius
       du(:)=du(:)-sign(trmax,gk)*bb(:,kk)
     end if
   end do
   dnrm=sqrt(du(1)**2+du(2)**2)
   if(dnrm>trmax) du(:)=du(:)*trmax/dnrm
   uu(:)=uu(:)+du(:)
   if(dnrm<tol) then
     write(6,'(a,i4,a)') ' Newton refinement converged in',it,' steps'
     xm(:)=exp(uu(:))
     return
   end if
 end do
 write(6,'(a,1p,e10.2)') ' WARNING: Newton refinement not converged, &
&      last step',dnrm
 xm(:)=exp(uu(:))
 end subroutine newton_refine

 end subroutine modcore3

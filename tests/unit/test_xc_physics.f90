! Unit tests of the physics of the exchange-correlation functionals: exact
! properties of the exchange and correlation energies, each checked for
! the built-in functionals over a range of densities (step 3 of
! doc/physics_tests_plan.md). rho is 4 pi n, as exc_dispatch takes it, and
! r_s = (3/rho)**(1/3). The LDAs' exchange is Slater's,
! eps_x = -(3/4)(3 n/pi)**(1/3), so their correlation is eps_xc - eps_x;
! PBE's two are apart in pbe_x and pbe_c (exc_off_pbe_m).
!   gmb        (iexc 3, 6, PBE's PW92) the high-density limit of the
!              correlation, eps_c = A ln r_s + B + O(r_s ln r_s), with
!              Gell-Mann and Brueckner's A = (1 - ln 2)/pi**2: the slope
!              in ln r_s between r_s 1e-5 and 1e-4, to 1e-3 relative (PZ's
!              A is rounded to 0.0311, 3e-4 off; the next order, 1e-4).
!              Wigner (1) and Hedin-Lundqvist (2) are not built for it;
!              for them (and all) eps_c is bounded at r_s 1e-5, |eps_c| <
!              1 Ha against |eps_x| = 5e4: the exchange is Slater's
!   derivative (iexc 1, 2, 3, 6, 7) the potential is the derivative of the
!              energy, v = d(rho eps)/d rho, by 4th-order differences in
!              rho (steps 1e-3 rho, error ~1e-12) at 60 r_s from 0.1 to
!              20 (PZ not within 1e-3 of r_s = 1, where it steps); 1e-10
!              of |v|
!   functional (iexc 1, 2, 3, 4, 6, 7) on a model atom, the potential is
!              the functional derivative of E_xc = int eps rho r**2 dr:
!              dE[rho + t drho]/dt = int v drho r**2 dr, the difference
!              4th order in t (steps 1e-3), the integrals the same
!              quadrature, on the program's mesh (amesh 1.006); for the
!              GGA its divergence form leaves the mesh's error, which falls
!              by 17 from amesh 1.012 to 1.006, to 5e-10: 2e-9 (the LDAs
!              ~1e-10, round-off)
!   eps_c      (iexc 1, 2, 3, 6, PW92) at 200 r_s from 1e-3 to 1e3:
!              eps_c <= 0; eps_c rises with r_s (the correlation weakens as
!              the density falls); r_s eps_c falls with r_s, which is the
!              scaling inequality E_c[n_g] > g E_c[n] for g > 1, n_g =
!              g**3 n(g r), of Levy and Perdew, for a uniform density
!   pz_branch  (iexc 3) eps_xc and v_xc are continuous at r_s = 1, where PZ
!              joins its two fits, only to the rounding of the fits'
!              coefficients in their 4th decimal: 1e-4 Ha
!   lieb_oxford (iexc 1, 2, 3, 4, 6, 7) E_xc >= -1.679 int n**(4/3) d3r, on
!              three model atoms; and PBE's exchange enhancement, the local
!              form of the bound, 1 <= F_x(s) <= 1 + kappa = 1.804, rising,
!              for s from 0 to 100
!   scaling_x  (PBE) uniform scaling, E_x[n_g] = g E_x[n]: pointwise,
!              eps_x(g**3 n, g**4 n') = g eps_x(n, n') and the same for
!              its potential, for g 0.5, 2, 10, over n and s; 1e-13
!              relative (round-off)
!   uniform    (PBE) at zero gradient PBE is the LDA it is built on:
!              eps_x is Slater's, H is 0 and eps_c is PW92's G(r_s), here
!              from its formula (with PBE's 0.0310907 for A); 1e-13
!              relative
!   linear_response (PBE) the second-order gradient terms of exchange,
!              mu s**2, and of correlation, beta t**2, cancel (mu =
!              beta pi**2/3), so PBE keeps the uniform gas' linear
!              response: eps_xc(n, n') - eps_xc(n, 0) over
!              eps_x(n, n') - eps_x(n, 0) is of the next order, C s**2,
!              from H = beta t**2 - (beta**2/(2 gamma)) t**4 and F_x =
!              1 + mu s**2 - (mu**2/kappa) s**4 with t**2 = (pi k_F/4) s**2:
!              C = beta**2 pi**3 k_F/(24 gamma mu) - mu/kappa; within 5% of
!              the larger of the two terms, at s 1e-3 and 1e-2 (the order
!              after, ~t**2 <= 0.02; at s 1e-3 the round-off, ~1e-2 at
!              r_s 10), r_s 0.1 to 10
!   relativistic (iexc 7 against 6) the relativistic exchange is weaker,
!              eps_xc(7) >= eps_xc(6), and the difference falls by 4 when
!              c doubles (O((k_F/c)**2)), within 4 (2 k_F/c)**2 of 4, at
!              c 1 and 2 times the physical value, r_s 0.1 to 20
! Stops with a non-zero code if any check failed.
 program test_xc_physics
 use constants_m, only: dp, pi
 use, intrinsic :: iso_fortran_env, only: stderr => error_unit
 use unit_check_m, only: check, check_done
 use xc_dispatch_m, only: exc_dispatch
 use exc_off_pbe_m, only: pbe_x, pbe_c
 use quad_m, only: quad_gregory_sum
 use relativity_m, only: c_light => c_light_default
 implicit none
 integer, parameter :: iout=6

!Local variables
 integer, parameter :: mmax=2176
 real(dp), parameter :: rr1=1.0e-4_dp
!> the LDAs and all the built-in functionals
 integer, parameter :: ildas(5)=[1,2,3,6,7], iall(6)=[1,2,3,4,6,7]
!> Gell-Mann and Brueckner's coefficient, and PBE's kappa and mu
 real(dp), parameter :: agmb=(1.0_dp-log(2.0_dp))/pi**2
 real(dp), parameter :: uk=0.8040_dp, um=0.2195149727645171_dp
!> PBE's beta and gamma
 real(dp), parameter :: bet=0.06672455060314922_dp
 real(dp), parameter :: gam=0.03109069086965489503494086371273_dp
 real(dp) :: al,rr(mmax)
 integer :: ii

 al=log(1.006_dp)
 do ii=1,mmax
   rr(ii)=rr1*exp(al*(ii-1))
 end do

 call check_gmb()
 call check_derivative()
 call check_functional()
 call check_eps_c()
 call check_pz_branch()
 call check_lieb_oxford()
 call check_scaling_x()
 call check_uniform()
 call check_linear_response()
 call check_relativistic()

 call check_done('xc_physics')

 contains

 elemental real(dp) function rho_rs(rs)
!! rho = 4 pi n of the uniform gas of radius r_s
 real(dp), intent(in) :: rs
 rho_rs=3.0_dp/rs**3
 end function rho_rs

 elemental real(dp) function eps_slater(rho)
!! Slater's exchange energy per electron of rho = 4 pi n
 real(dp), intent(in) :: rho
 eps_slater=-0.75_dp*(3.0_dp*rho/(4.0_dp*pi**2))**(1.0_dp/3.0_dp)
 end function eps_slater

 function rs_grid(n,rs1,rs2) result(rs)
!! n values of r_s from rs1 to rs2, evenly in ln r_s
 integer, intent(in) :: n
 real(dp), intent(in) :: rs1,rs2
 real(dp) :: rs(n)
 integer :: ii
 do ii=1,n
   rs(ii)=rs1*(rs2/rs1)**(real(ii-1,dp)/(n-1))
 end do
 end function rs_grid

 subroutine lda(iexc,rho,vxc,exc,cl)
!! exc_dispatch of the LDA iexc at the points rho, at the speed of light cl
!! (by default the physical one)
 integer, intent(in) :: iexc
 real(dp), intent(in) :: rho(:)
 real(dp), intent(out) :: vxc(:),exc(:)
 real(dp), intent(in), optional :: cl
 real(dp) :: cc
 integer :: ifail
 cc=c_light
 if(present(cl)) cc=cl
 vxc(:)=0.0_dp
 exc(:)=0.0_dp
 call exc_dispatch(iexc,cc,al,rho,vxc,exc,rr,size(rho),'test_xc_physics', &
&                  ifail)
 if(ifail/=0) error stop 1
 end subroutine lda

 real(dp) function eps_c_lda(iexc,rs)
!! the correlation energy per electron of the uniform gas: eps_xc less
!! Slater's exchange for an LDA, PBE's PW92 for iexc 4
 integer, intent(in) :: iexc
 real(dp), intent(in) :: rs
 real(dp) :: vv(1),ee(1),hh,dv,ds
 if(iexc==4) then
   dv=0.0_dp
   ds=0.0_dp
   call pbe_c(rho_rs(rs)/(4.0_dp*pi),0.0_dp,ee(1),hh,dv,ds)
   eps_c_lda=ee(1)+hh
 else
   call lda(iexc,[rho_rs(rs)],vv,ee)
   eps_c_lda=ee(1)-eps_slater(rho_rs(rs))
 end if
 end function eps_c_lda

 function tag(what,iexc) result(lab)
!! a check's label: what and iexc
 character(len=*), intent(in) :: what
 integer, intent(in) :: iexc
 character(len=64) :: lab
 write(lab,'(2a,i0)') what,' iexc ',iexc
 end function tag

 subroutine check_gmb()
!! gmb, and the bounded correlation of every LDA at high density
 integer, parameter :: igmb(3)=[3,6,4]
 real(dp) :: slope
 integer :: ii
 do ii=1,size(igmb)
   slope=(eps_c_lda(igmb(ii),1.0e-4_dp)-eps_c_lda(igmb(ii),1.0e-5_dp)) &
&        /log(10.0_dp)
   call check(tag('gmb',igmb(ii)),abs(slope/agmb-1.0_dp),1.0e-3_dp)
 end do
 do ii=1,4
   call check(tag('exchange limit |eps_c(1e-5)|',ildas(ii)), &
&             abs(eps_c_lda(ildas(ii),1.0e-5_dp)),1.0_dp)
 end do
 end subroutine check_gmb

 subroutine check_derivative()
!! derivative: v = d(rho eps)/d rho for the LDAs
 integer, parameter :: nrs=60
 real(dp), parameter :: hh=1.0e-3_dp
 real(dp) :: rs(nrs),rho(nrs),vv(nrs),ee(nrs),ff(nrs,-2:2),vt(nrs),err
 integer :: ix,kk,ii
 rs(:)=rs_grid(nrs,0.1_dp,20.0_dp)
 rho(:)=rho_rs(rs)
 do ix=1,size(ildas)
   call lda(ildas(ix),rho,vv,ee)
   do kk=-2,2
     call lda(ildas(ix),rho*(1.0_dp+kk*hh),vt,ee)
     ff(:,kk)=rho*(1.0_dp+kk*hh)*ee
   end do
   err=0.0_dp
   do ii=1,nrs
     if(ildas(ix)==3 .and. abs(rs(ii)-1.0_dp)<1.0e-3_dp) cycle
     err=max(err,abs((ff(ii,-2)-8.0_dp*ff(ii,-1)+8.0_dp*ff(ii,1)-ff(ii,2)) &
&            /(12.0_dp*hh*rho(ii))-vv(ii))/abs(vv(ii)))
   end do
   call check(tag('derivative',ildas(ix)),err,1.0e-10_dp)
 end do
 end subroutine check_derivative

 function shell(zeta,nr) result(sh)
!! the density 4 pi n of a hydrogen-like s shell of exponent zeta with
!! the power r**nr of its radial function, normalized to 1 electron
 real(dp), intent(in) :: zeta
 integer, intent(in) :: nr
 real(dp) :: sh(mmax)
 sh(:)=rr**(2*nr)*exp(-2.0_dp*zeta*rr) &
&      /(gamma(real(2*nr+3,dp))/(2.0_dp*zeta)**(2*nr+3))
 end function shell

 real(dp) function integral(ff)
!! int ff r**2 dr over the mesh (ff falls fast enough at both ends)
 real(dp), intent(in) :: ff(mmax)
 integral=al*quad_gregory_sum(ff*rr**3,mmax)
 end function integral

 subroutine check_functional()
!! functional: dE_xc/dt = int v drho r**2 dr on the model atom of test_fxc
 real(dp), parameter :: tt=1.0e-3_dp
 real(dp) :: rho(mmax),drho(mmax),vv(mmax),ee(mmax),es(-2:2),dedt,vint
 integer :: ix,kk,ifail
 rho(:)=2.0_dp*shell(14.0_dp,0)+8.0_dp*shell(5.0_dp,1) &
&      +2.0_dp*shell(1.5_dp,2)
 drho(:)=shell(1.5_dp,2)
 do ix=1,size(iall)
   do kk=-2,2
     call exc_dispatch(iall(ix),c_light,al,rho+kk*tt*drho,vv,ee,rr,mmax, &
&                      'test_xc_physics',ifail)
     es(kk)=integral(ee*(rho+kk*tt*drho))
   end do
   dedt=(es(-2)-8.0_dp*es(-1)+8.0_dp*es(1)-es(2))/(12.0_dp*tt)
   call exc_dispatch(iall(ix),c_light,al,rho,vv,ee,rr,mmax, &
&                    'test_xc_physics',ifail)
   vint=integral(vv*drho)
   call check(tag('functional',iall(ix)),abs(dedt/vint-1.0_dp),2.0e-9_dp)
 end do
 end subroutine check_functional

 subroutine check_eps_c()
!! eps_c: sign, monotonicity, and the scaling inequality
 integer, parameter :: nrs=200
 integer, parameter :: ic(5)=[1,2,3,6,4]
 real(dp) :: rs(nrs),ec(nrs),sgn,rise,scal
 integer :: ix,ii
 rs(:)=rs_grid(nrs,1.0e-3_dp,1.0e3_dp)
 do ix=1,size(ic)
   do ii=1,nrs
     ec(ii)=eps_c_lda(ic(ix),rs(ii))
   end do
! the largest violations, 0 when the property holds
   sgn=max(0.0_dp,maxval(ec))
   rise=max(0.0_dp,maxval(ec(1:nrs-1)-ec(2:nrs)))
   scal=max(0.0_dp,maxval(rs(2:nrs)*ec(2:nrs)-rs(1:nrs-1)*ec(1:nrs-1)))
   call check(tag('eps_c <= 0',ic(ix)),sgn,0.0_dp)
   call check(tag('eps_c rises with r_s',ic(ix)),rise,0.0_dp)
   call check(tag('r_s eps_c falls with r_s',ic(ix)),scal,0.0_dp)
 end do
 end subroutine check_eps_c

 subroutine check_pz_branch()
!! pz_branch: PZ's steps at r_s = 1
 real(dp) :: rho(2),vv(2),ee(2)
 rho(:)=rho_rs([1.0_dp-1.0e-12_dp,1.0_dp+1.0e-12_dp])
 call lda(3,rho,vv,ee)
 call check('pz_branch eps_xc step',abs(ee(2)-ee(1)),1.0e-4_dp)
 call check('pz_branch v_xc step',abs(vv(2)-vv(1)),1.0e-4_dp)
 end subroutine check_pz_branch

 subroutine check_lieb_oxford()
!! lieb_oxford: the integrated bound on model atoms, and PBE's F_x
 integer, parameter :: ns=1001
 real(dp) :: rho(mmax),vv(mmax),ee(mmax),ratio,fx(ns),ex,vx,dfds,dd
 integer :: im,ix,ifail,ii
 do im=1,3
   select case(im)
   case(1)
     rho(:)=2.0_dp*shell(14.0_dp,0)+8.0_dp*shell(5.0_dp,1) &
&          +2.0_dp*shell(1.5_dp,2)
   case(2)
     rho(:)=shell(1.0_dp,0)
   case default
     rho(:)=2.0_dp*shell(1.7_dp,0)
   end select
   do ix=1,size(iall)
     call exc_dispatch(iall(ix),c_light,al,rho,vv,ee,rr,mmax, &
&                      'test_xc_physics',ifail)
! E_xc over -int n**(4/3) d3r, at most 1.679
     ratio=integral(ee*rho)/(-integral((rho/(4.0_dp*pi))**(1.0_dp/3.0_dp)*rho))
     call check(tag('lieb_oxford',iall(ix)),ratio,1.679_dp)
   end do
 end do
! F_x(s) at n = 1 and s = 0..100, from pbe_x's eps_x over Slater's
 dd=1.0_dp
 do ii=1,ns
   call pbe_x(dd,grad_s(dd,0.1_dp*(ii-1)),ex,vx,dfds)
   fx(ii)=ex/eps_slater(4.0_dp*pi*dd)
 end do
 call check('lieb_oxford F_x(0) = 1',abs(fx(1)-1.0_dp),1.0e-15_dp)
 call check('lieb_oxford F_x <= 1 + kappa',max(0.0_dp,maxval(fx)-1.0_dp-uk), &
&           0.0_dp)
 call check('lieb_oxford F_x rises',max(0.0_dp,maxval(fx(1:ns-1)-fx(2:ns))), &
&           0.0_dp)
 end subroutine check_lieb_oxford

 real(dp) function grad_s(dd,ss)
!! |n'| of the density dd at the reduced gradient s = |n'|/(2 k_F n)
 real(dp), intent(in) :: dd,ss
 grad_s=ss*2.0_dp*(3.0_dp*pi**2*dd)**(1.0_dp/3.0_dp)*dd
 end function grad_s

 subroutine check_scaling_x()
!! scaling_x: PBE's exchange under n -> g**3 n(g r)
 real(dp), parameter :: gs(3)=[0.5_dp,2.0_dp,10.0_dp]
 real(dp), parameter :: ss(4)=[0.0_dp,0.3_dp,1.0_dp,10.0_dp]
 real(dp) :: dd,gg,ex,vx,dfds,exg,vxg,err
 integer :: ig,id,is
 err=0.0_dp
 do ig=1,size(gs)
   do id=-4,2
     dd=10.0_dp**id
     do is=1,size(ss)
       gg=grad_s(dd,ss(is))
       call pbe_x(dd,gg,ex,vx,dfds)
       call pbe_x(gs(ig)**3*dd,gs(ig)**4*gg,exg,vxg,dfds)
       err=max(err,abs(exg/(gs(ig)*ex)-1.0_dp),abs(vxg/(gs(ig)*vx)-1.0_dp))
     end do
   end do
 end do
 call check('scaling_x',err,1.0e-13_dp)
 end subroutine check_scaling_x

 real(dp) function pw92(rs)
!! Perdew and Wang's (1992) correlation of the unpolarized uniform gas,
!! G(r_s) = -2 A (1 + a1 r_s) ln(1 + 1/(2 A (b1 r_s**(1/2) + b2 r_s +
!! b3 r_s**(3/2) + b4 r_s**2))), with PBE's digits of A
 real(dp), intent(in) :: rs
 real(dp), parameter :: aa=0.0310907_dp,a1=0.21370_dp
 real(dp), parameter :: b1=7.5957_dp,b2=3.5876_dp,b3=1.6382_dp,b4=0.49294_dp
 pw92=-2.0_dp*aa*(1.0_dp+a1*rs)*log(1.0_dp+1.0_dp/(2.0_dp*aa &
&     *(b1*sqrt(rs)+b2*rs+b3*rs**1.5_dp+b4*rs**2)))
 end function pw92

 subroutine check_uniform()
!! uniform: PBE at zero gradient
 integer, parameter :: nrs=60
 real(dp) :: rs(nrs),dd,ex,vx,ds,ec,hh,errx,errc,errh
 integer :: ii
 rs(:)=rs_grid(nrs,0.01_dp,100.0_dp)
 errx=0.0_dp
 errc=0.0_dp
 errh=0.0_dp
 do ii=1,nrs
   dd=rho_rs(rs(ii))/(4.0_dp*pi)
   call pbe_x(dd,0.0_dp,ex,vx,ds)
   vx=0.0_dp
   ds=0.0_dp
   call pbe_c(dd,0.0_dp,ec,hh,vx,ds)
   errx=max(errx,abs(ex/eps_slater(4.0_dp*pi*dd)-1.0_dp))
   errc=max(errc,abs(ec/pw92(rs(ii))-1.0_dp))
   errh=max(errh,abs(hh))
 end do
 call check('uniform eps_x is Slater',errx,1.0e-13_dp)
 call check('uniform eps_c is PW92',errc,1.0e-13_dp)
 call check('uniform H = 0',errh,0.0_dp)
 end subroutine check_uniform

 subroutine check_linear_response()
!! linear_response: PBE's gradient terms cancel to second order
 integer, parameter :: nrs=9
 real(dp), parameter :: ss(2)=[1.0e-3_dp,1.0e-2_dp]
 real(dp) :: rs(nrs),dd,gg,ex0,ex,ec0,ec,h0,hh,vv,ds,dx,dxc,kf,t1,t2
 integer :: ii,is
 character(len=64) :: lab
 rs(:)=rs_grid(nrs,0.1_dp,10.0_dp)
 do is=1,size(ss)
   do ii=1,nrs
     dd=rho_rs(rs(ii))/(4.0_dp*pi)
     gg=grad_s(dd,ss(is))
     call pbe_x(dd,0.0_dp,ex0,vv,ds)
     call pbe_x(dd,gg,ex,vv,ds)
     call pbe_c(dd,0.0_dp,ec0,h0,vv,ds)
     call pbe_c(dd,gg,ec,hh,vv,ds)
     dx=ex-ex0
     dxc=(ex+ec+hh)-(ex0+ec0+h0)
     write(lab,'(a,es7.1,a,f0.3)') 'linear_response s ',ss(is),' r_s ', &
&      rs(ii)
     kf=(9.0_dp*pi/4.0_dp)**(1.0_dp/3.0_dp)/rs(ii)
     t1=bet**2*pi**3*kf/(24.0_dp*gam*um)
     t2=um/uk
     call check(lab,abs(dxc/dx/ss(is)**2-(t1-t2))/max(t1,t2),0.05_dp)
   end do
 end do
 end subroutine check_linear_response

 subroutine check_relativistic()
!! relativistic: iexc 7 against 6, at c and 2c
 integer, parameter :: nrs=60
 real(dp) :: rs(nrs),rho(nrs),vv(nrs),e6(nrs),e7(nrs),e72(nrs),kf(nrs)
 real(dp) :: weak,err
 integer :: ii
 rs(:)=rs_grid(nrs,0.1_dp,20.0_dp)
 rho(:)=rho_rs(rs)
 kf(:)=(9.0_dp*pi/4.0_dp)**(1.0_dp/3.0_dp)/rs
 call lda(6,rho,vv,e6)
 call lda(7,rho,vv,e7)
 call lda(7,rho,vv,e72,2.0_dp*c_light)
 weak=max(0.0_dp,maxval(e6-e7))
 call check('relativistic exchange is weaker',weak,0.0_dp)
 err=0.0_dp
 do ii=1,nrs
   err=max(err,abs((e7(ii)-e6(ii))/(e72(ii)-e6(ii))-4.0_dp) &
&      /(4.0_dp*(2.0_dp*kf(ii)/c_light)**2))
 end do
 call check('relativistic O(1/c**2)',err,1.0_dp)
 end subroutine check_relativistic

 end program test_xc_physics

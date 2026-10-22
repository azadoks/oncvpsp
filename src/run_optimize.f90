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
 subroutine run_optimize(eig,ll,kap,irel,mmax,mxprj,rr,uua,upa,qq,&
&                      freetl,rc,irc,qcut,qmsbf,ncon_in,nbas_in,npr, &
&                      psopt,vpsp,vkb,vae,cvgplt,ncurve,ecurve,pswfx,chix, &
&                      ektail)

!! calls routines to generate optimized pseudo wave function and semi-local
!! pseudopotential and prints diagnostic information on process and
!! convergence performance
!!
!! @param eig  energy at which pseudo wave function is computed
!! @param ll  angular momentum
!! @param kap  =l, -(l+1) for j=l -/+ 1/2 (irel=2 only)
!! @param irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!! @param mmax  dimension of log grid
!! @param mxprj  dimension of number of projectors
!! @param rr  log radial grid
!! @param uua  all-electron wave functions (*rr) for the projectors
!! @param upa  their derivatives in the mesh index, as the solvers return them
!! @param qq  matrix of norms and overlaps of the uua inside rc
!! @param freetl  for each projector, whether E_resid takes the free tail
!!                of its well state (free_tail_m) instead of the uua
!! @param rc  core radius (bohr)
!! @param irc  last mesh point at or inside rc
!! @param qcut  q cutoff defining residual energy
!! @param qmsbf  maximum q in sbfs for this ll
!! @param ncon_in  number of constraints for matching psuedo and AE wave
!!                 functions
!! @param nbas_in  number of basis functions for pseudo wave function
!! @param npr  number of projectors for this ll
!! @param psopt  optimized pseudo wave function(s)
!! @param vpsp  corresponding pseudopotential
!! @param vkb  Vanderbilt-Kleinman-Bylander projectors (without local v
!!             correction)
!! @param vae  all-electron potential; the well potentials of the scattering
!!             states equal it out to rr(irc+5), so the derivatives of the uua
!!             at rc are those of the radial equation in vae
!! @param cvgplt  energy error vs. cutoff energy for plotting
!! @param ncurve  dimension of ecurve
!! @param ecurve  cutoff energy and energy error per electron on the E_resid(q)
!!                mesh, zero beyond its end
!! @param pswfx  the inside pieces of the pseudo wave functions (their Bessel
!!               sums) on the nodes to past rc, 0 beyond (pspot)
!! @param chix  (eig - T) pswfx, the projectors' seeds, the same
!! @param ektail  the r-space kinetic energy of each optimized wave less the
!!                Fourier one to q_infinity (Ha): the energy above
!!                q_infinity, which E_resid leaves out, plus the r-space
!!                integral's error

 use constants_m, only: dp, nerrlv, nerrdg, pi
 use, intrinsic :: iso_fortran_env, only: qp => real128
 use ugrid_m, only: ug_cross_down
 use free_tail_m, only: mesh_tail, free_tail, tail_kinetic
 use logmesh_m, only: lm_integral
 implicit none

!Input variables
 integer, intent(in) :: ll,kap,irel,mmax,mxprj,ncon_in,nbas_in,npr,ncurve
 real(dp), intent(in) :: rr(mmax)
 real(dp), intent(in) :: uua(mmax,mxprj),upa(mmax,mxprj),vae(mmax)
 real(dp), intent(in) :: qq(mxprj,mxprj)
 real(dp), intent(in) :: eig(mxprj),qcut,rc
 logical, intent(in) :: freetl(mxprj)

!Input/output variables
 integer, intent(inout) :: irc
 real(dp), intent(inout) :: vpsp(mmax),vkb(mmax,mxprj)
 real(dp), intent(inout) :: cvgplt(2,nerrdg,mxprj)

!Output variables
 real(dp), intent(out) :: qmsbf
 real(dp), intent(out) :: psopt(mmax,mxprj)
 real(dp), intent(out) :: ecurve(2,ncurve,mxprj)
 real(dp), intent(out) :: pswfx(mmax,mxprj),chix(mmax,mxprj)
 real(dp), intent(out) :: ektail(mxprj)

!Local variables
 real(dp), parameter :: Ha_eV=27.21138386d0 ! 1 Hartree, in eV
 real(dp) :: uord(5)
 real(dp) :: al,ulgd,tt
 real(dp) :: sn,ps0norm,amesh,ro
 real(dp) :: err,qerr,ehaerr,eeverr,dq,q0
 real(dp) :: sbf(mmax)
 real(dp) :: cons(10)
 real(qp), allocatable :: orbasis(:,:),orbasis_der(:,:)
 real(qp), allocatable :: pswf0_sbq(:),pswfnull_sbq(:,:)
 real(dp), allocatable :: pswf0_sb(:),pswfnull_sb(:,:)
 real(dp), allocatable :: pswf0_or(:),pswfnull_or(:,:)
 real(dp), allocatable :: work(:),ulin(:),dlin(:),tke(:)
 integer :: ii,iprj,jj,l1,lmax,nbas,ncon,nconmx
 logical :: found_root

 integer :: nnull,nqout,nq,iq,kk,nr

 real(dp), allocatable :: qroot(:),sbasis(:)
 real(dp), allocatable :: eresid0(:),eresiddot(:,:)
 real(dp), allocatable :: eresidmat(:,:,:),pswfresid_sb(:,:)
 real(dp), allocatable :: qout(:),eresq(:),eresg(:)
 real(dp), allocatable :: e0(:),edot(:,:),emat(:,:,:),pswfopt_nu(:)

 real(dp), allocatable :: pswfopt_sb(:),pswfopt_or(:),psb(:,:)
 real(dp) :: eresidmin
 real(dp) :: ekin_anal,ekin_num,qmax,hh,ekin_tail


 al = 0.01d0 * log(rr(101) / rr(1))
 l1=ll+1

 nconmx=ncon_in+npr-1

 allocate(work(mmax))

!maximum basis size
 nbas=nbas_in+npr-1
 allocate(qroot(nbas))
 allocate(psb(nbas,mxprj))
 psb(:,:)=0.0d0

!calculate derivatives of all-electron wave function at r_c from the radial
!equation for the first projector to set sbf wave vectors for all projectors

 call ode_rc_der(irel,ll,kap,eig(1),rr,vae,uua(1,1),upa(1,1),mmax,rc,uord)
 ulgd=uord(2)/uord(1)

 qmax=200.0d0  !very large value should always allow nbas q's to be found

!select wave vectors for spherical Bessel functions, scanned in steps of
!0.02 (doc/smoothness_plan.md, Q-5, D-3)
 call qroots(ll,rc,ulgd,nbas,0.02d0,qmax,qroot)

 qmsbf=qroot(nbas)

 psopt(:,:)=0.0d0
 ektail(:)=0.0d0

! loop over projectors
 do iprj=1,npr

   write(6,'(/a,i4)') 'Calculating optimized projector #',iprj, &
&        ' for l=',ll

!basis size and number of constraints for this projector
   nbas=nbas_in+iprj-1
   ncon=ncon_in+iprj-1
   nnull=nbas-ncon

   allocate(sbasis(nbas))
   allocate(orbasis(nbas,nbas))
   allocate(orbasis_der(ncon,nbas))
   allocate(pswf0_sb(nbas),pswf0_or(nbas))
   allocate(pswfopt_sb(nbas),pswfopt_or(nbas))
   allocate(pswfresid_sb(nbas,nnull))
   allocate(pswfnull_sb(nbas,nnull))
   allocate(pswf0_sbq(nbas),pswfnull_sbq(nbas,nnull))
   allocate(pswfnull_or(nbas,nnull))

!calculate derivatives of all-electron wave function at r_c from the radial
!equation to define basis set for current projector and derivative
!constraints

   call ode_rc_der(irel,ll,kap,eig(iprj),rr,vae,uua(1,iprj),upa(1,iprj), &
&                  mmax,rc,uord)

!q considered infinity for E_resid calculation: twice the largest basis
!wave vector, at least 20 and at least twice qcut, by smooth maxima of
!width 1 bohr^-1 (doc/smoothness_plan.md, Q-3); the q grid of E_resid has
!150 points per 2 pi/rc, down from qmax to below dq (eresid)
   qmax=smax(smax(2.0d0*qroot(nbas),20.0d0),2.0d0*qcut)
   dq=2.0d0*pi/(150.0d0*rc)
   nq=1+int(qmax/dq)
   q0=qmax-dble(nq-1)*dq

!E_resid is wanted at the exact qcut and at q=0 (the kinetic energy)
   nqout=2

   allocate(qout(nqout),eresq(nqout),eresid0(nqout))
   allocate(eresiddot(nnull,nqout),eresidmat(nnull,nnull,nqout))
   allocate(e0(nq),edot(nq,nnull),emat(nq,nnull,nnull),eresg(nq))
   allocate(pswfopt_nu(nnull))

   qout(1)=qcut
   qout(2)=0.0d0

!calculate orthogonal basis and constraint matrix
   call sbf_basis_con(ll,rc,nbas,qroot,psb,nbas_in+npr-1,orbasis, &
&                     orbasis_der,iprj,mxprj,ncon,ncon_in)


!load constraint vector for value/derivative matching
   cons(:)=0.0d0
   do jj=1,ncon_in
     cons(jj)=uord(jj)
   end do
!load overlap constraints if present
   if(iprj>=2) then
     do jj=1,iprj-1
       cons(ncon_in+jj)=qq(jj,iprj)
     end do
   end if

!calculate constrained basis for residual energy minimization
!satisfying off-diagonal norm conservation
   call const_basis(nbas,ncon,cons,orbasis,orbasis_der, &
&                   pswf0_or,pswfnull_or, &
&                   pswf0_sbq,pswfnull_sbq,ps0norm)
   pswf0_sb(:)=real(pswf0_sbq(:),dp)
   pswfnull_sb(:,:)=real(pswfnull_sbq(:,:),dp)

!calculate eigenvectors, eigenvalues, and inhomogeneous terms for
!the residual energy for a set of q lower cutoffs

   write(6,'(/a,f10.6)') '    Fraction of norm inside rc',qq(iprj,iprj)

!the all-electron tail on the linear mesh of E_resid, 75 points per
!period of qmax: the wave function itself, or the free tail of a well
!state (doc/smoothness_plan.md, phase 8)
   hh=2.0d0*pi/(75.0d0*qmax)
   if(freetl(iprj)) then
     call free_tail(irel,ll,kap,eig(iprj),rr,vae,uua(1,iprj),upa(1,iprj), &
&                   mmax,rc,qcut,hh,nr,ulin,dlin)
     ekin_tail=tail_kinetic(ll,rc,hh,nr,ulin,dlin)
     deallocate(dlin)
   else
     call mesh_tail(rr,uua(1,iprj),mmax,rc,hh,nr,ulin)
   end if

!E_resid on the q grid, and at the exact qcut and the points of the
!convergence curve (Q-1, Q-4)
   call eresid(ll,rc,nnull,nbas,nq,dq,nqout,qout,qmax,qroot,hh,nr,ulin, &
&              pswf0_sbq,pswfnull_sbq,e0,edot,emat, &
&              eresid0,eresiddot,eresidmat)
   deallocate(ulin)

   write(6,'(a,f7.2,a,f7.2,a)') '    Optimizing pswf for qcut=',qout(1), &
&   ' a_B^-1,  ecut=',0.5d0*qout(1)**2,' Ha'
   write(6,'(a,f6.2,a,f7.1,a)') '    q_infinity defining residual KE=',&
&   qmax,'  (E_inf=',0.5d0*qmax**2,' Ha)'

!find the null-basis coefficients which minimize the eresid while
!satisfying diagonal norm conservation

   call optimize(nnull,nbas,pswf0_sb,pswf0_or,nqout,qout,&
&                eresid0,eresiddot,eresidmat,&
&                pswfnull_sb,pswfnull_or,qq(iprj,iprj),ps0norm,eresidmin, &
&                pswfopt_sb,pswfopt_or,ekin_anal,eresq,pswfopt_nu)

   write(6,'(a,1p,d10.2,a)') '    Residual kinetic energy error=', &
&        eresidmin,' Ha'
!the optimized coefficients from the quad-precision basis vectors (those
!of optimize, from the rounded ones, carry their 1e4 cancellation)
   do ii=1,nbas
     pswfopt_sb(ii)=real(pswf0_sbq(ii) &
&                        +sum(pswfnull_sbq(ii,:)*real(pswfopt_nu(:),qp)),dp)
   end do
   psb(1:nbas,iprj)=pswfopt_sb(1:nbas)

!the convergence curve of the optimized pseudo wave function on the q
!grid, for the cutoffs at each error level (below) and the model core's
!budget (modcore_budget_table): all of it, as a table of every nstep-th
!node changed its interpolant when nstep did (M-13)
   do iq=1,nq
     eresg(iq)=e0(iq)
     do jj=1,nnull
       eresg(iq)=eresg(iq)+2.0d0*edot(iq,jj)*pswfopt_nu(jj)
       do kk=1,nnull
         eresg(iq)=eresg(iq)+pswfopt_nu(kk)*emat(iq,kk,jj)*pswfopt_nu(jj)
       end do
     end do
   end do
   if(nq>ncurve) then
     write(6,'(/2a,i0,a,i0)') 'run_optimize: ERROR q grid longer than', &
&          ' the convergence curve, nq=',nq,', ncurve=',ncurve
     error stop 1
   end if
   ecurve(:,:,iprj)=0.0d0
   do iq=1,nq
     ecurve(1,iq,iprj)=0.5d0*(q0+dble(iq-1)*dq)**2
     ecurve(2,iq,iprj)=eresg(iq)
   end do

! find the Vanderbilt projectors and optimized wave functions

   call pspot(iprj,ll,rr,irc,rc,mmax,al,nbas,qroot,eig(iprj),uua(1,iprj), &
&             pswfopt_sb,pswfx(1,iprj),chix(1,iprj),psopt(1,iprj),vae,work, &
&             vkb(1,iprj),ekin_num)

!semi-local potential
   if(iprj==1) then
     do ii=1,irc
       vpsp(ii)=work(ii)
     end do
     do ii=irc+1,mmax
       vpsp(ii)=vae(ii)
     end do
   end if

!with the free tail, the r-space kinetic energy is that of the same
!function: the Bessel sum's to rc (T pswfx = eig pswfx - chix, as u T u
!goes as r**(2l+2) at the origin) and the tail's
   if(freetl(iprj)) then
     allocate(tke(mmax))
     tke(:)=pswfx(:,iprj)*(eig(iprj)*pswfx(:,iprj)-chix(:,iprj))
     call lm_integral(rr,tke,mmax,2*ll+2,0.0d0,rc,ekin_num)
     ekin_num=ekin_num+ekin_tail
     deallocate(tke)
   end if

   ektail(iprj)=ekin_num-ekin_anal
   write(6,'(/a)') '    Total kinetic energy consistency test'
   write(6,'(a)') '      Fourier integration compared to d^2/dr^2 integral'
   write(6,'(a,f12.8,a,f12.8,a,f12.8)') '      Fourier',ekin_anal, &
&    '  r-space=',ekin_num,'  ratio=',ekin_anal/ekin_num
   write(6,'(a)') '    Potential consistency test at r_c'
   write(6,'(a,f12.8,a,f12.8,a,1p,d10.2)') '    "vpsp"=',work(irc),'  vae=', &
&    vae(irc),'  difference=',vae(irc)-work(irc)

!interpolate the convergence behavior of the optimized pseudo wave function

   write(6,'(/a)') '    Energy error per electron        Cutoff'
   write(6,'(a)') '         Ha          eV             Ha'
!the cutoffs at which the curve falls to each error level, by
!root-finding on it; the table shows the first nerrlv, the DIAGNOSTICS
!block all of them
   err=0.01d0
   do jj=1,nerrdg
     call ug_cross_down(eresg,nq,dq,q0,err,qerr,found_root)
     if(found_root) then
       cvgplt(1,jj,iprj)= 0.5d0*qerr**2
       cvgplt(2,jj,iprj)= err
       if(mod(jj,2)/=0 .and. jj<=nerrlv) then
         write(6,'(4x,2f12.5,f12.2)') err,err*Ha_eV,0.5d0*qerr**2
       end if
     end if
     err=sqrt(0.1d0)*err
   end do

   deallocate(sbasis)
   deallocate(orbasis)
   deallocate(orbasis_der)
   deallocate(pswf0_sb,pswf0_or)
   deallocate(pswfopt_sb,pswfopt_or)
   deallocate(pswfresid_sb)
   deallocate(pswfnull_sb)
   deallocate(pswf0_sbq,pswfnull_sbq)
   deallocate(pswfnull_or)
   deallocate(eresiddot,eresidmat)
   deallocate(qout,eresq,eresid0)
   deallocate(e0,edot,emat,eresg,pswfopt_nu)

 end do !iprj

 deallocate(work,qroot,psb)

 return

 contains

 real(dp) function smax(aa,bb)

!! a C-infinity maximum of two wave vectors, within 0.5 bohr^-1 of the
!! larger: (a + b + sqrt((a-b)**2 + 1))/2

 real(dp), intent(in) :: aa,bb

 smax=0.5d0*(aa+bb+sqrt((aa-bb)**2+1.0d0))
 end function smax

 end subroutine run_optimize

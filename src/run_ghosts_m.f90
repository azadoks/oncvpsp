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

 module run_ghosts_m

 implicit none
 private
 public :: run_ghosts

 contains

 subroutine run_ghosts(lmax,la,ea,nc,nv,lloc,rc,qmsbf, &
&                    vlocpw,vkbpw,evkb,nproj,rr,vp,mmax,adams, &
&                    ngh,ghn,ghe,ghp,iout,ifail)

!! Two tests for ghosts;  The local potential is terminated by a hard wall
!! barrier at a radius of mbfact*rc, presently 3*rc.  A basis set is formed
!! from eigenstates this combination and used to compute a hamiltonian matrix
!! whose off-diagonal terms come exclusively from the non-local pojectors.
!! A cutoff is chosen to simulate a reasonably well converged plane-wave
!! calculation.  The hamiltoian is diagonalized, and the results used
!! as follows:
!!
!! The wall is at 3*rc exactly: the barrier states are found on the log mesh
!! scaled by s in [1, amesh) to put a node on it, the local potential and
!! the projectors interpolated there from their pieces. The cutoff is
!! smooth: each barrier state's overlaps with the projectors are weighted
!! by cinf_taper of its energy, 1 up to ebcut and 0 from ebfact*ebcut, so
!! a state enters or leaves the basis with weight 0. The eigenvalues are
!! then smooth in the inputs (doc/smoothness_plan.md, X-3), and so are the
!! margins returned for the diagnostics block.
!!
!! Test 1) Negative eigenvalues are compared with those of the radial
!! Schroedinger equation for the full-range pseudopotential.  Those
!! lying below the Schroedinger Eq. results and which presumably have
!! more nodes are reported and tagged as GHOST(-).
!!
!! Test2) For barrier states with  positive eigenvalues, the average radius
!! is computed.  Such states could be ghost resonances, and highly localized
!! ones with <r>/rc<1 are reported as GHOST(+).  Real fairly localized
!! resonances can occur, particularly at lower energies, and the GHOST(+)
!! report is an indication that the log-derivative plot should be examined
!! carefully around the reported energy to see if the all-electron plot
!! has a corresponding drop of pi  (arctan(log der), of course).
!!
!! @param lmax  maximum angular momentum
!! @param la  angular momenta
!! @param ea  all-electron bound-state energies (Ha)
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param lloc  l of the local potential, 4 for the polynomial
!! @param rc  core radius (bohr) of each l
!! @param qmsbf  maximum q in sbf basis for each l
!! @param vlocpw  the screened local potential as pieces
!! @param vkbpw  the projectors as pieces
!! @param evkb  coefficients of the projectors (Ha)
!! @param nproj  number of projectors for each l
!! @param rr  log radial mesh (bohr)
!! @param vp  semi-local pseudopotentials (vp(:,5) is local potential if linear
!!            comb.)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param ngh  rows of the bound-state test per l
!! @param ghn  their principal quantum numbers
!! @param ghe  their energies: (1) of the nonlocal Schroedinger equation
!!             (0 for the vacuum level), (2) of the basis diagonalization
!! @param ghp  per l, the basis state of least <r>/rc above the reference,
!!             <r>/rc divided by a weight rising from 0 at 0.025 Ha above
!!             it to 1 at 0.05 Ha: (1) its energy, (2) that ratio; 0 if none
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_numerical
 use constants_m, only: dp, pi, mxstate, mxl, mxlin, mxprj
 use pieces_m, only: piecewise_t, pw_integral, pw_value, pw_resample, &
&                    pw_free, pw_node_values
 use logmesh_m, only: lm_interp, lm_step, lm_index
 use quad_m, only: quad_log_norm
 use cinf_taper_m, only: cinf_taper
 use lschpsbar_m, only: lschpsbar
 use lschvkbb_m, only: lschvkbb
 use lapack_m, only: dsyev
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: nc,nv,lmax,lloc,mmax
 type(adams_t), intent(in) :: adams
 integer, intent(in) :: la(mxstate),nproj(mxlin+1)
 real(dp), intent(in) :: rr(mmax),vp(mmax,mxl+2),rc(mxlin+1)
 real(dp), intent(in) :: qmsbf(mxlin+1),ea(mxstate),evkb(mxprj,mxl+1)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1)
 integer, intent(in) :: iout

!Output variables
 integer, intent(out) :: ngh(mxl+1),ghn(mxstate,mxl+1)
 real(dp), intent(out) :: ghe(2,mxstate,mxl+1),ghp(2,mxl+1)
 integer, intent(out) :: ifail

!Local variables
 integer :: ll,l1,jj,ierr,mbar,npgh
 integer :: nprj,nbas,nbmax,lwork,info
 integer :: nodes(mxl+1),lpgh(20)
 real(dp) :: al,amesh,ee,ecut,eps,emin,emax,tht
 real(dp) :: rbar,vmin,ebtop,ww
 real(dp) :: cpgh(20),ebcut(mxlin+1),epgh(20),rpgh(20)
 real(dp), allocatable :: uu(:),up(:),rb(:),vb(:),vsw(:)
 real(dp), allocatable :: gg0(:,:),hmat(:,:),hmev(:),work(:)
 real(dp), allocatable :: uua(:,:)
 type(piecewise_t) :: vkbb(mxprj)
! the wall at mbfact*rc; the basis weights fall from 1 at basfact times
! the largest kinetic energy of the sbf basis to 0 at ebfact times that
 real(dp), parameter :: mbfact=3.0d0,basfact=1.25d0,ebfact=1.5d0

 ifail=0
 allocate(uu(mmax),up(mmax),rb(mmax),vb(mmax),vsw(mmax))

! the local potential is vp(:,lloc+1), the pieces vlocpw plus the switch
! to the ionic tail beyond rcmax (ionic_screening), a smooth function on
! the nodes that is interpolated, as in outwf_setup; the reference
! energies of lschvkbb are those of the same potential
! (doc/numerics_audit.md, C-7)
 call pw_node_values(rr,mmax,vlocpw,vsw)
 vsw(:)=vp(:,lloc+1)-vsw(:)

 al=lm_step(rr,mmax)
 amesh = exp(al)
 eps=1.0d-4 !threshold for negative-energy ghost report
 tht=0.0d0 !sets barrier boundary conditions to zero

! loop for diagnostic output using Vanderbilt Kleinman-Bylander projectors

 write(iout,'(/a)') ' Testing for bound ghosts'
 write(iout,'(5a)') '   n','   l','    E NL Schr. Eq', &
&      '   E Basis Diag.','   E Cutoff'

 nodes(:)=0
 npgh=0
 ngh(:)=0
 ghn(:,:)=0
 ghe(:,:,:)=0.0d0
 ghp(:,:)=0.0d0
 do l1 = 1, lmax + 1
   nprj=nproj(l1)
   if(nprj==0) cycle

   write(iout,'()')
   ll = l1 - 1

   call barrier_mesh(l1,ifail)
   if(ifail/=0) return

! Basis cutoff energy is based on largest q in sbf optimization basis
   ebcut(l1)= basfact*0.5d0*qmsbf(l1)**2
   ebtop=ebfact*ebcut(l1)

! Estimate maximum basis set size nbmax: the states of the wall's sphere
! below ebtop over the deepest local potential, and two to spare

   nbmax=int(sqrt(2.0d0*(ebtop-vmin))*rbar/pi)+2

   lwork=3*nbmax

   allocate(gg0(nprj,nbmax),hmat(nbmax,nbmax),hmev(nbmax),work(lwork))
   allocate(uua(mmax,nbmax))

   gg0(:,:)=0.0d0
   hmat(:,:)=0.0d0
   hmev(:)=0.0d0


   call barrier_basis(l1)

! add non-local terms to hamiltonian matrix
   call add_nonlocal(l1)

! find eigenvalues of the hamiltonian matrix

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

   call dsyev( 'V', 'U', nbas, hmat, nbmax, hmev, work, lwork, info )
   if(info /= 0) then
     write(iout,'(a,i4)') 'run_ghosts: hamiltonian matrix eigenvalue ERROR, &
&          info=',info
     ifail=fail_numerical
     return
   end if
   call bound_tests(l1,ifail)
   if(ifail/=0) return

   call localized_states(l1)


   deallocate(gg0,hmat,hmev,work)
   deallocate(uua)
 end do !l1

 write(iout,'(/a)') ' Testing for highly-localized positive-energy ghosts'
 write(iout,'(5a)') '    ','   l','    <radius>/rc  ', &
&      '   E Basis Diag.','   E Cutoff'

 if(npgh>0) then
   do jj=1,npgh
     write(iout,'(/4x,i4,2f16.6,f10.2,a)') lpgh(jj),rpgh(jj), &
&          epgh(jj),cpgh(jj),'  WARNING - GHOST(+)'
   end do
 end if

 deallocate(uu,up,rb,vb,vsw)
 return

 contains

 subroutine barrier_mesh(l1,ifail)
!! the basis of local-pseudopotential plus hard-wall barrier
!! eigenfunctions of l1-1: the barrier at mbfact*rc, node mbar of the mesh
!! rb = s*rr, the last node of rr at or inside it scaled up onto it; the
!! local potential vb on rb (vmin its least inside the barrier) and the
!! projectors vkbb
 use failure_m, only: fail_input
 integer, intent(in) :: l1
 integer :: ii,iprj
 real(dp) :: vs
 integer, intent(out) :: ifail
 ifail=0
 rbar=mbfact*rc(l1)
 mbar=lm_index(rr,mmax,rbar)
! lschpsbar starts the solution at mbar:mbar+4, inside rr(mmax)
 if(mbar+5>mmax) then
   write(iout,'(/2a,i0,a,i0)') 'run_ghosts: ERROR barrier radius beyond', &
&          ' radial mesh, mbar=',mbar,', mmax=',mmax
   ifail=fail_input
   return
 end if
 rb(:)=(rbar/rr(mbar))*rr(:)
 rb(mbar)=rbar
! the local potential on rb, its last node beyond rr(mmax) given the last
! value
 vmin=0.0d0
 do ii=1,mmax
   if(rb(ii)<=rr(mmax)) then
     call lm_interp(rr,vsw,mmax,rb(ii),vs)
     vb(ii)=pw_value(rr,mmax,vlocpw,rb(ii))+vs
   else
     vb(ii)=vp(mmax,lloc+1)
   end if
   if(ii<=mbar) vmin=min(vmin,vb(ii))
 end do
 do iprj=1,nprj
   call pw_resample(rr,mmax,vkbpw(iprj,l1),rb,mbar+4,vkbb(iprj))
 end do
 end subroutine barrier_mesh

 subroutine barrier_basis(l1)
!! the barrier states of l1-1 below ebtop (nbas of them), the diagonal
!! of the hamiltonian, and their overlaps with the projectors, weighted by
!! their energies; ecut the cutoff of the basis
 integer, intent(in) :: l1

 integer :: jj,iprj,nn
! nbas stays nbmax unless lschpsbar runs out of states below ebtop
 nbas=nbmax
 ecut=ebcut(l1)
 do jj=1,nbmax
   nn=ll+jj

   emin=0.0d0
   emax=ebtop
   ee=0.0d0

   call lschpsbar(nn,ll,ierr,ee,emin,emax,rb,vb,uu,up, &
&                   mmax,adams,mbar,tht)

   if(ierr/=0) then
     if(emax/=ebtop) then
! a barrier state below the cutoff that does not converge, a high one on
! a coarse mesh whose wavelength at the barrier is a few mesh steps: the
! basis stops below it, and the ghost test covers energies up to there
       write(iout,'(/a,i0,a,i0,a,f12.6,a/a,f12.6,a)') &
&              'run_ghosts: WARNING barrier state n=',nn,' l=',ll, &
&              ' near ',ee,' Ha not converged;', &
&              'run_ghosts: WARNING basis cut there, ghost test to ', &
&              emin,' Ha'
       ecut=emin
     end if
     nbas=jj-1
     exit
   end if

! barrier potential eigenvalues are contirbution to hamiltonian matrix diagonal
   hmat(jj,jj)=ee

   uua(:,jj)=uu(:)

! overlap of eigenfunction jj and projectors, per piece (the barrier
! state is smooth across the projectors' kinks), weighted by its energy
   ww=cinf_taper(ee,ebcut(l1),ebtop)
   do iprj=1,nprj
     call pw_integral(rb,mmax,2*ll+2,uu,vkbb(iprj),gg0(iprj,jj))
     gg0(iprj,jj)=ww*gg0(iprj,jj)
   end do

 end do !jj
 if(nbas==nbmax) then
   write(iout,'(/2a,i0,a,i0)') 'run_ghosts: WARNING basis full below the', &
&          ' cutoff, l=',ll,', nbmax=',nbmax
 end if
 do iprj=1,nprj
   call pw_free(vkbb(iprj))
 end do
 end subroutine barrier_basis

 subroutine add_nonlocal(l1)
!! the non-local terms of the hamiltonian matrix of l1-1
 integer, intent(in) :: l1
 integer :: ii,jj,iprj
 do jj=1,nbas
   do ii=1,jj
     do iprj=1,nprj
       hmat(ii,jj)=hmat(ii,jj)+evkb(iprj,l1)*gg0(iprj,ii)*gg0(iprj,jj)
     end do
     hmat(jj,ii)=hmat(ii,jj)
   end  do
 end do
 end subroutine add_nonlocal

 subroutine bound_tests(l1,ifail)
!! the bound-state test of l1-1: each valence state of l against the
!! basis eigenvalues, a lower one a GHOST(-); without valence states, any
!! nodeless bound state, or the vacuum level (0.0); ee the energy of the
!! last
 integer, intent(in) :: l1
 integer :: ii,jj,kk,nn,mch
 real(dp) :: et
!treat valence states for this l
 integer, intent(out) :: ifail
 ifail=0
 nn=ll
 jj=0
 do kk=1,nv
   if(la(nc+kk)/=ll) cycle
   jj=jj+1
   nn=nn+1
   ee=ea(nc+kk)
   emax=0.9d0*ea(nc+kk)
   emin=1.1d0*ea(nc+kk)
   call lschvkbb(nn,ll,nprj,ierr,ee,emin,emax, &
&                  rr,vp(1,lloc+1),vlocpw,vkbpw(1:nprj,l1),evkb(1,l1), &
&                  uu,up,mmax,adams,mch,iout,ifail)
   if(ifail/=0) return
   if(ierr/=0) then
     write(iout,'(a,3i4,2f12.6)') 'run_ghosts: lschvkbb ERROR', &
&            nn,ll,ierr,ea(nc+kk),ee
   end if

   do ii=1,10
     et=min(ee-eps,ee*(1.0d0+eps))
     if(hmev(jj)<et) then
       write(iout,'(4x,i4,16x,f16.6,f10.2,a)') ll,hmev(jj),ecut, &
&              '  WARNING - GHOST(-)'
       jj=jj+1
     else
       exit
     end if
   end do !ii

   write(iout,'(2i4,2f16.6,f10.2)') nn,ll,ee,hmev(jj),ecut
   ngh(l1)=ngh(l1)+1
   ghn(ngh(l1),l1)=nn
   ghe(:,ngh(l1),l1)=[ee,hmev(jj)]
 end do !kk

!if there were no valence states, look for any nodeless bound state
!if this fails, use vacuum level (0.0) to compare to basis-set states

 if(jj==0) then
   ee=-1.0d5
   do kk=1,nv
     ee=max(ee,ea(nc+kk))
   end do
   emax=0.0d0
   emin=2.0d0*ee
   call lschvkbb(ll+1,ll,nprj,ierr,ee,emin,emax, &
&                  rr,vp(1,lloc+1),vlocpw,vkbpw(1:nprj,l1),evkb(1,l1), &
&                  uu,up,mmax,adams,mch,iout,ifail)
   if(ifail/=0) return
   if(ierr/=0) ee=0.0d0
   do jj=1,10
     et=min(ee-eps,ee*(1.0d0+eps))
     if(hmev(jj)<et) then
       write(iout,'(4x,i4,16x,f16.6,f10.2,a)') ll,hmev(jj),ecut, &
&              '  WARNING - GHOST(-)'
     else
       if(ee<0.0d0) then
         write(iout,'(2i4,2f16.6,f10.2)') ll+1,ll,ee,hmev(jj),ecut
       end if
       ngh(l1)=1
       ghn(1,l1)=ll+1
       ghe(:,1,l1)=[ee,hmev(jj)]
       exit
     end if
   end do
 end if
 end subroutine bound_tests

 subroutine localized_states(l1)
!! the mean radius of each eigenstate of l1-1, to look for highly
!! localized states, which indicate positive-energy ghosts
 integer, intent(in) :: l1
 integer :: ii,jj
 real(dp) :: sn,ww
 do jj=1,nbas
   uu(:)=0.0d0
   do ii=1,nbas
     uu(:)=uu(:)+hmat(ii,jj)*uua(:,ii)
   end do

! <r> on the barrier's mesh, by the integral that normalized the states
! (quad_log_norm), r u**2 going as r**(2l+3)
   sn=quad_log_norm(rb,rb(1:mbar)*uu(1:mbar)**2,mbar,ll,dble(ll)+1.5d0)

! the margin: the least <r>/rc over the states above the reference, each
! divided by a C-infinity weight in its energy, 1 from 0.05 Ha above the
! reference (the states the test took before) and 0 at 0.025 Ha, so that
! a state entering does not make the least jump; a GHOST(+) is a state
! whose ratio is below 1, the margin below 1 (doc/numerics_audit.md, C-8)
   ww=1.0d0-cinf_taper(hmev(jj),ee+0.025d0,ee+0.05d0)
   if(ww>0.0d0) then
     if(ghp(2,l1)==0.0d0 .or. sn/(rc(l1)*ww)<ghp(2,l1)) then
       ghp(:,l1)=[hmev(jj),sn/(rc(l1)*ww)]
     end if
     if(sn<rc(l1)*ww .and. npgh<20) then
       npgh=npgh+1
       lpgh(npgh)=ll
       epgh(npgh)=hmev(jj)
       rpgh(npgh)=sn/(rc(l1)*ww)
       cpgh(npgh)=ecut
     end if
   end if
 end do
 end subroutine localized_states

 end subroutine run_ghosts

 end module run_ghosts_m

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
 subroutine run_ghosts(lmax,la,ea,nc,nv,lloc,irc,qmsbf, &
&                    vlocpw,vkbpw,evkb,nproj,rr,vp,mmax,mxprj)

!! Two tests for ghosts;  The local potential is terminated by a hard wall
!! barrier at a radius of mbfact*rc, presently 3*rc.  A basis set is formed
!! from eigenstates this combination and used to compute a hamiltonian matrix
!! whose off-diagonal terms come exclusively from the non-local pojectors.
!! A cutoff is chosen to simulate a reasonably well converged plane-wave
!! calculation.  The hamiltoian is diagonalized, and the results used
!! as follows:
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
!! @param la  angular momenta of all-electron states
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param ea  all-electron bound-state energies
!! @param lloc  l for local potential
!! @param irc  core radii indices
!! @param vlocpw  the screened local potential's pieces
!! @param vkbpw  VKB projectors as functions of pieces
!! @param qmsbf  maximum q in sbf basis for each l
!! @param evkb  coefficients of VKB projectors
!! @param nproj  number of vkb projectors for each l
!! @param rr  log radial grid
!! @param vp  semi-local pseudopotentials (vp(:,5) is local potential if linear
!!            comb.)
!! @param mmax  size of radial grid
!! @param mxprj  dimension of number of projectors

 use constants_m, only: dp, pi, mxstate, mxl, mxlin
 use pieces_m, only: piecewise_t, pw_integral
 implicit none

!Input variables
 integer, intent(in) :: nc,nv,lmax,lloc,mmax,mxprj
 integer, intent(in) :: la(mxstate),irc(mxlin+1),nproj(mxlin+1)
 real(dp), intent(in) :: rr(mmax),vp(mmax,mxl+2)
 real(dp), intent(in) :: qmsbf(mxlin+1),ea(mxstate),evkb(mxprj,mxl+1)
 type(piecewise_t), intent(in) :: vlocpw,vkbpw(mxprj,mxl+1)

!Local variables
 integer :: ll,l1,ii,jj,kk,ierr,mbfact,mbar,mch,nn,npgh
 integer :: iprj,nprj,nbas,nbmax,lwork,info
 integer :: nodes(mxl+1),lpgh(20)
 real(dp) :: al,amesh,basfact,ee,ecut,eps,et,emin,emax,sn,ro,tht
 real(dp) :: cpgh(20),ebcut(mxlin+1),epgh(20),rpgh(20)
 real(dp), allocatable :: uu(:),up(:)
 real(dp), allocatable :: gg0(:,:),hmat(:,:),hmev(:),work(:)
 real(dp), allocatable :: uua(:,:)

 allocate(uu(mmax),up(mmax))

 al = 0.01d0 * log(rr(101)/rr(1))
 amesh = exp(al)
 eps=1.0d-4 !threshold for negative-energy ghost report
 tht=0.0d0 !sets barrier boundary conditions to zero

! loop for diagnostic output using Vanderbilt Kleinman-Bylander projectors

 write(6,'(/a)') ' Testing for bound ghosts'
 write(6,'(5a)') '   n','   l','    E NL Schr. Eq', &
&      '   E Basis Diag.','   E Cutoff'

 nodes(:)=0
 npgh=0
 do l1 = 1, lmax + 1
   nprj=nproj(l1)
   if(nprj==0) cycle

   write(6,'()')
   ll = l1 - 1

!Construct basis of local-pseidupotential plus hard-wall barrier eigenfunctions
!Barrier is set at mbfact*rc

   mbfact=3
   mbar = int(log(mbfact*rr(irc(l1))/rr(1))/al+1.0d0)
! lschpsbar starts the solution at mbar:mbar+4
   if(mbar+4>mmax) then
     write(6,'(/2a,i0,a,i0)') 'run_ghosts: ERROR barrier radius beyond', &
&          ' radial mesh, mbar=',mbar,', mmax=',mmax
     error stop 1
   end if

! Basis cutoff energy is based on largest q in sbf optimization basis
   basfact=1.25d0
   ebcut(l1)= basfact*0.5d0*qmsbf(l1)**2

! Estimate maximum basis set size nbmax based on ll=0 sbf
! with 20% allowance for binding by local potential

   nbmax=int(2.0d0*sqrt(1.2d0*ebcut(l1))*rr(mbar)/pi)

   lwork=3*nbmax

   allocate(gg0(nprj,nbmax),hmat(nbmax,nbmax),hmev(nbmax),work(lwork))
   allocate(uua(mmax,nbmax))

   gg0(:,:)=0.0d0
   hmat(:,:)=0.0d0
   hmev(:)=0.0d0


! nbas stays nbmax unless lschpsbar runs out of states below ebcut
   nbas=nbmax
   ecut=ebcut(l1)
   do jj=1,nbmax
     nn=ll+jj

     emin=0.0d0
     emax=ebcut(l1)
     ee=0.0d0

     call lschpsbar(nn,ll,ierr,ee,emin,emax,rr,vp(1,lloc+1),uu,up, &
&                   mmax,mbar,tht)

     if(ierr/=0) then
       if(emax/=ebcut(l1)) then
! a barrier state below the cutoff that does not converge, a high one on
! a coarse mesh whose wavelength at the barrier is a few mesh steps: the
! basis stops below it, and the ghost test covers energies up to there
         write(6,'(/a,i0,a,i0,a,f12.6,a/a,f10.2,a)') &
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
! state is smooth across the projectors' kinks)
     do iprj=1,nprj
       call pw_integral(rr,mmax,2*ll+2,uu,vkbpw(iprj,l1),gg0(iprj,jj))
     end do

   end do !jj

! add non-local terms to hamiltonian matrix
   do jj=1,nbas
     do ii=1,jj
       do iprj=1,nprj
         hmat(ii,jj)=hmat(ii,jj)+evkb(iprj,l1)*gg0(iprj,ii)*gg0(iprj,jj)
       end do
       hmat(jj,ii)=hmat(ii,jj)
     end  do
   end do

! find eigenvalues of the hamiltonian matrix

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

   call dsyev( 'V', 'U', nbas, hmat, nbmax, hmev, work, lwork, info )
   if(info /= 0) then
     write(6,'(a,i4)') 'run_ghosts: hamiltonian matrix eigenvalue ERROR, &
&          info=',info
     error stop 1
   end if
!treat valence states for this l
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
&                  uu,up,mmax,mch)
     if(ierr/=0) then
       write(6,'(a,3i4,2f12.6)') 'run_ghosts: lschvkbb ERROR', &
&            nn,ll,ierr,ea(nc+kk),ee
     end if

     do ii=1,10
       et=min(ee-eps,ee*(1.0d0+eps))
       if(hmev(jj)<et) then
         write(6,'(4x,i4,16x,f16.6,f10.2,a)') ll,hmev(jj),ecut, &
&              '  WARNING - GHOST(-)'
         jj=jj+1
       else
         exit
       end if
     end do !ii

     write(6,'(2i4,2f16.6,f10.2)') nn,ll,ee,hmev(jj),ecut
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
&                  uu,up,mmax,mch)
     if(ierr/=0) ee=0.0d0
     do jj=1,10
       et=min(ee-eps,ee*(1.0d0+eps))
       if(hmev(jj)<et) then
         write(6,'(4x,i4,16x,f16.6,f10.2,a)') ll,hmev(jj),ecut, &
&              '  WARNING - GHOST(-)'
       else if(ee<0.0d0) then
         write(6,'(2i4,2f16.6,f10.2)') ll+1,ll,ee,hmev(jj),ecut
         exit
       else
         exit
       end if
     end do
   end if

!calculate mean radius of all eigenstates to look for highly localized
!states which indicate postive-energy ghosts
   do jj=1,nbas
     uu(:)=0.0d0
     do ii=1,nbas
       uu(:)=uu(:)+hmat(ii,jj)*uua(:,ii)
     end do

     ro=rr(1)/sqrt(amesh)
     sn=(ro*uu(1))**2/real(2*ll+4,dp)

     do ii=1,mbar-3
       sn=sn+al*(rr(ii)*uu(ii))**2
     end do

     sn=sn + al*(23.0d0*(rr(mbar-2)*uu(mbar-2))**2 &
&              + 28.0d0*(rr(mbar-1)*uu(mbar-1))**2 &
&              +  9.0d0*(rr(mbar  )*uu(mbar  ))**2)/24.0d0

     if(sn<rr(irc(l1)) .and. hmev(jj)>ee+0.05d0 .and. npgh<20) then
       npgh=npgh+1
       lpgh(npgh)=ll
       epgh(npgh)=hmev(jj)
       rpgh(npgh)=sn/rr(irc(l1))
       cpgh(npgh)=ecut
     end if
   end do


   deallocate(gg0,hmat,hmev,work)
   deallocate(uua)
 end do !l1

 write(6,'(/a)') ' Testing for highly-localized positive-energy ghosts'
 write(6,'(5a)') '    ','   l','    <radius>/rc  ', &
&      '   E Basis Diag.','   E Cutoff'

 if(npgh>0) then
   do jj=1,npgh
     write(6,'(/4x,i4,2f16.6,f10.2,a)') lpgh(jj),rpgh(jj), &
&          epgh(jj),cpgh(jj),'  WARNING - GHOST(+)'
   end do
 end if

 deallocate(uu,up)
 return
 end subroutine run_ghosts

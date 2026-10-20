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
 subroutine sbf_basis_con(ll,rc,nbas,qroot,psb,nbasmx,orbasis,orbasis_der, &
&                         iprj,mxprj,ncon,ncon_in)

!! orthonormalize basis functions and derivatives at rc and form constraint
!! matrix based on derivative to be matched and overlaps with prior
!! optimized wave functions. Every overlap is the closed-form integral of
!! two spherical Bessel functions inside rc (sbf_overlap), so the
!! orthonormal basis is exact to round-off and smooth in rc: the basis is
!! nearly dependent (the extra wave vectors of qroots; the overlap matrix
!! has an eigenvalue of 5e-9), and the mesh integrals used before were
!! accurate to about 1e-8, so its last orthonormal function was mostly
!! their error (doc/smoothness_plan.md, R-5)
!!
!! @param ll  angujlar momentum
!! @param rc  core radius (bohr)
!! @param nbas  number of basis functions
!! @param qroot  q values for j_l(q*r)
!! @param psb  basis coefficients of the optimized pseudo wave functions of
!!             the lower iprj (projector kk uses the first nbas-iprj+kk)
!! @param nbasmx  first dimension of psb
!! @param orbasis  matrix for orthonormal basis (rows j_l, columns basis
!!                 vectors)
!! @param orbasis_der  values and derivatives of orthonormal basis set at rc,
!!                     and norm constraint vectors from already-computed psopt
!! @param iprj  index of the projector being optimized
!! @param mxprj  second dimension of psb, maximum number of projectors
!! @param ncon  total number of constraints, ncon_in+iprj-1
!! @param ncon_in  number of value and derivative constraints at rc

 use constants_m, only: dp
 implicit none
 external :: sbf_overlap, sbf_rc_der, dsyev

!Input variables
 integer, intent(in) :: ll,nbas,nbasmx,iprj,mxprj,ncon,ncon_in
 real(dp), intent(in) :: rc
 real(dp), intent(in) :: qroot(nbas),psb(nbasmx,mxprj)

!Output variables
 real(dp), intent(out) :: orbasis(nbas,nbas),orbasis_der(ncon,nbas)

!Local variables
 integer :: ii,jj,kk,ibas,info,nbk
 real(dp) :: sn,tt
 real(dp) :: sbfder(5)
 real(dp), allocatable :: sev(:),work(:),sovlp(:,:),ovl(:,:),vv(:)

 allocate(sev(nbas),work(5*nbas),sovlp(nbas,nbas),ovl(nbas,nbas),vv(nbas))

!overlap matrix of the basis functions inside rc, in closed form
 do ibas=1,nbas
   do jj=ibas,nbas
     call sbf_overlap(ll,qroot(ibas),qroot(jj),rc,sn)
     ovl(ibas,jj)=sn
     ovl(jj,ibas)=sn
   end do
 end do
 sovlp(:,:)=ovl(:,:)

!find eigenvalues and eigenvectors of the overlap matrix

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

 call dsyev( 'V', 'U', nbas, sovlp, nbas, sev, work, 5*nbas, info )
 if(info /= 0) then
   write(6,'(a,i4)') 'sbf_basis: overlap matrix eigenvalue error, info=',info
   error stop 1
 end if

!scale eigenvectors to form orthonormal basis coefficients for sbf's
!note that we are reversing order so that the leading eigenvector is the
!most linearly independent
 do ibas=1,nbas
   if(sev(ibas)>0.0d0) then
     tt=1.0d0/sqrt(sev(ibas))
   else
     write(6,'(a,f12.6)') 'sbfbasis: negative eigenvalue of overlap matrix'
     error stop 1
   end if
   do jj=1,nbas
     orbasis(jj,nbas-ibas+1)=tt*sovlp(jj,ibas)
   end do
 end do

!find rc derivatives of basis funtction
 orbasis_der(:,:)=0.0d0
 do jj=1,nbas !sbf loop
   call sbf_rc_der(ll,qroot(jj),rc,sbfder)
   do ibas=1,nbas !orbasis loop
     do ii=1,ncon_in
       orbasis_der(ii,ibas)=orbasis_der(ii,ibas)+orbasis(jj,ibas)*sbfder(ii)
     end do
   end do
 end do

!fill in last rows of orbasis_der constraint matrix with the overlaps of
!the orthonormal functions with the previous optimized wave functions,
!<phi_j|psopt_kk> = sum_i orbasis(i,j) sum_m ovl(i,m) psb(m,kk)

 if(iprj>=2) then
   do kk=1,iprj-1
     nbk=nbas-iprj+kk
     do ii=1,nbas
       vv(ii)=0.0d0
       do jj=1,nbk
         vv(ii)=vv(ii)+ovl(ii,jj)*psb(jj,kk)
       end do
     end do
     do jj=1,nbas
       sn=0.0d0
       do ii=1,nbas
         sn=sn+orbasis(ii,jj)*vv(ii)
       end do
       orbasis_der(ncon_in+kk,jj)=sn
     end do !jj
   end do !kk
 end if !iprj>=2

 deallocate(sev,work,sovlp,ovl,vv)

 return
 end subroutine sbf_basis_con

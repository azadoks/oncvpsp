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

 module sbf_basis_con_m

 implicit none
 private
 public :: sbf_basis_con

 contains

 subroutine sbf_basis_con(ll,rc,nbas,qroot,psb,nbasmx,orbasis,orbasis_der, &
&                         iprj,ncon,ncon_in,iout,ifail)

!! orthonormalize basis functions and derivatives at rc and form constraint
!! matrix based on derivative to be matched and overlaps with prior
!! optimized wave functions, in quad precision: the basis is nearly
!! dependent (the extra wave vectors of qroots; its overlap matrix has
!! eigenvalues down to 5e-9), so an orthonormal function is a combination
!! with coefficients up to 1e4, and in double precision its derivatives at
!! rc, the constraint matrix, were good to 1e-12 only. Every overlap is the
!! closed-form integral of two spherical Bessel functions inside rc
!! (sbfq_overlap), so the basis is exact to the round-off of quad precision
!! and smooth in rc (doc/smoothness_plan.md, R-5)
!!
!! @param ll  angular momentum
!! @param rc  core radius (bohr)
!! @param nbas  number of basis functions
!! @param qroot  q values for j_l(q*r)
!! @param psb  basis coefficients of the optimized pseudo wave functions of
!!             the lower iprj (projector kk uses the first nbas-iprj+kk)
!! @param nbasmx  first dimension of psb
!! @param orbasis  matrix for orthonormal basis (rows j_l, columns basis
!!                 vectors), quad precision
!! @param orbasis_der  values and derivatives of orthonormal basis set at rc,
!!                     and norm constraint vectors from already-computed psopt,
!!                     quad precision
!! @param iprj  index of the projector being optimized
!! @param ncon  total number of constraints, ncon_in+iprj-1
!! @param ncon_in  number of value and derivative constraints at rc
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: fail_numerical
 use constants_m, only: dp, qp, mxprj
 use sbf_q_m, only: sbfq_overlap, sbfq_rc_der
 use qlinalg_m, only: qjacobi

!Input variables
 integer, intent(in) :: ll,nbas,nbasmx,iprj,ncon,ncon_in
 real(dp), intent(in) :: rc
 real(dp), intent(in) :: qroot(nbas),psb(nbasmx,mxprj)
 integer, intent(in) :: iout

!Output variables
 real(qp), intent(out) :: orbasis(nbas,nbas),orbasis_der(ncon,nbas)
 integer, intent(out) :: ifail

!Local variables
 integer :: ii,jj,kk,ibas,nbk
 real(qp) :: sn,tt,rcq
 real(qp) :: sbfder(5)
 real(qp), allocatable :: sev(:),svec(:,:),ovl(:,:),vv(:)

 ifail=0
 allocate(sev(nbas),svec(nbas,nbas),ovl(nbas,nbas),vv(nbas))
 rcq=real(rc,qp)

!overlap matrix of the basis functions inside rc, in closed form
 do ibas=1,nbas
   do jj=ibas,nbas
     call sbfq_overlap(ll,real(qroot(ibas),qp),real(qroot(jj),qp),rcq,sn)
     ovl(ibas,jj)=sn
     ovl(jj,ibas)=sn
   end do
 end do

!find eigenvalues and eigenvectors of the overlap matrix
 call qjacobi(nbas,ovl,sev,svec)

!scale eigenvectors to form orthonormal basis coefficients for sbf's
!note that we are reversing order so that the leading eigenvector is the
!most linearly independent
 do ibas=1,nbas
   if(sev(ibas)>0.0_qp) then
     tt=1.0_qp/sqrt(sev(ibas))
   else
     write(iout,'(a,es12.4)') &
&          'sbfbasis: negative eigenvalue of overlap matrix', &
&         real(sev(ibas),dp)
     ifail=fail_numerical
     return
   end if
   do jj=1,nbas
     orbasis(jj,nbas-ibas+1)=tt*svec(jj,ibas)
   end do
 end do

!find rc derivatives of basis funtction
 orbasis_der(:,:)=0.0_qp
 do jj=1,nbas !sbf loop
   call sbfq_rc_der(ll,real(qroot(jj),qp),rcq,sbfder)
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
       vv(ii)=0.0_qp
       do jj=1,nbk
         vv(ii)=vv(ii)+ovl(ii,jj)*real(psb(jj,kk),qp)
       end do
     end do
     do jj=1,nbas
       sn=0.0_qp
       do ii=1,nbas
         sn=sn+orbasis(ii,jj)*vv(ii)
       end do
       orbasis_der(ncon_in+kk,jj)=sn
     end do !jj
   end do !kk
 end if !iprj>=2

 deallocate(sev,svec,ovl,vv)

 return
 end subroutine sbf_basis_con

 end module sbf_basis_con_m

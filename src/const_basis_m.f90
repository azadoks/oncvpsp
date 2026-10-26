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

 module const_basis_m

 implicit none
 private
 public :: const_basis

 contains

!calculates basis for residual energy minimization

 subroutine const_basis(nbas,ncon,cons,orbasis,orbasis_der,&
&                       pswf0_or,pswfnull_or, &
&                       pswf0_sb,pswfnull_sb,ps0norm,ifail)

!! the particular solution of the matching constraints of least norm in the
!! orthonormal basis, and an orthonormal basis of their null space, in
!! quad precision (doc/smoothness_plan.md, Q-8, R-5): each row of the
!! constraint matrix (a derivative of one order at rc, or an overlap) is
!! scaled to unit norm with its constraint, which changes neither the
!! solution set nor the null space; then from the eigen-decomposition of
!! A^T A, the null space its eigenvectors of eigenvalue zero and the
!! solution sum_k v_k (v_k . A^T b)/lambda_k over the others
!!
!! @param ncon  number of constraints
!! @param nbas  number of orthogonalized sbf basis functions
!! @param cons  constraints (value at rc and derivatives up to 4th)
!! @param orbasis  coefficients of sbfs in orthogonal basis (quad)
!! @param orbasis_der  values and derivatives of orthogonalized sbf basis
!!                     functions at rc (constraint matrix, quad)
!! @param pswf0_or  pswf0_sb in the orthogonalized basis
!! @param pswfnull_or  pswfnull_sb in the orthogonalized basis
!! @param pswf0_sb  linear combination of sbfs that matches cons at rc (quad)
!! @param pswfnull_sb  nbas-ncon linear combinations of sbfs satisfying value
!!                     and ncon-1 derivatives at rc == 0 (null space of
!!                     constraint matrix) These form an orthonormal set, and are
!!                     all orthogonal to the pswf0 combination (which is not
!!                     normalized) (quad)
!! @param ps0norm  charge contained in ps0 inside rc
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: qp => real128
 use qlinalg_m, only: qjacobi
 implicit none

!Input variables
 integer, intent(in) :: ncon,nbas
 real(qp), intent(in) :: orbasis(nbas,nbas),orbasis_der(ncon,nbas)
 real(dp), intent(in) :: cons(ncon)

!Output variables
 real(dp), intent(out) :: pswf0_or(nbas),pswfnull_or(nbas,nbas-ncon)
 real(qp), intent(out) :: pswf0_sb(nbas),pswfnull_sb(nbas,nbas-ncon)
 real(dp), intent(out) :: ps0norm
 integer, intent(out) :: ifail

!Local variables
 real(qp) :: amat(ncon,nbas),bq(ncon),gmat(nbas,nbas),lam(nbas),vv(nbas,nbas)
 real(qp) :: atb(nbas),xor(nbas),con_tst(ncon),rscale,sn
! the constraints are rank-deficient below this ratio of the smallest to
! the largest singular value of the scaled constraint matrix
 real(qp), parameter :: rankeps=1.0e-10_qp
 integer :: ii,jj,kk,nnull

 ifail=0
 nnull=nbas-ncon

!the scaled constraint matrix and vector
 do ii=1,ncon
   rscale=1.0_qp/sqrt(sum(orbasis_der(ii,:)**2))
   amat(ii,:)=rscale*orbasis_der(ii,:)
   bq(ii)=rscale*real(cons(ii),qp)
 end do

!A^T A and its eigen-decomposition, eigenvalues ascending
 do ii=1,nbas
   do jj=1,nbas
     gmat(ii,jj)=sum(amat(:,ii)*amat(:,jj))
   end do
 end do
 call qjacobi(nbas,gmat,lam,vv)
 if(lam(nnull+1)<rankeps**2*lam(nbas)) then
   write(6,'(/a,1p,e10.2)') &
&        'const_basis: ERROR the constraints are linearly dependent, ', &
&        'singular value ratio', &
&        real(sqrt(max(lam(nnull+1),0.0_qp)/lam(nbas)),dp)
   ifail=1
   return
 end if

!the least-norm particular solution
 do ii=1,nbas
   atb(ii)=sum(amat(:,ii)*bq(:))
 end do
 xor(:)=0.0_qp
 do kk=nnull+1,nbas
   xor(:)=xor(:)+vv(:,kk)*sum(vv(:,kk)*atb(:))/lam(kk)
 end do
 pswf0_or(:)=real(xor(:),dp)

!the null space
 do kk=1,nnull
   pswfnull_or(:,kk)=real(vv(:,kk),dp)
 end do

!now, test result
 do ii=1,ncon
   con_tst(ii)=sum(orbasis_der(ii,:)*xor(:))
 end do

 write(6,'(/a)') '    Constraint consistency test'
 write(6,'(a)') '      pswf0 val/der aewf val/der     diff'
 do ii=1,ncon
   write(6,'(4x,2f14.8,1p,d12.2)') real(con_tst(ii),dp),cons(ii), &
&        real(con_tst(ii),dp)-cons(ii)
 end do

!calculate charge contained in ps0 inside rc
 ps0norm=real(sum(xor(:)**2),dp)

!convert ps0 and psnull in orthogonal representation to sbf represention
 do jj=1,nbas
   sn=0.0_qp
   do kk=1,nbas
     sn=sn+orbasis(jj,kk)*xor(kk)
   end do
   pswf0_sb(jj)=sn
   do ii=1,nnull
     sn=0.0_qp
     do kk=1,nbas
       sn=sn+orbasis(jj,kk)*vv(kk,ii)
     end do
     pswfnull_sb(jj,ii)=sn
   end do
 end do

 return
 end subroutine const_basis

 end module const_basis_m

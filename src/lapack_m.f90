!
! Copyright (c) 2026 by Austin Zadoks
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
!> Interfaces of the LAPACK and BLAS routines that ONCVPSP calls, as the
!> reference implementation declares them, so that the compiler checks
!> every call. The libraries are external; only the interfaces are here.

 module lapack_m

 use constants_m, only: dp
 implicit none
 private
 public :: dgels,dgemm,dgesv,dgglse,dhseqr,dsyev

 interface

   subroutine dgels(trans,mm,nn,nrhs,aa,lda,bb,ldb,work,lwork,info)
!! least squares of a full-rank system, by QR or LQ
   import :: dp
   implicit none
   character(len=1), intent(in) :: trans
   integer, intent(in) :: mm,nn,nrhs,lda,ldb,lwork
   real(dp), intent(inout) :: aa(lda,*),bb(ldb,*)
   real(dp), intent(out) :: work(*)
   integer, intent(out) :: info
   end subroutine dgels

   subroutine dgemm(transa,transb,mm,nn,kk,alpha,aa,lda,bb,ldb,beta, &
&                   cc,ldc)
!! cc = alpha op(aa) op(bb) + beta cc
   import :: dp
   implicit none
   character(len=1), intent(in) :: transa,transb
   integer, intent(in) :: mm,nn,kk,lda,ldb,ldc
   real(dp), intent(in) :: alpha,beta,aa(lda,*),bb(ldb,*)
   real(dp), intent(inout) :: cc(ldc,*)
   end subroutine dgemm

   subroutine dgesv(nn,nrhs,aa,lda,ipiv,bb,ldb,info)
!! solves aa x = bb by LU factorization
   import :: dp
   implicit none
   integer, intent(in) :: nn,nrhs,lda,ldb
   real(dp), intent(inout) :: aa(lda,*),bb(ldb,*)
   integer, intent(out) :: ipiv(*),info
   end subroutine dgesv

   subroutine dgglse(mm,nn,pp,aa,lda,bb,ldb,cc,dd,xx,work,lwork,info)
!! least squares of aa x - cc subject to bb x = dd
   import :: dp
   implicit none
   integer, intent(in) :: mm,nn,pp,lda,ldb,lwork
   real(dp), intent(inout) :: aa(lda,*),bb(ldb,*),cc(*),dd(*)
   real(dp), intent(out) :: xx(*),work(*)
   integer, intent(out) :: info
   end subroutine dgglse

   subroutine dhseqr(job,compz,nn,ilo,ihi,hh,ldh,wr,wi,zz,ldz,work, &
&                    lwork,info)
!! eigenvalues (and Schur form) of an upper Hessenberg matrix
   import :: dp
   implicit none
   character(len=1), intent(in) :: job,compz
   integer, intent(in) :: nn,ilo,ihi,ldh,ldz,lwork
   real(dp), intent(inout) :: hh(ldh,*),zz(ldz,*)
   real(dp), intent(out) :: wr(*),wi(*),work(*)
   integer, intent(out) :: info
   end subroutine dhseqr

   subroutine dsyev(jobz,uplo,nn,aa,lda,ww,work,lwork,info)
!! eigenvalues (and eigenvectors) of a real symmetric matrix
   import :: dp
   implicit none
   character(len=1), intent(in) :: jobz,uplo
   integer, intent(in) :: nn,lda,lwork
   real(dp), intent(inout) :: aa(lda,*)
   real(dp), intent(out) :: ww(*),work(*)
   integer, intent(out) :: info
   end subroutine dsyev

 end interface

 end module lapack_m

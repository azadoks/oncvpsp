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
!calculates basis for residual energy minimization

 subroutine const_basis(nbas,ncon,cons,orbasis,orbasis_der,&
&                       pswf0_or,pswfnull_or, &
&                       pswf0_sb,pswfnull_sb,ps0norm)

!! @param ncon  number of constraints
!! @param nbas  number of orthogonalized sbf basis functions
!! @param cons  constraints (value at rc and derivatives up to 4th)
!! @param orbasis  coefficients of sbfs in orthogonal basis
!! @param orbasis_der  values and derivatives of orthogonalized sbf basis
!!                     functions at rc (constraint matrix)
!! @param pswf0_or  pswf0_sb in the orthogonalized basis
!! @param pswfnull_or  pswfnull_sb in the orthogonalized basis
!! @param pswf0_sb  linear combination of sbfs that matches cons at rc
!! @param pswfnull_sb  nbas-ncon linear combinations of sbfs satisfying value
!!                     and ncon-1 derivatives at rc == 0 (null space of
!!                     constraint matrix) These form an orthonormal set, and are
!!                     all orthogonal to the pswf0 combination (which is not
!!                     normalized)
!! @param ps0norm  charge contained in ps0 inside rc

 use constants_m, only: dp
 implicit none

!Input variables
 integer, intent(in) :: ncon,nbas
 real(dp), intent(in) :: orbasis(nbas,nbas),orbasis_der(ncon,nbas),cons(ncon)

!Output variables
 real(dp), intent(out) :: pswf0_or(nbas),pswfnull_or(nbas,nbas-ncon)
 real(dp), intent(out) :: pswf0_sb(nbas),pswfnull_sb(nbas,nbas-ncon)
 real(dp), intent(out) :: ps0norm

!Local variables
 real(dp) :: usvd(10,10),ss(10),work(100),con_tst(10),rscale(10)
 real(dp), allocatable :: amat(:,:),vt(:,:)
 real(dp) :: sn
! the constraints are rank-deficient below this ratio of the smallest to
! the largest singular value of the scaled constraint matrix
 real(dp), parameter :: rankeps=1.d-10
 integer :: ii,jj,kk,info


!DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
!     $                   WORK, LWORK, INFO )
 allocate(amat(10,nbas),vt(nbas,nbas))

!save orbasis_der from overwrite; each row of the constraint matrix (a
!derivative of one order at rc, or an overlap) is scaled to unit norm, with
!its constraint, which changes neither the solution set nor the null space
!but the accuracy of the decomposition: the rows differ in size by the
!derivative orders, (q rc)**4 and more, which the particular solution
!returned as round-off at eps times that ratio (1e-11, doc/smoothness_plan.md
!Q-8)
 amat(:,:)=0.0d0
 do ii=1,ncon
   rscale(ii)=1.0d0/sqrt(sum(orbasis_der(ii,:)**2))
   amat(ii,:)=rscale(ii)*orbasis_der(ii,:)
 end do

!singular value decomposition of constraint matrix
 call dgesvd('A','A',ncon,nbas,amat,10,ss,usvd,10,vt,nbas, &
&             work,64,info)

 if(info /= 0) then
   write(6,'(/a,i4)') 'const_basis: ERTOR const_basis: dgesvd info = ',info
   error stop 1
 end if

 if(ss(ncon)<rankeps*ss(1)) then
   write(6,'(/a,1p,e10.2)') &
&        'const_basis: ERROR the constraints are linearly dependent, ', &
&        'singular value ratio',ss(ncon)/ss(1)
   error stop 1
 end if

!operate on constraint vector with usvd^transpose
 work(:)=0.0d0
 do ii=1,ncon
   work(ii)=0.0d0
   do jj=1,ncon
     work(ii)=work(ii)+usvd(jj,ii)*rscale(jj)*cons(jj)
   end do
 end do

!scale work with inverse transpose of singular value "matrix" ss
 do ii=1,ncon
   work(ii)=work(ii)/ss(ii)
 end do

!form coefficient vector for basic matching pseudo wave function
!null vectors (with ss=0) can be added to improve kinetic residual
!convergence without changing match at rc
 do ii=1,nbas
   pswf0_or(ii)=0.0d0
   do jj=1,ncon
     pswf0_or(ii)=pswf0_or(ii)+work(jj)*vt(jj,ii)
   end do
 end do

!copy coefficient vectors spanning null space of constraint matrix
!note that vt (transpose) rows become columns of pswfnull
 do jj=1,nbas-ncon
   sn=0.0d0
   do ii=1,nbas
     pswfnull_or(ii,jj)=vt(ncon+jj,ii)
     sn=sn+pswfnull_or(ii,jj)**2
   end do
 end do


!now, test result
 do ii=1,ncon
   con_tst(ii)=0.0d0
   do jj=1,nbas
     con_tst(ii)=con_tst(ii)+orbasis_der(ii,jj)*pswf0_or(jj)
   end do
 end do

 write(6,'(/a)') '    Constraint consistency test'
 write(6,'(a)') '      pswf0 val/der aewf val/der     diff'
 do ii=1,ncon
   write(6,'(4x,2f14.8,1p,d12.2)') con_tst(ii),cons(ii),con_tst(ii)-cons(ii)
 end do

!calculate charge contained in ps0 inside rc
 ps0norm=0.0d0
 do ii=1,nbas
   ps0norm=ps0norm+pswf0_or(ii)**2
 end do

!convert ps0 and psnull in orthogonal representation to sbf represention

 pswf0_sb(:)=0.0d0
 pswfnull_sb(:,:)=0.0d0
 do jj=1,nbas
   do kk=1,nbas
     pswf0_sb(jj)     =pswf0_sb(jj)     +orbasis(jj,kk)*pswf0_or(kk)
     do ii=1,nbas-ncon
       pswfnull_sb(jj,ii)=pswfnull_sb(jj,ii)+orbasis(jj,kk)*pswfnull_or(kk,ii)
     end do
   end do
 end do

 deallocate(amat,vt)
 return
 end subroutine const_basis

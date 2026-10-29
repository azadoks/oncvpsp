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
!> The energy iteration the bound-state solvers share (lschfb, ldiracfb,
!> lschvkbb, lschvkbbe, lschpsbar): the convergence test on the energy
!> correction, and the steps of the energy inside its bracket, by the
!> correction or, where the node count is wrong, by bisection; and the
!> energy of the box state that stands in where there is no bound state;
!> and the whole iteration of the all-electron solvers (lschfb, ldiracfb),
!> bs_solve, over the steps that depend on the equation.
 module boundsolve_m
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: error_unit
 implicit none
 private
 public :: bs_converged, bs_correct, bs_bisect, bs_box_energy, bs_solve, &
&          bs_turning_point

 abstract interface
   subroutine node_count(ee,nnode,ifail)
!! the node count nnode of the outward solution over the whole mesh at the
!! energy ee; ifail, 0 or the kind of failure (failure_m), whose message is
!! written
   import :: dp
   implicit none
   real(dp), intent(in) :: ee
   integer, intent(out) :: nnode,ifail
   end subroutine node_count

   subroutine bs_setup(ee,mch)
!! the equation's coefficients at the energy ee, and the classical turning
!! point mch from them (bs_turning_point)
   import :: dp
   implicit none
   real(dp), intent(in) :: ee
   integer, intent(out) :: mch
   end subroutine bs_setup

   subroutine bs_outward(mch,mlast,nnode,ifail)
!! the outward solution at the energy of the last bs_setup, up to mlast;
!! nnode, its nodes inside the turning point mch; ifail, 0 or the kind of
!! failure (failure_m), whose message is written
   implicit none
   integer, intent(in) :: mch,mlast
   integer, intent(out) :: nnode,ifail
   end subroutine bs_outward

   subroutine bs_match(ee,mch,nin,lo,hi,de,ifail)
!! the inward solution from nin, matched to the outward one across the
!! window lo..hi around mch and normalized; de, the energy correction;
!! ifail, 0 or the kind of failure (failure_m), whose message is written
   import :: dp
   implicit none
   real(dp), intent(in) :: ee
   integer, intent(in) :: mch,nin,lo,hi
   real(dp), intent(out) :: de
   integer, intent(out) :: ifail
   end subroutine bs_match

   subroutine bs_box(elow,etop,ee,ifail)
!! the box state (bs_box_energy) above elow, widened above etop, and its
!! energy ee; ifail, 0 or the kind of failure (failure_m), whose message
!! is written
   import :: dp
   implicit none
   real(dp), intent(in) :: elow,etop
   real(dp), intent(out) :: ee
   integer, intent(out) :: ifail
   end subroutine bs_box
 end interface

 contains

 logical function bs_converged(de,deprev,ee,eps)

!! whether the energy ee has converged: its correction de is below
!! eps*max(|ee|,eig_floor) and either at round-off (root_tol) or no longer
!! halving from the one before, deprev
!!
!! @param de  the energy correction
!! @param deprev  the correction before
!! @param ee  the energy
!! @param eps  the relative convergence threshold

 use constants_m, only: eig_floor
 use root_m, only: root_tol

!Input variables
 real(dp), intent(in) :: de,deprev,ee,eps

 bs_converged=abs(de)<max(abs(ee),eig_floor)*eps .and. &
&             (abs(de)<=root_tol(ee,eig_floor) .or. &
&              abs(de)>=0.5d0*abs(deprev))
 return
 end function bs_converged

 subroutine bs_correct(de,ee,emin,emax,half)

!! the energy moved by its correction de, after the bracket's end on its
!! side moves to ee (half: halfway to ee); bisection of the bracket if the
!! step leaves it
!!
!! @param de  the energy correction
!! @param ee  the energy
!! @param emin,emax  the bracket
!! @param half  whether the end moves only halfway to ee

!Input variables
 real(dp), intent(in) :: de
 logical, intent(in) :: half

!Input/output variables
 real(dp), intent(inout) :: ee,emin,emax

 if(de>0.0d0) then
   if(half) then
     emin=0.5d0*(emin+ee)
   else
     emin=ee
   end if
 else
   if(half) then
     emax=0.5d0*(emax+ee)
   else
     emax=ee
   end if
 end if
 ee=ee+de
 if(ee>emax .or. ee<emin) ee=0.5d0*(emax+emin)
 return
 end subroutine bs_correct

 subroutine bs_bisect(ee,emin,emax,too_low,half)

!! the bracket's lower end (too_low: too few nodes) or upper end moved to
!! ee (half: halfway to ee), and ee to the middle of the bracket
!!
!! @param ee  the energy
!! @param emin,emax  the bracket
!! @param too_low  whether ee is below the state
!! @param half  whether the end moves only halfway to ee

!Input variables
 logical, intent(in) :: too_low,half

!Input/output variables
 real(dp), intent(inout) :: ee,emin,emax

 if(too_low) then
   if(half) then
     emin=0.5d0*(emin+ee)
   else
     emin=ee
   end if
 else
   if(half) then
     emax=0.5d0*(emax+ee)
   else
     emax=ee
   end if
 end if
 ee=0.5d0*(emin+emax)
 return
 end subroutine bs_bisect

 subroutine bs_box_energy(elow,etop,nnode,eps,nodes,ee,ifail)

!! the energy of the box state (lschfb, ldiracfb): the state with nnode
!! nodes and u=0 at the mesh end, above elow (the box state of
!! RadialAtomicSolvers.jl). The outward solution over the whole mesh gains
!! a node each time u at the mesh end crosses zero, so the energy is
!! bisected on its node count: the bracket from elow is widened above etop
!! by doubling until the count exceeds nnode, then bisected to a relative
!! width eps (of max(|E|,0.2)); ee is its lower end, the highest energy
!! found with nnode nodes or fewer
!!
!! @param elow  lower end of the bracket
!! @param etop  energy above which the bracket's top is widened
!! @param nnode  number of nodes of the box state
!! @param eps  relative width of the final bracket
!! @param nodes  the node count of the outward solution over the whole
!!               mesh at an energy
!! @param ee  the box state's energy
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: eig_floor
 use failure_m, only: fail_convergence, fail_not_found, failed, &
&                    msg_box_energy_bisection, msg_box_energy_nodes

!Input variables
 real(dp), intent(in) :: elow,etop,eps
 integer, intent(in) :: nnode
 procedure(node_count) :: nodes

!Output variables
 real(dp), intent(out) :: ee
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: eblo,ebhi
 integer :: kk,nb

 ee=elow
 eblo=elow
 ebhi=etop+1.0d0
 do kk=1,60
   call nodes(ebhi,nb,ifail)
   if(ifail/=0) return
   if(nb>nnode) exit
   eblo=ebhi
   ebhi=etop+2.0d0*(ebhi-etop)
 end do
 if(kk>60) then
   write(error_unit,'(/a,i0,a,es12.4)') &
&       'bs_box_energy: ERROR no energy with more than ',nnode, &
&       ' nodes up to ',ebhi
   ifail=failed(fail_not_found,msg_box_energy_nodes)
   return
 end if

 do kk=1,200
   ee=0.5d0*(eblo+ebhi)
   call nodes(ee,nb,ifail)
   if(ifail/=0) return
   if(nb>nnode) then
     ebhi=ee
   else
     eblo=ee
   end if
   if(ebhi-eblo<max(abs(ee),eig_floor)*eps) exit
 end do
 if(kk>200) then
   write(error_unit,'(/a,es12.4)') &
&       'bs_box_energy: ERROR bisection not converged in 200 steps, e=',ee
   ifail=failed(fail_convergence,msg_box_energy_bisection)
   return
 end if

 ee=eblo

 return
 end subroutine bs_box_energy

 subroutine bs_turning_point(cf,mmax,mch)

!! the classical turning point: the outermost mch with cf(mch-1)<=0 and
!! cf(mch)>0, cf the coefficient of u in the second-order equation;
!! mch=0 where none is found and cf(mmax)<=0 (classically allowed at the
!! mesh end, the energy too high), mch=-1 where none is found and
!! cf(mmax)>0 (forbidden everywhere, too low)
!!
!! @param cf  the coefficient, positive where classically forbidden
!! @param mmax  size of the log radial mesh
!! @param mch  the turning point, or 0 or -1

!Input variables
 integer, intent(in) :: mmax
 real(dp), intent(in) :: cf(mmax)

!Output variables
 integer, intent(out) :: mch

!Local variables
 integer :: ii

 mch=0
 do ii=mmax,2,-1
   if(cf(ii-1)<=0.0d0 .and. cf(ii)>0.0d0) then
     mch=ii
     exit
   end if
 end do
 if(mch==0 .and. cf(mmax)>0.0d0) mch=-1
 return
 end subroutine bs_turning_point

 subroutine bs_solve(nn,ll,rr,vv,mmax,adams,eps,setup,outward,match,box, &
&                    ee,emin,emax,mch,ierr,ifail)

!! the energy iteration of the all-electron bound-state solvers: from ee
!! in the bracket emin..emax, the turning point, the outward solution and
!! its node count; with the right count the inward solution matched to it
!! and the energy moved by its correction, otherwise the bracket bisected.
!! Once a correction is below eps the passes go on at the corrected
!! energy until it is converged (bs_converged), and the energy and
!! function are those of the last pass. A bracket that closes before that
!! has no bound state, and the box state stands in for it
!!
!! @param nn  principal quantum number
!! @param ll  angular momentum
!! @param rr  log radial mesh (bohr)
!! @param vv  potential (Ha)
!! @param mmax  size of the log radial mesh
!! @param adams  order and coefficients of the Adams-Moulton steps
!! @param eps  relative convergence threshold of the energy
!! @param setup  the coefficients and turning point at an energy
!! @param outward  the outward solution and its node count
!! @param match  the inward solution, the match, normalization and the
!!               energy correction
!! @param box  the box state
!! @param ee  the start energy on input, the energy on output (Ha)
!! @param emin,emax  the bracket
!! @param mch  the matching point, mmax for the box state
!! @param ierr  0 on success; the non-zero values differ by routine:
!!              ierr_box_state (failure_m) for the box state,
!!              ierr_iterations if not converged in 60 passes
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use failure_m, only: ierr_box_state, ierr_iterations
 use match_m, only: match_k
 use inward_start_m, only: inward_start
 use constants_m, only: eig_floor
 use adams_m, only: adams_t

!Input variables
 integer, intent(in) :: nn,ll,mmax
 real(dp), intent(in) :: rr(mmax),vv(mmax),eps
 type(adams_t), intent(in) :: adams
 procedure(bs_setup) :: setup
 procedure(bs_outward) :: outward
 procedure(bs_match) :: match
 procedure(bs_box) :: box

!Input/output variables
 real(dp), intent(inout) :: ee,emin,emax

!Output variables
 integer, intent(out) :: mch,ierr
 integer, intent(out) :: ifail

!Local variables
 real(dp) :: de,deprev,etop,elow
 integer :: nint,node,nin,lo,hi
 logical :: polish

 ifail=0
 deprev=huge(1.0d0)
 polish=.false.
 ierr=ierr_iterations
 mch=0

! top of the bracket: no bound state lies above it on the mesh;
! highest energy found to have too few nodes
 etop=emax
 elow=etop

! return point for bound state convergence
 do nint=1,60

! a bracket closed without a converged energy shift: there is no bound
! state, which an occupied state can lack in an early, overshooting SCF
! potential; the box state stands in for it. Once a shift is below eps
! the passes polish a converging energy, whose alternating shifts close
! the bracket too
   if(.not. polish .and. emax-emin<max(abs(ee),eig_floor)*eps) then
     call box(elow,etop,ee,ifail)
     if(ifail/=0) return
     mch=mmax
     ierr=ierr_box_state
     exit
   end if

   call setup(ee,mch)

   if(mch==-1) then
! classically forbidden everywhere: too low
     mch=0
     emin=ee
     elow=ee
     ee=0.5d0*(emin+emax)
     cycle
   else if(mch==0) then
! classically allowed at rr(mmax): too high
     emax=ee
     ee=0.5d0*(emin+emax)
     cycle
   end if

! the outward solution, carried through the matching window
! mch-1-match_k..mch+match_k (match_m), and its nodes inside the turning
! point
   lo=mch-1-match_k
   hi=mch+match_k
   call outward(mch,min(hi,mmax),node,ifail)
   if(ifail/=0) return

! start inward integration with simple exponential where the
! wave function has decayed to machine precision
   call inward_start(ll,ee,rr,vv,mch,mmax,adams,nin)

   if(node-nn+ll+1==0 .and. nin<=hi+4) then
! no room for the inward integration: too high
     emax=ee
     ee=0.5d0*(emin+emax)

   else if(node-nn+ll+1==0) then

     call match(ee,mch,nin,lo,hi,de,ifail)
     if(ifail/=0) return

! convergence test and possible exit
     if(bs_converged(de,deprev,ee,eps)) then
       ierr=0
       exit
     end if
! below eps but not converged: the energy is being polished
     if(abs(de)<max(abs(ee),eig_floor)*eps) polish=.true.
     deprev=de

     call bs_correct(de,ee,emin,emax,.false.)

   else if(node-nn+ll+1<0) then
! too few nodes
     elow=ee
     call bs_bisect(ee,emin,emax,.true.,.false.)

   else
! too many nodes
     call bs_bisect(ee,emin,emax,.false.,.false.)
   end if

 end do

 return
 end subroutine bs_solve

 end module boundsolve_m

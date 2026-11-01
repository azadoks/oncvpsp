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
!> Unit test of the transferability of the pseudopotential: the excitation
!> energies of the pseudo-atom follow the all-electron ones to second order
!> in the change of the occupations. By Janak's theorem the first order of
!> an excitation energy is sum_i eps_i df_i, for the all-electron atom and
!> for the pseudo-atom alike, so the excitation error de_ae - de_ps is
!>   sum_i (eps_ae,i - eps_ps,i) df_i + O(df**2)
!> with the reference configuration's eigenvalues; what is left, R, is
!> second order (the two atoms' hardnesses differ: the frozen core, the
!> model core and the pseudization), and R(d) = 4 R(d/2) + O(d**3).
!> Checked for two moves of d = 1/64 and d/2 (the test configurations of
!> the run, set here): the highest valence state ionized by d, and d moved
!> to it from the highest valence state of another l (if there is one):
!>   |R(d) - 4 R(d/2)| <= 0.1 (|R(d)| + 3e-9 Ha)
!> The cubic term is (d/2) times the relative change of the hardness per
!> electron, ~1 for valence shells (3 at most in the cases: R/d**2 changes
!> by 1-12% from d = 0.01 to 0.04), so 0.1 is 5-20 times it at this d; the
!> 3e-10 Ha it allows besides is twice the noise of R(d) - 4 R(d/2), up to
!> 1.4e-10 Ha (Nd), from the total energies (where the hardnesses nearly
!> agree, R(d/2) is as small as 1.6e-10 Ha: Li fr, VWN). A first-order
!> term left over would make the ratio 2.
!> The functional is the test's (VWN or PBE): PZ's (iexc 3) energy steps
!> at r_s = 1 (its coefficients are rounded), so its total energies on the
!> mesh jump by ~1e-6 Ha wherever the radius of r_s = 1 crosses a node (O
!> at d = 0.0107 on amesh 1.006; VWN is smooth there). irel 0 and 2 only:
!> the sr all-electron energy is not the variational energy of its
!> eigenvalues (test_ae_physics), so its excitation error has a first-order
!> term (Si's R is linear in d)
!> (doc/physics_tests_plan.md, Step 5).
!> Stops with a non-zero code if any check failed.
!>
!> usage: test_transfer_physics IREL INPUT.dat IEXC

 program test_transfer_physics
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: stderr => error_unit
 use unit_check_m, only: check, check_done
 use input_m, only: input_t
 use input_text_m, only: read_input_text
 use oncvpsp_run_m, only: oncvpsp_run
 use failure_m, only: failure_t
 use results_m, only: results_t
 implicit none

! 1/64, so that the occupations of a configuration sum exactly (check_data
! refuses a sum above Z by round-off as a negative ion)
 real(dp), parameter :: delta=0.015625d0,tol_ratio=0.1d0,noise=3.0d-9
 type(results_t) :: psp
 type(failure_t) :: at
 type(input_t) :: inp
 integer :: irel,iout,ifail,unit,inline,nc,nv,itop,ifrom,jj,kk,nmove
 real(dp) :: rr(5)
 character(len=1024) :: arg,input_filename

 call get_command_argument(1,arg)
 read(arg,*) irel
 call get_command_argument(2,input_filename)
 open(newunit=iout,file='/dev/null',action='write')
 open(newunit=unit,file=trim(input_filename),status='old',action='read')
 call read_input_text(unit,inline,inp,iout,ifail)
 close(unit)
 if(ifail/=0) error stop 'test_transfer_physics: read_input_text failed'
 call get_command_argument(3,arg)
 read(arg,*) inp%iexc

! the highest valence state, and the highest of another l below it
 nc=inp%nc
 nv=inp%nv
 itop=nc+nv
 ifrom=0
 do kk=nc+nv-1,nc+1,-1
   if(inp%la(kk)/=inp%la(itop) .and. inp%fa(kk)>=delta) then
     ifrom=kk
     exit
   end if
 end do
 nmove=1
 if(ifrom>0) nmove=2

! the test configurations: ionized by delta and delta/2, then (if there
! is a state to move from) delta and delta/2 moved to the highest state
 inp%ncnf=2*nmove
 do jj=2,inp%ncnf+1
   inp%cnf(jj)%nv=nv
   inp%cnf(jj)%na(nc+1:nc+nv)=inp%na(nc+1:nc+nv)
   inp%cnf(jj)%la(nc+1:nc+nv)=inp%la(nc+1:nc+nv)
   inp%cnf(jj)%fa(nc+1:nc+nv)=inp%fa(nc+1:nc+nv)
 end do
 inp%cnf(2)%fa(itop)=inp%fa(itop)-delta
 inp%cnf(3)%fa(itop)=inp%fa(itop)-0.5d0*delta
 if(nmove==2) then
   inp%cnf(4)%fa(ifrom)=inp%fa(ifrom)-delta
   inp%cnf(4)%fa(itop)=inp%fa(itop)+delta
   inp%cnf(5)%fa(ifrom)=inp%fa(ifrom)-0.5d0*delta
   inp%cnf(5)%fa(itop)=inp%fa(itop)+0.5d0*delta
 end if

 call oncvpsp_run(irel,0,inp,'','',.false.,psp,at,iout,ifail)
 if(ifail/=0) error stop 'test_transfer_physics: oncvpsp_run failed'
 close(iout)

 do jj=2,inp%ncnf+1
   rr(jj)=remainder(jj)
 end do
 call check_move('ionization',rr(2),rr(3))
 if(nmove==2) call check_move('transfer',rr(4),rr(5))
 call check_done('transfer_physics')

 contains

 real(dp) function remainder(jj)
!! the excitation error of configuration jj less its first order, the
!! occupations' changes times the reference's eigenvalue errors (by state
!! and kappa, as the run's eigenvalue table)
 integer, intent(in) :: jj
 integer :: kk,mm
 associate(cf=>psp%tests%configs)
   if(cf%nrow(jj)<0) error stop 'test_transfer_physics: a configuration failed'
   remainder=cf%de(3,jj)
   do kk=1,cf%nrow(jj)
     do mm=1,cf%nrow(1)
       if(all(cf%state(:,kk,jj)==cf%state(:,mm,1))) then
         remainder=remainder-(cf%ee(1,kk,jj)-cf%ee(1,mm,1)) &
&                  *(cf%ee(2,mm,1)-cf%ee(3,mm,1))
       end if
     end do
   end do
 end associate
 end function remainder

 subroutine check_move(what,r1,r2)
!! the remainders r1 of delta and r2 of delta/2 in the ratio of second
!! order: |r1 - 4 r2| relative to |r1|, with the energies' noise
 character(len=*), intent(in) :: what
 real(dp), intent(in) :: r1,r2
 write(stderr,'(2a,2es14.5)') what,': R(d), R(d/2) (Ha)',r1,r2
 call check(what//': |R(d) - 4 R(d/2)|/(|R(d)| + 3e-9)', &
&           abs(r1-4.0d0*r2)/(abs(r1)+noise),tol_ratio)
 end subroutine check_move

 end program test_transfer_physics

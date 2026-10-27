! Unit test of oncvpsp_run called from several threads at once (OpenMP):
! every case run alone, then all of them twice over concurrently, whose
! results must be those of the run alone bit for bit, so that the library
! keeps no state between or across runs. The cases differ in the
! relativistic treatment and the speed of light (one with c_light set),
! which the runs must not share. Arguments: the cases, irel:input:c_light
! (c_light 0 for the input's own). Stops with a non-zero code at the end if
! any check failed.
 program test_run_threads
 use, intrinsic :: iso_fortran_env, only: dp => real64, stderr => error_unit
 use input_m, only: input_t
 use input_text_m, only: read_input_text
 use oncvpsp_run_m, only: oncvpsp_run
 use results_m, only: results_t
 implicit none

!> the copies of each case run concurrently
 integer, parameter :: ncopy=2

!Local variables
 type(results_t), allocatable :: alone(:),together(:)
 integer, allocatable :: irel(:)
 real(dp), allocatable :: c_light(:)
 character(len=1024), allocatable :: input_filename(:)
 integer :: ncase,icase,itask,nfail,i1,i2
 character(len=1024) :: arg

 nfail=0
 ncase=command_argument_count()
 allocate(irel(ncase),c_light(ncase),input_filename(ncase))
 do icase=1,ncase
   call get_command_argument(icase,arg)
   i1=index(arg,':')
   i2=index(arg,':',back=.true.)
   read(arg(1:i1-1),*) irel(icase)
   input_filename(icase)=arg(i1+1:i2-1)
   read(arg(i2+1:),*) c_light(icase)
 end do

!each case alone
 allocate(alone(ncase),together(ncase*ncopy))
 do icase=1,ncase
   call run(icase,alone(icase))
 end do

!all of them at once, each ncopy times
!$omp parallel do schedule(static,1) num_threads(ncase*ncopy)
 do itask=1,ncase*ncopy
   call run(mod(itask-1,ncase)+1,together(itask))
 end do
!$omp end parallel do

 do itask=1,ncase*ncopy
   call check_same(alone(mod(itask-1,ncase)+1),together(itask), &
&                  trim(arg_of(mod(itask-1,ncase)+1)))
 end do

 if(nfail/=0) then
   write(stderr,'(i0,a)') nfail,' checks failed'
   error stop 1
 end if
 write(*,'(a,i0,a,i0,a)') 'test_run_threads: ',ncase*ncopy, &
&     ' concurrent runs of ',ncase,' cases match the runs alone'

 contains

 subroutine run(icase,psp)
!! reads the input of a case and runs its calculation and tests, the text
!! to /dev/null
!!
!! @param icase  the case
!! @param psp  the results
 integer, intent(in) :: icase
 type(results_t), intent(out) :: psp
 type(input_t) :: inp
 integer :: unit,iout,inline,ifail
 open(newunit=iout,file='/dev/null',action='write')
 open(newunit=unit,file=trim(input_filename(icase)),status='old', &
&     action='read')
 call read_input_text(unit,inline,inp,iout,ifail)
 close(unit)
 if(ifail/=0) error stop 'test_run_threads: read_input_text failed'
 if(c_light(icase)>0.0d0) inp%c_light=c_light(icase)
 call oncvpsp_run(irel(icase),0,inp,'','',.false.,psp,iout,ifail)
 if(ifail/=0) error stop 'test_run_threads: oncvpsp_run failed'
 close(iout)
 end subroutine run

 function arg_of(icase) result(text)
!! the argument of a case, for the messages
!!
!! @param icase  the case
 integer, intent(in) :: icase
 character(len=1024) :: text
 call get_command_argument(icase,text)
 end function arg_of

 subroutine check_same(pa,pb,what)
!! the results of a run alone and of a concurrent one, bit for bit
!!
!! @param pa  those of the run alone
!! @param pb  those of the concurrent run
!! @param what  the case, for the messages
 type(results_t), intent(in) :: pa,pb
 character(len=*), intent(in) :: what
 call differ(pa%pot%c_light/=pb%pot%c_light,what,'c_light')
 call differ(any(pa%pot%evkb/=pb%pot%evkb),what,'evkb')
 call differ(any(pa%pot%vkb/=pb%pot%vkb),what,'vkb')
 call differ(any(pa%pot%vpuns/=pb%pot%vpuns),what,'vpuns')
 call differ(any(pa%pot%rhomod/=pb%pot%rhomod),what,'rhomod')
 call differ(any(pa%tests%logder%ae/=pb%tests%logder%ae),what, &
&            'log derivatives ae')
 call differ(any(pa%tests%logder%ps/=pb%tests%logder%ps),what, &
&            'log derivatives ps')
 call differ(any(pa%tests%ghosts%e/=pb%tests%ghosts%e),what,'ghosts')
 call differ(any(pa%tests%configs%e/=pb%tests%configs%e),what, &
&            'test eigenvalues')
 call differ(any(pa%tests%configs%de/=pb%tests%configs%de),what, &
&            'test energies')
 call differ(any(pa%tests%convergence%curve/= &
&                pb%tests%convergence%curve),what,'curves')
 call differ(any(pa%tests%vloc_cut%cut/=pb%tests%vloc_cut%cut),what, &
&            'vloc cutoffs')
 end subroutine check_same

 subroutine differ(dd,what,result)
!! counts a failed check, and writes its message, if the runs differ
!!
!! @param dd  whether they differ
!! @param what  the case
!! @param result  the result
 logical, intent(in) :: dd
 character(len=*), intent(in) :: what,result
 if(.not. dd) return
 nfail=nfail+1
 write(stderr,'(4a)') 'FAIL: ',what,': concurrent run differs in ',result
 end subroutine differ

 end program test_run_threads

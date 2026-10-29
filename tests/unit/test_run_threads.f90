! Unit test of oncvpsp_run called from several threads at once (OpenMP):
! every case run alone, then all of them twice over concurrently, whose
! results and text output must be those of the run alone bit for bit, so
! that the library keeps no state between or across runs and the runs'
! text does not mix. The cases differ in the relativistic treatment, the
! speed of light (one with c_light set) and the functional (Libxc's
! among them), which the runs must not share. The inputs are read before
! the threads start, and each run writes to a scratch file of its own
! (gfortran does not open one file, /dev/null included, on two units at
! once). Arguments: the cases, irel:input:c_light (c_light 0 for the
! input's own). Stops with a non-zero code at the end if any check
! failed.
 program test_run_threads
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: stderr => error_unit
 use unit_check_m, only: check_done, nfail
 use input_m, only: input_t
 use input_text_m, only: read_input_text
 use oncvpsp_run_m, only: oncvpsp_run
 use failure_m, only: failure_t
 use results_m, only: results_t
 implicit none

!> the copies of each case run concurrently
 integer, parameter :: ncopy=2

!> the text output of a run, line by line
 type text_t
   character(len=2048), allocatable :: lines(:)
 end type text_t

!Local variables
 type(results_t), allocatable :: alone(:),together(:)
 type(text_t), allocatable :: text_alone(:),text_together(:)
 type(input_t), allocatable :: inps(:)
 integer, allocatable :: irel(:)
 real(dp), allocatable :: c_light(:)
 character(len=1024), allocatable :: input_filename(:)
 integer :: ncase,icase,itask,i1,i2,unit,iout,inline,ifail
 character(len=1024) :: arg

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

!the inputs, read once
 allocate(inps(ncase))
 open(newunit=iout,status='scratch',action='readwrite')
 do icase=1,ncase
   open(newunit=unit,file=trim(input_filename(icase)),status='old', &
&       action='read')
   call read_input_text(unit,inline,inps(icase),iout,ifail)
   close(unit)
   if(ifail/=0) error stop 'test_run_threads: read_input_text failed'
   if(c_light(icase)>0.0d0) inps(icase)%c_light=c_light(icase)
 end do
 close(iout)

!each case alone
 allocate(alone(ncase),together(ncase*ncopy))
 allocate(text_alone(ncase),text_together(ncase*ncopy))
 do icase=1,ncase
   call run(icase,alone(icase),text_alone(icase))
 end do

!all of them at once, each ncopy times
!$omp parallel do schedule(static,1) num_threads(ncase*ncopy)
 do itask=1,ncase*ncopy
   call run(mod(itask-1,ncase)+1,together(itask),text_together(itask))
 end do
!$omp end parallel do

 do itask=1,ncase*ncopy
   call check_same(alone(mod(itask-1,ncase)+1),together(itask), &
&                  trim(arg_of(mod(itask-1,ncase)+1)))
   call check_text(text_alone(mod(itask-1,ncase)+1),text_together(itask), &
&                  trim(arg_of(mod(itask-1,ncase)+1)))
 end do

 call check_done('run_threads')
 write(*,'(a,i0,a,i0,a)') 'test_run_threads: ',ncase*ncopy, &
&     ' concurrent runs of ',ncase,' cases match the runs alone'

 contains

 subroutine run(icase,psp,text)
!! runs the calculation and tests of a case on a copy of its input, the
!! text to a scratch file of its own, read back
!!
!! @param icase  the case
!! @param psp  the results
!! @param text  the text output
 integer, intent(in) :: icase
 type(results_t), intent(out) :: psp
 type(text_t), intent(out) :: text
 type(input_t) :: inp
 type(failure_t) :: at
 integer :: iout,ifail,nline,ios
 character(len=2048) :: line
 inp=inps(icase)
 open(newunit=iout,status='scratch',action='readwrite')
 call oncvpsp_run(irel(icase),0,inp,'','',.false.,psp,at,iout,ifail)
 if(ifail/=0) error stop 'test_run_threads: oncvpsp_run failed'
 rewind(iout)
 nline=0
 do
   read(iout,'(a)',iostat=ios) line
   if(ios/=0) exit
   nline=nline+1
 end do
 allocate(text%lines(nline))
 rewind(iout)
 do nline=1,size(text%lines)
   read(iout,'(a)') text%lines(nline)
 end do
 close(iout)
 end subroutine run

 subroutine check_text(ta,tb,what)
!! the text of a run alone and of a concurrent one, line by line
!!
!! @param ta  that of the run alone
!! @param tb  that of the concurrent run
!! @param what  the case, for the messages
 type(text_t), intent(in) :: ta,tb
 character(len=*), intent(in) :: what
 integer :: ii
 if(size(ta%lines)/=size(tb%lines)) then
   call differ(.true.,what,'the number of text lines')
   return
 end if
 do ii=1,size(ta%lines)
   if(ta%lines(ii)/=tb%lines(ii)) then
     call differ(.true.,what,'text line '//trim(ta%lines(ii)))
     return
   end if
 end do
 end subroutine check_text

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
 call differ(any(pa%tests%ghosts%ee/=pb%tests%ghosts%ee),what,'ghosts')
 call differ(any(pa%tests%configs%ee/=pb%tests%configs%ee),what, &
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

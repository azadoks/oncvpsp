! Unit tests of oncvpsp_run called as from the library: the test results it
! returns (test_results_m) against the lines of the DIAGNOSTICS block it
! writes, which print them to 17 digits, so exactly; and a second run in
! the same process, writing to /dev/null, whose results must be the same
! bit for bit. The ghost lists, the ghost test's basis and the model
! core's parameters against the lines of the text that print them, to
! their printed digits; and the warnings, as many of each as the text
! writes. Arguments: irel (0, 1 or 2) and the input file (text).
! Stops with a non-zero code at the end if any check failed.
 program test_run
 use constants_m, only: dp
 use, intrinsic :: iso_fortran_env, only: stderr => error_unit
 use unit_check_m, only: check_done, nfail
 use input_m, only: input_t
 use input_text_m, only: read_input_text
 use input_echo_m, only: fe
 use kappa_m, only: kap_count
 use oncvpsp_run_m, only: oncvpsp_run
 use failure_m, only: failure_t, warnings_t, mxwarn
 use results_m, only: results_t
 use test_results_m, only: tests_t
 implicit none

!Local variables
 type(results_t) :: psp1,psp2
 integer :: irel,iout,ifail,nline
 character(len=1024) :: arg,input_filename

 call get_command_argument(1,arg)
 read(arg,*) irel
 call get_command_argument(2,input_filename)

!the first run, its text to a scratch file, read back below
 open(newunit=iout,status='scratch',action='readwrite')
 call run(psp1,iout)
 rewind(iout)
 call check_diagnostics(psp1%tests,psp1%pot%nkap,iout)
 rewind(iout)
 call check_text(psp1,iout)
 rewind(iout)
 call check_warnings(psp1%warnings,iout)
 close(iout)
 if(nline==0) call fail('no DIAGNOSTICS lines checked')

!the second, its text to /dev/null
 open(newunit=iout,file='/dev/null',action='write')
 call run(psp2,iout)
 close(iout)
 call check_same(psp1,psp2)

 call check_done('run')
 write(*,'(a,i0,a)') 'test_run: ',nline,' DIAGNOSTICS lines match'

 contains

 subroutine run(psp,iout)
!! reads the input and runs the calculation and its tests
!!
!! @param psp  the results
!! @param iout  unit of the text output
 type(results_t), intent(out) :: psp
 integer, intent(in) :: iout
 type(input_t) :: inp
 type(failure_t) :: at
 integer :: unit,inline
 open(newunit=unit,file=trim(input_filename),status='old',action='read')
 call read_input_text(unit,inline,inp,iout,ifail)
 close(unit)
 if(ifail/=0) error stop 'test_run: read_input_text failed'
 call oncvpsp_run(irel,0,inp,'','',.false.,psp,at,iout,ifail)
 if(ifail/=0) error stop 'test_run: oncvpsp_run failed'
 end subroutine run

 subroutine check_diagnostics(tt,nkap,iu)
!! each line of the DIAGNOSTICS block that a test result prints against
!! that result
!!
!! @param tt  the results of the tests
!! @param nkap  number of kappa per l
!! @param iu  unit of the text output, rewound
 type(tests_t), intent(in) :: tt
 integer, intent(in) :: nkap
 integer, intent(in) :: iu
 character(len=512) :: line
 character(len=32) :: tag
 integer :: ios,ll,ikap,nn,iprj,icfg,jj,irow,ild,ighost(4,2),ivl
 real(dp) :: xx(4)
 logical :: inblock
 inblock=.false.
 nline=0
 ild=0
 ighost(:,:)=0
 irow=0
 ivl=0
 do
   read(iu,'(a)',iostat=ios) line
   if(ios/=0) exit
   if(line=='Begin DIAGNOSTICS') inblock=.true.
   if(line=='End DIAGNOSTICS') inblock=.false.
   if(.not. inblock .or. line(1:1)=='#' .or. len_trim(line)==0) cycle
   read(line,*) tag
   select case(trim(tag))
   case('ld_r')
     read(line,*) tag,ll,ikap,xx(1)
     call same(xx(1),tt%logder%rl(ll+1),line)
     ild=0
   case('ld')
     read(line,*) tag,ll,ikap,xx(1:3)
     ild=ild+1
     call same(xx(1),tt%logder%ee(ild),line)
     call same(xx(2),tt%logder%ae(ild,ll+1,ikap),line)
     call same(xx(3),tt%logder%ps(ild,ll+1,ikap),line)
   case('ghost')
     read(line,*) tag,ll,ikap,nn,xx(1:2)
     ighost(ll+1,ikap)=ighost(ll+1,ikap)+1
     jj=ighost(ll+1,ikap)
     if(nn/=tt%ghosts%nn(jj,ll+1,ikap)) call fail(line)
     call same(xx(1),fe(tt%ghosts%ee(1,jj,ll+1,ikap)),line)
     call same(xx(2),fe(tt%ghosts%ee(2,jj,ll+1,ikap)),line)
   case('ghost_r')
     read(line,*) tag,ll,ikap,xx(1:2)
     call same(xx(1),tt%ghosts%least(1,ll+1,ikap),line)
     call same(xx(2),tt%ghosts%least(2,ll+1,ikap),line)
   case('test_config')
     read(line,*) tag,icfg,xx(1:3)
     irow=0
     do jj=1,3
       call same(xx(jj),fe(tt%configs%de(jj,icfg+1)),line)
     end do
   case('test_eig')
     read(line,*) tag,icfg,nn,ll,ikap,xx(1:3)
     irow=irow+1
     if(any([nn,ll,ikap]/=tt%configs%state(:,irow,icfg+1))) call fail(line)
     do jj=1,3
       call same(xx(jj),fe(tt%configs%ee(jj,irow,icfg+1)),line)
     end do
   case('ecut_vloc')
     read(line,*) tag,xx(1:2)
     do
       ivl=ivl+1
       if(tt%vloc_cut%cut(1,ivl)/=0.0d0) exit
     end do
     call same(xx(1),fe(tt%vloc_cut%cut(2,ivl)),line)
     call same(xx(2),tt%vloc_cut%cut(1,ivl),line)
   case('vloc_check')
     read(line,*) tag,xx(1)
     call same(xx(1),fe(tt%vloc_cut%check),line)
   case('ekin_tail','ekin_check')
     read(line,*) tag,ll,ikap,iprj,xx(1)
     jj=merge(1,2,trim(tag)=='ekin_tail')
     call same(xx(1),fe(tt%convergence%tail(jj,iprj,ll+1,ikap)),line)
   case('opt')
     read(line,*) tag,ll,ikap,iprj,xx(1:4)
     do jj=1,4
       call same(xx(jj),fe(tt%convergence%opt(jj,iprj,ll+1,ikap)),line)
     end do
   case('ecut')
     read(line,*) tag,ll,ikap,iprj,xx(1:2)
     if(.not. any(tt%convergence%cross(1,:,iprj,ll+1,ikap)==xx(2) .and. &
&       fe(tt%convergence%cross(2,:,iprj,ll+1,ikap))==xx(1))) then
       call fail(line)
     end if
   case('ecut_w')
     read(line,*) tag,xx(1:2)
     if(.not. tt%weighted_cut%used) call fail(line)
     call same(xx(1),fe(tt%weighted_cut%errw),line)
     call same(xx(2),tt%weighted_cut%ecutw,line)
   case('ecut_w_weight')
     read(line,*) tag,ll,ikap,iprj,xx(1)
     call same(xx(1),fe(tt%weighted_cut%weight(iprj,ll+1,ikap)),line)
   case('vkb_b')
     read(line,*) tag,ll,ikap,iprj,xx(1)
     call same(xx(1),fe(tt%bbeig(iprj,ll+1,ikap)),line)
   case default
     cycle
   end select
   nline=nline+1
 end do
!every row of the tables
 do jj=1,tt%configs%ncfg
   if(tt%configs%nrow(jj)<=0) call fail('a test configuration failed')
 end do
 do ll=0,3
   do ikap=1,kap_count(nkap,ll)
     if(ighost(ll+1,ikap)/=tt%ghosts%nrow(ll+1,ikap)) call fail('ghost rows')
   end do
 end do
 if(size(tt%convergence%curve,2)==0) call fail('no convergence curve')
 end subroutine check_diagnostics

 subroutine check_text(psp,iu)
!! the lines of the text outside the DIAGNOSTICS block that print the
!! GHOST(-) and GHOST(+) states, the ghost test's basis notes and the
!! model core's parameters, against the results to the printed digits
!!
!! @param psp  the results
!! @param iu  unit of the text output, rewound
 type(results_t), intent(in) :: psp
 integer, intent(in) :: iu
 character(len=512) :: line
 integer :: ios,ll,lbar,ikap,jj,nminus(4,2),nplus(4,2),nbasis(4,2)
 real(dp) :: xx(3)
 ikap=1
 lbar=0
 nminus(:,:)=0
 nplus(:,:)=0
 nbasis(:,:)=0
 associate(gh=>psp%tests%ghosts,bb=>psp%core%budget,tt=>psp%core%teter)
   do
     read(iu,'(a)',iostat=ios) line
     if(ios/=0) exit
     if(index(line,'Ghost test for J = L - 1/2')>0) ikap=2
     if(index(line,'WARNING - GHOST(-)')>0) then
       read(line,*) ll,xx(1:2)
       nminus(ll+1,ikap)=nminus(ll+1,ikap)+1
       jj=nminus(ll+1,ikap)
       call near(xx(1),gh%minus(jj,ll+1,ikap),1.0d-6,line)
       call near(xx(2),gh%ecut(ll+1,ikap),1.0d-2,line)
     else if(index(line,'WARNING - GHOST(+)')>0) then
       read(line,*) ll,xx(1:3)
       nplus(ll+1,ikap)=nplus(ll+1,ikap)+1
       jj=nplus(ll+1,ikap)
       call near(xx(1),gh%plus(2,jj,ll+1,ikap),1.0d-6,line)
       call near(xx(2),gh%plus(1,jj,ll+1,ikap),1.0d-6,line)
       call near(xx(3),gh%ecut(ll+1,ikap),1.0d-2,line)
     else if(index(line,'run_ghosts: WARNING barrier state')>0) then
       lbar=nint(num_after(line,' l='))
     else if(index(line,'run_ghosts: WARNING basis cut there')>0) then
       nbasis(lbar+1,ikap)=1
       call near(num_after(line,'ghost test to'),gh%ecut(lbar+1,ikap), &
&                1.0d-6,line)
     else if(index(line,'run_ghosts: WARNING basis full')>0) then
       nbasis(nint(num_after(line,' l='))+1,ikap)=2
     else if(index(line,'rmatch, rhocmatch')==1) then
       read(line(18:),*) xx(1:2)
       call near(xx(1),tt(3),1.0d-4,line)
       call near(xx(2),tt(4),1.0d-4,line)
     else if(line=='amplitude prefactor, scale prefactor') then
       read(iu,'(a)') line
       read(line,*) xx(1:2)
       call near(xx(1),tt(1),1.0d-4,line)
       call near(xx(2),tt(2),1.0d-4,line)
     else if(index(line,'  cutoff budget')==1) then
       call near(num_after(line,'budget'),bb%ecut,1.0d-4,line)
       call near(num_after(line,'error of'),bb%error,1.0d-2*bb%error,line)
     else if(index(line,'  slack')==1) then
       call near(num_after(line,'slack'),bb%slack,1.0d-4,line)
     else if(index(line,'  model radius r_s')==1) then
       call near(num_after(line,'r_s'),bb%r_s,1.0d-6,line)
       call near(num_after(line,'smallest rc'),bb%rc_min,1.0d-6,line)
     else if(index(line,'  last radius with core charge')==1) then
       call near(num_after(line,'charge'),bb%r_core,1.0d-6,line)
     else if(line=='  fill minimum tail') then
       if(bb%fill/=1) call fail(line)
     else if(line=='  fill Martins') then
       if(bb%fill/=2) call fail(line)
     else if(index(line,'  cutoff the model core needs')==1) then
       call near(num_after(line,'needs'),bb%ecut_needed,1.0d-4,line)
     else if(index(line,'  smallest fitting radius')==1) then
       call near(num_after(line,'(icmod 6)'),bb%r_s6,1.0d-6,line)
     else if(index(line,'  hardness mismatch')==1) then
       call near(num_after(line,'mismatch'),bb%mismatch, &
&                1.0d-4*bb%mismatch,line)
     else
       cycle
     end if
     nline=nline+1
   end do
!every state of the lists printed (the GHOST(+) lines stop at 20 in all)
   if(any(nminus/=gh%nminus(:,1:2))) call fail('GHOST(-) count')
   if(sum(gh%nplus)<20 .and. any(nplus/=gh%nplus(:,1:2))) then
     call fail('GHOST(+) count')
   end if
   if(any(nbasis/=gh%basis(:,1:2))) call fail('ghost basis')
 end associate
 end subroutine check_text

 subroutine check_warnings(warns,iu)
!! the warnings of the run against the text: of each wrn_*, as many
!! records as lines of the text with its message, every one with a stage
!!
!! @param warns  the warnings
!! @param iu  unit of the text output, rewound
 type(warnings_t), intent(in) :: warns
 integer, intent(in) :: iu
!the line of each warning (wrn_* in order) that the text writes once
 character(len=40), parameter :: pat(29)=[character(len=40) :: &
&  'ae_converged: WARNING scf stopped', 'sratom: WARNING failed to converge', &
&  'relatom: WARNING failed to converge', &
&  'WARNING self-consistency failed to', 'runconfig: WARNING lschfb', &
&  'runconfig: WARNING ldiracfb', 'WARNING for fully non-local PS atom', &
&  'WARNING lschvkbb convergence error', 'potential failed to converge', &
&  'optimize: WARNING E_resid flat', 'WARNING the well''s tail changes it by', &
&  'WARNING qcut is not above', 'WARNING this well state''s free tail', &
&  'WARNING wellstate: localized bound', &
&  'WARNING wellstate: moderately localized', &
&  'WARNING wellstate: negative energy', 'lschvkbb ERROR', &
&  'lschvkbbe ERROR', 'run_ghosts: lschvkbb ERROR', 'WARNING - GHOST(-)', &
&  'WARNING - GHOST(+)', 'WARNING basis cut there', 'WARNING basis full', &
&  'write_waves: WARNING lschvkbb', 'sr_so_r: WARNING', &
&  'does not fit the budget', 'no minimum of the mismatch', &
&  'root at a jump', 'WARNING icmod 6 and 7 are deprecated']
 character(len=512) :: line
 integer :: ios,id,ntext(29),nrec(29),ii
 character(len=8) :: tag
 ntext(:)=0
 do
   read(iu,'(a)',iostat=ios) line
   if(ios/=0) exit
   do id=1,size(pat)
     if(index(line,trim(pat(id)))==0) cycle
!run_diag's lschvkbb ERROR, not run_ghosts'
     if(id==17 .and. index(line,'run_ghosts')>0) cycle
     ntext(id)=ntext(id)+1
   end do
 end do
 if(warns%nwarn>mxwarn) return
 nrec(:)=0
 do ii=1,warns%nwarn
   id=warns%rec(1,ii)
   if(id<1 .or. id>size(pat) .or. warns%rec(2,ii)==0) then
     call fail('a warning record without its warning or stage')
     cycle
   end if
   nrec(id)=nrec(id)+1
 end do
 do id=1,size(pat)
   if(nrec(id)/=ntext(id)) then
     write(tag,'(i0)') id
     call fail('warning '//trim(tag)//': '//trim(pat(id)))
   end if
   nline=nline+ntext(id)
 end do
 end subroutine check_warnings

 function num_after(line,key) result(xx)
!! the number after key in line
!!
!! @param line  the line
!! @param key  the text before it
 character(len=*), intent(in) :: line,key
 real(dp) :: xx
 integer :: ii,ios
 xx=huge(1.0d0)
 ii=index(line,key)
 if(ii==0) return
 read(line(ii+len(key):),*,iostat=ios) xx
 if(ios/=0) xx=huge(1.0d0)
 end function num_after

 subroutine near(xx,yy,dd,line)
!! fails unless the value read is the result to the printed digits
!!
!! @param xx  the value read
!! @param yy  the result
!! @param dd  the last printed digit's unit
!! @param line  the line, for the message
 real(dp), intent(in) :: xx,yy,dd
 character(len=*), intent(in) :: line
 if(abs(xx-yy)>0.5000001d0*dd) call fail(line)
 end subroutine near

 subroutine check_same(pa,pb)
!! the results of two runs, bit for bit
!!
!! @param pa  those of the first
!! @param pb  those of the second
 type(results_t), intent(in) :: pa,pb
 call differ(any(pa%pot%evkb/=pb%pot%evkb),'evkb')
 call differ(any(pa%pot%vpuns/=pb%pot%vpuns),'vpuns')
 call differ(any(pa%pot%rhomod/=pb%pot%rhomod),'rhomod')
 call differ(any(pa%tests%logder%ae/=pb%tests%logder%ae),'log derivatives ae')
 call differ(any(pa%tests%logder%ps/=pb%tests%logder%ps),'log derivatives ps')
 call differ(any(pa%tests%ghosts%ee/=pb%tests%ghosts%ee),'ghosts')
 call differ(any(pa%tests%configs%ee/=pb%tests%configs%ee),'test eigenvalues')
 call differ(any(pa%tests%configs%de/=pb%tests%configs%de),'test energies')
 call differ(any(pa%tests%convergence%curve/= &
&                pb%tests%convergence%curve),'curves')
 call differ(any(pa%tests%vloc_cut%cut/=pb%tests%vloc_cut%cut),'vloc cutoffs')
 call differ(any(pa%tests%ghosts%minus/=pb%tests%ghosts%minus),'GHOST(-)')
 call differ(any(pa%tests%ghosts%plus/=pb%tests%ghosts%plus),'GHOST(+)')
 call differ(any(pa%core%teter/=pb%core%teter),'Teter parameters')
 call differ(pa%core%budget%r_s/=pb%core%budget%r_s,'budget r_s')
 end subroutine check_same

 subroutine differ(dd,what)
!! fails if the two runs differ in what
!!
!! @param dd  whether they differ
!! @param what  the result, for the message
 logical, intent(in) :: dd
 character(len=*), intent(in) :: what
 if(dd) call fail('second run: '//what)
 end subroutine differ

 subroutine same(xx,yy,line)
!! fails unless the value read is the result
!!
!! @param xx  the value read
!! @param yy  the result
!! @param line  the line, for the message
 real(dp), intent(in) :: xx,yy
 character(len=*), intent(in) :: line
 if(xx/=yy) call fail(line)
 end subroutine same

 subroutine fail(msg)
!! counts a failed check and writes its message
!!
!! @param msg  the message
 character(len=*), intent(in) :: msg
 nfail=nfail+1
 write(stderr,'(2a)') 'FAIL: ',trim(msg)
 end subroutine fail

 end program test_run

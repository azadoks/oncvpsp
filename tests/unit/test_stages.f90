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
!> Unit test of the stages of oncvpsp_run called on their own: the
!> all-electron atom solved once (solve_ae), and every channel optimized
!> twice from it (optimize_channel), each time bit for bit the same, and
!> the same as the convergence data of the whole run; then the local
!> potential (local_potential) and every channel's projectors
!> (channel_vkb) from those channels alone, the same as the whole run's.
!>
!> usage: test_stages IREL INPUT.dat

 program test_stages
 use constants_m, only: dp
 use unit_check_m, only: check_done, nfail
 use input_m, only: input_t, config_t
 use input_text_m, only: read_input_text
 use oncvpsp_run_m, only: oncvpsp_run, complete_input
 use failure_m, only: failure_t
 use ae_state_m, only: ae_state_t, solve_ae
 use channel_m, only: channel_t, optimize_channel
 use run_vkb_m, only: local_t, local_potential, channel_vkb
 use pieces_m, only: piecewise_t
 use main_setup_m, only: core_radius_points
 use results_m, only: results_t
 use constants_m, only: mxcnf, mxlin, mxl, mxprj
 implicit none

 type(results_t) :: psp
 type(failure_t) :: at
 type(input_t) :: inp
 type(config_t) :: cnf(mxcnf)
 type(ae_state_t) :: ae
 type(channel_t) :: ch1,ch2,chs(mxlin+1)
 type(local_t) :: lp
 type(piecewise_t), allocatable :: vkbpw(:,:)
 real(dp), allocatable :: evkb(:,:),vkb(:,:,:),bbeig(:,:)
 integer :: irel,nkap,iout,ifail,ierr,l1,lc,nprj,nq,ncheck
 integer :: irc(mxlin+1)
 real(dp) :: rcmax
 character(len=1024) :: arg,input_filename

 call get_command_argument(1,arg)
 read(arg,*) irel
 call get_command_argument(2,input_filename)
 nkap=1
 if(irel==2) nkap=2

 open(newunit=iout,file='/dev/null',action='write')
 call read_input(inp)
 call oncvpsp_run(irel,0,inp,'','',.false.,psp,at,iout,ifail)
 if(ifail/=0) error stop 'test_stages: oncvpsp_run failed'

 call read_input(inp)
 call complete_input(inp,cnf,iout,ifail)
 if(ifail/=0) error stop 'test_stages: complete_input failed'
 call solve_ae(irel,nkap,0,inp,cnf,'',ae,ierr,iout,ifail)
 if(ifail/=0 .or. ierr/=0) error stop 'test_stages: solve_ae failed'
 call core_radius_points(inp%lmax,inp%lloc,inp%rc,ae%rr,ae%mmax,irc,rcmax, &
&                        iout,ifail)
 if(ifail/=0) error stop 'test_stages: core_radius_points failed'

 nq=size(psp%tests%convergence%curve,2)
 ncheck=0
 do l1=1,inp%lmax+1
   nprj=inp%nproj(l1)
   if(l1==inp%lloc+1) nprj=1
   call optimize_channel(ae,inp,l1,nprj,irc(l1),ch1,iout,ifail)
   if(ifail/=0) error stop 'test_stages: optimize_channel failed'
   call optimize_channel(ae,inp,l1,nprj,irc(l1),ch2,iout,ifail)
   if(ifail/=0) error stop 'test_stages: optimize_channel failed'
   call same('pswf',[ch1%pswf],[ch2%pswf])
   call same('vp',[ch1%vp],[ch2%vp])
   call same('vkb',[ch1%vkb],[ch2%vkb])
   call same('ecurve',[ch1%ecurve],[ch2%ecurve])
   call same('epa',[ch1%epa],[ch2%epa])
   if(ch1%irc/=ch2%irc .or. any(ch1%npa/=ch2%npa)) call fail('irc, npa')
! the whole run's
   call same('cross',[ch1%cvgplt],[psp%tests%convergence%cross(:,:,:,l1,:)])
   call same('curve',[ch1%ecurve(:,1:nq,:,:)], &
&            [             psp%tests%convergence%curve(:,:,:,l1,:)])
   call same('tail',[ch1%ektail],[psp%tests%convergence%tail(:,:,l1,:)])
   call same('opt',[ch1%optdg],[psp%tests%convergence%opt(:,:,l1,:)])
   call same('wshare',[ch1%wshare],[psp%tests%convergence%wshare(:,:,l1,:)])
   chs(l1)=ch1
   ncheck=ncheck+1
 end do

! the local potential, from the local channel for lloc<4, and the
! projectors of the others
 lc=min(inp%lloc+1,inp%lmax+1)
 call local_potential(inp%lloc,inp%lpopt,inp%dvloc0,inp%rc,irc,ae%rr, &
&                     ae%mmax,ae%vfull,nkap,chs(lc)%pswfx,chs(lc)%chix, &
&                     chs(lc)%vrloc,lp,iout,ifail)
 if(ifail/=0) error stop 'test_stages: local_potential failed'
 if(any(lp%vlocpw%vals/=psp%pot%vlocpw%vals)) then
   l1=0
   call fail('vlocpw')
 end if
 allocate(evkb(mxprj,nkap),vkb(ae%mmax,mxprj,nkap),bbeig(mxprj,nkap))
 allocate(vkbpw(mxprj,nkap))
 do l1=1,inp%lmax+1
   if(l1==inp%lloc+1 .or. l1>mxl+1) cycle
   evkb(:,:)=0.0d0
   vkb(:,:,:)=0.0d0
   call channel_vkb(l1,inp%lloc,lp,inp%rc(l1),chs(l1)%irc,inp%nproj(l1), &
&                   ae%rr,ae%mmax,ae%vfull,chs(l1)%pswfx,chs(l1)%chix, &
&                   chs(l1)%dlta,irel,nkap,evkb,vkb,vkbpw,bbeig,iout,ifail)
   if(ifail/=0) error stop 'test_stages: channel_vkb failed'
   call same('evkb',[evkb],[psp%pot%evkb(:,l1,:)])
   call same('vkb',[vkb],[psp%pot%vkb(:,:,l1,:)])
   call same('bbeig',[bbeig],[psp%tests%bbeig(:,l1,:)])
 end do
 close(iout)

 call check_done('stages')
 write(*,'(a,i0,a)') 'test_stages: ',ncheck,' channels match'

 contains

 subroutine read_input(inp)
 type(input_t), intent(out) :: inp
 integer :: unit,inline
 open(newunit=unit,file=trim(input_filename),status='old',action='read')
 call read_input_text(unit,inline,inp,iout,ifail)
 close(unit)
 if(ifail/=0) error stop 'test_stages: read_input_text failed'
 end subroutine read_input

 subroutine same(what,aa,bb)
 character(len=*), intent(in) :: what
 real(dp), intent(in) :: aa(:),bb(:)
 if(size(aa)/=size(bb)) then
   call fail(what//': size')
 else if(any(aa/=bb)) then
   call fail(what)
 end if
 end subroutine same

 subroutine fail(what)
 character(len=*), intent(in) :: what
 write(*,'(2a,i0)') 'test_stages: FAIL ',what,', l1=',l1
 nfail=nfail+1
 end subroutine fail

 end program test_stages

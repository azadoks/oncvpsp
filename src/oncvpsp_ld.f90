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
! Log derivatives, oncvpsp_ld.x: run_phsft's all-electron and
! pseudopotential phases atan(r u'/u - 1) at any radii and energies, from
! the all-electron data file of an _ae.x (ae_file_m) and the
! pseudopotential data file of a run with --ps-file (ps_file_m), for any
! relativistic treatment.  The phases at all the radii come from one
! solution per energy (phsft_m).  With the defaults (the radius and the
! energies of the run's input) it prints the log derivatives of the run's
! DIAGNOSTICS block; each radius adds a block of the same lines:
!   ld_r l ikap r(bohr)
!   ld l ikap e(Ha) phase_ae phase_ps
! for l=0..3 (or those of --channels) and each kappa of l, ikap 1 for
! j=l+1/2 and 2 for j=l-1/2.

 program oncvpsp_ld
 use constants_m, only: dp
 use kappa_m, only: kap_of
 use ps_file_m, only: ps_file_t, ps_file_free
 use phsft_m, only: phsft_pair
 use cli_m, only: cli_next_arg, cli_radii, cli_channels, cli_energies
 use ae_ps_pair_m, only: load_pair
 use build_info_m, only: write_build_info
 implicit none
 integer, parameter :: iout=6

!> the most radii of a run
 integer, parameter :: mxq=1000

 character(len=1024) :: ae_filename,ps_filename,arg
 integer :: iarg,irel,mmax,nq,npsh,ii,iq,l1,ll,ikap,kap,nvkb,lc
 real(dp) :: zz,epsh1,epsh2,depsh,amesh,rq(mxq),ee(3)
!> the speed of light of the files
 real(dp) :: c_light
 real(dp), allocatable :: rr(:),vfull(:),pshae(:,:),pshps(:,:),rql(:)
 logical :: given_radii,given_energies
!> the angular momenta whose log derivatives are computed
 logical :: chan(0:3)
 type(ps_file_t) :: pf
 integer :: ifail

 ae_filename=''
 ps_filename=''
 given_radii=.false.
 given_energies=.false.
 chan(:)=.true.
 nq=0
 iarg=1
 do while(iarg<=command_argument_count())
   call get_command_argument(iarg,arg)
   select case(arg)
   case('-h','--help')
     call write_help()
     stop
   case('-v','--version')
     call write_build_info(6,'oncvpsp_ld.x','any relativistic treatment')
     stop
   case('-a','--ae-file')
     call cli_next_arg('oncvpsp_ld',iarg,ae_filename)
   case('-p','--ps-file')
     call cli_next_arg('oncvpsp_ld',iarg,ps_filename)
   case('-r','--radii')
     call cli_next_arg('oncvpsp_ld',iarg,arg)
     call cli_radii('oncvpsp_ld',arg,mxq,nq,rq)
     given_radii=.true.
   case('-l','--channels')
     call cli_next_arg('oncvpsp_ld',iarg,arg)
     call cli_channels('oncvpsp_ld',arg,chan)
   case('-e','--energies')
     call cli_next_arg('oncvpsp_ld',iarg,arg)
     call cli_energies('oncvpsp_ld',arg,ee)
     given_energies=.true.
   case default
     write(iout,'(2a)') 'oncvpsp_ld: ERROR unknown argument ',trim(arg)
     error stop 1
   end select
   iarg=iarg+1
 end do
 if(len_trim(ae_filename)==0 .or. len_trim(ps_filename)==0) then
   write(iout,'(a)') 'oncvpsp_ld: ERROR requires --ae-file and --ps-file'
   error stop 1
 end if

 call load_pair('oncvpsp_ld',trim(ae_filename),trim(ps_filename),irel, &
&               c_light,zz,mmax,rr,vfull,pf,iout=iout,ifail=ifail)
 if(ifail/=0) error stop 1

!the energies, from epsh2 down by depsh, as oncvpsp_main
 if(given_energies) then
   epsh1=ee(1)
   epsh2=ee(2)
   depsh=ee(3)
 else
   epsh1=pf%epsh1
   epsh2=pf%epsh2
   depsh=pf%depsh
 end if
 npsh=int(((epsh2-epsh1)/depsh)-0.5d0)+1
 amesh=rr(2)/rr(1)

 write(iout,'(a)') '# ld_r l ikap r(bohr)'
 write(iout,'(a)') '# ld l ikap e(Ha) atan(r*u''/u)_ae atan(r*u''/u)_ps'
 do l1=1,4
   ll=l1-1
   if(.not. chan(ll)) cycle
!the input's radius of l, as run_phsft
   if(.not. given_radii) then
     nq=1
     if(pf%rxpsh>=0.0d0) then
       rq(1)=pf%rxpsh
     else if(ll<=pf%pot%lmax) then
       rq(1)=pf%rc(l1)*amesh**2
     else
       lc=pf%pot%lloc+1
       rq(1)=pf%rc(lc)
     end if
   end if
   rql=rq(1:nq)
   allocate(pshae(npsh,nq),pshps(npsh,nq))
   do ikap=1,merge(1,pf%pot%nkap,l1==1)
     kap=kap_of(ll,ikap)
     nvkb=pf%nvkb(l1,ikap)
     call phsft_pair(irel,c_light,ll,kap,nvkb,rr,vfull,zz,pf%vloc(:,ikap), &
&                    pf%pot%vlocpw,pf%pot%vkbpw(1:nvkb,l1,ikap), &
&                    pf%pot%evkb(1:nvkb,l1,ikap),mmax,nq,rql,npsh,epsh2,depsh, &
&                    pf%epa(l1,ikap),pshae,pshps,iout,ifail)
     if(ifail/=0) error stop 1
     do iq=1,nq
       write(iout,'(a,2i4,1p,e24.16)') 'ld_r',ll,ikap,rql(iq)
       do ii=1,npsh
         write(iout,'(a,2i4,1p,3e24.16)') 'ld',ll,ikap, &
&              epsh2-depsh*real(ii-1,dp),pshae(ii,iq),pshps(ii,iq)
       end do
     end do
   end do
   deallocate(pshae,pshps)
 end do

 call ps_file_free(pf)
 deallocate(rr,vfull)
 stop

 contains

 subroutine write_help()
!! the usage: what the program does and its options
 write(iout,'(a)') 'Usage: oncvpsp_ld.x -a AE_FILE -p PS_FILE [options]'
 write(iout,'(a)') ''
 write(iout,'(a)') &
&      'The log derivatives of a pseudopotential and its all-electron'
 write(iout,'(a)') &
&      'atom, atan(r u''/u - 1) at radii r over a grid of energies, as'
 write(iout,'(a)') 'continuous phases (the pseudopotential''s aligned with the'
 write(iout,'(a)') 'all-electron ones at the reference energy of each l), for'
 write(iout,'(a)') 'l = 0..3 (or those of -l) and, fully relativistic, each j.'
 write(iout,'(a)') 'All the radii come from one solution per energy. Without'
 write(iout,'(a)') '-r, -l or -e it prints the log derivatives of the run''s'
 write(iout,'(a)') 'DIAGNOSTICS block. The output, on stdout, is a block per'
 write(iout,'(a)') 'l, j and radius:'
 write(iout,'(a)') '  ld_r l ikap r(bohr)'
 write(iout,'(a)') '  ld l ikap e(Ha) phase_ae phase_ps'
 write(iout,'(a)') &
&      'ikap 1 for j = l+1/2 (and every l if not fully relativistic),'
 write(iout,'(a)') '2 for j = l-1/2.'
 write(iout,'(a)') ''
 write(iout,'(a)') 'Options:'
 write(iout,'(a)') '  -h, --help'
 write(iout,'(a)') '      Show this help and exit'
 write(iout,'(a)') '  -v, --version'
 write(iout,'(a)') '      Show the build and exit'
 write(iout,'(a)') '  -a, --ae-file FILE'
 write(iout,'(a)') '      The all-electron data file of an _ae.x (required)'
 write(iout,'(a)') '  -p, --ps-file FILE'
 write(iout,'(a)') '      The pseudopotential data file of a run of oncvpsp.x,'
 write(iout,'(a)') &
&      '      oncvpspnr.x, oncvpspr.x or their _ps.x with --ps-file,'
 write(iout,'(a)') '      for the same atom, treatment and mesh (required)'
 write(iout,'(a)') '  -r, --radii R1,R2,...'
 write(iout,'(a)') &
&      '      The radii (bohr), each for every l. The two phases are'
 write(iout,'(a)') '      comparable only beyond rc of that l and of the local'
 write(iout,'(a)') '      potential; inside, they differ by construction.'
 write(iout,'(a)') '      Default the input''s rxpsh, or for rxpsh < 0'
 write(iout,'(a)') '      rc(l)*amesh**2 (rc of the local potential above lmax)'
 write(iout,'(a)') '  -l, --channels L1,L2,...'
 write(iout,'(a)') '      The angular momenta, from 0 to 3, whose log'
 write(iout,'(a)') '      derivatives are computed (both j of each when fully'
 write(iout,'(a)') '      relativistic). Default all four'
 write(iout,'(a)') '  -e, --energies EMIN,EMAX,DE'
 write(iout,'(a)') &
&      '      The energies (Ha), from EMAX down by DE to above EMIN.'
 write(iout,'(a)') '      Default the input''s epsh1,epsh2,depsh'
 write(iout,'(a)') ''
 write(iout,'(a)') 'Errors are reported on stdout, with a nonzero exit status.'
 end subroutine write_help

 end program oncvpsp_ld

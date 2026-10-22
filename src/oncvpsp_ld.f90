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
 use ae_file_m, only: read_ae_potential
 use ps_file_m, only: ps_file_t, read_ps_file, ps_file_free
 use phsft_m, only: phsft_pair
 use relativity_m, only: set_relativity
 use build_info_m, only: write_build_info
 implicit none

!> relative tolerance of the comparison of the two files' meshes and
!> constants, as in ae_file_m
 real(dp), parameter :: rtol=1.0d-12
!> the most radii of a run
 integer, parameter :: mxq=1000

 character(len=1024) :: ae_filename,ps_filename,arg
 character(len=2) :: atsym
 character(len=19) :: treatment
 integer :: iarg,irel,mmax,nq,npsh,ii,iq,l1,ll,ikap,kap,nvkb,ios,lc,nl
 integer :: lsel(4)
 real(dp) :: zz,c_light_f,epsh1,epsh2,depsh,amesh,rq(mxq),ee(3)
 real(dp), allocatable :: rr(:),vfull(:),pshae(:,:),pshps(:,:),rql(:)
 logical :: given_radii,given_energies
!> the angular momenta whose log derivatives are computed
 logical :: chan(0:3)
 type(ps_file_t) :: pf

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
     call next_arg(ae_filename)
   case('-p','--ps-file')
     call next_arg(ps_filename)
   case('-r','--radii')
     call next_arg(arg)
     nq=count_list(arg)
     if(nq>mxq) then
       write(6,'(a,i0,a)') 'oncvpsp_ld: ERROR more than ',mxq,' radii'
       error stop 1
     end if
     read(arg,*,iostat=ios) rq(1:nq)
     if(ios/=0 .or. nq<1 .or. any(rq(1:nq)<=0.0d0)) then
       write(6,'(2a)') 'oncvpsp_ld: ERROR --radii must be positive ', &
&            'numbers separated by commas'
       error stop 1
     end if
     given_radii=.true.
   case('-l','--channels')
     call next_arg(arg)
     nl=min(count_list(arg),size(lsel)+1)
     lsel(:)=-1
     ios=1
     if(nl<=size(lsel)) read(arg,*,iostat=ios) lsel(1:nl)
     if(ios/=0) lsel(1)=-1
     if(any(lsel(1:min(nl,size(lsel)))<0) &
&       .or. any(lsel(1:min(nl,size(lsel)))>3)) then
       write(6,'(2a)') 'oncvpsp_ld: ERROR --channels must be values ', &
&            'of l from 0 to 3 separated by commas'
       error stop 1
     end if
     chan(:)=.false.
     chan(lsel(1:nl))=.true.
   case('-e','--energies')
     call next_arg(arg)
     read(arg,*,iostat=ios) ee(1:3)
     if(ios/=0 .or. count_list(arg)/=3 .or. ee(3)<=0.0d0 &
&       .or. ee(2)<=ee(1)) then
       write(6,'(2a)') 'oncvpsp_ld: ERROR --energies must be ', &
&            'emin,emax,de with emin < emax and de > 0'
       error stop 1
     end if
     given_energies=.true.
   case default
     write(6,'(2a)') 'oncvpsp_ld: ERROR unknown argument ',trim(arg)
     error stop 1
   end select
   iarg=iarg+1
 end do
 if(len_trim(ae_filename)==0 .or. len_trim(ps_filename)==0) then
   write(6,'(a)') 'oncvpsp_ld: ERROR requires --ae-file and --ps-file'
   error stop 1
 end if

 call read_ae_potential(trim(ae_filename),irel,atsym,zz,mmax,rr,c_light_f, &
&                       vfull)
 call read_ps_file(trim(ps_filename),pf)
 if(pf%irel/=irel .or. pf%atsym/=atsym .or. differ(pf%zz,zz) &
&   .or. pf%mmax/=mmax .or. differ(pf%c_light,c_light_f)) then
   write(6,'(2a)') 'oncvpsp_ld: ERROR the files differ in relativistic ', &
&       'treatment, atom, mesh size or c_light'
   error stop 1
 end if
 if(any(abs(pf%rr(:)-rr(:))>rtol*rr(:))) then
   write(6,'(a)') 'oncvpsp_ld: ERROR the files'' radial meshes differ'
   error stop 1
 end if
 call set_relativity(c_light_f)

 if(irel==0) then
   treatment='non-relativistic'
 else if(irel==1) then
   treatment='scalar-relativistic'
 else
   treatment='fully-relativistic'
 end if
 call write_build_info(6,'oncvpsp_ld.x',trim(treatment))

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

 write(6,'(a)') '# ld_r l ikap r(bohr)'
 write(6,'(a)') '# ld l ikap e(Ha) atan(r*u''/u)_ae atan(r*u''/u)_ps'
 do l1=1,4
   ll=l1-1
   if(.not. chan(ll)) cycle
!the input's radius of l, as run_phsft
   if(.not. given_radii) then
     nq=1
     if(pf%rxpsh>=0.0d0) then
       rq(1)=pf%rxpsh
     else if(ll<=pf%lmax) then
       rq(1)=pf%rc(l1)*amesh**2
     else
       lc=pf%lloc+1
       rq(1)=pf%rc(lc)
     end if
   end if
   rql=rq(1:nq)
   allocate(pshae(npsh,nq),pshps(npsh,nq))
   do ikap=1,merge(1,pf%nkap,l1==1)
     if(ikap==1) kap=-(ll+1)
     if(ikap==2) kap=ll
     nvkb=pf%nvkb(l1,ikap)
     call phsft_pair(irel,ll,kap,nvkb,rr,vfull,zz,pf%vloc(:,ikap), &
&                    pf%vlocpw,pf%vkbpw(1:nvkb,l1,ikap), &
&                    pf%evkb(1:nvkb,l1,ikap),mmax,nq,rql,npsh,epsh2,depsh, &
&                    pf%epa(l1,ikap),pshae,pshps)
     do iq=1,nq
       write(6,'(a,2i4,1p,e24.16)') 'ld_r',ll,ikap,rql(iq)
       do ii=1,npsh
         write(6,'(a,2i4,1p,3e24.16)') 'ld',ll,ikap, &
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
 write(6,'(a)') 'Usage: oncvpsp_ld.x -a AE_FILE -p PS_FILE [options]'
 write(6,'(a)') ''
 write(6,'(a)') 'The log derivatives of a pseudopotential and its all-electron'
 write(6,'(a)') 'atom, atan(r u''/u - 1) at radii r over a grid of energies, as'
 write(6,'(a)') 'continuous phases (the pseudopotential''s aligned with the'
 write(6,'(a)') 'all-electron ones at the reference energy of each l), for'
 write(6,'(a)') 'l = 0..3 (or those of -l) and, fully relativistic, each j.'
 write(6,'(a)') 'All the radii come from one solution per energy. Without'
 write(6,'(a)') '-r, -l or -e it prints the log derivatives of the run''s'
 write(6,'(a)') 'DIAGNOSTICS block. The output, on stdout, is a block per'
 write(6,'(a)') 'l, j and radius:'
 write(6,'(a)') '  ld_r l ikap r(bohr)'
 write(6,'(a)') '  ld l ikap e(Ha) phase_ae phase_ps'
 write(6,'(a)') 'ikap 1 for j = l+1/2 (and every l if not fully relativistic),'
 write(6,'(a)') '2 for j = l-1/2.'
 write(6,'(a)') ''
 write(6,'(a)') 'Options:'
 write(6,'(a)') '  -h, --help'
 write(6,'(a)') '      Show this help and exit'
 write(6,'(a)') '  -v, --version'
 write(6,'(a)') '      Show the build and exit'
 write(6,'(a)') '  -a, --ae-file FILE'
 write(6,'(a)') '      The all-electron data file of an _ae.x (required)'
 write(6,'(a)') '  -p, --ps-file FILE'
 write(6,'(a)') '      The pseudopotential data file of a run of oncvpsp.x,'
 write(6,'(a)') '      oncvpspnr.x, oncvpspr.x or their _ps.x with --ps-file,'
 write(6,'(a)') '      for the same atom, treatment and mesh (required)'
 write(6,'(a)') '  -r, --radii R1,R2,...'
 write(6,'(a)') '      The radii (bohr), each for every l; inside the core'
 write(6,'(a)') '      radii too. Default the input''s rxpsh, or for rxpsh < 0'
 write(6,'(a)') '      rc(l)*amesh**2 (rc of the local potential above lmax)'
 write(6,'(a)') '  -l, --channels L1,L2,...'
 write(6,'(a)') '      The angular momenta, from 0 to 3, whose log'
 write(6,'(a)') '      derivatives are computed (both j of each when fully'
 write(6,'(a)') '      relativistic). Default all four'
 write(6,'(a)') '  -e, --energies EMIN,EMAX,DE'
 write(6,'(a)') '      The energies (Ha), from EMAX down by DE to above EMIN.'
 write(6,'(a)') '      Default the input''s epsh1,epsh2,depsh'
 write(6,'(a)') ''
 write(6,'(a)') 'Errors are reported on stdout, with a nonzero exit status.'
 end subroutine write_help

 subroutine next_arg(val)
!! the value of the option at iarg, which it passes
 character(len=1024), intent(out) :: val
 if(iarg+1>command_argument_count()) then
   call get_command_argument(iarg,val)
   write(6,'(3a)') 'oncvpsp_ld: ERROR ',trim(val),' requires a value'
   error stop 1
 end if
 iarg=iarg+1
 call get_command_argument(iarg,val)
 end subroutine next_arg

 integer function count_list(str)
!! the number of comma-separated items of str
 character(len=*), intent(in) :: str
 integer :: kk
 count_list=1
 do kk=1,len_trim(str)
   if(str(kk:kk)==',') count_list=count_list+1
 end do
 end function count_list

 elemental logical function differ(aa,bb)
!! .true. if aa and bb differ by more than rtol relative to the larger
 real(dp), intent(in) :: aa,bb
 differ=abs(aa-bb)>rtol*max(1.0d0,abs(aa),abs(bb))
 end function differ

 end program oncvpsp_ld

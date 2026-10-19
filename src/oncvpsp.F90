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
#if defined(ONCVPSP_AE_ONLY) || defined(ONCVPSP_PS_ONLY)
#define ONCVPSP_SPLIT
#endif
#if !defined(ONCVPSP_EXE)
#define ONCVPSP_EXE 'oncvpsp.x'
#endif
#if defined(ONCVPSP_NONREL)
#define ONCVPSP_TREATMENT 'non-relativistic'
#else
#define ONCVPSP_TREATMENT 'scalar-relativistic'
#endif
 program oncvpsp
!
! Creates and tests optimized norm-conserving Vanderbilt or Kleinman-Bylander
! pseudopotentials based on D. R. Hamann, Phys. Rev. B 88, 085117 (2013)
! and references therein.
!
!   D. R. Hamann
!   Mat-Sim Research LLC
!   P.O. Box 742
!   Murray Hill, NJ 07974
!   USA
!
!   Developed from original "gncpp" code of March 8,1987
!
!   Output format for ABINIT pspcod=8 and upf format for quantumespresso
!
!   Built with ONCVPSP_AE_ONLY (oncvpsp_ae.x), this program solves only the
!   all-electron atom and the all-electron halves of the test
!   configurations, and writes them to the file given by --ae-file.  Built
!   with ONCVPSP_PS_ONLY (oncvpsp_ps.x), it reads that file instead of
!   solving the all-electron atom and otherwise runs as oncvpsp.x, with the
!   same output.  Both read the same input.  See ae_file_m.
!
 use, intrinsic :: iso_fortran_env, only: stdin => input_unit, &
&                                         stdout => output_unit, stderr => error_unit
 use m_psmlout, only: psmlout
 use build_info_m, only: write_build_info
 use input_text_m, only: read_input_text
 use config_ae_m, only: config_ae_t, run_config_ae
#if defined(ONCVPSP_AE_ONLY)
 use ae_file_m, only: write_ae_file
#elif defined(ONCVPSP_PS_ONLY)
 use ae_file_m, only: read_ae_file, ae_file_config
#endif
#if (defined WITH_TOML)
 use input_toml_m, only: read_input_toml
#endif
 use constants_m, only: dp, mxprj
 use input_echo_m, only: write_icmod, write_epsh
 implicit none

!
 integer :: ii,ierr,iexc,iexct,ios,iprint,irps,it,icmod,lpopt
 integer :: jj,kk,ll,l1,lloc,lmax,lt,inline
 integer :: mch,mchf,mmax,n1,n2,nc,nlim,nlloc,nlmax,irpsh,nrl
 integer :: nv,irct,ncnf,nvt
 integer :: iprj
 integer,allocatable :: npa(:,:)
!
 integer :: dtime(8),na(30),la(30),np(6)
 integer :: nacnf(30,5),lacnf(30,5),nvcnf(5)
 integer :: nat(30),lat(30),indxr(30),indxe(30)
 integer :: irc(6),nodes(4)
 integer :: nproj(6),npx(6),lpx(6)
 integer :: ncon(6),nbas(6)

 real(dp) :: al,amesh,csc,csc1,deblt,depsh,depsht,drl,eeel
 real(dp) :: eeig,eexc
 real(dp) :: emax,epsh1,epsh1t,epsh2,epsh2t,rxpsh
 real(dp) :: et,etest,emin,sls
 real(dp) :: fcfact,fcfact0,rcfact,dvloc0
 real(dp) :: fcfact_min,fcfact_max,fcfact_step
 real(dp) :: rcfact_min,rcfact_max,rcfact_step
 real(dp) :: fcgrid(3),rcgrid(3)
 real(dp) :: rr1,rcion,rciont,rcmax,rct,rlmax,rpkt
 real(dp) :: amesh_in,rr1_in,rmax_in,rmax,log_mesh(3)
 real(dp) :: sf,zz,zion,zval,etot
 real(dp) :: xdummy
!
 real(dp) :: cl(6),debl(6),ea(30),ep(6),fa(30),facnf(30,5)
 real(dp) :: fat(30,2)
 real(dp) :: fnp(6),fp(6)
 real(dp) :: qcut(6),qmsbf(6),rc(6),rc0(6)
 real(dp) :: rpk(30)
 real(dp) :: epx(6),fpx(6)
 real(dp) :: epstot
 real(dp), parameter :: eps=1.0d-8

 real(dp), allocatable :: evkb(:,:),cvgplt(:,:,:,:),qq(:,:)
 real(dp), allocatable :: ecurve(:,:,:,:),wcurve(:,:)
 real(dp) :: ecutw
 integer, parameter :: ncurve=1000
 real(dp), allocatable :: rr(:)
 real(dp), allocatable :: rho(:),rhoc(:),rhot(:)
 real(dp), allocatable :: uu(:),up(:)
 real(dp), allocatable :: vp(:,:),vfull(:),vkb(:,:,:),pswf(:,:,:)
 real(dp), allocatable :: vwell(:)
 real(dp), allocatable :: vpuns(:,:)
 real(dp), allocatable :: vo(:),vxc(:)
 real(dp), allocatable :: rhomod(:,:),rhoae(:,:),rhops(:,:),rhotae(:)
 real(dp), allocatable :: uupsa(:,:) !pseudo-atomic orbitals array
 real(dp), allocatable :: epa(:,:),fpa(:,:)
 real(dp), allocatable :: uua(:,:),upa(:,:)
 real(dp), allocatable :: vr(:,:,:)

 type(config_ae_t) :: cae(5)

 character(len=2) :: atsym
 character(len=4) :: psfile

 logical :: srel,cset

 integer, parameter :: INPUT_STDIN=1, INPUT_TEXT=2, INPUT_TOML=3
 integer :: input_mode
 integer :: unit
 character(len=1024) :: input_filename
 character(len=1024) :: ae_filename=''

 input_mode = INPUT_STDIN

 write(6,'(a/a//)') &
&      'ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&      ONCVPSP_TREATMENT//' version 4.0.1 03/01/2019'

 write(6,'(a/a/a//)') &
&      'While it is not required under the terms of the GNU GPL, it is',&
&      'suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)', &
&      'in any publication utilizing these pseudopotentials.'

 call write_build_info(6,ONCVPSP_EXE,ONCVPSP_TREATMENT)
 write(6,'(/)')

#if defined(ONCVPSP_NONREL)
 srel=.false.
#else
 srel=.true.
#endif

 parse_args: block
    integer :: i
    character(len=256) :: arg
    do i = 1, command_argument_count()
       call get_command_argument(i, arg)
       select case(arg)
       case('-h', '--help')
          write (stdout, '(a)') 'Usage: '//ONCVPSP_EXE//' [options]'
          write (stdout, '(a)') 'Options:'
          write (stdout, '(a)') '  -h, --help       Show this help message and exit'
          write (stdout, '(a)') '  -v, --version    Show version information and exit'
          write (stdout, '(a)') '  -i, --input      Read the input from a text file ' // &
&                               'instead of stdin'
          write (stdout, '(a)') '  -t, --toml-input Read the input from a TOML file'
#if defined(ONCVPSP_AE_ONLY)
          write (stdout, '(a)') '  -a, --ae-file    All-electron data file to write (required)'
#elif defined(ONCVPSP_PS_ONLY)
          write (stdout, '(a)') '  -a, --ae-file    All-electron data file written by ' // &
&                               'the matching _ae.x (required)'
#endif
          stop 0
       case('-v', '--version')
#if defined(ONCVPSP_NONREL)
          write (stdout, '(a)') 'ONCVPSP (non-relativistic) version 4.0.1'
#else
          write (stdout, '(a)') 'ONCVPSP (scalar-relativistic) version 4.0.1'
#endif
          call write_build_info(stdout,ONCVPSP_EXE,ONCVPSP_TREATMENT)
          stop 0
       case('-i', '--input')
         if (i + 1 > command_argument_count()) then
            write (stderr, '(a)') 'Error: --input requires a filename argument'
            stop 1
         end if
         call get_command_argument(i + 1, input_filename)
         input_mode = INPUT_TEXT
#if (defined WITH_TOML)
       case('-t', '--toml-input')
         if (i + 1 > command_argument_count()) then
            write (stderr, '(a)') 'Error: --toml-input requires a filename argument'
            stop 1
         end if
         call get_command_argument(i + 1, input_filename)
         input_mode = INPUT_TOML
#else
        case('-t', '--toml-input')
          error stop 'Error: TOML input support not enabled in this build.'
#endif
#if defined(ONCVPSP_SPLIT)
       case('-a', '--ae-file')
         if (i + 1 > command_argument_count()) then
            write (stderr, '(a)') 'Error: --ae-file requires a filename argument'
            stop 1
         end if
         call get_command_argument(i + 1, ae_filename)
#endif
       case default
       end select
    end do
 end block parse_args

#if defined(ONCVPSP_SPLIT)
 if (len_trim(ae_filename) == 0) then
    write (stderr, '(a)') 'Error: '//ONCVPSP_EXE//' requires --ae-file'
    stop 1
 end if
#endif

 ios = 0
 select case(input_mode)
    case(INPUT_STDIN)
       call read_input_text(stdin, inline, &
                            atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
                            lmax, rc, ep, ncon, nbas, qcut, &
                            lloc, lpopt, dvloc0, nproj, debl, &
                            icmod, fcfact, rcfact, &
                            fcfact_min, fcfact_max, fcfact_step, &
                            rcfact_min, rcfact_max, rcfact_step, &
                            epsh1, epsh2, depsh, rxpsh, &
                            rlmax, drl, &
                            ncnf, nvcnf, nacnf, lacnf, facnf, &
                            amesh_in, rr1_in, rmax_in)
    case(INPUT_TEXT)
       open(newunit=unit, file=input_filename, status='old', action='read', iostat=ios)
       if (ios /= 0) then
          write(6, '(/2a)') 'oncvpsp: ERROR cannot open input file, file=', &
             trim(input_filename)
          error stop 1
       end if
       call read_input_text(unit, inline, &
                            atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
                            lmax, rc, ep, ncon, nbas, qcut, &
                            lloc, lpopt, dvloc0, nproj, debl, &
                            icmod, fcfact, rcfact, &
                            fcfact_min, fcfact_max, fcfact_step, &
                            rcfact_min, rcfact_max, rcfact_step, &
                            epsh1, epsh2, depsh, rxpsh, &
                            rlmax, drl, &
                            ncnf, nvcnf, nacnf, lacnf, facnf, &
                            amesh_in, rr1_in, rmax_in)
       close(unit)
#if (defined WITH_TOML)
    case(INPUT_TOML)
       open(newunit=unit, file=input_filename, status='old', action='read', iostat=ios)
       if (ios /= 0) then
          write(6, '(/2a)') 'oncvpsp: ERROR cannot open input file, file=', &
             trim(input_filename)
          error stop 1
       end if
       call read_input_toml(unit, &
                            atsym, zz, nc, nv, iexc, psfile, na, la, fa, &
                            lmax, rc, ep, ncon, nbas, qcut, &
                            lloc, lpopt, dvloc0, nproj, debl, &
                            icmod, fcfact, rcfact, &
                            fcfact_min, fcfact_max, fcfact_step, &
                            rcfact_min, rcfact_max, rcfact_step, &
                            epsh1, epsh2, depsh, rxpsh, &
                            rlmax, drl, &
                            ncnf, nvcnf, nacnf, lacnf, facnf, &
                            amesh_in, rr1_in, rmax_in)
        close(unit)
#endif
    case default
        write (stderr, '(a,i2)') 'Error: Unknown input mode =', input_mode
        stop 1
 end select

 nvcnf(1)=nv
 do ii=1,nc+nv
   nacnf(ii,1)=na(ii)
   lacnf(ii,1)=la(ii)
   facnf(ii,1)=fa(ii)
 end do

 do jj=2,ncnf+1
   do ii=1,nc
     nacnf(ii,jj)=na(ii)
     lacnf(ii,jj)=la(ii)
     facnf(ii,jj)=fa(ii)
   end do
 end do

 if(icmod==0) then
   fcfact=0.0d0
 end if

!icmod=4 coarse-search grids: min, max, step
 fcgrid=[fcfact_min,fcfact_max,fcfact_step]
 rcgrid=[rcfact_min,rcfact_max,rcfact_step]

 call check_data(atsym,zz,fcfact,rcfact,fcgrid,rcgrid,epsh1,epsh2,depsh,rxpsh,rlmax,drl, &
&                fa,facnf, &
&                rc,ep,qcut,debl,nc,nv,iexc,lmax,lloc,lpopt,icmod, &
&                ncnf,na,la,nvcnf,nacnf,lacnf,ncon,nbas,nproj,psfile, &
&                amesh_in,rr1_in,rmax_in)

 nrl=int((rlmax/drl)-0.5d0)+1

!PWSCF wants an even number of mesh pointe
!if(trim(psfile)=='upf') then
  if(mod(nrl,2)/=0) nrl=nrl+1
!end if

!log mesh: the input values, or the built-in ones for those 0 or less
 if(amesh_in>0.0d0) then
  amesh=amesh_in
 else
!amesh=1.012d0
  amesh=1.006d0
!amesh=1.003d0
!amesh=1.0015d0
 end if

 al=log(amesh)
 if(rr1_in>0.0d0) then
  rr1=rr1_in
 else
  rr1=0.0005d0/zz
  rr1=min(rr1,0.0005d0/10)
 end if
 if(rmax_in>0.0d0) then
  rmax=rmax_in
 else
  rmax=45.0d0
 end if
 mmax=log(rmax/rr1)/al
!the log mesh for the echoes of the input, only when the input sets it, so that the
!output of an input without it is unchanged
 log_mesh(:)=0.0d0
 if(amesh_in>0.0d0 .or. rr1_in>0.0d0 .or. rmax_in>0.0d0) log_mesh(:)=[amesh,rr1,rmax]

!calculate zion for output
 zion=zz
 do ii=1,nc
  zion=zion-fa(ii)
 end do


 allocate(rr(mmax))
 allocate(rho(mmax),rhoc(mmax),rhot(mmax))
 allocate(uu(mmax),up(mmax),uupsa(mmax,30))
 allocate(evkb(mxprj,4), cvgplt(2,7,mxprj,4),qq(mxprj,mxprj))
 allocate(ecurve(2,ncurve,mxprj,4),wcurve(mxprj,4))
 allocate(vp(mmax,5),vfull(mmax),vkb(mmax,mxprj,4),pswf(mmax,mxprj,4))
 allocate(vwell(mmax))
 allocate(vpuns(mmax,5))
 allocate(vo(mmax),vxc(mmax))
 allocate(rhoae(mmax,nv),rhops(mmax,nv),rhotae(mmax))
 allocate(npa(mxprj,6))
 allocate(epa(mxprj,6),fpa(mxprj,6))
 allocate(uua(mmax,mxprj),upa(mmax,mxprj))
 allocate(vr(mmax,mxprj,6))

 vr(:,:,:)=0.0d0
 vp(:,:)=0.0d0
 vkb(:,:,:)=0.0d0
 epa(:,:)=0.0d0

 do ii=1,mmax
  rr(ii)=rr1*exp(al*(ii-1))
 end do

!
! full potential atom solution
!
#if defined(ONCVPSP_PS_ONLY)
   call read_ae_file(ae_filename,srel,atsym,zz,nc,nv,iexc,na,la,fa, &
&                    mmax,rr,it,etot,ea,rhoc,rho,vfull)
!sratom prints the libxc functional description on its first exc_libxc call
   if(iexc<0) then
     block
       real(dp) :: vxct(mmax),exct(mmax)
       call exc_libxc(iexc,al,rho,vxct,exct,rr,mmax)
     end block
   end if
#else
   call sratom(na,la,ea,fa,rpk,nc,nc+nv,it,rhoc,rho, &
&              rr,vfull,zz,mmax,iexc,etot,ierr,srel)
#endif
!
!

! Drop digits beyond 5 decimals for input rcs before making any use of them
 do l1=1,max(lmax+1,lloc+1)
    jj=int(rc(l1)*10.0d5)
    rc(l1)=jj/10.0d5
 end do

 rcmax=0.0d0
 do l1=1,lmax+1
   rcmax=max(rcmax,rc(l1))
 end do
 do l1=1,lmax+1
   if(rc(l1)==0.0d0) then
     rc(l1)=rcmax
   end if
 end do

 nproj(lloc+1)=0
 rc0(:)=rc(:)
!input value, for the echoes (modcore3 resets fcfact for icmod=4)
 fcfact0=fcfact

! output printing (echos input data, with all-electron eigenvalues added)

 write(6,'(a)') '# ATOM AND REFERENCE CONFIGURATION'
 write(6,'(a)') '# atsym  z   nc   nv     iexc    psfile'
 write(6,'(a,a,f6.2,2i5,i8,2a)') '  ',trim(atsym),zz,nc,nv,iexc, &
&      '      ',psfile
 write(6,'(a/a)') '#','#   n    l    f        energy (Ha)'
 do ii=1,nc+nv
   write(6,'(2i5,f8.2,1pe18.7)') na(ii),la(ii),fa(ii),ea(ii)
 end do

 write(6,'(a/a/a)') '#','# PSEUDOPOTENTIAL AND OPTIMIZATION','# lmax'
 write(6,'(i5)')  lmax
 write(6,'(a/a)') '#','#   l,   rc,      ep,       ncon, nbas, qcut'
 do l1=1,lmax+1
   write(6,'(i5,2f10.5,2i5,f10.5)') l1-1,rc(l1),ep(l1),ncon(l1),&
&        nbas(l1),qcut(l1)
 end do

 write(6,'(a/a/a,a)') '#','# LOCAL POTENTIAL','# lloc, lpopt,  rc(5),', &
&      '   dvloc0'
 write(6,'(2i5,f10.5,a,f10.5)') lloc,lpopt,rc(5),'   ',dvloc0

 write(6,'(a/a/a)') '#','# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs', &
&      '# l, nproj, debl'
 do l1=1,lmax+1
   write(6,'(2i5,f10.5)') l1-1,nproj(l1),debl(l1)
 end do

 call write_icmod(icmod,fcfact,rcfact,fcgrid,rcgrid)

 call write_epsh(epsh1,epsh2,depsh,rxpsh)

 write(6,'(a/a/a)') '#','# OUTPUT GRID','# rlmax, drl'
 write(6,'(2f8.4)') rlmax,drl

 write(6,'(a/a/a)') '#','# TEST CONFIGURATIONS','# ncnf'
 write(6,'(i5)') ncnf
 write(6,'(a/a)') '# nvcnf','#   n    l    f'
 do jj=2,ncnf+1
   write(6,'(i5)') nvcnf(jj)
   do ii=nc+1,nc+nvcnf(jj)
     write(6,'(2i5,f8.2)') nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
   end do
   write(6,'(a)') '#'
 end do

 if(log_mesh(1)>0.0d0) then
   write(6,'(a/a/a)') '#','# LOG MESH','# amesh, rr1, rmax'
   write(6,'(3es24.16)') log_mesh
 end if

 write(6,'(//a)') 'Reference configufation results'
 write(6,'(a,i6)') '  iterations',it
 if(it >= 100) then
  write(6,'(a)') 'oncvpsp: ERROR all-electron reference atom not converged'
  error stop 1
 end if
 write(6,'(a,1p,d18.8)') '  all-electron total energy (Ha)',etot

#if defined(ONCVPSP_AE_ONLY)
!all-electron halves of the test configurations; their pseudo-atom halves
!and all output about them come from oncvpsp_ps.x
 do jj=1,ncnf+1
   call run_config_ae(jj,nacnf,lacnf,facnf,nc,nvcnf,rr,zz, &
&                     mmax,iexc,ea,srel,cae(jj))
 end do

 call write_ae_file(ae_filename,srel,atsym,zz,nc,nv,iexc,na,la,fa, &
&                   mmax,rr,it,etot,ea,rhoc,rho,vfull, &
&                   ncnf,nvcnf,nacnf,lacnf,facnf,cae)

 write(6,'(/2a)') 'All-electron data written to ',trim(ae_filename)
 stop
#endif

!find log mesh point nearest input rc
 rcmax=0.0d0
 irc(:)=0
 do l1=1,max(lmax+1,lloc+1)
  rct=rc(l1)
  irc(l1)=0
  do ii=2,mmax
   if(rr(ii)>rct) then
    irc(l1)=ii
    rc(l1)=rr(ii)
    exit
   end if
  end do
  rcmax=max(rcmax,rc(l1))
 end do

!
 cvgplt(:,:,:,:)=0.0d0
 ecurve(:,:,:,:)=0.0d0
 wcurve(:,:)=0.0d0
!
! loop to construct pseudopotentials for all angular momenta
!
 write(6,'(/a/a)') 'Begin loop to  construct optimized pseudo wave functions',&
&      'and semi-local pseudopoentials for all angular momenta'

!temporarily set this to 1 so that the pseudo wave function needed for the
!local potential will be generated.  Reset after run_vkb.
 nproj(lloc+1)=1
 do l1=1,lmax+1
   ll=l1-1
   uu(:)=0.0d0; qq(:,:)=0.0d0
   iprj=0
!electrons filling each pseudo-wave in the cutoff of icmod 6: the shell's
!for a scattering or empty state, else the occupation
   wcurve(1:nproj(l1),l1)=dble(2*(2*ll+1))

!get principal quantum number for the highest core state for this l
   npa(1,l1)=l1
   do kk=1,nc
     if(la(kk)==l1-1) npa(1,l1)=na(kk)+1
   end do !kk

!get all-electron bound states for projectors
   if(nv/=0) then
     do kk=nc+1,nc+nv
       if(la(kk)==l1-1) then
         iprj=iprj+1
         et=ea(kk)
         call lschfb(na(kk),la(kk),ierr,et, &
&                      rr,vfull,uu,up,zz,mmax,mch,srel)
         if(ierr /= 0) then
           write(6,'(/a,3i4)') 'oncvpsp-387: lschfb convergence ERROR n,l,iter=', &
&           na(kk),la(kk),it
           error stop 1
         end if
         epa(iprj,l1)=ea(kk)
         if(fa(kk)>0.0d0) wcurve(iprj,l1)=fa(kk)
         npa(iprj,l1)=na(kk)
         uua(:,iprj)=uu(:)
         upa(:,iprj)=up(:)
       end if !la(kk)==l1-1
       if(iprj==nproj(l1)) exit
     end do !kk
   end if !nv/=0

!get all-electron well states for projectors
!if there were no valence states, use ep from input data for 1st well state
!otherwise shift up by input debl
   if(iprj==0) epa(1,l1)=ep(l1)
   if(iprj<nproj(l1))then
     do kk=1,nproj(l1)-iprj
       iprj=iprj+1
       if(iprj>1 .and. debl(l1)<=0.0d0) then
         write(6,'(a,f8.3,a,i2/a/a)') 'oncvpsp: ERROR debl =',debl(l1), &
&              ' for l=',ll, &
&              ' ERROR not allowed with 2 or more scattering states', &
&              'program will stop'
         error stop 1
       end if
       if(iprj>1) then
         epa(iprj,l1)=epa(iprj-1,l1)+debl(l1)
         npa(iprj,l1)=npa(iprj-1,l1)+1
       end if

       call wellstate(npa(iprj,l1),ll,irc(l1),epa(iprj,l1),rr, &
&                     vfull,uu,up,zz,mmax,mch,srel)
       uua(:,iprj)=uu(:)
       upa(:,iprj)=up(:)
     end do !kk
   end if !iprj<nproj(l1)

   do iprj=1,nproj(l1)

!calculate relativistic correction to potential to force projectors to 0 at rc
     call vrel(ll,epa(iprj,l1),rr,vfull,vr(1,iprj,l1),uua(1,iprj),upa(1,iprj), &
&              zz,mmax,irc(l1),srel)

   end do

!get all-electron overlap matrix
   do jj=1,nproj(l1)
     do ii=1,jj
       call fpovlp(uua(1,ii),uua(1,jj),irc(l1),ll,zz,qq(ii,jj),rr,srel)
       qq(jj,ii)=qq(ii,jj)
     end do
   end do

   call run_optimize(epa(1,l1),ll,mmax,mxprj,rr,uua,qq, &
&                    irc(l1),qcut(l1),qmsbf(l1),ncon(l1),nbas(l1),nproj(l1), &
&                    pswf(1,1,l1),vp(1,l1),vkb(1,1,l1),vfull,cvgplt(1,1,1,l1), &
&                    ncurve,ecurve(1,1,1,l1))

 end do !l1

! construct Vanderbilt / Kleinman-Bylander projectors

 write(6,'(/a,a)') 'Construct Vanderbilt / Kleinmman-Bylander projectors'

 call run_vkb(lmax,lloc,lpopt,dvloc0,irc,nproj,rr,mmax,mxprj,pswf,vfull,vp, &
&             evkb,vkb,nlim,vr)

!restore this to its proper value
 nproj(lloc+1)=0

 deallocate(uua,upa)

! accumulate charge and eigenvalues
! pseudo wave functions are calculated with VKB projectors for
! maximum consistency of unscreening
! get all-electron and pseudopotential valence-state by valence-state
! charge densities

! null charge and eigenvalue accumulators
 uupsa(:,:)=0.0d0
 eeig=0.0d0
 zval=0.0d0
 rho(:)=0.0d0
 nodes(:)=0
 rhotae(:)=0.0d0
 irps=0
 do kk=1,nv

   et=ea(nc+kk)
   ll=la(nc+kk)
   l1=ll+1
   call lschfb(na(nc+kk),ll,ierr,et, &
&                rr,vfull,uu,up,zz,mmax,mch,srel)
   if(ierr /= 0) then
     write(6,'(/a,3i4)') 'oncvpsp-461: lschfb convergence ERROR n,l,iter=', &
&     na(nc+kk),la(nc+kk),it
     error stop 1
   end if

   rhoae(:,kk)=(uu(:)/rr(:))**2

   rhotae(:)=rhotae(:) + fa(nc+kk)*rhoae(:,kk)

   emax=0.75d0*et
   emin=1.25d0*et

   call lschvkbb(ll+nodes(l1)+1,ll,nproj(l1),ierr,et,emin,emax, &
&                rr,vp(1,lloc+1),vkb(1,1,l1),evkb(1,l1), &
&                uu,up,mmax,mch)

   if(ierr/=0) then
     write(6,'(a,3i4)') 'oncvpsp: lschvkbb ERROR',ll+nodes(l1)+1,ll,ierr
     flush(6)
     error stop 1
   end if

! save valence pseudo wave functions for upfout
   uupsa(:,kk)=uu(:)

   rhops(:,kk)=(uu(:)/rr(:))**2
   rho(:)=rho(:)+fa(nc+kk)*rhops(:,kk)
   eeig=eeig+fa(nc+kk)*et

   zval=zval+fa(nc+kk)
   nodes(l1)=nodes(l1)+1
   irps=max(irps,irc(l1))
 end do !kk

 allocate(rhomod(mmax,5))

 rhomod(:,:)=0.0d0

! construct model core charge based on monotonic polynomial fit
! or Teter function fit

 if(icmod==1) then
   call modcore(icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
&               fcfact,rcfact,irps,mmax,rr,nc,nv,la,zion,iexc)

 else if(icmod==2) then
   call modcore2(icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
&               fcfact,rcfact,irps,mmax,rr,nc,nv,la,zion,iexc)

 else if(icmod==6) then
   ecutw=0.0d0
   if(fcfact<=0.0d0) then
     call modcore_budget_ecut(ecurve,wcurve,ncurve,mxprj*4,rcfact,ecutw)
   end if
   call modcore_budget(rho,rhoc,rhomod,fcfact,ecutw,rcfact, &
&               minval(irc(1:lmax+1)),mmax,rr,iexc)

 else if(icmod>=3) then
   call modcore3(icmod,rhops,rho,rhoc,rhoae,rhotae,rhomod, &
&               fcfact,rcfact,fcgrid,rcgrid,irps,mmax,rr,nc,nv,la,zion,iexc)

 end if

! screening potential for pseudocharge

 call vout(1,rho,rhomod(1,1),vo,vxc,zval,eeel,eexc,rr,mmax,iexc)

! total energy output

 epstot= eeig + eexc - 0.5d0*eeel
 write(6,'(/a,f12.6/)') 'Pseudoatom total energy', epstot

 call run_diag(lmax,npa,epa,lloc,irc, &
&                    vkb,evkb,nproj,rr,vfull,vp,zz,mmax,mxprj,srel)

 call run_ghosts(lmax,la,ea,nc,nv,lloc,irc,qmsbf, &
&                    vkb,evkb,nproj,rr,vp,mmax,mxprj)

! unscreen semi-local potentials

 do l1=1,max(lmax+1,lloc+1)
   vpuns(:,l1)=vp(:,l1)-vo(:)
 end do

!fix unscreening error due to greater range of all-electron charge
 do ii=mmax,1,-1
   if(rho(ii)==0.0d0) then
     do l1=1,max(lmax+1,lloc+1)
       vpuns(ii,l1)=-zion/rr(ii)
     end do
   else
     exit
   end if
 end do

! smoothly force vpuns to -zion/r beyond rcmax, over a window of 1 a.u.,
! with the 5th-order "smootherstep" 6x^5-15x^4+10x^3 (ported from METAPSP
! ba95225), so the psp correction no longer depends on where the
! unscreening data stops
 smooth_vloc: block
   real(dp) :: xs,ss
   do l1=1,max(lmax+1,lloc+1)
     do ii=1,mmax
       if(rr(ii)>rcmax) then
         xs=1.0d0+(rcmax-rr(ii))
         if(xs<=0.0d0) then
           ss=0.0d0
         else
           ss=xs**3*(10.0d0+xs*(-15.0d0+6.0d0*xs))
         end if
         vpuns(ii,l1)=-zion/rr(ii)+(vpuns(ii,l1)+zion/rr(ii))*ss
       end if
     end do
   end do
 end block smooth_vloc


! loop over reference plus test atom configurations
!if(.false.) then

!ncnf=0
 rhot(:)=rho(:)
 do jj=1,ncnf+1

   write(6,'(/a,i2)') 'Test configuration',jj-1

! charge density is initialized to that of reference configuration

   rhot(:)=rho(:)

#if defined(ONCVPSP_PS_ONLY)
   call ae_file_config(jj,nc,nvcnf,nacnf,lacnf,facnf,cae(jj))
#else
   call run_config_ae(jj,nacnf,lacnf,facnf,nc,nvcnf,rr,zz, &
&                     mmax,iexc,ea,srel,cae(jj))
#endif

   call run_config_ps(jj,cae(jj),nc,rhot,rhomod,rr, &
&                  rcmax,mmax,mxprj,iexc,etot,epstot,nproj,vpuns, &
&                  lloc,vkb,evkb)

 end do !jj

 call run_plot(lmax,npa,epa,lloc,irc, &
&                    vkb,evkb,nproj,rr,vfull,vp,vpuns,zz,mmax,mxprj,drl,nrl, &
&                    rho,rhoc,rhomod,srel,cvgplt)



 call run_phsft(lmax,lloc,nproj,epa,epsh1,epsh2,depsh,rxpsh,vkb,evkb, &
&               rr,vfull,vp,zz,mmax,mxprj,irc,srel)

 call gnu_script(epa,evkb,lmax,lloc,mxprj,nproj)

 if(trim(psfile)=='psp8' .or. trim(psfile)=='both') then


  call linout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
&             rhotae,rhoc,zz,zion,mmax,mxprj,iexc,icmod,nrl,drl,atsym, &
&             na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&             fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,fcfact0,rcfact,fcgrid,rcgrid, &
&             epsh1,epsh2,depsh,rxpsh,rlmax,psfile,log_mesh)
 end if

 if(trim(psfile)=='upf' .or. trim(psfile)=='both') then
  call upfout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
&             zz,zion,mmax,mxprj,iexc,icmod,nrl,drl,atsym,epstot, &
&             na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&             fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact0,rcfact,fcgrid,rcgrid, &
&             epsh1,epsh2,depsh,rxpsh,rlmax,psfile,uupsa,ea,log_mesh)
 end if


 if(trim(psfile)=='psml' .or. trim(psfile)=='both') then
!
! Write info for PSML format
!
!
! irct: first point beyond which the model core equals the all-electron core
!
   irct=mmax
   do ii=mmax,1,-1
     if(rhomod(ii,1)/=rhoc(ii)) exit
     irct=ii
   end do
   call psmlout(lmax,lloc,rc,vkb,evkb,mxprj,nproj,rr,vpuns,rho,rhomod, &
&             irct, srel, &
&             zz,zion,mmax,iexc,icmod,nrl,drl,atsym,epstot, &
&             na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&             fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact0,rcfact,fcgrid,rcgrid, &
&             epsh1,epsh2,depsh,rxpsh,rlmax,psfile,log_mesh)
 end if

 stop
 end program oncvpsp

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

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
! interpolates various arrays onto linear radial mesh to create file
! for PWSCF input using the UPF file format, with the spin-orbit data when
! fully relativistic, or takes them at every istr-th point of the log mesh
! (upf_mesh > 0), as generated, with nothing interpolated

 subroutine upfout(lmax,lloc,rc,vkbpw,evkb,nproj,rr,vpuns,vlocpw,rho,corepw, &
&                  zz,zion,mmax,mxprj,iexc,icmod,nrlin,drl,atsym,epstot, &
&                  na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&                  fa,rc0,ep,qcut,ep_as_given,well_onset,well_tail_free, &
&                  debl,facnf, &
&                  dvloc0,fcfact,rcfact, &
&                  slack,fcgrid,rcgrid, &
&                  epsh1,epsh2,depsh,rxpsh,rlmax,psfile,uupsa,ea,log_mesh, &
&                  upf_mesh,nkap)


!! @param lmax  maximum angular momentum
!! @param lloc  l for local potential
!! @param rc  core radii
!! @param vkbpw  VKB projectors' pieces, last index 1 for j=l+1/2, 2 for
!!               j=l-1/2 when nkap=2
!! @param evkb  coefficients of VKB projectors, last index as for vkbpw
!! @param nproj  number of vkb projectors for each l
!! @param rr  log radial grid
!! @param vpuns  unscreened semi-local pseudopotentials (vp(:,5) is local
!!               potential if linear combination is used)
!! @param vlocpw  the screened local potential's pieces; vpuns(:,lloc+1)
!!                less their node values is smooth
!! @param rho  valence pseudocharge
!! @param corepw  the model core charge and its first 4 derivatives, as
!!                pieces split at its joint (modcore_m, core_pieces)
!! @param zz  atomic number
!! @param zion  at this point, total valence charge (becomes psuedoion charge)
!! @param mmax  size of log radial grid
!! @param mxprj  dimension of number of projectors
!! @param iexc  type of exchange-correlation
!! @param icmod  1 if model core charge is used, otherwise 0
!! @param nrlin  size of linear radial grid (enlarged locally to fit the
!!               orbitals)
!! @param drl  spacing of linear radial grid
!! @param atsym  atomic symbol
!! @param epstot  pseudoatom total energy
!! @param na  principal quantum numbers of the reference configuration states
!! @param la  angular momenta of the reference configuration states
!! @param ncon  number of matching constraints for each l
!! @param nbas  number of basis functions for each l
!! @param nvcnf  numbers of valence states of the reference and test
!!               configurations
!! @param nacnf  principal quantum numbers of the reference and test
!!               configurations
!! @param lacnf  angular momenta of the reference and test configurations
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param lpopt  type of polynomial continuation for lloc=4
!! @param ncnf  number of test configurations
!! @param fa  occupations of the reference configuration states
!! @param rc0  core radii as read from the input (bohr), rc0(5) for lloc=4
!! @param ep  energies at which the pseudopotentials are generated (Ha), second
!!            index 1 for j=l+1/2, 2 for j=l-1/2 (ep(:,1) is echoed)
!! @param qcut  residual-energy cutoff wave vectors for each l
!! @param debl  energy shifts for the second projectors for each l (Ha)
!! @param facnf  occupancies of the reference and test configurations
!! @param dvloc0  r=0 shift of the lloc=4 local potential (Ha)
!! @param fcfact  prefactor for model core charge amplitude, as read from the
!!                input
!! @param rcfact  prefactor for model core charge scale (icmod=3), or error
!!                tolerance (icmod=6)
!! @param slack  share of the cutoff the model core must fit (icmod=6, 7)
!! @param fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!! @param rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!! @param epsh1  lower limit of the log derivative analysis energy range (Ha)
!! @param epsh2  upper limit of the log derivative analysis energy range (Ha)
!! @param depsh  energy step of the log derivative analysis (Ha)
!! @param rxpsh  radius of the log derivative analysis (bohr), negative for
!!               default
!! @param rlmax  maximum radius of the linear output mesh (bohr)
!! @param psfile  pseudopotential file format, 'upf' or 'both'
!! @param log_mesh  amesh, rr1, rmax of the log mesh, echoed when the input sets
!!                  it (else 0)
!! @param upf_mesh  0 for the linear output mesh, else the most points of the
!!                  log mesh of the file: every istr-th point of rr, out to
!!                  the extent of the linear mesh
!! @param uupsa  pseudo-atomic orbitals, second index 1 for j=l+1/2, 2 for
!!               j=l-1/2 when nkap=2
!! @param ea  psuedo-orbital eigenvalues (Ha), second index as for uupsa
!! @param nkap  number of kappa per l: 1, or 2 for the fully-relativistic code

 use constants_m, only: dp, pi, mxstate, mxcnf, mxl, mxlin, nrhomod
 use input_echo_m, only: write_input_echo
 use outmesh_m, only: om_values, om_piece_values, om_local_values
 use pieces_m, only: piecewise_t
 implicit none


!Input variables
 integer, intent(in) :: lmax,lloc,iexc,mmax,mxprj,icmod,nrlin,nkap,upf_mesh
 integer, intent(in) :: nproj(mxlin+1)
 real(dp), intent(in) :: drl,fcfact,rcfact,zz,zion,epstot
 real(dp), intent(in) :: slack
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)
 real(dp), intent(in) :: rr(mmax),vpuns(mmax,mxl+2),rho(mmax)
 type(piecewise_t), intent(in) :: vkbpw(mxprj,mxl+1,nkap),vlocpw
 type(piecewise_t), intent(in) :: corepw(nrhomod)
 real(dp), intent(in) :: rc(mxlin+1),evkb(mxprj,mxl+1,nkap)
 character(len=2), intent(in) :: atsym
 real(dp), intent(in) :: uupsa(mmax,nkap,nv)

!additional input for upf output to echo input file, all as defined
! in the main progam
 integer, intent(in) :: na(mxstate),la(mxstate),ncon(mxlin+1),nbas(mxlin+1)
 integer, intent(in) :: nvcnf(mxcnf),nacnf(mxstate,mxcnf),lacnf(mxstate,mxcnf)
 integer, intent(in) :: nc,nv,lpopt,ncnf
 real(dp), intent(in) :: fa(mxstate),rc0(mxlin+1),ep(mxlin+1,nkap),qcut(mxlin+1)
 logical, intent(in) :: ep_as_given
 real(dp), intent(in) :: well_onset
 logical, intent(in) :: well_tail_free
 real(dp), intent(in) :: debl(mxlin+1)
 real(dp), intent(in) :: facnf(mxstate,mxcnf)
 real(dp), intent(in) :: ea(mxstate,nkap)
 real(dp), intent(in) :: dvloc0,epsh1,epsh2,depsh,rxpsh,rlmax
 real(dp), intent(in) :: log_mesh(3)
 character(len=4), intent(in) :: psfile

!Output variables - printing only

!Local variables
 integer :: nrl,istr,nlog
 integer :: ii,jj,ll,l1,iproj,ntotproj,nrlproj,nwfc
 integer :: ikap,mkap,jstart
 integer :: dtime(8)
 real(dp) :: djjj
 real(dp) :: al,nrmsum,uurcut,fj,rcmx
 real(dp), allocatable :: rhomodl(:,:),dmat(:,:),rabl(:)
 real(dp),allocatable :: rhol(:),rl(:),vkbl(:,:,:,:),vpl(:,:),uual(:,:,:)
 character(len=5) :: lnames
 character(len=2) :: pspd(3)
 character(len=19) :: version

 lnames = "SPDFG"
 if(nkap==1) then
   version='scalar-relativistic'
 else
   version='fully-relativistic'
 end if
 nrl = nrlin

! adjust nrl to properly accomodate atomic orbitals
 al = log(rr(2)/rr(1))
 uurcut = 0.d0
 do ii=1,nv
   l1 = la(nc+ii)
   if(nkap==1 .or. l1==0) then
     mkap=1
   else
     mkap=2
   end if
   do ikap=1,mkap
     nrmsum = 0.d0
!the search starts at mmax-1 in the sr code and at mmax in the fr code
     if(nkap==1) then
       jstart=mmax-1
     else
       jstart=mmax
     end if
     do jj=jstart,1,-1
       nrmsum = nrmsum + (uupsa(jj,ikap,ii)**2) * rr(jj)*al
       if(nrmsum > 1.d-6) then
         exit  ! Cutoff radius such that uu norm accurate to 10^-6
       end if
     end do
     jj=max(jj,1)
     if(rr(jj) > uurcut) uurcut = rr(jj)
   end do
 end do
 if(uurcut > drl*dble(nrl-1)) then
   nrl = 1 + int(uurcut/drl)
   if(mod(nrl,2)/=0) nrl=nrl+1
   write(6,'(a,i5,a,f10.5)') "Updating nrl = ", nrl, " for uurcut = ", uurcut
 end if

 istr=0
 if(upf_mesh>0) then
! log output mesh: the points of rr out to the extent of the linear mesh,
! every istr-th of them, istr the smallest stride leaving at most upf_mesh
   nlog=1
   do ii=1,mmax
     if(rr(ii)>drl*dble(nrl-1)) exit
     nlog=ii
   end do
   istr=(nlog-1)/upf_mesh+1
   nrl=(nlog-1)/istr+1
 end if

 allocate(rhol(nrl),rl(nrl),vkbl(nrl,mxprj,mxl+1,nkap),vpl(nrl,mxl+2), &
&          rhomodl(nrl,nrhomod), &
&          uual(nrl,nkap,nv),rabl(nrl))

! everything on the linear output mesh or the log mesh points, the
! projectors and the local potential from their pieces (outmesh_m)

 if(istr==0) then
   do  ii=1,nrl
     rl(ii)=drl*dble(ii-1)
     rabl(ii)=drl
   end do
 else
   do  ii=1,nrl
     rl(ii)=rr(1+istr*(ii-1))
     rabl(ii)=dble(istr)*al*rl(ii)
   end do
 end if
!
 vpl(:,:)=0.0d0
 call om_local_values(rr,vpuns(1,lloc+1),mmax,vlocpw,rl,vpl(1,lloc+1),nrl)
! the value at the origin that of the first point of the log mesh (in the
! sr code only)
 if(nkap==1 .and. istr==0) vpl(1,lloc+1)=vpuns(1,lloc+1)

 do l1=1,lmax+1
   if(l1==lloc+1) cycle
   ll=l1-1
   if(nkap==1 .or. ll==0) then
     mkap=1
   else
     mkap=2
   end if
! loop on J = ll +/- 1/2
   do ikap=1,mkap
     do jj=1,nproj(l1)

       call om_piece_values(rr,mmax,vkbpw(jj,l1,ikap),rl, &
&                           vkbl(1,jj,l1,ikap),nrl)

     end do !jj
   end do !ikap
 end do !l1

 call om_values(rr,rho,mmax,rl,rhol,nrl)

 do jj=1,5
   call om_piece_values(rr,mmax,corepw(jj),rl,rhomodl(1,jj),nrl)
 end do

 nwfc = 0
 do ii=1,nv
   l1 = la(nc+ii)
   if(nkap==1 .or. l1==0) then
     mkap=1
   else
     mkap=2
   end if
   nwfc = nwfc + mkap
   do ikap=1,mkap
     call om_values(rr,uupsa(1,ikap,ii),mmax,rl,uual(1,ikap,ii),nrl)
   end do
 end do

 call date_and_time(VALUES=dtime)
 ii=dtime(1)-2000
 if(ii<10) then
   write(pspd(1),'(a,i1)') '0',ii
 else
   write(pspd(1),'(i2)') ii
 end if
 ii=dtime(2)
 if(ii<10) then
   write(pspd(2),'(a,i1)') '0',ii
 else
   write(pspd(2),'(i2)') ii
 end if
 ii=dtime(3)
 if(ii<10) then
   write(pspd(3),'(a,i1)') '0',ii
 else
   write(pspd(3),'(i2)') ii
 end if

! section for upf output for pwscf

 write(6,'(/a)') 'Begin PSP_UPF'

 write(6,'(a,/a)') &
&      '<UPF version="2.0.1">', &
&      '  <PP_INFO>'
 write(6,'(/t2,a/t2,a/t2,a/t2,a/t2,a/t2,a//)') &
&      'This pseudopotential file has been produced using the code', &
&      'ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&      trim(version)//' version 4.0.1 06/20/2107 by D. R. Hamann', &
&      'The code is available through a link at URL www.mat-simresearch.com.', &
&      'Documentation with the package provides a full discription of the', &
&      'input data below.'

 write(6,'(t2,a/t2,a/t2,a/t2,a//)') &
&      'While it is not required under the terms of the GNU GPL, it is',&
&      'suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)',&
&      'in any publication using these pseudopotentials.'

 write(6,'(a)') &
&      '    <PP_INPUTFILE>'

! output printing (echos input data)

 call write_input_echo(nkap,2,atsym,zz,nc,nv,iexc,psfile, &
&                       na,la,fa,lmax,rc0,ep(:,1),ncon,nbas,qcut,ep_as_given, &
&                       well_onset,well_tail_free,lloc,lpopt,dvloc0, &
&                       nproj,debl, &
&                       icmod,fcfact,rcfact,slack,fcgrid,rcgrid, &
&                       epsh1,epsh2,depsh,rxpsh,rlmax,drl, &
&                       ncnf,nvcnf,nacnf,lacnf,facnf,log_mesh,upf_mesh)

 write(6,'(a)') &
&      '    </PP_INPUTFILE>'

 write(6,'(a,3(/a))') &
&      '  </PP_INFO>', &
&      '  <!--                               -->', &
&      '  <!-- END OF HUMAN READABLE SECTION -->', &
&      '  <!--                               -->'

 write(6,'(a)') &
&      '    <PP_HEADER'

 write(6,'(t8,a)') &
&      'generated="Generated using ONCVPSP code by D. R. Hamann"'
 write(6,'(t8,a)') &
&      'author="anonymous"'
 write(6,'(t8,5a)') &
&      'date="',pspd(:),'"'
 write(6,'(t8,a)') &
&      'comment=""'
 write(6,'(t8,a,a2,a)') &
&      'element="',atsym,'"'
 write(6,'(t8,a)') &
&      'pseudo_type="NC"'
 if(nkap==1) then
   write(6,'(t8,a)') &
&        'relativistic="scalar"'
 else
   write(6,'(t8,a)') &
&        'relativistic="full"'
 end if
 write(6,'(t8,a)') &
&      'is_ultrasoft="F"'
 write(6,'(t8,a)') &
&      'is_paw="F"'
 write(6,'(t8,a)') &
&      'is_coulomb="F"'
 if(nkap==1) then
   write(6,'(t8,a)') &
&        'has_so="F"'
 else
   write(6,'(t8,a)') &
&        'has_so="T"'
 end if
 write(6,'(t8,a)') &
&      'has_wfc="F"'
 write(6,'(t8,a)') &
&      'has_gipaw="F"'

 if(icmod>=1) then
   write(6,'(t8,a)') &
&      'core_correction="T"'
 else
   write(6,'(t8,a)') &
&      'core_correction="F"'
 end if

 if(iexc==2) then
   write(6,'(t8,a)') &
&      'functional="SLA  HL   NOGX NOGC"'

 else if(iexc==-000004) then
   write(6,'(t8,a)') &
&      'functional="NOX+HL"'

 else if(iexc==3 .or. iexc==-001009) then
   write(6,'(t8,a)') &
&      'functional="PZ"'

 else if(iexc==4 .or. iexc==-101130) then
   write(6,'(t8,a)') &
&      'functional="PBE"'

 else if(iexc==6 .or. iexc==-001007) then
   write(6,'(t8,a)') &
&      'functional="SLA  VWN  NOGX NOGC"'

!QE's relativistic Slater exchange (RXC)
 else if(iexc==7 .or. iexc==-532007) then
   write(6,'(t8,a)') &
&      'functional="RXC  VWN  NOGX NOGC"'

 else if(iexc==-001012) then
   write(6,'(t8,a)') &
&      'functional="SLA  PW   NOGX NOGC"'

 else if(iexc==-109134) then
   write(6,'(t8,a)') &
&      'functional="PW91"'

 else if(iexc==-116133) then
   write(6,'(t8,a)') &
&      'functional="PBESOL"'

 else if(iexc==-102130) then
   write(6,'(t8,a)') &
&      'functional="REVPBE"'

 else if(iexc==-106132) then
   write(6,'(t8,a)') &
&      'functional="BP"'

 else if(iexc==-106131) then
   write(6,'(t8,a)') &
&      'functional="BLYP"'

 else if(iexc==-118130) then
   write(6,'(t8,a)') &
&      'functional="WC"'

 else
   write(6,'(t8,a,i0,a)') &
&      'upfout: ERROR iexc = ',iexc,' is presently unsupported for UPF output'
   error stop 1
 end if

 write(6,'(t8,a,f8.2,a)') &
&      'z_valence="',zion,'"'

 write(6,'(t8,a,1p,e24.16,a)') &
&      'total_psenergy="',2*epstot,'"'

 write(6,'(t8,a,1p,e24.16,a)') &
&      'rho_cutoff="',rl(nrl),'"'

 write(6,'(t8,a,i1,a)') &
&      'l_max="',lmax,'"'

 if(lloc==4) then
   write(6,'(t8,a)') &
&      'l_local="-1"'
 else
   write(6,'(t8,a,i1,a)') &
&      'l_local="',lloc,'"'
 end if

! calculate total number of projectors and maximum linear mesh point for
! projectors

 ntotproj=0
 nrlproj=0
 do l1=1,lmax+1
   if(l1/=lloc+1) then
     if(nkap==1 .or. l1==1) then
       ntotproj=ntotproj+nproj(l1)
     else
       ntotproj=ntotproj+2*nproj(l1)
     end if
   end if
   nrlproj=max(nrlproj,4+int(rc(l1)/drl))
 end do
 allocate(dmat(ntotproj,ntotproj))

 if(mod(nrlproj,2)/=0) nrlproj=nrlproj+1
!make divisible by 4 to use for 0 compression below
 if(mod(nrlproj,4)/=0) nrlproj=nrlproj+2

!on the log mesh, the first point beyond every core radius (where the
!projectors are zero), made divisible by 4 as above
 if(istr>0) then
   rcmx=0.0d0
   do l1=1,lmax+1
     rcmx=max(rcmx,rc(l1))
   end do
   nrlproj=nrl
   do ii=nrl,1,-1
     if(rl(ii)<=rcmx) exit
     nrlproj=ii
   end do
   nrlproj=min(4*((nrlproj+3)/4),nrl)
 end if

 write(6,'(t8,a,i6,a)') &
&      'mesh_size="',nrl,'"'

 write(6,'(t8,a,i0,a)') &
&      'number_of_wfc="',nwfc,'"'

 write(6,'(t8,a,i0,a)') &
&      'number_of_proj="',ntotproj,'"/>' !end of PP_HEADER

 write(6,'(t2,a)') &
&      '<PP_MESH>'

 write(6,'(t4,a,i4,a)') &
&      '<PP_R type="real"  size="',nrl,'" columns="4">'

 write(6,'(1p,4e24.16)') (rl(ii),ii=1,nrl)

 write(6,'(t4,a)') &
&      '</PP_R>'

 write(6,'(t4,a,i4,a)') &
&      '<PP_RAB type="real"  size="',nrl,'" columns="4">'

 if(istr==0) then
   write(6,'(1p,4e24.16)') (drl,ii=1,nrl)
 else
   write(6,'(1p,4e24.16)') (rabl(ii),ii=1,nrl)
 end if

 write(6,'(t4,a)') &
&      '</PP_RAB>'

 write(6,'(t2,a)') &
&      '</PP_MESH>'

! write local potential with factor of 2 for Rydberg units
 write(6,'(a,i4,a)') &
&      '  <PP_LOCAL type="real"  size="',nrl,'" columns="4">'

 write(6,'(1p,4e24.16)') (2.0d0*vpl(ii,lloc+1),ii=1,nrl)

 write(6,'(a)') &
&      '  </PP_LOCAL>'

! loop on angular mommentum for projector outputs

 write(6,'(t2,a)') &
&      '<PP_NONLOCAL>'

 dmat(:,:)=0.0d0
 iproj=0

 do l1=1,lmax+1
   if(l1==lloc+1) cycle
   ll=l1-1
   if(nkap==1 .or. ll==0) then
     mkap=1
   else
     mkap=2
   end if
   do jj=1,nproj(l1)
! loop on J = ll +/- 1/2
     do ikap=mkap,1,-1
       if(ikap==1) then
         djjj=ll+0.5d0
       else
         djjj=ll-0.5d0
       end if
       iproj=iproj+1
       dmat(iproj,iproj)=2.0d0*evkb(jj,l1,ikap) !2 for Rydbergs
       if(iproj<=9) then
         write(6,'(t4,a,i1)') &
&              '<PP_BETA.',iproj
       else
         write(6,'(t4,a,i2)') &
&              '<PP_BETA.',iproj
       end if
       write(6,'(t8,a)') &
&            'type="real"'
       write(6,'(t8,a,i4,a)') &
&            'size="',nrl,'"'
       write(6,'(t8,a)') &
&            'columns="4"'
       write(6,'(t8,a,i0,a)') &
&           'index="',iproj,'"'
       write(6,'(t8,a,i1,a)') &
&            'angular_momentum="',l1-1,'"'
       write(6,'(t8,a,i4,a)') &
&            'cutoff_radius_index="',nrlproj,'"'
       if(istr==0) then
         write(6,'(t8,a,1p,e24.16,a)') &
&              'cutoff_radius="',(nrlproj-1)*drl,'" >'
       else
         write(6,'(t8,a,1p,e24.16,a)') &
&              'cutoff_radius="',rl(nrlproj),'" >'
       end if

       write(6,'(1p,4e24.16)') (vkbl(ii,jj,l1,ikap),ii=1,nrlproj)
       write(6,'(1p,4f3.0)') (0.d0,ii=nrlproj+1,nrl)

       if(iproj<=9) then
         write(6,'(t4,a,i1,a)') &
&              '</PP_BETA.',iproj,'>'
       else
         write(6,'(t4,a,i2,a)') &
&              '</PP_BETA.',iproj,'>'
       end if
     end do !ikap
   end do !jj (1,nproj(l1))
 end do !l1

 write(6,'(t4,a,i4,a)') &
&      '<PP_DIJ type="real"  size="',ntotproj**2,'" columns="4">'

 write(6,'(1p,4e24.16)') ((dmat(ii,jj),ii=1,ntotproj),jj=1,ntotproj)

 write(6,'(t4,a)') &
&      '</PP_DIJ>'

 write(6,'(t2,a)') &
&      '</PP_NONLOCAL>'

 write(6,'(t2,a)') &
&      '<PP_PSWFC>'
 nwfc = 0
 do ii=1,nv
   l1 = la(nc+ii)
   if(nkap==1 .or. l1==0) then
     mkap=1
   else
     mkap=2
   end if
   do ikap=1,mkap
     nwfc = nwfc + 1
     if(nkap==1) then
       fj=fa(nc+ii)
     else if(ikap==1) then
       fj=dble(l1+1)*fa(nc+ii)/dble(2*l1+1)
     else
       fj=dble(l1)*fa(nc+ii)/dble(2*l1+1)
     end if
     if(nwfc <= 9) then
       write(6,'(t4,a,i1)') &
&            '<PP_CHI.',nwfc
     else
       write(6,'(t4,a,i2)') &
&            '<PP_CHI.',nwfc
     end if
     write(6,'(t8,a)') &
&          'type="real"'
     write(6,'(t8,a,i4,a)') &
&          'size="',nrl,'"'
     write(6,'(t8,a)') &
&          'columns="4"'
     if(nwfc <= 9) then
       write(6,'(t8,a,i1,a)') &
&            'index="',nwfc,'"'
     else
       write(6,'(t8,a,i2,a)') &
&            'index="',nwfc,'"'
     end if
     write(6,'(t8,a,f6.3,a)') &
&          'occupation="',fj,'"'
     write(6,'(t8,a,1p,e24.16,a)') &
&          'pseudo_energy="',2*ea(nc+ii,ikap),'"'
     write(6,'(t8,a,i1,a,a)') &
&          'label="',na(nc+ii),lnames(l1+1:l1+1),'"'
     write(6,'(t8,a,i1,a)') &
&             'l="',l1,'" >'

     write(6,'(1p,4e24.16)') (uual(jj,ikap,ii),jj=1,nrl)

     if(nwfc <= 9) then
       write(6,'(t4,a,i1,a)') &
&            '</PP_CHI.',nwfc,'>'
     else
       write(6,'(t4,a,i2,a)') &
&            '</PP_CHI.',nwfc,'>'
     end if
   end do
 end do
 write(6,'(t2,a)') &
&      '</PP_PSWFC>'

 if(icmod>=1) then
   write(6,'(t2,a,i4,a)') &
&        '<PP_NLCC type="real"  size="',nrl,'" columns="4">'

   write(6,'(1p,4e24.16)') (rhomodl(ii,1)/(4.0d0*pi),ii=1,nrl)

   write(6,'(t2,a)') &
&        '</PP_NLCC>'
 end if

 write(6,'(t2,a,i4,a)') &
&      '<PP_RHOATOM type="real"  size="',nrl,'" columns="4">'

 write(6,'(1p,4e24.16)') ((rl(ii)**2)*rhol(ii),ii=1,nrl)

 write(6,'(t2,a)') &
&      '</PP_RHOATOM>'

! spin-orbit data, fully-relativistic only
 if(nkap==2) then
   write(6,'(t2,a)') &
&      '<PP_SPIN_ORB>'

   iproj=0
   do l1=1,lmax+1
     if(l1==lloc+1) cycle
     ll=l1-1
     if(ll==0) then
       mkap=1
     else
       mkap=2
     end if
! loop on J = ll +/- 1/2
     do jj=1,nproj(l1)
       do ikap=mkap,1,-1
         if(ikap==1) then
           djjj=ll+0.5d0
         else
           djjj=ll-0.5d0
         end if
         iproj=iproj+1
         if(iproj<=9) then
           write(6,'(t4,a,i1,a,i1,a,i1,a,f3.1,a)') &
&              '<PP_RELBETA.',iproj,'  index="',iproj,'"  lll="',ll, &
&              '" jjj="',djjj,'"/>'
         else
           write(6,'(t4,a,i2,a,i2,a,i1,a,f3.1,a)') &
&              '<PP_RELBETA.',iproj,' index="',iproj,'" lll="',ll, &
&              '" jjj="',djjj,'"/>'
         end if
       end do !ikap
     end do !jj (1,nproj(l1))
   end do !l1

   nwfc = 0
   do ii=1,nv
     l1 = la(nc+ii)
     if(l1==0) then
       mkap=1
     else
       mkap=2
     end if
     do ikap=1,mkap
       nwfc = nwfc + 1
       if(ikap==1) then
         djjj=l1+0.5d0
       else
         djjj=l1-0.5d0
       end if
       if(nwfc <= 9) then
         write(6,'(t4,a,i1,a,i1,a,i1,a,f3.1,a,i1,a)') &
&            '<PP_RELWFC.',nwfc,'  index="',nwfc,'"  lchi="',l1, &
&            '" jchi="',djjj,'" nn="',ii,'"/>'
       else
         write(6,'(t4,a,i2,a,i2,a,i1,a,f3.1,a,i1,a)') &
&            '<PP_RELWFC.',nwfc,'  index="',nwfc,'"  lchi="',l1, &
&            '" jchi="',djjj,'" nn="',ii,'"/>'
       end if
     end do
   end do

   write(6,'(t2,a)') &
&      '</PP_SPIN_ORB>'
 end if

 write(6,'(a)') &
&      '</UPF>'

! write termination signal
 write(6,'(/a)') 'END_PSP'

 deallocate(rhol,rl,vkbl,vpl,rhomodl,dmat,rabl)

 return
 end subroutine upfout

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

 module upfout_m

 implicit none
 private
 public :: upfout

 contains

 subroutine upfout(inp,c_light,vkbpw,evkb,rr,vpuns,vlocpw,rho,corepw,zion, &
&                  mmax,nrl_in,epstot,rc_in,fcfact_in,fcgrid,rcgrid,uupsa,ea, &
&                  log_mesh,nkap,iout,ifail)

!! the pseudopotential in the UPF file format of PWSCF, written to iout,
!! with the spin-orbit data when fully relativistic: interpolated to the
!! linear output mesh, or taken at every istr-th point of the log mesh
!! (upf_mesh > 0), as generated, with nothing interpolated
!!
!! @param inp  the input values (input_m)
!! @param c_light  speed of light (Hartree atomic units)
!! @param vkbpw  the projectors as pieces, last index 1 for j=l+1/2, 2 for
!!               j=l-1/2 when nkap=2
!! @param evkb  coefficients of the projectors (Ha), last index as for vkbpw
!! @param rr  log radial mesh (bohr)
!! @param vpuns  unscreened semi-local pseudopotentials (vp(:,5) is local
!!               potential if linear combination is used)
!! @param vlocpw  the screened local potential as pieces; vpuns(:,lloc+1)
!!                less their node values is smooth
!! @param rho  valence pseudocharge
!! @param corepw  the model core charge and its first 4 derivatives, as
!!                pieces split at its joint (modcore_m, core_pieces)
!! @param zion  ion charge
!! @param mmax  size of the log radial mesh
!! @param nrl_in  size of the linear output mesh (enlarged locally to fit the
!!                orbitals)
!! @param epstot  pseudoatom total energy (Ha)
!! @param rc_in  core radii as read from the input (bohr), rc_in(5) for lloc=4
!! @param fcfact_in  prefactor for model core charge amplitude, as read from the
!!                   input
!! @param fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!! @param rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!! @param uupsa  pseudo-atomic orbitals, second index 1 for j=l+1/2, 2 for
!!               j=l-1/2 when nkap=2
!! @param ea  all-electron eigenvalues (Ha) of the reference configuration,
!!            the orbitals' pseudo_energy; second index as for uupsa
!! @param log_mesh  amesh, rr1, rmax and the Adams order of the log mesh
!!                  line, echoed when the input sets it (else 0)
!! @param nkap  number of kappas per l: 1, or 2 for irel=2
!! @param iout  unit of the text output
!! @param ifail  0, or the kind of failure (failure_m), whose message is
!!               written

 use constants_m, only: dp, pi, mxstate, mxl, mxlin, nrhomod, mxprj
 use kappa_m, only: kap_count, kap_j_share
 use input_echo_m, only: write_input_echo, fe, date_stamp
 use outmesh_m, only: om_values, om_piece_values, om_local_values
 use pieces_m, only: piecewise_t
 use input_m, only: input_t
 use logmesh_m, only: lm_step


!Input variables
 type(input_t), intent(in) :: inp
 real(dp), intent(in) :: c_light
 integer, intent(in) :: mmax,nrl_in,nkap
 real(dp), intent(in) :: fcfact_in,zion,epstot
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)
 real(dp), intent(in) :: rr(mmax),vpuns(mmax,mxl+2),rho(mmax)
 type(piecewise_t), intent(in) :: vkbpw(mxprj,mxl+1,nkap),vlocpw
 type(piecewise_t), intent(in) :: corepw(nrhomod)
 real(dp), intent(in) :: evkb(mxprj,mxl+1,nkap)
 real(dp), intent(in) :: uupsa(mmax,nkap,inp%nv)

!the core radii as read (rc_in), for the echo
 real(dp), intent(in) :: rc_in(mxlin+1)
 real(dp), intent(in) :: ea(mxstate,nkap)
 real(dp), intent(in) :: log_mesh(4)
 integer, intent(in) :: iout

!Output variables - printing only

!Output variables
 integer, intent(out) :: ifail

!Local variables
 integer :: nrl,istr,ntotproj,nrlproj,nwfc
 real(dp) :: al
 real(dp), allocatable :: rhomodl(:,:),rabl(:)
 real(dp),allocatable :: rhol(:),rl(:),vkbl(:,:,:,:),vpl(:,:),uual(:,:,:)
!the projectors in the order of the file, each l but lloc's, each of
!its projectors, and for each one j=l-1/2 before j=l+1/2: l+1, projector
!and ikap
 integer, allocatable :: l1p(:),jjp(:),ikp(:)
 character(len=5) :: lnames
 character(len=2) :: pspd(3)
 character(len=19) :: version

 ifail=0
 lnames = "SPDFG"
 if(nkap==1) then
   version='scalar-relativistic'
 else
   version='fully-relativistic'
 end if
 nrl = nrl_in
 al=lm_step(rr,mmax)

 call orbital_extent()
 call log_stride()
 call output_values(ifail)
 if(ifail/=0) return
 call projector_order()
 call date_stamp(pspd)

! section for upf output for pwscf

 write(iout,'(/a)') 'Begin PSP_UPF'

 call write_info()
 call write_header(ifail)
 if(ifail/=0) return
 call write_mesh_local()
 call write_nonlocal()
 call write_pswfc()
 call write_densities()
! spin-orbit data, fully-relativistic only
 if(nkap==2) call write_spin_orb()

 write(iout,'(a)') &
&      '</UPF>'

! write termination signal
 write(iout,'(/a)') 'END_PSP'

 deallocate(rhol,rl,vkbl,vpl,rhomodl,uual,rabl,l1p,jjp,ikp)

 return

 contains

 subroutine orbital_extent()
!! nrl enlarged to accomodate the atomic orbitals: to the radius beyond
!! which each one's norm is below 1e-6
 integer :: ii,jj,ll,ikap,mkap,jstart
 real(dp) :: nrmsum,uurcut
 uurcut = 0.0d0
 do ii=1,inp%nv
   ll = inp%la(inp%nc+ii)
   mkap=kap_count(nkap,ll)
   do ikap=1,mkap
     nrmsum = 0.0d0
!the search starts at mmax-1 in the sr code and at mmax in the fr code
     if(nkap==1) then
       jstart=mmax-1
     else
       jstart=mmax
     end if
     do jj=jstart,1,-1
       nrmsum = nrmsum + (uupsa(jj,ikap,ii)**2) * rr(jj)*al
       if(nrmsum > 1.0d-6) then
         exit  ! Cutoff radius such that uu norm accurate to 10^-6
       end if
     end do
     jj=max(jj,1)
     if(rr(jj) > uurcut) uurcut = rr(jj)
   end do
 end do
 if(uurcut > inp%drl*dble(nrl-1)) then
   nrl = 1 + int(uurcut/inp%drl)
   if(mod(nrl,2)/=0) then
     if(inp%drl*dble(nrl)<=rr(mmax)) then
       nrl=nrl+1
     else
       nrl=nrl-1
     end if
   end if
   write(iout,'(a,i5,a,f10.5)') "Updating nrl = ", nrl, " for uurcut = ", uurcut
 end if
 end subroutine orbital_extent

 subroutine log_stride()
!! istr 0 for the linear output mesh, or with upf_mesh > 0 the log output
!! mesh: the points of rr out to the extent of the linear mesh, every
!! istr-th of them, istr the smallest stride leaving at most upf_mesh;
!! and nrl its points
 integer :: ii,nlog
 istr=0
 if(inp%upf_mesh>0) then
   nlog=1
   do ii=1,mmax
     if(rr(ii)>inp%drl*dble(nrl-1)) exit
     nlog=ii
   end do
   istr=(nlog-1)/inp%upf_mesh+1
   nrl=(nlog-1)/istr+1
 end if
 end subroutine log_stride

 subroutine output_values(ifail)
!! everything on the linear output mesh or the log mesh points, the
!! projectors and the local potential from their pieces (outmesh_m)
 integer :: ii,jj,ll,l1,ikap,mkap
 integer, intent(out) :: ifail
 ifail=0
 allocate(rhol(nrl),rl(nrl),vkbl(nrl,mxprj,mxl+1,nkap),vpl(nrl,mxl+2), &
&          rhomodl(nrl,nrhomod), &
&          uual(nrl,nkap,inp%nv),rabl(nrl))

 if(istr==0) then
   do  ii=1,nrl
     rl(ii)=inp%drl*dble(ii-1)
     rabl(ii)=inp%drl
   end do
 else
   do  ii=1,nrl
     rl(ii)=rr(1+istr*(ii-1))
     rabl(ii)=dble(istr)*al*rl(ii)
   end do
 end if
!
 vpl(:,:)=0.0d0
 call om_local_values(rr,vpuns(1,inp%lloc+1),mmax,vlocpw,rl,vpl(1, &
&                     inp%lloc+1),nrl,iout,ifail)
 if(ifail/=0) return

 do l1=1,inp%lmax+1
   if(l1==inp%lloc+1) cycle
   ll=l1-1
   mkap=kap_count(nkap,ll)
! loop on J = ll +/- 1/2
   do ikap=1,mkap
     do jj=1,inp%nproj(l1)

       call om_piece_values(rr,mmax,vkbpw(jj,l1,ikap),rl, &
&                           vkbl(1,jj,l1,ikap),nrl,iout,ifail)
       if(ifail/=0) return

     end do !jj
   end do !ikap
 end do !l1

 call om_values(rr,rho,mmax,rl,rhol,nrl,iout,ifail)
 if(ifail/=0) return

 do jj=1,5
   call om_piece_values(rr,mmax,corepw(jj),rl,rhomodl(1,jj),nrl,iout,ifail)
   if(ifail/=0) return
 end do

 nwfc = 0
 do ii=1,inp%nv
   ll = inp%la(inp%nc+ii)
   mkap=kap_count(nkap,ll)
   nwfc = nwfc + mkap
   do ikap=1,mkap
     call om_values(rr,uupsa(1,ikap,ii),mmax,rl,uual(1,ikap,ii),nrl,iout,ifail)
     if(ifail/=0) return
   end do
 end do
 end subroutine output_values

 subroutine projector_order()
!! the projectors in the order of the file (l1p, jjp, ikp), ntotproj of
!! them
 integer :: l1,jj,ikap,mkap
 ntotproj=0
 do l1=1,inp%lmax+1
   if(l1/=inp%lloc+1) then
     if(nkap==1 .or. l1==1) then
       ntotproj=ntotproj+inp%nproj(l1)
     else
       ntotproj=ntotproj+2*inp%nproj(l1)
     end if
   end if
 end do
 allocate(l1p(ntotproj),jjp(ntotproj),ikp(ntotproj))
 ntotproj=0
 do l1=1,inp%lmax+1
   if(l1==inp%lloc+1) cycle
   mkap=kap_count(nkap,l1-1)
   do jj=1,inp%nproj(l1)
! loop on J = ll +/- 1/2
     do ikap=mkap,1,-1
       ntotproj=ntotproj+1
       l1p(ntotproj)=l1
       jjp(ntotproj)=jj
       ikp(ntotproj)=ikap
     end do
   end do
 end do
 end subroutine projector_order

 subroutine write_info()
!! PP_INFO: the provenance and the input file
 write(iout,'(a,/a)') &
&      '<UPF version="2.0.1">', &
&      '  <PP_INFO>'
 write(iout,'(/t2,a/t2,a/t2,a/t2,a/t2,a/t2,a//)') &
&      'This pseudopotential file has been produced using the code', &
&      'ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&      trim(version)//' version 4.0.1 06/20/2107 by D. R. Hamann', &
&      'The code is available through a link at URL www.mat-simresearch.com.', &
&      'Documentation with the package provides a full discription of the', &
&      'input data below.'

 write(iout,'(t2,a/t2,a/t2,a/t2,a//)') &
&      'While it is not required under the terms of the GNU GPL, it is',&
&      'suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)',&
&      'in any publication using these pseudopotentials.'

 write(iout,'(a)') &
&      '    <PP_INPUTFILE>'

! output printing (echos input data)

 call write_input_echo(nkap,2,inp%atsym,inp%zz,inp%nc,inp%nv,inp%iexc,c_light, &
&                      inp%psfile, &
&                       inp%na,inp%la,inp%fa,inp%lmax,rc_in,inp%ep(:,1), &
&                       inp%ncon,inp%nbas,inp%qcut,inp%ep_as_given, &
&                       inp%well_onset,inp%well_tail,inp%wave_regularization, &
&                       inp%lloc, &
&                       inp%lpopt,inp%dvloc0, &
&                       inp%nproj,inp%debl,inp%well_n,inp%well_norm, &
&                       inp%ecut_weight, &
&                       inp%icmod,fcfact_in,inp%rcfact,inp%slack,inp%r_s, &
&                       inp%budget_cutoff, &
&                       fcgrid,rcgrid, &
&                       inp%epsh1,inp%epsh2,inp%depsh,inp%rxpsh,inp%rlmax, &
&                       inp%drl, &
&                       inp%ncnf,inp%cnf, &
&                       log_mesh,inp%upf_mesh,iout)

 write(iout,'(a)') &
&      '    </PP_INPUTFILE>'

 write(iout,'(a,3(/a))') &
&      '  </PP_INFO>', &
&      '  <!--                               -->', &
&      '  <!-- END OF HUMAN READABLE SECTION -->', &
&      '  <!--                               -->'
 end subroutine write_info

 subroutine write_header(ifail)
!! PP_HEADER, and nrlproj, the points of the projectors
 integer :: ii,l1
 real(dp) :: rcmx
 character(len=20) :: fname
 integer, intent(out) :: ifail
 ifail=0
 call upf_functional(fname,ifail)
 if(ifail/=0) return
 write(iout,'(a)') &
&      '    <PP_HEADER'

 write(iout,'(t8,a)') &
&      'generated="Generated using ONCVPSP code by D. R. Hamann"'
 write(iout,'(t8,a)') &
&      'author="anonymous"'
 write(iout,'(t8,5a)') &
&      'date="',pspd(:),'"'
 write(iout,'(t8,a)') &
&      'comment=""'
 write(iout,'(t8,a,a2,a)') &
&      'element="',inp%atsym,'"'
 write(iout,'(t8,a)') &
&      'pseudo_type="NC"'
 if(nkap==1) then
   write(iout,'(t8,a)') &
&        'relativistic="scalar"'
 else
   write(iout,'(t8,a)') &
&        'relativistic="full"'
 end if
 write(iout,'(t8,a)') &
&      'is_ultrasoft="F"'
 write(iout,'(t8,a)') &
&      'is_paw="F"'
 write(iout,'(t8,a)') &
&      'is_coulomb="F"'
 if(nkap==1) then
   write(iout,'(t8,a)') &
&        'has_so="F"'
 else
   write(iout,'(t8,a)') &
&        'has_so="T"'
 end if
 write(iout,'(t8,a)') &
&      'has_wfc="F"'
 write(iout,'(t8,a)') &
&      'has_gipaw="F"'

 if(inp%icmod>=1) then
   write(iout,'(t8,a)') &
&      'core_correction="T"'
 else
   write(iout,'(t8,a)') &
&      'core_correction="F"'
 end if

 write(iout,'(t8,3a)') &
&      'functional="',trim(fname),'"'

 write(iout,'(t8,a,f8.2,a)') &
&      'z_valence="',zion,'"'

 write(iout,'(t8,a,1p,e24.16,a)') &
&      'total_psenergy="',2*epstot,'"'

 write(iout,'(t8,a,1p,e24.16,a)') &
&      'rho_cutoff="',rl(nrl),'"'

 write(iout,'(t8,a,i1,a)') &
&      'l_max="',inp%lmax,'"'

 if(inp%lloc==4) then
   write(iout,'(t8,a)') &
&      'l_local="-1"'
 else
   write(iout,'(t8,a,i1,a)') &
&      'l_local="',inp%lloc,'"'
 end if

! the maximum linear mesh point for projectors

 nrlproj=0
 do l1=1,inp%lmax+1
   nrlproj=max(nrlproj,4+int(inp%rc(l1)/inp%drl))
 end do

 if(mod(nrlproj,2)/=0) nrlproj=nrlproj+1
!make divisible by 4 to use for 0 compression below, within the mesh
 if(mod(nrlproj,4)/=0) nrlproj=nrlproj+2
 nrlproj=min(nrlproj,nrl)

!on the log mesh, the first point beyond every core radius (where the
!projectors are zero), made divisible by 4 as above
 if(istr>0) then
   rcmx=0.0d0
   do l1=1,inp%lmax+1
     rcmx=max(rcmx,inp%rc(l1))
   end do
   nrlproj=nrl
   do ii=nrl,1,-1
     if(rl(ii)<=rcmx) exit
     nrlproj=ii
   end do
   nrlproj=min(4*((nrlproj+3)/4),nrl)
 end if

 write(iout,'(t8,a,i6,a)') &
&      'mesh_size="',nrl,'"'

 write(iout,'(t8,a,i0,a)') &
&      'number_of_wfc="',nwfc,'"'

 write(iout,'(t8,a,i0,a)') &
&      'number_of_proj="',ntotproj,'"/>' !end of PP_HEADER
 end subroutine write_header

 subroutine write_mesh_local()
!! PP_MESH, and PP_LOCAL with a factor of 2 for Rydberg units
 integer :: ii
 write(iout,'(t2,a)') &
&      '<PP_MESH>'

 write(iout,'(t4,a,i4,a)') &
&      '<PP_R type="real"  size="',nrl,'" columns="4">'

 write(iout,'(1p,4e24.16)') (rl(ii),ii=1,nrl)

 write(iout,'(t4,a)') &
&      '</PP_R>'

 write(iout,'(t4,a,i4,a)') &
&      '<PP_RAB type="real"  size="',nrl,'" columns="4">'

 if(istr==0) then
   write(iout,'(1p,4e24.16)') (inp%drl,ii=1,nrl)
 else
   write(iout,'(1p,4e24.16)') (fe(rabl(ii)),ii=1,nrl)
 end if

 write(iout,'(t4,a)') &
&      '</PP_RAB>'

 write(iout,'(t2,a)') &
&      '</PP_MESH>'

 write(iout,'(a,i4,a)') &
&      '  <PP_LOCAL type="real"  size="',nrl,'" columns="4">'

 write(iout,'(1p,4e24.16)') (fe(2.0d0*vpl(ii,inp%lloc+1)),ii=1,nrl)

 write(iout,'(a)') &
&      '  </PP_LOCAL>'
 end subroutine write_mesh_local

 subroutine write_nonlocal()
!! PP_NONLOCAL: the projectors and their coefficients, in Rydberg
 integer :: ii,jj,l1,ikap,iproj
 real(dp), allocatable :: dmat(:,:)
 write(iout,'(t2,a)') &
&      '<PP_NONLOCAL>'

 allocate(dmat(ntotproj,ntotproj))
 dmat(:,:)=0.0d0

 do iproj=1,ntotproj
   l1=l1p(iproj)
   jj=jjp(iproj)
   ikap=ikp(iproj)
   dmat(iproj,iproj)=2.0d0*evkb(jj,l1,ikap) !2 for Rydbergs
   write(iout,'(t4,a,i0)') &
&        '<PP_BETA.',iproj
   write(iout,'(t8,a)') &
&        'type="real"'
   write(iout,'(t8,a,i4,a)') &
&        'size="',nrl,'"'
   write(iout,'(t8,a)') &
&        'columns="4"'
   write(iout,'(t8,a,i0,a)') &
&       'index="',iproj,'"'
   write(iout,'(t8,a,i1,a)') &
&        'angular_momentum="',l1-1,'"'
   write(iout,'(t8,a,i4,a)') &
&        'cutoff_radius_index="',nrlproj,'"'
   if(istr==0) then
     write(iout,'(t8,a,1p,e24.16,a)') &
&          'cutoff_radius="',(nrlproj-1)*inp%drl,'" >'
   else
     write(iout,'(t8,a,1p,e24.16,a)') &
&          'cutoff_radius="',rl(nrlproj),'" >'
   end if

   write(iout,'(1p,4e24.16)') (fe(vkbl(ii,jj,l1,ikap)),ii=1,nrlproj)
   write(iout,'(1p,4f3.0)') (0.0d0,ii=nrlproj+1,nrl)

   write(iout,'(t4,a,i0,a)') &
&        '</PP_BETA.',iproj,'>'
 end do

 write(iout,'(t4,a,i4,a)') &
&      '<PP_DIJ type="real"  size="',ntotproj**2,'" columns="4">'

 write(iout,'(1p,4e24.16)') ((fe(dmat(ii,jj)),ii=1,ntotproj),jj=1,ntotproj)

 write(iout,'(t4,a)') &
&      '</PP_DIJ>'

 write(iout,'(t2,a)') &
&      '</PP_NONLOCAL>'
 deallocate(dmat)
 end subroutine write_nonlocal

 subroutine write_pswfc()
!! PP_PSWFC: the pseudo-atomic orbitals, of each J when fully relativistic
 integer :: ii,jj,ll,ikap,mkap,iwfc
 real(dp) :: fj
 write(iout,'(t2,a)') &
&      '<PP_PSWFC>'
 iwfc = 0
 do ii=1,inp%nv
   ll = inp%la(inp%nc+ii)
   mkap=kap_count(nkap,ll)
   do ikap=1,mkap
     iwfc = iwfc + 1
     fj=inp%fa(inp%nc+ii)
     if(nkap==2) fj=kap_j_share(ll,ikap,inp%fa(inp%nc+ii))
     write(iout,'(t4,a,i0)') &
&          '<PP_CHI.',iwfc
     write(iout,'(t8,a)') &
&          'type="real"'
     write(iout,'(t8,a,i4,a)') &
&          'size="',nrl,'"'
     write(iout,'(t8,a)') &
&          'columns="4"'
     write(iout,'(t8,a,i0,a)') &
&          'index="',iwfc,'"'
     write(iout,'(t8,a,f6.3,a)') &
&          'occupation="',fj,'"'
     write(iout,'(t8,a,1p,e24.16,a)') &
&          'pseudo_energy="',2*ea(inp%nc+ii,ikap),'"'
     write(iout,'(t8,a,i1,a,a)') &
&          'label="',inp%na(inp%nc+ii),lnames(ll+1:ll+1),'"'
     write(iout,'(t8,a,i1,a)') &
&             'l="',ll,'" >'

     write(iout,'(1p,4e24.16)') (fe(uual(jj,ikap,ii)),jj=1,nrl)

     write(iout,'(t4,a,i0,a)') &
&          '</PP_CHI.',iwfc,'>'
   end do
 end do
 write(iout,'(t2,a)') &
&      '</PP_PSWFC>'
 end subroutine write_pswfc

 subroutine write_densities()
!! PP_NLCC, the model core with icmod >= 1, and PP_RHOATOM
 integer :: ii
 if(inp%icmod>=1) then
   write(iout,'(t2,a,i4,a)') &
&        '<PP_NLCC type="real"  size="',nrl,'" columns="4">'

   write(iout,'(1p,4e24.16)') (fe(rhomodl(ii,1)/(4.0d0*pi)),ii=1,nrl)

   write(iout,'(t2,a)') &
&        '</PP_NLCC>'
 end if

 write(iout,'(t2,a,i4,a)') &
&      '<PP_RHOATOM type="real"  size="',nrl,'" columns="4">'

 write(iout,'(1p,4e24.16)') (fe((rl(ii)**2)*rhol(ii)),ii=1,nrl)

 write(iout,'(t2,a)') &
&      '</PP_RHOATOM>'
 end subroutine write_densities

 subroutine write_spin_orb()
!! PP_SPIN_ORB: l and j of each projector and each orbital
 integer :: ii,ll,ikap,mkap,iproj,iwfc
 real(dp) :: djjj
 write(iout,'(t2,a)') &
&    '<PP_SPIN_ORB>'

 do iproj=1,ntotproj
   ll=l1p(iproj)-1
   if(ikp(iproj)==1) then
     djjj=ll+0.5d0
   else
     djjj=ll-0.5d0
   end if
   if(iproj<=9) then
     write(iout,'(t4,a,i1,a,i1,a,i1,a,f3.1,a)') &
&        '<PP_RELBETA.',iproj,'  index="',iproj,'"  lll="',ll, &
&        '" jjj="',djjj,'"/>'
   else
     write(iout,'(t4,a,i2,a,i2,a,i1,a,f3.1,a)') &
&        '<PP_RELBETA.',iproj,' index="',iproj,'" lll="',ll, &
&        '" jjj="',djjj,'"/>'
   end if
 end do

 iwfc = 0
 do ii=1,inp%nv
   ll = inp%la(inp%nc+ii)
   mkap=kap_count(2,ll)
   do ikap=1,mkap
     iwfc = iwfc + 1
     if(ikap==1) then
       djjj=ll+0.5d0
     else
       djjj=ll-0.5d0
     end if
     if(iwfc <= 9) then
       write(iout,'(t4,a,i1,a,i1,a,i1,a,f3.1,a,i1,a)') &
&          '<PP_RELWFC.',iwfc,'  index="',iwfc,'"  lchi="',ll, &
&          '" jchi="',djjj,'" nn="',ii,'"/>'
     else
       write(iout,'(t4,a,i2,a,i2,a,i1,a,f3.1,a,i1,a)') &
&          '<PP_RELWFC.',iwfc,'  index="',iwfc,'"  lchi="',ll, &
&          '" jchi="',djjj,'" nn="',ii,'"/>'
     end if
   end do
 end do

 write(iout,'(t2,a)') &
&    '</PP_SPIN_ORB>'
 end subroutine write_spin_orb


 subroutine upf_functional(name,ifail)
!! the functional attribute of the UPF header for iexc; an error for one
!! UPF has no name for
 use failure_m, only: fail_input, failed, msg_upf_iexc
 character(len=20), intent(out) :: name
 integer, intent(out) :: ifail

 ifail=0
 select case(inp%iexc)
 case(2)
   name='SLA  HL   NOGX NOGC'
 case(-000004)
   name='NOX+HL'
 case(3,-001009)
   name='PZ'
 case(4,-101130)
   name='PBE'
 case(6,-001007)
   name='SLA  VWN  NOGX NOGC'
!QE's relativistic Slater exchange (RXC)
 case(7,-532007)
   name='RXC  VWN  NOGX NOGC'
 case(-001012)
   name='SLA  PW   NOGX NOGC'
 case(-109134)
   name='PW91'
 case(-116133)
   name='PBESOL'
 case(-102130)
   name='REVPBE'
 case(-106132)
   name='BP'
 case(-106131)
   name='BLYP'
 case(-118130)
   name='WC'
 case default
   write(iout,'(t8,a,i0,a)') &
&      'upfout: ERROR iexc = ',inp%iexc, &
&      ' is presently unsupported for UPF output'
   ifail=failed(fail_input,msg_upf_iexc)
   return
 end select
 end subroutine upf_functional

 end subroutine upfout

 end module upfout_m

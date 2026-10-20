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
! for Abinit input using pspcod=8, with the spin-orbit projectors when
! fully relativistic

 subroutine linout(lmax,lloc,rc,vkb,evkb,vsr,esr,vso,eso,nproj,rr,vpuns, &
&                  rho,rhomod, &
&                  rhotae,rhoc,zz,zion,mmax,mxprj,iexc,icmod,nrl,drl, &
&                  atsym,soscale, &
&                  na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&                  fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,fcfact0,rcfact, &
&                  slack,fcgrid,rcgrid, &
&                  epsh1,epsh2,depsh,rxpsh,rlmax,psfile,log_mesh,upf_mesh, &
&                  nkap)


!! @param lmax  maximum angular momentum
!! @param lloc  l for local potential
!! @param rc  core radii
!! @param vkb  VKB projectors, written when nkap=1
!! @param evkb  coefficients of VKB projectors
!! @param vsr  normalized scalar-relativistic projectors, written with vso when
!!             nkap=2
!! @param esr  energy coefficients of vsr
!! @param vso  normalized spin-orbit projectors
!! @param eso  energy coefficients of vso
!! @param nproj  number of vkb projectors for each l
!! @param rr  log radial grid
!! @param vpuns  unscreened semi-local pseudopotentials (vpuns(:,5) is local
!!               potential if linear combination is used)
!! @param rho  valence pseudocharge
!! @param rhotae  all-electron valence charge
!! @param rhoc  all-electron core charge
!! @param rhomod  model core charge
!! @param zz  atomic number
!! @param zion  at this point, total valence charge (becomes psuedoion charge)
!! @param mmax  size of log radial grid
!! @param mxprj  dimension of number of projectors
!! @param iexc  type of exchange-correlation
!! @param icmod  1 if model core charge is used, otherwise 0
!! @param fcfact  model core amplitude prefactor (1 for icmod=4, cutoff for
!!                icmod=6)
!! @param fcfact0  input fcfact, for the echo
!! @param rcfact  prefactor for the model core charge scale
!! @param nrl  size of linear radial grid
!! @param drl  spacing of linear radial grid
!! @param atsym  atomic symbol
!! @param soscale  possible factor to boost or cut spin-orbit strength, not in
!!                 use yet
!! @param nkap  number of kappa per l: 1, or 2 for the fully-relativistic code
!!
!! (the arguments from na on are input variables echoed in the psp8 file)
!!
!! @param na  principal quantum numbers, reference configuration
!! @param la  angular momenta, reference configuration
!! @param ncon  number of matching constraints for each l
!! @param nbas  number of basis functions for each l
!! @param nvcnf  number of valence states, all configurations
!! @param nacnf  principal quantum numbers, all configurations
!! @param lacnf  angular momenta, all configurations
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param lpopt  choice of polynomial for lloc==4
!! @param ncnf  number of test configurations
!! @param fa  occupancies, reference configuration
!! @param rc0  input core radii
!! @param ep  energies at which the pseudopotentials are generated (Ha), each
!!            kap
!! @param qcut  residual-energy cutoff wave vectors for each l
!! @param debl  energy shifts for the second projectors for each l (Ha)
!! @param facnf  occupancies, all configurations
!! @param dvloc0  amplitude at rr==0 of the function added to the local
!!                potential for lloc==4
!! @param slack  share of the cutoff the model core must fit (icmod=6, 7)
!! @param fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!! @param rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!! @param epsh1  lower limit of the log derivative analysis energy range (Ha)
!! @param epsh2  upper limit of the log derivative analysis energy range (Ha)
!! @param depsh  energy step of the log derivative analysis (Ha)
!! @param rxpsh  radius (bohr) of the log derivative analysis, negative for the
!!               default
!! @param rlmax  maximum radius of the linear output mesh (bohr)
!! @param psfile  output format, should be 'psp8' or 'both'
!! @param log_mesh  amesh, rr1, rmax of the log mesh, echoed when the input sets
!!                  it (else 0)
!! @param upf_mesh  most points of the log mesh of the UPF file, 0 for the
!!                  linear output mesh (echoed when not 0)


 use constants_m, only: dp
 use input_echo_m, only: write_input_echo
 implicit none

!Input variables
 integer, intent(in) :: lmax,lloc,iexc,mmax,mxprj,nrl,icmod,nkap
 integer, intent(in) :: nproj(6)
 real(dp), intent(in) :: drl,fcfact,fcfact0,rcfact,zz,zion,soscale
 real(dp), intent(in) :: slack
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)
 real(dp), intent(in) :: rr(mmax),vpuns(mmax,5),rho(mmax)
 real(dp), intent(in) :: rhotae(mmax),rhoc(mmax)
 real(dp), intent(in) :: vkb(mmax,mxprj,4,nkap),evkb(mxprj,4,nkap)
 real(dp), intent(in) :: vsr(mmax,2*mxprj,4),vso(mmax,2*mxprj,4)
 real(dp), intent(in) :: esr(2*mxprj,4),eso(2*mxprj,4)
 real(dp), intent(in) :: rhomod(mmax,5)
 real(dp), intent(in) :: rc(6)
 character(len=2), intent(in) :: atsym

!additional input for psp8 output to echo input file, all as defined
! in the main progam
 integer, intent(in) :: na(30),la(30),ncon(6),nbas(6)
 integer, intent(in) :: nvcnf(5),nacnf(30,5),lacnf(30,5)
 integer, intent(in) :: nc,nv,lpopt,ncnf,upf_mesh
 real(dp), intent(in) :: fa(30),rc0(6),ep(6,nkap),qcut(6),debl(6),facnf(30,5)
 real(dp), intent(in) :: dvloc0,epsh1,epsh2,depsh,rxpsh,rlmax
 real(dp), intent(in) :: log_mesh(3)
 character(len=4), intent(in) :: psfile

!Local variables
 integer :: ii,jj,ll,l1,ixc_abinit
 integer :: dtime(8),npr_sr(5),npr_so(5)
 real(dp) :: zero(20)
 real(dp) :: fchrg
 real(dp), allocatable :: rhomodl(:,:)
 real(dp),allocatable :: rhol(:),rl(:),vpl(:,:)
 real(dp),allocatable :: rhotael(:),rhocl(:)
 real(dp), allocatable :: vkbl(:,:,:),vsrl(:,:,:),vsol(:,:,:)
 character(len=2) :: pspd(3)
 character(len=48) :: version

 allocate(rhol(nrl),rl(nrl),vpl(nrl,5),rhomodl(nrl,5))
 allocate(rhotael(nrl),rhocl(nrl))
 allocate(vkbl(nrl,mxprj,4))
 allocate(vsrl(nrl,2*mxprj,4),vsol(nrl,2*mxprj,4))

! set up projector number for sr_so calculations based on non-zero
! coefficients
 npr_sr(:)=0 ; npr_so(:)=0
 if(nkap==2) then
   do l1=1,lmax+1
     do ii=1,2*nproj(l1)
       if(abs(esr(ii,l1))>0.0d0) npr_sr(l1)=npr_sr(l1)+1
       if(abs(eso(ii,l1))>0.0d0) npr_so(l1)=npr_so(l1)+1
     end do
     write(6,'(a,3i4)') 'l1-1,npr_sr,npr_so',l1-1,npr_sr(l1),npr_so(l1)
   end do
 end if

! interpolation of everything onto linear output mesh
!
 do  ii=1,nrl
   rl(ii)=drl*dble(ii-1)
 end do
!
 vpl(:,:)=0.0d0
 l1=lloc+1
 call dpnint(rr,vpuns(1,l1),mmax,rl,vpl(1,l1),nrl)

! override dpnint extrapolation to zero for vpl
 vpl(1,l1)=vpuns(1,l1)

 do l1=1,lmax+1
   if(l1 /= lloc+1) then
     do jj=1,nproj(l1)
       if(nkap==1) call dpnint(rr,vkb(1,jj,l1,1),mmax,rl,vkbl(1,jj,l1),nrl)
     end do
     do jj=1,npr_sr(l1)
       call dpnint(rr,vsr(1,jj,l1),mmax,rl,vsrl(1,jj,l1),nrl)
     end do
     do jj=1,npr_so(l1)
       call dpnint(rr,vso(1,jj,l1),mmax,rl,vsol(1,jj,l1),nrl)
     end do
   end if
 end do


 call dpnint(rr,rho,mmax,rl,rhol,nrl)
 call dpnint(rr,rhotae,mmax,rl,rhotael,nrl)
 call dpnint(rr,rhoc,mmax,rl,rhocl,nrl)

 do jj=1,5
   call dpnint(rr,rhomod(1,jj),mmax,rl,rhomodl(1,jj),nrl)
 end do

! Output for Abinit input using pspcod=8

 if(iexc==1) then
   ixc_abinit=4
 else if(iexc==2) then
   ixc_abinit=5
 else if(iexc==3) then
   ixc_abinit=2
 else if(iexc==4) then
   ixc_abinit=11
 else if(iexc<0) then
   ixc_abinit=iexc
 end if

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

!fchrg flags the model core for Abinit: fcfact, or 1 when fcfact<=0
!(icmod=6 with the cutoff taken from the pseudo-waves)
 fchrg=0.0d0
 if(icmod>=1) then
   fchrg=fcfact
   if(fchrg<=0.0d0) fchrg=1.0d0
 end if

 write(6,'(/a)') 'Begin PSPCODE8'
 write(6,'(3a,4f10.5)') atsym,'    ONCVPSP-4.0.1' &
&  ,'  r_core=',(rc(l1),l1=1,lmax+1)
 write(6,'(2f12.4, 5a)') zz,zion, '      ', pspd,  &
&  '    zatom,zion,pspd'
 write(6,'(i6,i8,i4,3i6, a)') 8, ixc_abinit,lmax,lloc, &
&  nrl, 0, '    pspcod,pspxc,lmax,lloc,mmax,r2well'
 write(6,'(3f12.8, a)') rl(nrl),fchrg, 0.0,  &
&  '    rchrg fchrg qchrg'
 if(nkap==1) then
   write(6,'(5i6, a)') (nproj(l1),l1=1,5),  &
&    '    nproj'
   write(6,'(2i6, a)') 1,1, &
&    '           extension_switch'
 else
   write(6,'(4i6, a)') npr_sr(1:4),'    nproj'

! this is where abinit is informed that we have spin-orbit
   write(6,'(2i6, a)') 3,1, &
&    '           extension_switch'
   if(soscale==1.d0) then
     write(6,'(3i6, a)') npr_so(2:4),'    nprojso'
   else
     write(6,'(3i6, a,f6.2)') npr_so(2:4),'    nprojso  X,soscale', soscale
   end if
 end if

! write the VKB projectors, or the scalar-relativistic ones, and the local
! potential if it is one of the semi-local potentials
 zero(:)=0.0d0

 do l1=1,lmax+1
   ll=l1-1
   if(ll==lloc) then
     write(6,'(i4)') ll
     do ii = 1,nrl
       write(6,'(i6,1p,2e21.13)') ii,rl(ii),vpl(ii,l1)
     end do
   else if(nkap==1 .and. nproj(l1)>0) then
     write(6,'(i4,23x,1p,6e21.13)') ll,(evkb(jj,l1,1),jj=1,nproj(l1))
     do ii = 1,nrl
       if(rl(ii)<=rc(l1)+2.0d0*drl) then
         write(6,'(i6,1p,6e21.13)') ii,rl(ii),(vkbl(ii,jj,l1),jj=1,nproj(l1))
       else
         write(6,'(i6,f6.2,6f4.0)') ii,rl(ii),(vkbl(ii,jj,l1),jj=1,nproj(l1))
       end if
     end do
   else if(nkap==2 .and. npr_sr(l1)>0) then
     write(6,'(i4,23x,1p,10e21.13)') ll,(esr(jj,l1),jj=1,npr_sr(l1))
     do ii = 1,nrl
       if(rl(ii)<=rc(l1)+2.0d0*drl) then
         write(6,'(i6,1p,11e21.13)') ii,rl(ii),(vsrl(ii,jj,l1),jj=1,npr_sr(l1))
       else
         write(6,'(i6,f6.2,10f4.0)') ii,rl(ii),(zero(jj),jj=1,npr_sr(l1))
       end if
     end do
   end if
 end do

! write general local potential if called for
 if(lloc>lmax) then
   write(6,'(i4)') lloc
   do ii = 1,nrl
     write(6,'(i6,1p,2e21.13)') ii,rl(ii),vpl(ii,lloc+1)
   end do
 end if

! write spin-orbit projectors (none when nkap=1)
 do l1=2,lmax+1
   ll=l1-1
   if(npr_so(l1)>0) then
     write(6,'(i4,23x,1p,20e21.13)') ll,(soscale*eso(jj,l1),jj=1,npr_so(l1))
     do ii = 1,nrl
       if(rl(ii)<=rc(l1)+2.0d0*drl) then
         write(6,'(i6,1p,11e21.13)') ii,rl(ii),(vsol(ii,jj,l1),jj=1,npr_so(l1))
       else
         write(6,'(i6,f6.2,10f4.0)') ii,rl(ii),(zero(jj),jj=1,npr_so(l1))
       end if
     end do
   end if
 end do

! write the model core charge if called for
 if(fchrg>0.0d0) then
   do ii=1,nrl
     write(6,'(i6,1p,6e21.13)') ii,rl(ii),  &
&      (rhomodl(ii,jj),jj=1,5)
   end do
 end if

! write pseudo valence, all-electron valence and all-electron core charges
 do ii=1,nrl
   write(6,'(i6,1p,4e21.13)') ii,rl(ii),rhol(ii),rhotael(ii),rhocl(ii)
 end do

 if(nkap==1) then
   version='#scalar-relativistic version 4.0.1 0r/27/2018'
 else
   version='#relativistic version 4.0.1 03/03/2109'
 end if
 write(6,'(a)') '<INPUT>'
 write(6,'(a/a/a/a)') &
&    '#', &
&    '#ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&    trim(version), &
&    '#'

 write(6,'(a/a/a/a)') &
&    '#While it is not required under the terms of the GNU GPL, it is',&
&    '#suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)', &
&    '#in any publication utilizing these pseudopotentials.', &
&    '#'

 if(nkap==2) write(6,'(a)') '#Echo of input data for oncvpsp-4.0.1'
 call write_input_echo(nkap,1,atsym,zz,nc,nv,iexc,psfile, &
&                       na,la,fa,lmax,rc0,ep(:,1),ncon,nbas,qcut, &
&                       lloc,lpopt,dvloc0,nproj,debl, &
&                       icmod,fcfact0,rcfact,slack,fcgrid,rcgrid, &
&                       epsh1,epsh2,depsh,rxpsh,rlmax,drl, &
&                       ncnf,nvcnf,nacnf,lacnf,facnf,log_mesh,upf_mesh)

! write termination signal
 write(6,'(a)') '</INPUT>'
 write(6,'(/a)') 'END_PSP'

 deallocate(rhol,rl,vpl,rhomodl,rhotael,rhocl)
 deallocate(vkbl,vsrl,vsol)

 return
 end subroutine linout

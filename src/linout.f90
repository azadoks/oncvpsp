!
! Copyright (c) 1989-2017 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
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
! for Abinit input using pspcod=8

 subroutine linout(lmax,lloc,rc,vkb,evkb,nproj,rr,vpuns,rho,rhomod, &
&                  rhotae,rhoc,zz,zion,mmax,mxprj,iexc,icmod,nrl,drl, &
&                  atsym, &
&                  na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&                  fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,fcfact0,rcfact, &
&                  fcgrid,rcgrid, &
&                  epsh1,epsh2,depsh,rxpsh,rlmax,psfile,log_mesh)


!lmax  maximum angular momentum
!lloc  l for local potential
!rc  core radii
!vkb  VKB projectors
!evkb  coefficients of VKB projectors
!nproj  number of vkb projectors for each l
!rr  log radial grid
!vpuns  unscreened semi-local pseudopotentials (vpuns(:,5) is local
!       potential if linear combination is used)
!rho  valence pseudocharge
!rhomod  model core charge
!rhotae  all-electron valence charge
!rhoc  all-electron core charge
!zz  atomic number
!zion  at this point, total valence charge (becomes psuedoion charge)
!mmax  size of log radial grid
!mxprj  dimension of number of projectors
!iexc  type of exchange-correlation
!icmod  1 if model core charge is used, otherwise 0
!fcfact  model core amplitude prefactor (1 for icmod=4, cutoff for icmod=6)
!fcfact0  input fcfact, for the echo
!rcfact  prefactor for the model core charge scale
!nrl  size of linear radial grid
!drl  spacing of linear radial grid
!atsym  atomic symbol
!(the arguments from na on are input variables echoed in the psp8 file)
!na  principal quantum numbers, reference configuration
!la  angular momenta, reference configuration
!ncon  number of matching constraints for each l
!nbas  number of basis functions for each l
!nvcnf  number of valence states, all configurations
!nacnf  principal quantum numbers, all configurations
!lacnf  angular momenta, all configurations
!nc  number of core states
!nv  number of valence states
!lpopt  choice of polynomial for lloc==4
!ncnf  number of test configurations
!fa  occupancies, reference configuration
!rc0  input core radii
!ep  energies at which the pseudopotentials are generated (Ha)
!qcut  residual-energy cutoff wave vectors for each l
!debl  energy shifts for the second projectors for each l (Ha)
!facnf  occupancies, all configurations
!dvloc0  amplitude at rr==0 of the function added to the local potential
!        for lloc==4
!fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!epsh1  lower limit of the log derivative analysis energy range (Ha)
!epsh2  upper limit of the log derivative analysis energy range (Ha)
!depsh  energy step of the log derivative analysis (Ha)
!rxpsh  radius (bohr) of the log derivative analysis, negative for the
!       default
!rlmax  maximum radius of the linear output mesh (bohr)
!psfile  output format, should be 'psp8' or 'both'
!log_mesh  amesh, rr1, rmax of the log mesh, echoed when the input sets
!          it (else 0)

 use constants_m, only: dp
 use input_echo_m, only: write_input_echo
 implicit none

!Input variables
 integer, intent(in) :: lmax,lloc,iexc,mmax,mxprj,nrl,icmod
 integer, intent(in) :: nproj(6)
 real(dp), intent(in) :: drl,fcfact,fcfact0,rcfact,zz,zion
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)
 real(dp), intent(in) :: rr(mmax),vpuns(mmax,5),rho(mmax),vkb(mmax,mxprj,4)
 real(dp), intent(in) :: rhotae(mmax),rhoc(mmax)
 real(dp), intent(in) :: rhomod(mmax,5)
 real(dp), intent(in) :: rc(6),evkb(mxprj,4)
 character(len=2), intent(in) :: atsym

!additional input for psp8 output to echo input file, all as defined
! in the main progam
 integer, intent(in) :: na(30),la(30),ncon(6),nbas(6)
 integer, intent(in) :: nvcnf(5),nacnf(30,5),lacnf(30,5)
 integer, intent(in) :: nc,nv,lpopt,ncnf
 real(dp), intent(in) :: fa(30),rc0(6),ep(6),qcut(6),debl(6),facnf(30,5)
 real(dp), intent(in) :: dvloc0,epsh1,epsh2,depsh,rxpsh,rlmax
 real(dp), intent(in) :: log_mesh(3)
 character(len=4), intent(in) :: psfile

!Local variables
 integer :: ii,iprj,jj,ll,l1,ixc_abinit
 integer :: dtime(8)
 real(dp) :: fchrg
 real(dp), allocatable :: rhomodl(:,:)
 real(dp),allocatable :: rhol(:),rl(:),vkbl(:,:,:),vpl(:,:)
 real(dp),allocatable :: rhotael(:),rhocl(:)
 character(len=2) :: pspd(3)


 allocate(rhol(nrl),rl(nrl),vkbl(nrl,mxprj,4),vpl(nrl,5),rhomodl(nrl,5))
 allocate(rhotael(nrl),rhocl(nrl))

!
! interpolation of everything onto linear output mesh

 do  ii=1,nrl
   rl(ii)=drl*dble(ii-1)
 end do
!
 do l1=1,max(lmax+1,lloc+1)
   call dpnint(rr,vpuns(1,l1),mmax,rl,vpl(1,l1),nrl)

! override dpnint extrapolation to zero for vpl
   vpl(1,l1)=vpuns(1,l1)

   if(l1 /= lloc + 1) then
     do iprj=1,nproj(l1)

       call dpnint(rr,vkb(1,iprj,l1),mmax,rl,vkbl(1,iprj,l1),nrl)

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
 write(6,'(5i6, a)') (nproj(l1),l1=1,5),  &
&  '    nproj'
 write(6,'(2i6, a)') 1,1, &
&  '           extension_switch'

! write the VKB projectors and the local potential
 do l1=1,max(lmax+1,lloc+1)
   ll=l1-1
   if(ll==lloc) then
     write(6,'(i4)') ll
     do ii = 1,nrl
       write(6,'(i6,1p,2e21.13)') ii,rl(ii),vpl(ii,l1)
     end do
   else if(nproj(l1)>0) then
     write(6,'(i4,23x,1p,6e21.13)') ll,(evkb(jj,l1),jj=1,nproj(l1))
     do ii = 1,nrl
       if(rl(ii)<=rc(l1)+2.0d0*drl) then
         write(6,'(i6,1p,6e21.13)') ii,rl(ii),(vkbl(ii,jj,l1),jj=1,nproj(l1))
       else
         write(6,'(i6,f6.2,6f4.0)') ii,rl(ii),(vkbl(ii,jj,l1),jj=1,nproj(l1))
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

 write(6,'(a)') '<INPUT>'
 write(6,'(a/a/a/a)') &
&    '#', &
&    '#ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)', &
&    '#scalar-relativistic version 4.0.1 0r/27/2018', &
&    '#'

 write(6,'(a/a/a/a)') &
&    '#While it is not required under the terms of the GNU GPL, it is',&
&    '#suggested that you cite D. R. Hamann, Phys. Rev. B 88, 085117 (2013)', &
&    '#in any publication utilizing these pseudopotentials.', &
&    '#'

 call write_input_echo(1,1,atsym,zz,nc,nv,iexc,psfile, &
&                       na,la,fa,lmax,rc0,ep,ncon,nbas,qcut, &
&                       lloc,lpopt,dvloc0,nproj,debl, &
&                       icmod,fcfact0,rcfact,fcgrid,rcgrid, &
&                       epsh1,epsh2,depsh,rxpsh,rlmax,drl, &
&                       ncnf,nvcnf,nacnf,lacnf,facnf,log_mesh)

! write termination signal
 write(6,'(a)')'</INPUT>'
 write(6,'(/a)') 'END_PSP'

 deallocate(rhol,rl,vkbl,vpl,rhomodl,rhotael,rhocl)

 return
 end subroutine linout

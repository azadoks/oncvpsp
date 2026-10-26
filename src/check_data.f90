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
 subroutine check_data(atsym,zz,fcfact,rcfact,slack,fcgrid,rcgrid,epsh1,epsh2, &
&                      depsh,rxpsh,rlmax,drl,upf_mesh,fa,facnf, &
&                      rc,ep,qcut,debl,nc,nv,iexc,lmax,lloc,lpopt,icmod, &
&                      ncnf,na,la,nvcnf,nacnf,lacnf,ncon,nbas,nproj,psfile, &
&                      amesh,rr1,rmax,c_light_in)

!! check input data for various violations of size or sign
!! see annotated input document for all definitions
!! *cnf variables are arrays of the basis atomic configuration variables
!! to be used for the test configurations
!!
!! @param atsym  atomic symbol
!! @param zz  atomic number
!! @param fcfact  prefactor for model core charge amplitude
!! @param rcfact  prefactor for model core charge scale (icmod=3), or error
!!                tolerance (icmod=6)
!! @param slack  share of the cutoff the model core must fit (icmod=6, 7)
!! @param fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!! @param rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!! @param epsh1  lower limit of the log derivative analysis energy range (Ha)
!! @param epsh2  upper limit of the log derivative analysis energy range (Ha)
!! @param depsh  energy step of the log derivative analysis (Ha)
!! @param rxpsh  log derivative radius (bohr), negative for the default
!! @param rlmax  maximum radius of the linear output mesh (bohr)
!! @param drl  spacing of the linear output mesh (bohr)
!! @param upf_mesh  most points of the log mesh of the UPF file, 0 for the
!!                  linear output mesh
!! @param fa  reference configuration occupancies
!! @param facnf  occupancies of the reference (1) and test configurations
!! @param rc  core radii for each l, rc(5) for the lloc=4 local potential (bohr)
!! @param ep  energies at which the pseudopotentials are generated (Ha;
!!            not checked)
!! @param qcut  residual-energy cutoff wave vectors for each l
!! @param debl  energy shifts for the second projectors for each l (Ha)
!! @param nc  number of core states
!! @param nv  number of valence states
!! @param iexc  exchange-correlation functional, negative for libxc
!! @param lmax  maximum angular momentum of the pseudopotential
!! @param lloc  angular momentum of the local potential, 4 for a polynomial
!! @param lpopt  type of polynomial continuation for lloc=4
!! @param icmod  model core charge option
!! @param ncnf  number of test configurations
!! @param na  reference configuration principal quantum numbers
!! @param la  reference configuration angular momenta
!! @param nvcnf  numbers of valence states of the reference and test
!!               configurations
!! @param nacnf  principal quantum numbers of the reference and test
!!               configurations
!! @param lacnf  angular momenta of the reference and test configurations
!! @param ncon  number of matching constraints for each l
!! @param nbas  number of basis functions for each l
!! @param nproj  number of projectors for each l
!! @param psfile  format of the pseudopotential file (psp8, upf, psml, both,
!!                none)
!! @param amesh  log mesh ratio, 0 or less for the default
!! @param rr1  first log mesh point (bohr), 0 or less for the default
!! @param rmax  outer radius of the log mesh (bohr), 0 or less for the default
!! @param c_light_in  speed of light, 0 or less for the default

 use constants_m, only: dp, mxstate, mxcnf, mxl, mxlin
 implicit none

!Input variables
 character(len=2), intent(in) :: atsym
 character(len=4), intent(in) :: psfile

 real(dp), intent(in) :: zz,fcfact,rcfact,epsh1,epsh2,depsh,rxpsh,rlmax,drl
 real(dp), intent(in) :: slack
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)
 real(dp), intent(in) :: amesh,rr1,rmax,c_light_in
 real(dp), intent(in) :: fa(mxstate),facnf(mxstate,mxcnf),rc(mxlin+1)
 real(dp), intent(in) :: ep(mxlin+1),qcut(mxlin+1)
 real(dp), intent(in) :: debl(mxlin+1)

 integer, intent(in) :: nc,nv,iexc,lmax,lloc,lpopt,icmod,ncnf,upf_mesh
 integer, intent(in) :: na(mxstate),la(mxstate),nvcnf(mxcnf)
 integer, intent(in) :: nacnf(mxstate,mxcnf),lacnf(mxstate,mxcnf)
 integer, intent(in) :: ncon(mxlin+1),nbas(mxlin+1),nproj(mxlin+1)

!Local variables
!the number of errors found
 integer :: ierr

 ierr=0

 call check_atom()

 call check_channels()

 call check_model_core()

 call check_meshes()

 call check_tests()

 if(ierr>0) then
   write(6,'(a,i4,a)') 'ERROR: test_data found',ierr,' ERROR; stopping'
   error stop 1
 end if

 return

 contains

 subroutine check_atom()
!! the atom, the file format, the functional and the reference
!! configuration
 integer :: ii
 real(dp) :: sf
 if(trim(atsym)=='') then
   write(6,'(a)') 'test_data: must have non-blank atsym'
   ierr=ierr+1
 end if

 if(trim(psfile)/='psp8' .and. trim(psfile)/='upf' .and. trim(psfile)/='psml' &
&   .and. trim(psfile)/='both' .and. trim(psfile)/='none') then
   write(6,'(a)') 'test_data: psfile must == psp8 or upf or psml'// &
&        ' or both or none'
   ierr=ierr+1
 end if

 if(zz<=0.0d0) then
   write(6,'(a)') 'test_data: must have positive atomic number z'
   ierr=ierr+1
 end if

 if(nc+nv<=0 .or. nc+nv>mxstate) then
   write(6,'(a,i0)') 'test_data: must have 1 <= nc+nv <=',mxstate
   ierr=ierr+1
 end if

!iexc=5 is r2SCAN in METAPSP, so it is kept free here
 if(iexc==5) then
   write(6,'(a)') 'test_data: iexc=5 is r2SCAN, METAPSP only'
   ierr=ierr+1
 else if(iexc==0 .or. iexc>7) then
   write(6,'(a)') 'test_data: must have iexc <=7, iexc /= 0, 5.'
   write(6,'(a)') '  negative values for libxc'
   ierr=ierr+1
 end if

 sf=0.0d0
 do ii=1,nc+nv
   sf=sf+fa(ii)
   if(la(ii)<0 .or. la(ii)>3) then
     write(6,'(a,a,i4)') 'test_data: must have 0<= l <=3,', &
&     ' reference configuration line',ii
     ierr=ierr+1
   end if
   if(na(ii)<=la(ii)) then
     write(6,'(a,i4)') 'test_data: l > n, reference configuration line',ii
     ierr=ierr+1
   end if
   if(fa(ii)<=0.0d0) then
     write(6,'(a,i4)') 'test_data: f <= 0.0, reference configuration line',ii
     ierr=ierr+1
   end if
 end do

 if(sf>zz) then
   write(6,'(a)') 'test_data: reference configuration is negative ion'
   ierr=ierr+1
 end if

 end subroutine check_atom

 subroutine check_channels()
!! lmax, and the core radius, constraints, basis, cutoff and projectors
!! of each l and of the local potential
 integer :: l1
 if(lmax<0 .or. lmax>mxl) then
   write(6,'(a,i0)') 'test_data: 0 <= lmax <= ',mxl
   ierr=ierr+1
 end if

 do l1=1,lmax+1
   if(rc(l1)<=0.0d0) then
     write(6,'(a,i4)') 'test_data: must have rc>0.0, l=',l1-1
     ierr=ierr+1
   end if
   if(ncon(l1)<3 .or. ncon(l1)>5) then
     write(6,'(a,i4)') 'test_data: must have 3 <= ncon <= 5, l=',l1-1
     ierr=ierr+1
   end if
   if(nproj(l1)==2) then
     if(nbas(l1)<ncon(l1)+3 .or. nbas(l1)>ncon(l1)+5) then
       write(6,'(a,i4/a)') 'test_data: must have ncon+3 <= nbas <= ncon+5, l=',&
&       l1-1,' for nproj=2'
       ierr=ierr+1
     end if
   else
     if(nbas(l1)<ncon(l1)+2 .or. nbas(l1)>ncon(l1)+5) then
       write(6,'(a,i4)') 'test_data: must have ncon+2 <= nbas <= ncon+5, l=', &
&        l1-1
       ierr=ierr+1
     end if
   end if
   if(qcut(l1)<=0.0d0) then
     write(6,'(a,i4)') 'test_data: must have qcut>0.0, l=',l1-1
     ierr=ierr+1
   end if
 end do

 if(lloc==4) then
   if(lpopt<1 .or. lpopt>5) then
     write(6,'(a)') 'test_data: must have 1 <= lpopt <= 5'
     ierr=ierr+1
   end if
   if(rc(5)<=0.0d0) then
     write(6,'(a)') 'test_data: for lloc==4, must have rc > 0.0'
     ierr=ierr+1
   end if
 else if(lloc<0 .or. lloc>lmax) then
   write(6,'(a)') 'test_data: must have 0 <= lloc <= lmax or lloc == 4'
   ierr=ierr+1
 end if

 if(lloc>=0 .and. lloc<=4) then
   do l1=1,lmax+1
     if(rc(l1)<rc(lloc+1)) then
       write(6,'(a,i2,a)') 'test_data: rc < rc(lloc) for l =',l1-1,&
&            '.  Not allowed.'
       ierr=ierr+1
     end if
   end do
 end if

 do l1=1,lmax+1
   if(nproj(l1)==0 .and. (lloc==4 .or. lloc==l1-1)) then
     cycle
   else if(nproj(l1)<1 .or. nproj(l1)>5) then
     write(6,'(a,i4,a)') 'test_data: must have nproj in [1,5]'// &
&          ' (0 OK lloc=4 or', &
&          l1-1,')'
     ierr=ierr+1
   end if
   if(debl(l1)<0.0d0) then
     write(6,'(a,i4)') 'test_data: must have debl>0.0)',l1-1
     ierr=ierr+1
   end if
 end do

 end subroutine check_channels

 subroutine check_model_core()
!! icmod and the parameters of its model core
 if((icmod<0 .or. icmod>4) .and. icmod/=6 .and. icmod/=7) then
   write(6,'(a)') 'test_data: must have 0<= icmod <=4 or icmod=6 or 7'
   ierr=ierr+1
 end if

 if((icmod==6 .or. icmod==7) .and. rcfact<=0.0d0) then
   write(6,'(a)') 'test_data: must have rcfact>0.0 for icmod= 6 or 7'
   ierr=ierr+1
 end if

 if((icmod==6 .or. icmod==7) .and. (slack<=0.0d0 .or. slack>1.0d0)) then
   write(6,'(a)') 'test_data: must have 0.0<slack<=1.0 for icmod= 6 or 7'
   ierr=ierr+1
 end if

 if((icmod==1 .or. icmod==3) .and. fcfact<=0.0d0) then
   write(6,'(a)') 'test_data: must have fcfact>0.0 for icmod= 1 or 3'
   ierr=ierr+1
 end if

 if((icmod==3) .and. rcfact<=0.0d0) then
   write(6,'(a)') 'test_data: must have rcfact>0.0 for icmod= 3'
   ierr=ierr+1
 end if

 if(icmod==4) then
   if(fcgrid(1)<=0.0d0 .or. fcgrid(2)<fcgrid(1) .or. fcgrid(3)<=0.0d0) then
     write(6,'(/a,3(a,f12.6))') 'check_data: ERROR icmod=4 needs', &
&          ' 0<fcfact_min<=fcfact_max and fcfact_step>0, fcfact_min=', &
&          fcgrid(1), &
&          ', fcfact_max=',fcgrid(2),', fcfact_step=',fcgrid(3)
     ierr=ierr+1
   end if
   if(rcgrid(1)<=0.0d0 .or. rcgrid(2)<rcgrid(1) .or. rcgrid(3)<=0.0d0) then
     write(6,'(/a,3(a,f12.6))') 'check_data: ERROR icmod=4 needs', &
&          ' 0<rcfact_min<=rcfact_max and rcfact_step>0, rcfact_min=', &
&          rcgrid(1), &
&          ', rcfact_max=',rcgrid(2),', rcfact_step=',rcgrid(3)
     ierr=ierr+1
   end if
 end if

 end subroutine check_model_core

 subroutine check_meshes()
!! the energies and radius of the log derivatives, and the output and log
!! meshes, against the largest core radius
 integer :: l1
 real(dp) :: rcmax
 if(epsh2<epsh1) then
   write(6,'(a)') 'test_data: must have epsh2 > epsh1'
   ierr=ierr+1
 end if
 if(depsh<=0.0d0) then
   write(6,'(a)') 'test_data: must have depsh>0.0'
   ierr=ierr+1
 end if

 rcmax=0.0d0
 do l1=1,lmax+1
   rcmax=max(rcmax,rc(l1))
 end do
 if(lloc==4) then
   rcmax=max(rcmax,rc(5))
 end if

 if(rlmax<rcmax) then
   write(6,'(a)') 'test_data: must have rlmax>rcmax'
   ierr=ierr+1
 end if

!a set rxpsh must lie beyond every core radius, so that the log derivatives
!are evaluated outside the pseudized region (run_phsft)
 if(rxpsh>=0.0d0 .and. rxpsh<=rcmax) then
   write(6,'(/2a,f12.6,a,f12.6)') 'check_data: ERROR rxpsh not beyond'// &
&        ' every core radius', &
&        ', rxpsh=',rxpsh,', rcmax=',rcmax
   ierr=ierr+1
 end if

 if(drl<=0.0d0) then
   write(6,'(a)') 'test_data: must have drl > 0.0'
   ierr=ierr+1
 end if

 if(upf_mesh<0) then
   write(6,'(a)') 'test_data: must have upf_mesh >= 0'
   ierr=ierr+1
 end if

 if(amesh>0.0d0 .and. amesh<=1.0d0) then
   write(6,'(a)') 'test_data: must have amesh > 1.0'
   ierr=ierr+1
 end if

 if(rr1>0.0d0 .and. rr1>=rcmax) then
   write(6,'(a)') 'test_data: must have rr1 < rcmax'
   ierr=ierr+1
 end if

 if(rmax>0.0d0 .and. rmax<=rlmax) then
   write(6,'(a)') 'test_data: must have rmax > rlmax'
   ierr=ierr+1
 end if

 if(c_light_in>0.0d0 .and. c_light_in<1.0d0) then
   write(6,'(a)') 'test_data: must have c_light >= 1.0, or 0 or less'// &
&        ' for the default'
   ierr=ierr+1
 end if

 end subroutine check_meshes

 subroutine check_tests()
!! the test configurations
 integer :: ii,jj
 real(dp) :: sf
 if(ncnf<0 .or. ncnf>mxcnf-1)then
   write(6,'(a,i0,a)') 'test_data: must have 0 <= ncnf <= ',mxcnf-1, &
&        ', skipping config. checks'
   ierr=ierr+1
 else
   do jj=2,ncnf+1
     sf=0.0d0
     do ii=1,nc+nvcnf(jj)
       sf=sf+facnf(ii,jj)
       if(lacnf(ii,jj)<0 .or. lacnf(ii,jj)>3) then
         write(6,'(a,a,i4,a,i4)') 'test_data: must have 0<= l <=3,', &
&         ' test configuration',jj-1, ' line',ii-nc
         ierr=ierr+1
       end if
       if(nacnf(ii,jj)<=lacnf(ii,jj)) then
         write(6,'(a,a,i4,a,i4)') 'test_data: l > n,', &
&         ' test configuration',jj-1, ' line',ii-nc
         ierr=ierr+1
       end if
       if(facnf(ii,jj)<0.0d0) then
         write(6,'(a,a,i4,a,i4)') 'test_data: f < 0.0,', &
&         ' test configuration',jj-1, ' line',ii-nc
         ierr=ierr+1
       end if
     end do

     if(sf>zz) then
       write(6,'(a,i4)') 'test_data: negative ion, test configuration',jj-1
       ierr=ierr+1
     end if
   end do
 end if

 end subroutine check_tests

 end subroutine check_data

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

 module check_data_m

 implicit none
 private
 public :: check_data

 contains

 subroutine check_data(inp,fcgrid,rcgrid,iout,ifail)

!! check input data for various violations of size or sign
!! see annotated input document for all definitions
!! *cnf variables are arrays of the basis atomic configuration variables
!! to be used for the test configurations
!!
!! @param inp  the input values (input_m), the reference configuration
!!             copied into the first test configuration
!! @param fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!! @param rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!! @param iout  unit of the text output
!! @param ifail  0, or 1 after a failure, whose message is written

 use constants_m, only: dp, mxstate, mxcnf, mxl, mxlin
 use input_m, only: input_t
 use exc_libxc_m, only: have_libxc
 implicit none

!Input variables
 type(input_t), intent(in) :: inp
 real(dp), intent(in) :: fcgrid(3),rcgrid(3)

!Output variables
 integer, intent(in) :: iout
 integer, intent(out) :: ifail

!Local variables
!the number of errors found
 integer :: ierr

 ifail=0
 ierr=0

 call check_atom()

 call check_channels()

 call check_model_core()

 call check_meshes()

 call check_tests()

 if(ierr>0) then
   write(iout,'(a,i4,a)') 'ERROR: test_data found',ierr,' ERROR; stopping'
   ifail=1
   return
 end if

 return

 contains

 subroutine check_atom()
!! the atom, the file format, the functional and the reference
!! configuration
 integer :: ii
 real(dp) :: sf
 if(trim(inp%atsym)=='') then
   write(iout,'(a)') 'test_data: must have non-blank atsym'
   ierr=ierr+1
 end if

 if(trim(inp%psfile)/='psp8' .and. trim(inp%psfile)/='upf' .and. &
&   trim(inp%psfile)/='psml' .and. trim(inp%psfile)/='both' .and. &
&   trim(inp%psfile)/='none') then
   write(iout,'(a)') 'test_data: psfile must == psp8 or upf or psml'// &
&        ' or both or none'
   ierr=ierr+1
 end if

 if(inp%zz<=0.0d0) then
   write(iout,'(a)') 'test_data: must have positive atomic number z'
   ierr=ierr+1
 end if

 if(inp%nc+inp%nv<=0 .or. inp%nc+inp%nv>mxstate) then
   write(iout,'(a,i0)') 'test_data: must have 1 <= nc+nv <=',mxstate
   ierr=ierr+1
 end if

!iexc=5 is r2SCAN in METAPSP, so it is kept free here
 if(inp%iexc==5) then
   write(iout,'(a)') 'test_data: iexc=5 is r2SCAN, METAPSP only'
   ierr=ierr+1
 else if(inp%iexc==0 .or. inp%iexc>7) then
   write(iout,'(a)') 'test_data: must have iexc <=7, iexc /= 0, 5.'
   write(iout,'(a)') '  negative values for libxc'
   ierr=ierr+1
 else if(inp%iexc<0 .and. .not. have_libxc) then
   write(iout,'(a,i0,a)') 'test_data: iexc = ',inp%iexc, &
&        ' needs libxc, which this build lacks'
   ierr=ierr+1
 end if

 sf=0.0d0
 do ii=1,inp%nc+inp%nv
   sf=sf+inp%fa(ii)
   if(inp%la(ii)<0 .or. inp%la(ii)>3) then
     write(iout,'(a,a,i4)') 'test_data: must have 0<= l <=3,', &
&     ' reference configuration line',ii
     ierr=ierr+1
   end if
   if(inp%na(ii)<=inp%la(ii)) then
     write(iout,'(a,i4)') 'test_data: l > n, reference configuration line',ii
     ierr=ierr+1
   end if
   if(inp%fa(ii)<=0.0d0) then
     write(iout,'(a,i4)') 'test_data: f <= 0.0, reference configuration line',ii
     ierr=ierr+1
   end if
 end do

 if(sf>inp%zz) then
   write(iout,'(a)') 'test_data: reference configuration is negative ion'
   ierr=ierr+1
 end if

 end subroutine check_atom

 subroutine check_channels()
!! lmax, and the core radius, constraints, basis, cutoff and projectors
!! of each l and of the local potential
 integer :: l1
 if(inp%lmax<0 .or. inp%lmax>mxl) then
   write(iout,'(a,i0)') 'test_data: 0 <= lmax <= ',mxl
   ierr=ierr+1
 end if

 do l1=1,inp%lmax+1
   if(inp%rc(l1)<=0.0d0) then
     write(iout,'(a,i4)') 'test_data: must have rc>0.0, l=',l1-1
     ierr=ierr+1
   end if
   if(inp%ncon(l1)<3 .or. inp%ncon(l1)>5) then
     write(iout,'(a,i4)') 'test_data: must have 3 <= ncon <= 5, l=',l1-1
     ierr=ierr+1
   end if
   if(inp%nproj(l1)==2) then
     if(inp%nbas(l1)<inp%ncon(l1)+3 .or. inp%nbas(l1)>inp%ncon(l1)+5) then
       write(iout,'(a,i4/a)') &
&            'test_data: must have ncon+3 <= nbas <= ncon+5, l=',&
&       l1-1,' for nproj=2'
       ierr=ierr+1
     end if
   else
     if(inp%nbas(l1)<inp%ncon(l1)+2 .or. inp%nbas(l1)>inp%ncon(l1)+5) then
       write(iout,'(a,i4)') &
&            'test_data: must have ncon+2 <= nbas <= ncon+5, l=', &
&        l1-1
       ierr=ierr+1
     end if
   end if
   if(inp%qcut(l1)<=0.0d0) then
     write(iout,'(a,i4)') 'test_data: must have qcut>0.0, l=',l1-1
     ierr=ierr+1
   end if
 end do

 if(inp%lloc==4) then
   if(inp%lpopt<1 .or. inp%lpopt>5) then
     write(iout,'(a)') 'test_data: must have 1 <= lpopt <= 5'
     ierr=ierr+1
   end if
   if(inp%rc(5)<=0.0d0) then
     write(iout,'(a)') 'test_data: for lloc==4, must have rc > 0.0'
     ierr=ierr+1
   end if
 else if(inp%lloc<0 .or. inp%lloc>inp%lmax) then
   write(iout,'(a)') 'test_data: must have 0 <= lloc <= lmax or lloc == 4'
   ierr=ierr+1
 end if

 if(inp%lloc>=0 .and. inp%lloc<=4) then
   do l1=1,inp%lmax+1
     if(inp%rc(l1)<inp%rc(inp%lloc+1)) then
       write(iout,'(a,i2,a)') 'test_data: rc < rc(lloc) for l =',l1-1,&
&            '.  Not allowed.'
       ierr=ierr+1
     end if
   end do
 end if

 do l1=1,inp%lmax+1
   if(inp%nproj(l1)==0 .and. (inp%lloc==4 .or. inp%lloc==l1-1)) then
     cycle
   else if(inp%nproj(l1)<1 .or. inp%nproj(l1)>5) then
     write(iout,'(a,i4,a)') 'test_data: must have nproj in [1,5]'// &
&          ' (0 OK lloc=4 or', &
&          l1-1,')'
     ierr=ierr+1
   end if
   if(inp%debl(l1)<0.0d0) then
     write(iout,'(a,i4)') 'test_data: must have debl>0.0)',l1-1
     ierr=ierr+1
   end if
 end do

 end subroutine check_channels

 subroutine check_model_core()
!! icmod and the parameters of its model core
 if((inp%icmod<0 .or. inp%icmod>4) .and. inp%icmod/=6 .and. inp%icmod/=7) then
   write(iout,'(a)') 'test_data: must have 0<= icmod <=4 or icmod=6 or 7'
   ierr=ierr+1
 end if

 if((inp%icmod==6 .or. inp%icmod==7) .and. inp%rcfact<=0.0d0) then
   write(iout,'(a)') 'test_data: must have rcfact>0.0 for icmod= 6 or 7'
   ierr=ierr+1
 end if

 if((inp%icmod==6 .or. inp%icmod==7) .and. &
&   (inp%slack<=0.0d0 .or. inp%slack>1.0d0)) then
   write(iout,'(a)') 'test_data: must have 0.0<slack<=1.0 for icmod= 6 or 7'
   ierr=ierr+1
 end if

 if((inp%icmod==1 .or. inp%icmod==3) .and. inp%fcfact<=0.0d0) then
   write(iout,'(a)') 'test_data: must have fcfact>0.0 for icmod= 1 or 3'
   ierr=ierr+1
 end if

 if((inp%icmod==3) .and. inp%rcfact<=0.0d0) then
   write(iout,'(a)') 'test_data: must have rcfact>0.0 for icmod= 3'
   ierr=ierr+1
 end if

 if(inp%icmod==4) then
   if(fcgrid(1)<=0.0d0 .or. fcgrid(2)<fcgrid(1) .or. fcgrid(3)<=0.0d0) then
     write(iout,'(/a,3(a,f12.6))') 'check_data: ERROR icmod=4 needs', &
&          ' 0<fcfact_min<=fcfact_max and fcfact_step>0, fcfact_min=', &
&          fcgrid(1), &
&          ', fcfact_max=',fcgrid(2),', fcfact_step=',fcgrid(3)
     ierr=ierr+1
   end if
   if(rcgrid(1)<=0.0d0 .or. rcgrid(2)<rcgrid(1) .or. rcgrid(3)<=0.0d0) then
     write(iout,'(/a,3(a,f12.6))') 'check_data: ERROR icmod=4 needs', &
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
 if(inp%epsh2<inp%epsh1) then
   write(iout,'(a)') 'test_data: must have epsh2 > epsh1'
   ierr=ierr+1
 end if
 if(inp%depsh<=0.0d0) then
   write(iout,'(a)') 'test_data: must have depsh>0.0'
   ierr=ierr+1
 end if

 rcmax=0.0d0
 do l1=1,inp%lmax+1
   rcmax=max(rcmax,inp%rc(l1))
 end do
 if(inp%lloc==4) then
   rcmax=max(rcmax,inp%rc(5))
 end if

 if(inp%rlmax<rcmax) then
   write(iout,'(a)') 'test_data: must have rlmax>rcmax'
   ierr=ierr+1
 end if

!a set rxpsh must lie beyond every core radius, so that the log derivatives
!are evaluated outside the pseudized region (run_phsft)
 if(inp%rxpsh>=0.0d0 .and. inp%rxpsh<=rcmax) then
   write(iout,'(/2a,f12.6,a,f12.6)') 'check_data: ERROR rxpsh not beyond'// &
&        ' every core radius', &
&        ', rxpsh=',inp%rxpsh,', rcmax=',rcmax
   ierr=ierr+1
 end if

 if(inp%drl<=0.0d0) then
   write(iout,'(a)') 'test_data: must have drl > 0.0'
   ierr=ierr+1
 end if

 if(inp%upf_mesh<0) then
   write(iout,'(a)') 'test_data: must have upf_mesh >= 0'
   ierr=ierr+1
 end if

 if(inp%amesh>0.0d0 .and. inp%amesh<=1.0d0) then
   write(iout,'(a)') 'test_data: must have amesh > 1.0'
   ierr=ierr+1
 end if

 if(inp%rr1>0.0d0 .and. inp%rr1>=rcmax) then
   write(iout,'(a)') 'test_data: must have rr1 < rcmax'
   ierr=ierr+1
 end if

 if(inp%rmax>0.0d0 .and. inp%rmax<=inp%rlmax) then
   write(iout,'(a)') 'test_data: must have rmax > rlmax'
   ierr=ierr+1
 end if

 if(inp%c_light>0.0d0 .and. inp%c_light<1.0d0) then
   write(iout,'(a)') 'test_data: must have c_light >= 1.0, or 0 or less'// &
&        ' for the default'
   ierr=ierr+1
 end if

 end subroutine check_meshes

 subroutine check_tests()
!! the test configurations
 integer :: ii,jj
 real(dp) :: sf
 if(inp%ncnf<0 .or. inp%ncnf>mxcnf-1)then
   write(iout,'(a,i0,a)') 'test_data: must have 0 <= ncnf <= ',mxcnf-1, &
&        ', skipping config. checks'
   ierr=ierr+1
 else
   do jj=2,inp%ncnf+1
     sf=0.0d0
     do ii=1,inp%nc+inp%cnf(jj)%nv
       sf=sf+inp%cnf(jj)%f(ii)
       if(inp%cnf(jj)%l(ii)<0 .or. inp%cnf(jj)%l(ii)>3) then
         write(iout,'(a,a,i4,a,i4)') 'test_data: must have 0<= l <=3,', &
&         ' test configuration',jj-1, ' line',ii-inp%nc
         ierr=ierr+1
       end if
       if(inp%cnf(jj)%n(ii)<=inp%cnf(jj)%l(ii)) then
         write(iout,'(a,a,i4,a,i4)') 'test_data: l > n,', &
&         ' test configuration',jj-1, ' line',ii-inp%nc
         ierr=ierr+1
       end if
       if(inp%cnf(jj)%f(ii)<0.0d0) then
         write(iout,'(a,a,i4,a,i4)') 'test_data: f < 0.0,', &
&         ' test configuration',jj-1, ' line',ii-inp%nc
         ierr=ierr+1
       end if
     end do

     if(sf>inp%zz) then
       write(iout,'(a,i4)') 'test_data: negative ion, test configuration',jj-1
       ierr=ierr+1
     end if
   end do
 end if

 end subroutine check_tests

 end subroutine check_data

 end module check_data_m

!
! Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
! Copyright (c) 2015 by Alberto Garcia, ICMAB-CSIC for xmlf90-wxml calls
! in support of the PSML format
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
! Support for PSML file creation
!
 module psmlout_m
 use constants_m, only: dp
 implicit none

 public :: psmlout
 private

 character(len=1), dimension(0:4), parameter :: lsymb = ['s','p','d','f','g']

 contains

 subroutine psmlout(lmax,lloc,rc,vkb,evkb,mxprj,nproj,rr,vpuns,rho,rhomod, &
&                  irct, &
&                  vsr,esr,vso,eso, &
&                  zz,zion,mmax,iexc,icmod,nrl,drl,atsym,epstot, &
&                  na,la,ncon,nbas,nvcnf,nacnf,lacnf,nc,nv,lpopt,ncnf, &
&                  fa,rc0,ep,qcut,debl,facnf,dvloc0,fcfact,rcfact, &
&                  fcgrid,rcgrid, &
&                  epsh1,epsh2,depsh,rxpsh,rlmax,psfile,log_mesh, &
&                  irel,nkap)

!writes the pseudopotential to the output in PSML format, between
!'Begin PSML' and 'END_PSP' markers

!lmax  maximum angular momentum
!lloc  l for local potential
!rc  core radii
!vkb  VKB projectors, last index 1 for j=l+1/2, 2 for j=l-1/2 when nkap=2
!evkb  coefficients of VKB projectors, last index as for vkb
!mxprj  dimension of number of projectors
!nproj  number of vkb projectors for each l
!rr  log radial grid
!vpuns  unscreened semi-local pseudopotentials (vp(:,5,1) is local
!  potential if linear combination is used), last index as for vkb
!rho  valence pseudocharge
!rhomod  model core charge and its first four derivatives
!irct  index of the point at which the model core charge is matched
!vsr  scalar-relativistic projectors (irel=2 only)
!esr  energy coefficients of vsr
!vso  spin-orbit projectors
!eso  energy coefficients of vso
!zz  atomic number
!zion  at this point, total valence charge (becomes psuedoion charge)
!mmax  size of log radial grid
!iexc  type of exchange-correlation
!icmod  model core charge option, 0 if no model core charge is used
!nrl  size of linear radial grid (not used)
!drl  spacing of linear radial grid
!atsym  atomic symbol
!epstot  pseudoatom total energy (not used)
!na  principal quantum numbers of the reference configuration states
!la  angular momenta of the reference configuration states
!ncon  number of matching constraints for each l
!nbas  number of basis functions for each l
!nvcnf  numbers of valence states of the reference and test configurations
!nacnf  principal quantum numbers of the reference and test configurations
!lacnf  angular momenta of the reference and test configurations
!nc  number of core states
!nv  number of valence states
!lpopt  type of polynomial continuation for lloc=4
!ncnf  number of test configurations
!fa  occupations of the reference configuration states
!rc0  core radii as read from the input (bohr), rc0(5) for lloc=4
!ep  energies at which the pseudopotentials are generated (Ha)
!qcut  residual-energy cutoff wave vectors for each l
!debl  energy shifts for the second projectors for each l (Ha)
!facnf  occupancies of the reference and test configurations
!dvloc0  r=0 shift of the lloc=4 local potential (Ha)
!fcfact  prefactor for model core charge amplitude, as read from the input
!rcfact  prefactor for model core charge scale (icmod=3), or error
!  tolerance (icmod=6)
!fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!epsh1  lower limit of the log derivative analysis energy range (Ha)
!epsh2  upper limit of the log derivative analysis energy range (Ha)
!depsh  energy step of the log derivative analysis (Ha)
!rxpsh  radius of the log derivative analysis (bohr), negative for default
!rlmax  maximum radius of the linear output mesh (bohr)
!psfile  pseudopotential file format as read from the input; for irel=2,
!  'psp8' writes scalar-relativistic and spin-orbit sets, otherwise l,j
!  sets are written
!log_mesh  amesh, rr1, rmax of the log mesh, echoed when the input sets it
!  (else 0)
!irel  0 non-relativistic, 1 scalar-relativistic, 2 fully-relativistic
!nkap  number of kappa per l: 1, or 2 for irel=2

! Alberto Garcia, February 1, 2015

 use xmlf90_wxml, only: xmlf_t, str, xml_OpenUnit, xml_Close, &
&                       xml_AddXMLDeclaration, xml_NewElement, &
&                       xml_EndElement, xml_AddArray
 use m_libxc_list, only: libxc_t, XC_NOT_IMPL, XC_LDA_X, XC_LDA_C_PZ, &
&                        XC_LDA_C_WIGNER, XC_LDA_C_HL, XC_GGA_X_PBE, &
&                        XC_GGA_C_PBE

 use constants_m, only: pi
 implicit none


!Input variables
 integer, intent(in) :: lmax
 integer, intent(in) :: lloc,iexc,mmax,nrl,icmod,irel,nkap
 integer, intent(in) :: nproj(6)
 integer, intent(in) :: irct ! index of point at which rho_core is matched
 integer, intent(in) :: mxprj ! maximum number of projectors for any l
 real(dp), intent(in) :: drl,fcfact,zz,zion,epstot
 real(dp), target, intent(in) :: rr(mmax),vpuns(mmax,5,nkap),rho(mmax)
 real(dp), target, intent(in) :: vkb(mmax,mxprj,4,nkap)
 real(dp), target, intent(in) :: rhomod(mmax,5)
 real(dp), intent(in) :: rc(6),evkb(mxprj,4,nkap)
 real(dp), target, intent(in) :: vsr(mmax,2*mxprj,4),vso(mmax,2*mxprj,4)
 real(dp), intent(in) :: esr(2*mxprj,4),eso(2*mxprj,4)
 character(len=2), intent(in) :: atsym

!additional input to echo the input file, all as defined
! in the main progam
 integer, intent(in) :: na(30),la(30),ncon(6),nbas(6)
 integer, intent(in) :: nvcnf(5),nacnf(30,5),lacnf(30,5)
 integer, intent(in) :: nc,nv,lpopt,ncnf
 real(dp), intent(in) :: fa(30),rc0(6),ep(6),qcut(6),debl(6),facnf(30,5)
 real(dp), intent(in) :: dvloc0,epsh1,epsh2,depsh,rxpsh,rlmax
 real(dp), intent(in) :: rcfact,fcgrid(3),rcgrid(3),log_mesh(3)
 character(len=4), intent(in) :: psfile

!Local variables
 integer :: dtime(8)
 integer :: npr_so(4), npr_sr(4)

!-------
!psml stuff

 type(xmlf_t) :: xf
 type(libxc_t) :: libxc_id(2)

 integer   :: npts
 integer   :: ii, jj, l1, jk
 real(dp)  :: rcore, jval
 real(dp)  :: total_valence_charge
 real(dp), pointer :: rgrid(:), chval(:), chcore(:)
 real(dp), pointer :: vlocal(:)
 real(dp), allocatable :: r0(:), f0(:), vps(:)

 character(len=100)    :: line
 character(len=4)      :: polattrib, coreattrib
 character(len=1)      :: pscode, char_dummy
 character(len=2) :: nameat
 character(len=40):: psflavor


 character(len=30) :: xcfuntype, xcfunparam
 integer          :: ncore, nval, ncp, norbs, npots

!nsh, lsh, fsh  n, l and occupation of the reference configuration states
!nn, ll, ff  n, l and occupation that label each semi-local potential
 integer, allocatable  :: nsh(:), lsh(:)
 integer, allocatable  :: nn(:), ll(:)
 real(dp), allocatable :: fsh(:), ff(:)
 real(dp), allocatable :: fdown(:), fup(:)
 integer :: nr
 character(len=10)     :: datestr

 logical :: tdopsp, nonrel, polarized, there_is_core, found

!---

 call date_and_time(VALUES=dtime)
 write(datestr,"(i4,'-',i2.2,'-',i2.2)") dtime(1:3)


! the document goes to the output, between markers, as do the psp8 and
! upf files
 write(6,'(/a)') 'Begin PSML'
 call xml_OpenUnit(6,xf, indent=.false.)

 call xml_AddXMLDeclaration(xf,"UTF-8")

 call xml_NewElement(xf,"psml")
 call my_add_attribute(xf,"version","0.8")
 call my_add_attribute(xf,"energy_unit","hartree")
 call my_add_attribute(xf,"length_unit","bohr")

 call xml_NewElement(xf,"provenance")
 call my_add_attribute(xf,"creator","ONCVPSP-2.1.2+psml")
 call my_add_attribute(xf,"date",datestr)
 call xml_NewElement(xf,"input-file")
 call my_add_attribute(xf,"name","oncvpsp-input")
!
 call echo_input(xf,atsym,zz,nc,nv,iexc,psfile,na,la,fa,lmax,rc0,ep,ncon, &
&                nbas,qcut,lloc,lpopt,dvloc0,nproj,debl,icmod,fcfact, &
&                rcfact,fcgrid,rcgrid,epsh1,epsh2,depsh,rxpsh,rlmax,drl, &
&                ncnf,nvcnf,nacnf,lacnf,facnf,log_mesh)

 call xml_EndElement(xf,"input-file")
!
 call xml_EndElement(xf,"provenance")

 nameat = atsym
 ncore  = nc
 nval   = nv

 norbs = ncore + nval
 allocate(nsh(norbs), lsh(norbs), fsh(norbs))

 total_valence_charge = 0.0_dp
 ncp = ncore + 1
 do ii = 1, norbs
   nsh(ii) = na(ii)
   lsh(ii) = la(ii)
   fsh(ii) = fa(ii)
   if(ii > ncore) then
     total_valence_charge =   total_valence_charge + fsh(ii)
   end if
 end do

 npots = lmax + 1
 allocate(ll(npots), nn(npots), ff(npots))
 do ii = 1, npots
   ll(ii) = ii - 1
   found = .false.
! look for the appropriate shell in the valence
   do jj = ncp, norbs
     if(lsh(jj) == ll(ii)) then
       found = .true.
       nn(ii) = nsh(jj)
       ff(ii) = fsh(jj)
       exit
     end if
   end do
   if(.not. found) then
! generate the appropriate effective n
     nn(ii) = ll(ii) + 1
     do jj = 1, ncore
       if(lsh(jj) == ll(ii)) then
         nn(ii) = nn(ii) + 1
       end if
       ff(ii) = 0.0_dp
     end do
   end if
 end do

 psflavor ="Hamann's oncvpsp"

 polarized = .false.
 polattrib = "no"
 there_is_core = (icmod >= 1)
 if(there_is_core) then
   coreattrib = "yes"
 else
   coreattrib = "no"
 end if

! XC name handling

 select case(iexc)

 case(1)
   xcfuntype    = 'LDA'
   xcfunparam   = 'Wigner'
   libxc_id = [ XC_LDA_X, XC_LDA_C_WIGNER ]
 case(2)
   xcfuntype    = 'LDA'
   xcfunparam   = 'Hedin-Lundqvist'
   libxc_id = [ XC_LDA_X, XC_LDA_C_HL ]
 case(3)
   xcfuntype    = 'LDA'
   xcfunparam   = 'Ceperley-Alder PZ'
   libxc_id = [ XC_LDA_X, XC_LDA_C_PZ ]

 case(4)
   xcfuntype    = 'GGA'
   xcfunparam   = 'Perdew-Burke-Ernzerhof'
   libxc_id = [ XC_GGA_X_PBE, XC_GGA_C_PBE ]

 case  default

   xcfuntype    = '-----'
   xcfunparam   = '-----'
   libxc_id = [  XC_NOT_IMPL, XC_NOT_IMPL ]   ! ???

 end select
!
!
 call xml_NewElement(xf,"header")
 call my_add_attribute(xf,"atomic-label",nameat)
 call my_add_attribute(xf,"atomic-number",str(zz))
 call my_add_attribute(xf,"z-pseudo",str(zion))
 call my_add_attribute(xf,"flavor",psflavor)
 if(irel==2) then
   call my_add_attribute(xf,"relativity","dirac")
 else if(irel==1) then
   call my_add_attribute(xf,"relativity","scalar")
 else
   call my_add_attribute(xf,"relativity","no")
 end if
 call my_add_attribute(xf,"polarized",polattrib)
 call my_add_attribute(xf,"core-corrections",coreattrib)

 call xml_NewElement(xf,"exchange-correlation")
 call xml_NewElement(xf,"annotation")
 call my_add_attribute(xf,"oncvpsp-xc-code",str(iexc))
 call my_add_attribute(xf,"oncvpsp-xc-type",trim(xcfuntype))
 call my_add_attribute(xf,"oncvpsp-xc-authors",trim(xcfunparam))
 call xml_EndElement(xf,"annotation")

 call xml_NewElement(xf,"libxc-info")
 call my_add_attribute(xf,"number-of-functionals","2")
 do ii = 1, 2
   call xml_NewElement(xf,"functional")
   call my_add_attribute(xf,"name",trim(libxc_id(ii)%name))
   call my_add_attribute(xf,"type",trim(libxc_id(ii)%xc_kind%str))
   call my_add_attribute(xf,"id",str(libxc_id(ii)%code))
   call xml_EndElement(xf,"functional")
 end do
 call xml_EndElement(xf,"libxc-info")
 call xml_EndElement(xf,"exchange-correlation")
!

!
 call do_configuration()
 call xml_EndElement(xf,"header")

!AG -- decide how to handle the case of Vlocal as one of the sl pots.


!AG: decide whether to use a single mesh (Hamann's own) or have the
! projectors use another one (linear, shorter)

 npts = mmax
 allocate(vps(npts))

 rgrid => rr(:)
 chval => rho(:)
 chcore => rhomod(:,1)

 nr = npts + 1
 allocate(r0(nr))
 call add_zero_r(rgrid,rgrid,r0)
 r0(1) = 0.0_dp

 call xml_NewElement(xf,"grid")
 call my_add_attribute(xf,"npts",str(nr))

 call xml_NewElement(xf,"annotation")
 call my_add_attribute(xf,"type","oncvpsp plus r=0")
 call my_add_attribute(xf,"oncvpsp-npts",str(npts))
 call my_add_attribute(xf,"oncvpsp-r1",trim(str(rgrid(1))))
 call my_add_attribute(xf,"oncvpsp-factor",trim(str(rgrid(2)/rgrid(1))))
 call xml_EndElement(xf,"annotation")

 call xml_NewElement(xf,"grid-data")
 call xml_AddArray(xf,r0(1:nr))
 call xml_EndElement(xf,"grid-data")

 call xml_EndElement(xf,"grid")

 allocate(f0(nr))
 if(irel==2) then
   if(trim(psfile)=="psp8") then

     call xml_NewElement(xf,"semilocal-potentials")
     call my_add_attribute(xf,"set","scalar_relativistic")
! sr components
     do ii = 1, npots
       l1 = ii
! last index:  1: j=l+1/2; 2: j=l-1/2; l=0,j=0 stored in index 1
       vps(:) = ((ll(ii)+1)*vpuns(:,l1,1)+ ll(ii)*vpuns(:,l1,2)) &
&             / dble(2*ll(ii)+1)
       call add_zero_r(vps(1:npts),rgrid,f0)
       call write_psml_item(xf, class="slps", &
&                         nn=nn(ii), ll=ll(ii), &
&                         rc=rc(ii),flavor=psflavor,&
&                         ff=f0)
     end do
     call xml_EndElement(xf,"semilocal-potentials")
!
! so components
!
     call xml_NewElement(xf,"semilocal-potentials")
     call my_add_attribute(xf,"set","spin_orbit")

     do ii = 2, npots
       l1 = ii
       vps(:) = 2*(vpuns(:,l1,1) - vpuns(:,l1,2)) / dble(2*ll(ii)+1)
       call add_zero_r(vps(1:npts),rgrid,f0)
       call write_psml_item(xf, class="slps", &
&                         nn=nn(ii), ll=ll(ii), &
&                         rc=rc(ii),flavor=psflavor,&
&                         ff=f0)
     end do

     call xml_EndElement(xf,"semilocal-potentials")

   else   ! upf

     call xml_NewElement(xf,"semilocal-potentials")
     call my_add_attribute(xf,"set","lj")

     do ii = 1, npots
       l1 = ii
       if(ii == 1) then
         jval = 0.0
! last index:  1: j=l+1/2; 2: j=l-1/2; l=0,j=0 stored in index 1
         vps(:) = vpuns(:,l1,1)
         call add_zero_r(vps(1:npts),rgrid,f0)
         call write_psml_item(xf, class="slps", &
&                     nn=nn(ii), ll=ll(ii), jj=jval, &
&                     rc=rc(ii),flavor=psflavor,&
&                     ff=f0)
       else
         do jj=1,2
           jval = ll(ii) + (3-2*jj)*0.5  ! convert (1,2) to (1,-1)*1/2
! last index:  1: j=l+1/2; 2: j=l-1/2; l=0,j=0 stored in index 1
           vps(:) = vpuns(:,l1,jj)
           call add_zero_r(vps(1:npts),rgrid,f0)
           call write_psml_item(xf, class="slps", &
&                    nn=nn(ii), ll=ll(ii), jj=jval, &
&                    rc=rc(ii),flavor=psflavor,&
&                    ff=f0)
         end do
       end if
     end do
     call xml_EndElement(xf,"semilocal-potentials")

   end if ! upf vs psp8
 else
   call xml_NewElement(xf,"semilocal-potentials")
   if(irel==1) then
     call my_add_attribute(xf,"set","scalar_relativistic")
   else
     call my_add_attribute(xf,"set","non_relativistic")
   end if
   vpsd: do ii = 1, npots
     vps(:) = vpuns(:,ii,1)
     call add_zero_r(vps(1:npts),rgrid,f0)
     call write_psml_item(xf, class="slps", &
&                         nn=nn(ii), ll=ll(ii), &
&                         rc=rc(ii),flavor=psflavor,&
&                         ff=f0)
   end do vpsd
   call xml_EndElement(xf,"semilocal-potentials")
 end if
!
!--------
 call xml_NewElement(xf,"valence-charge")
 call my_add_attribute(xf,"total-charge",  &
&                    str(total_valence_charge))
 call xml_NewElement(xf,"radfunc")

 call xml_NewElement(xf,"data")
 call add_zero_r(chval(1:npts),rgrid,f0)
 call xml_AddArray(xf,f0(1:nr))
 call xml_EndElement(xf,"data")
 call xml_EndElement(xf,"radfunc")
 call xml_EndElement(xf,"valence-charge")


 if(there_is_core) then
   rcore = rr(irct)
   call xml_NewElement(xf,"pseudocore-charge")
   call my_add_attribute(xf,"matching-radius",str(rcore))
   call my_add_attribute(xf,"number-of-continuous-derivatives", &
&                                 str(4))
   call my_add_attribute(xf,"annotation",  &
&               "Monotonic 8th-order polynomial with no linear term")
   call xml_NewElement(xf,"radfunc")

   call xml_NewElement(xf,"data")
   call add_zero_r(chcore(1:npts),rgrid,f0)
   call xml_AddArray(xf,f0(1:nr))
   call xml_EndElement(xf,"data")
   call xml_EndElement(xf,"radfunc")
   call xml_EndElement(xf,"pseudocore-charge")
   nullify(chcore)
 end if

 call xml_NewElement(xf,"pseudopotential-operator")
 call my_add_attribute(xf,"version","0.1")
 call my_add_attribute(xf,"energy_unit","hartree")
 call my_add_attribute(xf,"length_unit","bohr")

 call xml_NewElement(xf,"provenance")
 call my_add_attribute(xf,"creator","oncvpsp 2.1.2+psml")
 call xml_EndElement(xf,"provenance")

 vlocal => vpuns(:,lloc+1,1)
 call xml_NewElement(xf,"local-potential")
 if(lloc > lmax) then
   call my_add_attribute(xf,"type","oncv-fit")
 else
   call my_add_attribute(xf,"type","l="//str(lloc))
 end if
 call xml_NewElement(xf,"radfunc")
 call xml_NewElement(xf,"data")
 call add_zero_r(vlocal(1:npts),rgrid,f0)
 call xml_AddArray(xf,f0(1:nr))
 call xml_EndElement(xf,"data")
 call xml_EndElement(xf,"radfunc")
 call xml_EndElement(xf,"local-potential")

 if(irel==2) then
   if(trim(psfile)=="psp8") then

! set up projector number for sr_so calculations based on non-zero coefficients
     npr_sr(:)=0
     npr_so(:)=0
     do l1=1,lmax+1
       do ii=1,2*nproj(l1)
         if(abs(esr(ii,l1))>0.0d0) npr_sr(l1)=npr_sr(l1)+1
         if(abs(eso(ii,l1))>0.0d0) npr_so(l1)=npr_so(l1)+1
       end do
     end do

     call xml_NewElement(xf,"projectors")
     call my_add_attribute(xf,"set","scalar_relativistic")

     do l1=1,lmax+1
       if(l1==lloc+1) cycle
       do jj=1,npr_sr(l1)
         call add_zero_r(vsr(:,jj,l1),rgrid,f0)
         call write_psml_item(xf, class="proj", &
&                           seq=jj, ll=l1-1, &
&                           ekb=esr(jj,l1), &
&                           type="oncv", ff=f0)
       end do
     end do
     call xml_EndElement(xf,"projectors")
!
!     Spin-orbit part
!
     call xml_NewElement(xf,"projectors")
     call my_add_attribute(xf,"set","spin_orbit")

     do l1=2,lmax+1
       if(l1==lloc+1) cycle
       do jj=1,npr_so(l1)
         call add_zero_r(vso(:,jj,l1),rgrid,f0)
         call write_psml_item(xf, class="proj", &
&                           seq=jj, ll=l1-1, &
&                           ekb=eso(jj,l1), &
&                           type="oncv", ff=f0)
       end do
     end do

     call xml_EndElement(xf,"projectors")

   else  ! upf

     call xml_NewElement(xf,"projectors")
     call my_add_attribute(xf,"set","lj")

     do l1=1,lmax+1
       if(l1==lloc+1) cycle

       if(l1 == 1) then
! l=0, only one j=0
         jval = 0.0

         do jj=1,nproj(l1)
           call add_zero_r(vkb(:,jj,l1,1),rgrid,f0)
           call write_psml_item(xf, class="proj", &
&                          seq=jj, ll=l1-1, jj=jval, &
&                          ekb=evkb(jj,l1,1), &
&                          type="oncv", ff=f0)
         end do

       else  ! l1 /=1  (l/=0)

! two j values
         do jk=1,2
           jval = l1-1 + (3-2*jk)*0.5  ! convert (1,2) to (1,-1)*1/2
           do jj=1,nproj(l1)
             call add_zero_r(vkb(:,jj,l1,jk),rgrid,f0)
             call write_psml_item(xf, class="proj", &
&                          seq=jj, ll=l1-1, jj=jval, &
&                          ekb=evkb(jj,l1,jk), &
&                          type="oncv", ff=f0)
           end do
         end do

       end if  ! l1 == 1

     end do ! over l shells
     call xml_EndElement(xf,"projectors")
!
   end if ! upf vs psp8
 else
   call xml_NewElement(xf,"projectors")
   if(irel==1) then
     call my_add_attribute(xf,"set","scalar_relativistic")
   else
     call my_add_attribute(xf,"set","non_relativistic")
   end if

   do l1=1,lmax+1
     if(l1==lloc+1) cycle
     do jj=1,nproj(l1)
       call add_zero_r(vkb(:,jj,l1,1),rgrid,f0)
       call write_psml_item(xf, class="proj", &
&                           seq=jj, ll=l1-1, &
&                           ekb=evkb(jj,l1,1), &
&                           type="oncv", ff=f0)
     end do
   end do

   call xml_EndElement(xf,"projectors")
 end if

 call xml_EndElement(xf,"pseudopotential-operator")

 call xml_EndElement(xf,"psml")


 call xml_Close(xf)
 write(6,'(/a)') 'END_PSP'

 deallocate(f0,r0)
 return

 contains
 subroutine do_configuration()

!writes the valence-configuration element of the header; it has no
!arguments and uses the host's xf, ncp, norbs, nsh, lsh, fsh, polarized,
!fdown and fup

 call xml_NewElement(xf,"valence-configuration")
 call my_add_attribute(xf,"total-valence-charge", str(total_valence_charge))
 do ii = ncp, norbs
   if(fsh(ii) < 1.0e-10_dp) cycle
   call xml_NewElement(xf,"shell")
   call my_add_attribute(xf,"n",str(nsh(ii)))
   call my_add_attribute(xf,"l",lsymb(lsh(ii)))
   call my_add_attribute(xf,"occupation",str(fsh(ii)))
   if(polarized) then
     call my_add_attribute(xf,"occupation-down",str(fdown(ii)))
     call my_add_attribute(xf,"occupation-up",str(fup(ii)))
   end if
   call xml_EndElement(xf,"shell")
 end do
 call xml_EndElement(xf,"valence-configuration")
 return
 end subroutine do_configuration

 end subroutine psmlout

 subroutine echo_input(xf,atsym,zz,nc,nv,iexc,psfile,na,la,fa,lmax,rc0,ep, &
&                      ncon,nbas,qcut,lloc,lpopt,dvloc0,nproj,debl,icmod, &
&                      fcfact,rcfact,fcgrid,rcgrid,epsh1,epsh2,depsh, &
&                      rxpsh,rlmax,drl, &
&                      ncnf,nvcnf,nacnf,lacnf,facnf,log_mesh)

!echoes the input data into the input-file element of the provenance, a
!line at a time, in the form of the psp8 echo (linout), so that it can be
!run again

!xf  xml file handle
!atsym  atomic symbol
!zz  atomic number
!iexc  type of exchange-correlation
!lmax  maximum angular momentum
!lloc  l for local potential
!nproj  number of vkb projectors for each l
!icmod  model core charge option
!drl  spacing of linear radial grid
!na  principal quantum numbers of the reference configuration states
!la  angular momenta of the reference configuration states
!ncon  number of matching constraints for each l
!nbas  number of basis functions for each l
!nvcnf  numbers of valence states of the reference and test configurations
!nacnf  principal quantum numbers of the reference and test configurations
!lacnf  angular momenta of the reference and test configurations
!nc  number of core states
!nv  number of valence states
!lpopt  type of polynomial continuation for lloc=4
!ncnf  number of test configurations
!fa  occupations of the reference configuration states
!rc0  core radii as read from the input (bohr), rc0(5) for lloc=4
!ep  energies at which the pseudopotentials are generated (Ha)
!qcut  residual-energy cutoff wave vectors for each l
!debl  energy shifts for the second projectors for each l (Ha)
!facnf  occupancies of the reference and test configurations
!dvloc0  r=0 shift of the lloc=4 local potential (Ha)
!fcfact  prefactor for model core charge amplitude, as read from the input
!rcfact  prefactor for model core charge scale (icmod=3), or error
!  tolerance (icmod=6)
!fcgrid  min, max and step of the fcfact coarse-search grid (icmod=4)
!rcgrid  min, max and step of the rcfact coarse-search grid (icmod=4)
!epsh1  lower limit of the log derivative analysis energy range (Ha)
!epsh2  upper limit of the log derivative analysis energy range (Ha)
!depsh  energy step of the log derivative analysis (Ha)
!rxpsh  radius of the log derivative analysis (bohr), negative for default
!rlmax  maximum radius of the linear output mesh (bohr)
!psfile  pseudopotential file format as read from the input
!log_mesh  amesh, rr1, rmax of the log mesh, echoed when the input sets it
!  (else 0)

 use xmlf90_wxml, only: xmlf_t, xml_AddPcData
 use input_echo_m, only: icmod_lines, epsh_lines, echo_len, echo_max
 implicit none

!Input variables
 character(len=2), intent(in) :: atsym
 character(len=4), intent(in) :: psfile
 integer, intent(in) :: nc,nv,iexc,lmax,lloc,lpopt,icmod,ncnf
 integer, intent(in) :: na(30),la(30),ncon(6),nbas(6),nproj(6)
 integer, intent(in) :: nvcnf(5),nacnf(30,5),lacnf(30,5)
 real(dp), intent(in) :: zz,fa(30),rc0(6),ep(6),qcut(6),debl(6),facnf(30,5)
 real(dp), intent(in) :: dvloc0,fcfact,epsh1,epsh2,depsh,rlmax,drl
 real(dp), intent(in) :: rcfact,fcgrid(3),rcgrid(3),rxpsh,log_mesh(3)

!Input/output variables
 type(xmlf_t), intent(inout) :: xf

!Local variables
 character(len=100) :: line
 character(len=echo_len) :: text(echo_max)
 integer :: ii,jj,l1,ntext

 call add('# ATOM AND REFERENCE CONFIGURATION')
 call add('# atsym  z   nc   nv     iexc    psfile')
 write(line,'(a,a,f6.2,2i5,i8,2a)') '  ',trim(atsym),zz,nc,nv,iexc, &
&     '      ',trim(psfile)
 call add(line)
 call add('#')
 call add('#   n    l    f')
 do ii=1,nc+nv
   write(line,'(2i5,f8.2)') na(ii),la(ii),fa(ii)
   call add(line)
 end do

 call add('#')
 call add('# PSEUDOPOTENTIAL AND OPTIMIZATION')
 call add('# lmax')
 write(line,'(i5)')  lmax
 call add(line)
 call add('#')
 call add('#   l,   rc,      ep,       ncon, nbas, qcut')
 do l1=1,lmax+1
   write(line,'(i5,2f10.5,2i5,f10.5)') l1-1,rc0(l1),ep(l1),ncon(l1), &
&          nbas(l1),qcut(l1)
   call add(line)
 end do

 call add('#')
 call add('# LOCAL POTENTIAL')
 call add('# lloc, lpopt,  rc(5),   dvloc0')
 write(line,'(2i5,f10.5,a,f10.5)') lloc,lpopt,rc0(5),'   ',dvloc0
 call add(line)

 call add('#')
 call add('# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs')
 call add('# l, nproj, debl')
 do l1=1,lmax+1
   write(line,'(2i5,f10.5)') l1-1,nproj(l1),debl(l1)
   call add(line)
 end do

 call icmod_lines(icmod,fcfact,rcfact,fcgrid,rcgrid,text,ntext)
 do ii=1,ntext
   call add(text(ii))
 end do

 call epsh_lines(epsh1,epsh2,depsh,rxpsh,text,ntext)
 do ii=1,ntext
   call add(text(ii))
 end do

 call add('#')
 call add('# OUTPUT GRID')
 call add('# rlmax, drl')
 write(line,'(2f8.2)') rlmax,drl
 call add(line)

 call add('#')
 call add('# TEST CONFIGURATIONS')
 call add('# ncnf')
 write(line,'(i5)') ncnf
 call add(line)
 call add('# nvcnf')
 call add('#   n    l    f')
 do jj=2,ncnf+1
   write(line,'(i5)') nvcnf(jj)
   call add(line)
   do ii=nc+1,nc+nvcnf(jj)
     write(line,'(2i5,f8.2)') nacnf(ii,jj),lacnf(ii,jj),facnf(ii,jj)
     call add(line)
   end do
   call add('#')
 end do

 if(log_mesh(1)>0.0d0) then
   call add('#')
   call add('# LOG MESH')
   call add('# amesh, rr1, rmax')
   write(line,'(3es24.16)') log_mesh
   call add(line)
 end if
 return

 contains

 subroutine add(text)

!adds a line of text to the input-file element

!text  the line

!Input variables
 character(len=*), intent(in) :: text
 call xml_AddPcData(xf,trim(text),line_feed=.true.)
 return
 end subroutine add

 end subroutine echo_input

 subroutine my_add_attribute(xf,name,value)

!adds an attribute to the current xml element, with trailing blanks
!removed from its value

!xf  xml file handle
!name  attribute name
!value  attribute value

 use xmlf90_wxml, only: xmlf_t, xml_AddAttribute

!Input variables
 character(len=*), intent(in)  :: name
 character(len=*), intent(in)  :: value

!Input/output variables
 type(xmlf_t), intent(inout)   :: xf

 call xml_AddAttribute(xf,name,trim(value))
 return
 end subroutine my_add_attribute

 subroutine add_zero_r(ff,rr,f0)

! Adds an r=0 element to a grid function f0, extrapolating

!ff  function on the log radial grid
!rr  log radial grid
!f0  ff with an extrapolated value at r=0 prepended

!Input variables
 real(dp), intent(in)  :: ff(:), rr(:)

!Output variables
 real(dp), intent(out) :: f0(:)

!Local variables
 integer :: ii, npts
 real(dp) :: r2

 npts = size(ff)
 if(size(f0) /= npts +1) error stop "nr /= npts + 1 in add_zero_r"
 do ii = 1, npts
   f0(ii+1) = ff(ii)
 end do
 r2 = rr(1)/(rr(2)-rr(1))
 f0(1) = ff(2) - (ff(3)-ff(2))*r2
 return

 end subroutine add_zero_r


 subroutine write_psml_item(xf,class, &
&                           nn, ll, jj, ss, &
&                           seq, &
&                           rc, ekb, &
&                           flavor, type, set, &
&                           ff)

!writes one radial function (a semi-local potential or a projector) as a
!PSML element with its attributes

!xf  xml file handle
!class  element name, "slps" or "proj"
!nn  principal quantum number (optional)
!ll  angular momentum (optional)
!jj  total angular momentum (optional)
!ss  spin, +0.5 or -0.5 (optional)
!seq  index of the projector for this l (optional)
!rc  core radius (bohr, optional)
!ekb  projector energy coefficient (Ha, optional)
!flavor  generation method (optional)
!type  projector type (optional)
!set  set the function belongs to (optional)
!ff  the function on the grid, including r=0

 use xmlf90_wxml, only: xmlf_t, str, xml_NewElement, xml_EndElement, &
&                       xml_AddArray

!Input variables
 character(len=*), intent(in)  :: class

 integer, intent(in), optional  :: nn
 integer, intent(in), optional  :: ll
 real(dp), intent(in), optional  :: jj
 real(dp), intent(in), optional  :: ss

! for sl potentials
 real(dp), intent(in), optional  :: rc
 character(len=*), intent(in), optional  :: flavor

! for projectors
 integer, intent(in), optional   :: seq
 real(dp), intent(in), optional  :: ekb
 character(len=*), intent(in), optional  :: type

 character(len=*), intent(in), optional  :: set
 real(dp), intent(in)  :: ff(:)

!Input/output variables
 type(xmlf_t), intent(inout)   :: xf

 call xml_NewElement(xf,trim(class))
 if(present(set))  call my_add_attribute(xf,"set",set)

! we might want to check input values
 if(present(nn))  call my_add_attribute(xf,"n",str(nn))
 if(present(ll))  call my_add_attribute(xf,"l",lsymb(ll))
 if(present(jj))  call my_add_attribute(xf,"j", &
&                                        str(jj,format="(f3.1)"))
! spin: +0.5 or -0.5
 if(present(ss))  call my_add_attribute(xf,"s", &
&                                        str(ss,format="(f4.1)"))

 if(present(seq)) call my_add_attribute(xf,"seq",str(seq))

 if(present(rc))  call my_add_attribute(xf,"rc",str(rc))
 if(present(ekb))  call my_add_attribute(xf,"ekb",str(ekb))

 if(present(flavor))  call my_add_attribute(xf,"flavor",flavor)
 if(present(type))  call my_add_attribute(xf,"type",type)

 call xml_NewElement(xf,"radfunc")
 call xml_NewElement(xf,"data")
 call xml_AddArray(xf,ff(:))
 call xml_EndElement(xf,"data")
 call xml_EndElement(xf,"radfunc")

 call xml_EndElement(xf,trim(class))
 return

 end subroutine write_psml_item

 end module psmlout_m
